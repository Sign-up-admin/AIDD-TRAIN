import torch
import os
import json
import gc
import time
from itertools import product
import argparse
import sys

# --- PATH CORRECTION ---
# Add project root to the Python path to resolve 'src' module not found
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
# --- END PATH CORRECTION ---

# --- Real Model and Data Processing Imports ---
from torch_geometric.data import Batch
from src.model import ViSNetPDB
from src.data_processing import process_item
from config import CONFIG

# --- Global variable for the processed 1jmf data ---
# We process it once and reuse it to avoid I/O overhead in each probe.
PROCESSED_1JMF_DATA = None

def prepare_real_data():
    """
    Loads and processes the 1jmf PDB data to be used for all stress tests.
    This function will run only once.
    """
    global PROCESSED_1JMF_DATA
    if PROCESSED_1JMF_DATA is not None:
        return

    print("\n--- Preparing 1JMF Real-World Test Case ---")
    # We need to create a mock 'item' dictionary that process_item expects.
    # The paths must be absolute for robustness.
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    item = {
        'pdb_code': '1jmf',
        'protein_path': os.path.join(base_dir, '1jmf', '1jmf_protein.pdb'),
        'ligand_path': os.path.join(base_dir, '1jmf', '1jmf_ligand.sdf'), # Assuming SDF, will need to verify
        'binding_data': 'Kd=1.0nM' # Dummy value, not used for hardware test
    }

    # Override the max_atoms check in the config for this specific case
    # to ensure our large molecule can be processed.
    print(f"Temporarily setting max_atoms to 20000 to process {item['pdb_code']}...")
    original_max_atoms = CONFIG.get('max_atoms')
    CONFIG['max_atoms'] = 20000 

    PROCESSED_1JMF_DATA = process_item(item)

    # Restore the original config value
    if original_max_atoms is not None:
        CONFIG['max_atoms'] = original_max_atoms
    else:
        del CONFIG['max_atoms'] # Clean up if it wasn't there before

    if PROCESSED_1JMF_DATA is None:
        print("CRITICAL: Failed to process 1jmf data. Aborting.")
        exit()
    
    print(f"--- Successfully processed 1jmf. Atom count: {PROCESSED_1JMF_DATA.num_nodes} ---")

def probe_config(config, stress_iterations=20, prefix=""):
    """
    Tests a given configuration using the real 1jmf data.
    """
    base_text = f"{prefix}--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations..."
    bar_length = 20
    try:
        # --- Use the Real ViSNetPDB Model ---
        model = ViSNetPDB(
            hidden_channels=config['visnet_hidden_channels'],
            num_layers=config['visnet_num_layers'],
            lmax=2, # Keep other params as needed
            vecnorm_type='max_min'
        ).cuda()
        
        # --- Create a batch from the real 1jmf data ---
        batch_size = config['batch_size']
        data_list = [PROCESSED_1JMF_DATA] * batch_size
        batch = Batch.from_data_list(data_list).cuda()

        optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(batch)
            # We need a scalar for the loss
            loss = output.sum() * 100
            loss.backward()
            optimizer.step()
            
            progress = (i + 1) / stress_iterations
            filled_length = int(bar_length * progress)
            bar = '█' * filled_length + '-' * (bar_length - filled_length)
            print(f'\r{base_text} [{bar}] {progress: >4.0%}', end='', flush=True)

        print(f'\r{base_text} SUCCESS' + ' ' * (bar_length + 10))
        return True
    except torch.cuda.OutOfMemoryError:
        print(f'\r{base_text} FAILED (OOM)' + ' ' * (bar_length + 10))
        return False
    finally:
        if 'model' in locals():
            del model, batch, optimizer
            gc.collect()
            torch.cuda.empty_cache()

def get_search_space(vram_gb, mode_to_optimize):
    """
    Generates a mode-specific search space for model architecture and batch size.
    Each mode has a completely independent set of architectures to try, ensuring
    that, for example, 'prototyping' finds a genuinely small and fast model,
    while 'production' searches through powerful, larger models.
    """
    print(f"\nGenerating dedicated search space for '{mode_to_optimize}' mode...")

    # Base definitions for different levels of model complexity
    arch_definitions = {
        'production': {
            'layers': [8, 7, 6, 5, 4],
            'channels': [256, 192, 128, 96]
        },
        'validation': {
            'layers': [6, 5, 4, 3],
            'channels': [128, 96, 64]
        },
        'prototyping': {
            'layers': [4, 3, 2],
            'channels': [64, 48, 32]
        },
        'smoke_test': {
            'layers': [2, 1],
            'channels': [32, 24, 16]
        }
    }

    # Base definitions for batch sizes and stress test intensity
    mode_params = {
        'production':  {'bs': 16, 'stress': 35},
        'validation':  {'bs': 32, 'stress': 35},
        'prototyping': {'bs': 64, 'stress': 25},
        'smoke_test':  {'bs': 128, 'stress': 15}
    }

    base_num_layers = arch_definitions[mode_to_optimize]['layers']
    base_hidden_channels = arch_definitions[mode_to_optimize]['channels']
    
    start_batch_size = mode_params[mode_to_optimize]['bs']
    stress_iterations = mode_params[mode_to_optimize]['stress']

    # --- VRAM Scaling: Adjust parameters based on available VRAM ---
    # This keeps the relative differences between modes but scales them to the GPU's capacity.
    if vram_gb > 20: vram_factor = 1.8
    elif vram_gb > 10: vram_factor = 1.4
    elif vram_gb > 7: vram_factor = 1.1
    else: vram_factor = 1.0
    
    print(f"--- Applying VRAM scaling factor of {vram_factor} (based on {vram_gb:.2f} GB) ---")

    num_attention_heads = 8  # ViSNet default

    # Scale and clean the architecture search space
    hidden_channels_list = [
        (int(c * vram_factor) // num_attention_heads) * num_attention_heads
        for c in base_hidden_channels
    ]
    hidden_channels_list = sorted([c for c in list(set(hidden_channels_list)) if c > 0], reverse=True)
    
    num_layers_list = sorted(list(set(int(l * vram_factor) for l in base_num_layers)), reverse=True)
    
    # Scale the starting batch size
    start_batch_size = int(start_batch_size * vram_factor)

    print(f"Search space: Layers={num_layers_list}, Channels={hidden_channels_list}")
    print(f"Starting BS: {start_batch_size}, Stress Iterations: {stress_iterations}")

    return start_batch_size, hidden_channels_list, num_layers_list, stress_iterations

def find_max_batch_size_by_stressing(base_config, start_batch_size, stress_iters):
    """
    Finds the maximum stable batch size by first forcing an OOM error to find the ceiling,
    then using binary search to find the highest stable value, and finally confirming stability.
    """
    oom_ceiling = -1

    # Phase 1: Find the OOM ceiling by increasing batch size until it fails.
    print("--- Strategy: Forcing OOM to find VRAM ceiling ---")
    probe_bs = start_batch_size
    while True:
        print(f"[1/3] Finding ceiling...", end='')
        config = {**base_config, 'batch_size': probe_bs}
        if probe_config(config, stress_iterations=5, prefix="\t"):
            print(f"\t> OK. Trying batch_size={probe_bs * 2}")
            if probe_bs == 1: probe_bs = 2
            else: probe_bs *= 2
        else:
            oom_ceiling = probe_bs
            print(f"\t> OOM at batch_size={oom_ceiling}. Ceiling found.")
            break

    # Phase 2: Use binary search to find the highest stable batch size.
    print(f"\n--- Strategy: Binary searching from ceiling ({oom_ceiling}) to find stable edge ---")
    bs_candidate = 0
    low, high = 1, oom_ceiling - 1

    while low <= high:
        mid = (low + high) // 2
        if mid == 0: break
        print(f"[2/3] Probing edge...", end='')
        config = {**base_config, 'batch_size': mid}
        if probe_config(config, stress_iterations=stress_iters, prefix="\t"):
            bs_candidate = mid
            low = mid + 1
            print(f"\t> Stable at batch_size={mid}. Trying higher.")
        else:
            high = mid - 1
            print(f"\t> OOM at batch_size={mid}. Trying lower.")

    if bs_candidate == 0:
        print("Could not find a stable configuration. This model may be too large for VRAM.")
        return None
    
    print(f"\t> Stable edge found at batch_size={bs_candidate}.")

    # Phase 3: Confirm the found batch size is stable.
    print(f"\n--- Strategy: Final stability check for batch_size={bs_candidate} ---")
    print(f"[3/3] Stability confirmation...", end='')
    config = {**base_config, 'batch_size': bs_candidate}
    # PRACTICAL ADJUSTMENT: The final stability check is important, but doesn't need to be excessively long.
    if probe_config(config, stress_iterations=35, prefix="\t"): # CORRECTED: Fixed to 35 iterations
        print(f"\t> Final configuration is stable.")
        return bs_candidate
    else:
        print(f"\t> Stability check failed. Reducing batch size as a precaution.")
        return max(1, bs_candidate - 1)

def save_results(final_profile):
    """
    Saves the final, optimized hardware profile to a JSON file.
    """
    output_path = 'hardware_profile.json'
    print(f"\n--- Saving Final Results to '{output_path}' ---")
    with open(output_path, 'w') as f:
        json.dump(final_profile, f, indent=4)
    print("Update complete.")

def find_optimal_configs(modes_to_optimize):
    """
    Finds the best hardware configuration for each specified mode independently.
    It iterates through each mode, gets a dedicated search space for it, and
    then performs a full search for the best-performing, stable model architecture
    and its corresponding maximum batch size.
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting.")
        return

    prepare_real_data()

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    output_path = 'hardware_profile.json'
    final_profile = {}
    if os.path.exists(output_path):
        try:
            with open(output_path, 'r') as f: final_profile = json.load(f)
        except (json.JSONDecodeError, IOError): pass

    optimization_hierarchy = ['production', 'validation', 'prototyping', 'smoke_test']
    modes_to_run = [m for m in optimization_hierarchy if m in modes_to_optimize]

    try:
        for mode in modes_to_run:
            print(f"\n================ Optimizing for: {mode.upper()} ================")
            
            # Each mode gets its own, independent search space. No more inheritance.
            start_bs, hidden_channels, num_layers, stress_iters = get_search_space(vram_gb, mode)
            
            best_config_for_this_mode = None

            # Perform a full, independent search for the best architecture for this mode.
            print(f"--- Strategy: Full search for best model and batch size for '{mode}' mode ---")
            param_combinations = list(product(num_layers, hidden_channels))
            
            for i, (layers, channels) in enumerate(param_combinations):
                print(f"\n--- Trying Model Architecture {(i+1)}/{len(param_combinations)} (L={layers}, C={channels}) ---")
                base_config = {
                    'visnet_hidden_channels': channels, 
                    'visnet_num_layers': layers
                }
                
                # First, check if the model architecture is viable at all (with batch size 1)
                if not probe_config({**base_config, 'batch_size': 1}, stress_iterations=5, prefix="\t"):
                    print("\t> Model architecture too large even for bs=1, skipping...")
                    continue

                # If it's viable, find the max batch size for this specific architecture
                optimal_bs = find_max_batch_size_by_stressing(base_config, start_bs, stress_iters)
                
                if optimal_bs:
                    # This is the largest, stable architecture found so far for this mode.
                    # Because we search from largest to smallest, we can stop here.
                    best_config_for_this_mode = {**base_config, 'batch_size': optimal_bs}
                    break 
            
            if best_config_for_this_mode:
                final_profile[mode] = best_config_for_this_mode
                print(f"\n>>> Found optimal configuration for '{mode}'!")
                print(json.dumps(best_config_for_this_mode, indent=4))
            else:
                print(f"\n--- Optimization FAILED for '{mode}'. No stable configuration found. ---")

    except KeyboardInterrupt:
        print("\n\n--- User interrupted optimization. Saving best configuration found so far. ---")

    if final_profile:
        save_results(final_profile)
    else:
        print("\nNo valid configurations were found.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find the optimal hardware configuration for each specified mode.")
    parser.add_argument(
        '--modes',
        nargs='+',
        default=['production', 'validation', 'prototyping', 'smoke_test'],
        choices=['production', 'validation', 'prototyping', 'smoke_test'],
        help="A list of development modes to optimize for. Runs in hierarchical order."
    )
    args = parser.parse_args()
    find_optimal_configs(args.modes)
