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
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    item = {
        'pdb_code': '1jmf',
        'protein_path': os.path.join(base_dir, '1jmf', '1jmf_protein.pdb'),
        'ligand_path': os.path.join(base_dir, '1jmf', '1jmf_ligand.sdf'),
        'binding_data': 'Kd=1.0nM'
    }

    print(f"Temporarily setting max_atoms to 10200 to process {item['pdb_code']}...")
    original_max_atoms = CONFIG.get('max_atoms')
    try:
        CONFIG['max_atoms'] = 10200
        PROCESSED_1JMF_DATA = process_item(item)
    finally:
        if original_max_atoms is not None:
            CONFIG['max_atoms'] = original_max_atoms
        elif 'max_atoms' in CONFIG:
            del CONFIG['max_atoms']

    if PROCESSED_1JMF_DATA is None:
        print("CRITICAL: Failed to process 1jmf data. Aborting.")
        exit()
    
    print(f"--- Successfully processed 1jmf. Atom count: {PROCESSED_1JMF_DATA.num_nodes} ---")

def _setup_probe_environment(config):
    """Helper to create model, batch, and optimizer for probing."""
    model = ViSNetPDB(
        hidden_channels=config['visnet_hidden_channels'],
        num_layers=config['visnet_num_layers'],
        lmax=2,
        vecnorm_type='max_min'
    ).cuda()
    
    batch_size = config['batch_size']
    data_list = [PROCESSED_1JMF_DATA] * batch_size
    batch = Batch.from_data_list(data_list).cuda()

    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    return model, batch, optimizer

def probe_config(config, stress_iterations=20, prefix=""):
    """
    Tests a given configuration using the real 1jmf data.
    """
    base_text = f"{prefix}--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations..."
    bar_length = 20
    model, batch, optimizer = None, None, None
    try:
        model, batch, optimizer = _setup_probe_environment(config)

        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(batch)
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
        if model is not None:
            del model, batch, optimizer
            gc.collect()
            torch.cuda.empty_cache()

def estimate_cycle_time(config, num_batches=450, warmup_iters=10, timed_iters=30):
    """
    Estimates the time for a training cycle of a given number of batches.
    """
    print(f"\t\t> Estimating cycle time for {num_batches} batches...")
    model, batch, optimizer = None, None, None
    try:
        model, batch, optimizer = _setup_probe_environment(config)

        # Warm-up runs
        for _ in range(warmup_iters):
            optimizer.zero_grad()
            output = model(batch)
            loss = output.sum()
            loss.backward()
            optimizer.step()
        
        torch.cuda.synchronize() # Wait for all kernels to finish

        # Timed runs
        start_time = time.time()
        for _ in range(timed_iters):
            optimizer.zero_grad()
            output = model(batch)
            loss = output.sum()
            loss.backward()
            optimizer.step()
        
        torch.cuda.synchronize() # Wait for all kernels to finish
        end_time = time.time()

        avg_time_per_iter = (end_time - start_time) / timed_iters
        estimated_total_seconds = avg_time_per_iter * num_batches
        estimated_total_minutes = estimated_total_seconds / 60

        print(f"\t\t> Avg time/batch: {avg_time_per_iter:.4f}s. Estimated cycle time: {estimated_total_minutes:.2f} mins.")
        return estimated_total_minutes

    except Exception as e:
        print(f"\n\t\t> Error during time estimation: {e}")
        return float('inf') # Return infinity on error
    finally:
        if model is not None:
            del model, batch, optimizer
            gc.collect()
            torch.cuda.empty_cache()

def get_search_space(vram_gb, mode_to_optimize):
    """
    Generates a VRAM-aware search space. For GPUs with very low VRAM (<= 6GB),
    it uses an 'ultra-conservative' set of architectures to ensure a baseline is found.
    """
    print(f"\nGenerating dedicated search space for '{mode_to_optimize}' mode...")

    # This dictionary defines the different search spaces for each mode.
    # This is a critical part of the multi-stage optimization philosophy.
    if vram_gb <= 6:
        print("--- Using ultra-conservative architecture for very low VRAM GPU (<=6GB) ---")
        arch_definitions = {
            'production':  {'layers': [4, 3, 2], 'channels': [48, 32, 24]},
            'validation':  {'layers': [3, 2], 'channels': [32, 24, 16]},
            'prototyping': {'layers': [2, 1], 'channels': [16, 8]},
            'smoke_test':  {'layers': [1], 'channels': [8]}
        }
    elif vram_gb <= 8:
        print("--- Using conservative architecture for low VRAM GPU (6-8GB) ---")
        arch_definitions = {
            'production':  {'layers': [6, 5, 4], 'channels': [96, 64, 48]},
            'validation':  {'layers': [4, 3, 2], 'channels': [64, 48, 32]},
            'prototyping': {'layers': [3, 2, 1], 'channels': [32, 24, 16]},
            'smoke_test':  {'layers': [1], 'channels': [16, 8]}
        }
    else:
        print("--- Using standard architecture for high VRAM GPU (>8GB) ---")
        arch_definitions = {
            'production':  {'layers': [8, 7, 6, 5], 'channels': [256, 192, 128]},
            'validation':  {'layers': [6, 5, 4, 3], 'channels': [128, 96, 64]},
            'prototyping': {'layers': [4, 3, 2], 'channels': [64, 48, 32]},
            'smoke_test':  {'layers': [2, 1], 'channels': [32, 24, 16]}
        }

    mode_params = {
        'production':  {'bs': 16, 'stress': 20},
        'validation':  {'bs': 32, 'stress': 20},
        'prototyping': {'bs': 64, 'stress': 20},
        'smoke_test':  {'bs': 128, 'stress': 1} # Minimal stress for smoke test
    }

    base_num_layers = arch_definitions[mode_to_optimize]['layers']
    base_hidden_channels = arch_definitions[mode_to_optimize]['channels']
    start_batch_size = mode_params[mode_to_optimize]['bs']
    stress_iterations = mode_params[mode_to_optimize]['stress']

    if vram_gb > 24: vram_factor = 1.5
    elif vram_gb > 12: vram_factor = 1.2
    elif vram_gb > 8: vram_factor = 1.1
    else: vram_factor = 1.0
    
    print(f"--- Applying VRAM scaling factor of {vram_factor} (based on {vram_gb:.2f} GB) ---")

    num_attention_heads = 8
    hidden_channels_list = sorted([c for c in list(set((int(c * vram_factor) // num_attention_heads) * num_attention_heads for c in base_hidden_channels)) if c > 0], reverse=True)
    num_layers_list = sorted([l for l in list(set(int(l * vram_factor) for l in base_num_layers)) if l > 0], reverse=True)
    start_batch_size = int(start_batch_size * vram_factor)

    if not num_layers_list or not hidden_channels_list:
        print("CRITICAL: Search space is empty. Check VRAM factor and base definitions.")
        exit()

    print(f"Search space: Layers={num_layers_list}, Channels={hidden_channels_list}")
    print(f"Starting BS: {start_batch_size}, Stress Iterations: {stress_iterations}")

    return start_batch_size, hidden_channels_list, num_layers_list, stress_iterations

def find_max_batch_size_by_stressing(base_config, start_batch_size, stress_iters):
    """
    Finds the maximum stable batch size using a robust search strategy.
    """
    oom_ceiling = -1

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

    print(f"\n--- Strategy: Final stability check for batch_size={bs_candidate} ---")
    print(f"[3/3] Stability confirmation...", end='')
    config = {**base_config, 'batch_size': bs_candidate}
    if probe_config(config, stress_iterations=20, prefix="\t"): 
        print(f"\t> Final configuration is stable.")
        return bs_candidate
    else:
        print(f"\t> Final stability check failed. This configuration is unstable.")
        # If the candidate fails, we cannot trust it. Return the next lower value, which will be 0 if candidate was 1.
        return bs_candidate - 1

def save_results(final_profile):
    """
    Saves the final, optimized hardware profile to a JSON file.
    """
    output_path = 'hardware_profile.json'
    print(f"\n--- Saving Final Results to '{output_path}' ---")
    with open(output_path, 'w') as f:
        json.dump(final_profile, f, indent=4)
    print("Update complete.")

def get_parameter_cap(dataset_size):
    """
    Returns a reasonable parameter cap based on dataset size to prevent overfitting.
    """
    if dataset_size < 1000:
        cap = 50_000  # 50k params
        level = "small"
    elif dataset_size < 10000:
        cap = 100_000 # 100k params
        level = "medium"
    else:
        cap = 250_000 # 250k params
        level = "large"
    print(f"--- Dataset size: {dataset_size} ({level}). Capping model parameters at ~{cap/1e3:.0f}k to prevent overfitting. ---")
    return cap

def find_optimal_configs(modes_to_optimize, dataset_size):
    """
    Finds the best hardware configuration for each specified mode independently.
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting.")
        return

    print("\n" + "="*80)
    print("!> WARNING: This optimizer uses a single data sample ('1jmf') for all tests.")
    print("!> The results are only reliable if '1jmf' is representative of your dataset.")
    print("!> For best results, consider replacing it with a sample of average or 95th-percentile size from your own data.")
    print("="*80 + "\n")
    prepare_real_data()

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    max_params = get_parameter_cap(dataset_size)

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
            
            start_bs, hidden_channels, num_layers, stress_iters = get_search_space(vram_gb, mode)
            
            best_config_for_this_mode = None
            param_combinations = list(product(num_layers, hidden_channels))

            if mode == 'production':
                # [CORE CONCEPT] Two-stage optimization for QUALITY. Goal: Find the highest-quality model that fits data and hardware, then maximize its throughput. This philosophy must not be removed.
                # 1. Find the highest-quality (largest) model architecture that is both data-appropriate (respects param cap) and fits the hardware (can run with bs=1).
                # 2. For that single best architecture, find the maximum possible batch size to fully utilize the remaining GPU performance.
                print(f"--- Strategy: Two-stage optimization for '{mode}' mode. ---")

                # --- Stage 1: Find the best possible architecture ---
                print("\n[1/2] Searching for the highest-quality model architecture...")
                best_architecture = None
                for i, (layers, channels) in enumerate(param_combinations):
                    print(f"\t> Probing Arch {(i+1)}/{len(param_combinations)} (L={layers}, C={channels})...", end='')
                    
                    temp_model = ViSNetPDB(hidden_channels=channels, num_layers=layers)
                    num_params = sum(p.numel() for p in temp_model.parameters() if p.requires_grad)
                    del temp_model
                    
                    if num_params > max_params:
                        print(f" SKIPPED (too complex: {num_params/1e3:.0f}k > {max_params/1e3:.0f}k cap)")
                        continue
                    
                    base_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers}
                    if not probe_config({**base_config, 'batch_size': 1}, stress_iterations=5, prefix=""):
                        print(" SKIPPED (too large for VRAM)")
                        continue

                    print(f" FOUND! (Model has {num_params/1e3:.0f}k params and fits VRAM)")
                    best_architecture = base_config
                    break

                # --- Stage 2: Find the max batch size for the chosen architecture ---
                if best_architecture:
                    print("\n[2/2] Architecture selected. Now finding max batch size to utilize remaining VRAM...")
                    optimal_bs = find_max_batch_size_by_stressing(best_architecture, start_bs, stress_iters)
                    if optimal_bs:
                        best_config_for_this_mode = {**best_architecture, 'batch_size': optimal_bs}
                else:
                    print("\nCould not find any suitable model architecture that fits both data and hardware constraints.")

            elif mode == 'prototyping' or mode == 'validation':
                # [CORE CONCEPT] Optimization for a FLEXIBLE TIME BUDGET. This philosophy must not be removed.
                # Prototyping Goal: A fast training cycle (~20 mins) for rapid iteration.
                # Validation Goal: A medium training cycle (~90 mins) to validate ideas with near-production quality.
                # This is achieved by finding the model with the highest throughput (max batch size) that completes a fixed-size training cycle (450 batches) 
                # WITHIN the mode's flexible time budget (soft target + 20 min flexibility).
                TIME_BUDGETS = {'prototyping': 20, 'validation': 90}
                TIME_FLEXIBILITY_MINUTES = 20
                CYCLE_BATCHES = 450
                time_budget = TIME_BUDGETS[mode]
                hard_time_limit = time_budget + TIME_FLEXIBILITY_MINUTES

                print(f"--- Strategy: Find MAX THROUGHPUT config with cycle time <= {time_budget} mins (+{TIME_FLEXIBILITY_MINUTES} min flex) for '{mode}' mode. ---")
                
                best_config_for_this_mode = None
                best_throughput_so_far = 0

                for i, (layers, channels) in enumerate(param_combinations):
                    print(f"\n--- Trying Model Architecture {(i+1)}/{len(param_combinations)} (L={layers}, C={channels}) ---")
                    
                    temp_model = ViSNetPDB(hidden_channels=channels, num_layers=layers)
                    num_params = sum(p.numel() for p in temp_model.parameters() if p.requires_grad)
                    del temp_model
                    
                    if num_params > max_params:
                        print(f"\t> SKIPPED: Model has {num_params/1e3:.1f}k params, exceeding the {max_params/1e3:.0f}k cap.")
                        continue
                    
                    base_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers}
                    if not probe_config({**base_config, 'batch_size': 1}, stress_iterations=5, prefix="\t"):
                        print("\t> Model architecture too large even for bs=1, skipping...")
                        continue

                    print("\t> Finding max stable batch size for this architecture...")
                    optimal_bs = find_max_batch_size_by_stressing(base_config, start_bs, stress_iters)

                    if not optimal_bs:
                        print("\t> Could not find a stable batch size. Skipping architecture.")
                        continue

                    print(f"\t> Max throughput found: BS={optimal_bs}. Checking if it fits the ~{time_budget} min time budget...")
                    candidate_config = {**base_config, 'batch_size': optimal_bs}
                    estimated_time = estimate_cycle_time(candidate_config, num_batches=CYCLE_BATCHES)

                    if estimated_time <= hard_time_limit:
                        print(f"\t> SUCCESS: Cycle time ({estimated_time:.2f}m) is within the flexible budget ({hard_time_limit}m). This is a valid candidate.")
                        if optimal_bs > best_throughput_so_far:
                            best_throughput_so_far = optimal_bs
                            best_config_for_this_mode = candidate_config
                            print(f"\t> NEW BEST for '{mode}': BS={optimal_bs}, L={layers}, C={channels} (Time: {estimated_time:.2f}m)")
                    else:
                        print(f"\t> FAILED: Max throughput cycle time ({estimated_time:.2f}m) exceeds flexible budget ({hard_time_limit}m). Discarding architecture.")

            elif mode == 'smoke_test':
                # [CORE CONCEPT] Smoke Test: Verify environment. Do NOT search. Use the absolute smallest model for a single run.
                print(f"--- Strategy: Minimal check for '{mode}' mode ---")
                layers, channels = param_combinations[-1] # Smallest is last in the list
                
                temp_model = ViSNetPDB(hidden_channels=channels, num_layers=layers)
                num_params = sum(p.numel() for p in temp_model.parameters() if p.requires_grad)
                del temp_model
                if num_params > max_params:
                     print(f"\nWARNING: The smallest model ({num_params/1e3:.1f}k params) already exceeds the data cap ({max_params/1e3:.0f}k).")

                base_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': 1}
                print(f"\n--- Trying smallest model (L={layers}, C={channels}, BS=1) for a single iteration ---")
                if probe_config(base_config, stress_iterations=1, prefix="\t"): # Override to 1 iteration
                    best_config_for_this_mode = base_config

            # --- REPORTING RESULTS ---
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
    parser = argparse.ArgumentParser(
        description="Find the optimal hardware configuration for each specified mode, respecting both hardware and data constraints.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--modes',
        nargs='+',
        default=['production', 'validation', 'prototyping', 'smoke_test'],
        choices=['production', 'validation', 'prototyping', 'smoke_test'],
        help="A list of development modes to optimize for. Runs in hierarchical order."
    )
    parser.add_argument(
        '--dataset-size',
        type=int,
        default=None,
        help="The total number of samples in the dataset.\nThis is crucial for preventing overfitting by pruning overly complex models.\nIf not provided, you will be prompted to enter it."
    )
    args = parser.parse_args()

    dataset_size = args.dataset_size
    if dataset_size is None:
        while True:
            try:
                raw_input = input("Please enter the dataset size (number of samples): ")
                dataset_size = int(raw_input)
                if dataset_size > 0:
                    break
                else:
                    print("Please enter a positive number.")
            except (ValueError, TypeError):
                print("Invalid input. Please enter a valid integer.")
    
    find_optimal_configs(args.modes, dataset_size)
