import torch
import os
import json
import gc
import time
from itertools import product
import argparse
import sys
import logging

# --- IMPORTS ---
try:
    from skopt import gp_minimize
    from skopt.space import Integer
except ImportError:
    gp_minimize, Integer = None, None
    print("ERROR: scikit-optimize library not found. Please install it with: pip install scikit-optimize")
    sys.exit(1)

from torch_geometric.data import Batch # type: ignore
from ..training.model import ViSNetPDB
from ..data.processing import process_item
from ..config import CONFIG
from .config import (
    ARCH_DEFINITIONS, VRAM_SCALING_FACTORS, MODE_PARAMS, TIME_RANGES,
    CYCLE_BATCHES, PARAMETER_CAPS, OPTIMIZATION_HIERARCHY, TEST_SAMPLE_CONFIG
)

# --- GLOBAL LOGGER SETUP ---
logger = logging.getLogger("HardwareOptimizer")
progress_logger = logging.getLogger("ProgressBar")

def setup_logging(log_level="INFO"):
    logger.setLevel(log_level)
    progress_logger.setLevel(log_level)
    logger.propagate = False
    progress_logger.propagate = False
    if not logger.handlers:
        # Ensure the main logs directory exists and place the log file inside it.
        log_dir = 'logs'
        os.makedirs(log_dir, exist_ok=True)
        log_file_path = os.path.join(log_dir, "hardware_optimizer.log")
 
        file_handler = logging.FileHandler(log_file_path, mode='w')
        file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter("%(levelname)s: %(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    if not progress_logger.handlers:
        progress_handler = logging.StreamHandler(sys.stdout)
        progress_handler.setFormatter(logging.Formatter('%(message)s'))
        progress_logger.addHandler(progress_handler)

PROCESSED_TEST_DATA = None

def prepare_real_data(project_root_path):
    global PROCESSED_TEST_DATA

    sample_name = TEST_SAMPLE_CONFIG['pdb_code']
    logger.info(f"--- Preparing Real-World Test Case: '{sample_name}' ---")

    # --- WORKAROUND for path issue ---
    protein_path = os.path.join(project_root_path, 'compass', 'optimizer', '1jmf', '1jmf_pocket.pdb')
    ligand_path = os.path.join(project_root_path, 'compass', 'optimizer', '1jmf', '1jmf_ligand.sdf')
    item = {
        'pdb_code': sample_name,
        'protein_path': protein_path,
        'ligand_path': ligand_path,
        'binding_data': TEST_SAMPLE_CONFIG['binding_data']
    }

    max_atoms_for_sample = TEST_SAMPLE_CONFIG['max_atoms']
    logger.info(f"Temporarily setting max_atoms to {max_atoms_for_sample} to process {sample_name}...")
    # TODO: Refactor `process_item` to accept a config dictionary directly,
    # which would remove this temporary and fragile patching of the global CONFIG.
    # For now, this maintains existing behavior.
    original_max_atoms = CONFIG.get('max_atoms')
    try:
        CONFIG['max_atoms'] = max_atoms_for_sample
        PROCESSED_TEST_DATA = process_item(item)
    finally:
        if original_max_atoms is not None:
            CONFIG['max_atoms'] = original_max_atoms
        elif 'max_atoms' in CONFIG:
            del CONFIG['max_atoms']

    if PROCESSED_TEST_DATA is None:
        logger.critical(f"Failed to process {sample_name} data. Aborting."); exit()

    logger.info(f"--- Successfully processed {sample_name}. Atom count: {PROCESSED_TEST_DATA.num_nodes} ---")

def _setup_probe_environment(config):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ViSNetPDB(
        hidden_channels=config['visnet_hidden_channels'],
        num_layers=config['visnet_num_layers'],
        lmax=2, vecnorm_type='max_min'
    ).to(device)
    data_list = [PROCESSED_TEST_DATA] * config['batch_size']
    batch = Batch.from_data_list(data_list).to(device)  # type: ignore
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    return model, batch, optimizer

def probe_config(config, stress_iterations=20, prefix=""):
    base_text = f"{prefix}--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations..."
    bar_length = 20
    model, batch, optimizer = None, None, None
    try:
        model, batch, optimizer = _setup_probe_environment(config)
        torch.cuda.synchronize()
        start_time = time.time()
        for i in range(stress_iterations):
            optimizer.zero_grad(); output = model(batch); loss = output.sum() * 100; loss.backward(); optimizer.step()
            progress = (i + 1) / stress_iterations
            filled_length = int(bar_length * progress)
            bar = '█' * filled_length + '-' * (bar_length - filled_length)
            progress_logger.info(f'\r{base_text} [{bar}] {progress: >4.0%}')
        torch.cuda.synchronize()
        end_time = time.time()
        avg_time_per_iter = (end_time - start_time) / stress_iterations if stress_iterations > 0 else float('inf')
        logger.info(f'\r{base_text} SUCCESS' + ' ' * (bar_length + 10))
        return True, avg_time_per_iter
    except torch.cuda.OutOfMemoryError:
        logger.warning(f'\r{base_text} FAILED (OOM)' + ' ' * (bar_length + 10))
        return False, float('inf')
    finally:
        if model is not None: del model, batch, optimizer; gc.collect(); torch.cuda.empty_cache()

def get_search_space(vram_gb, mode_to_optimize):
    logger.info(f"Generating dedicated search space for '{mode_to_optimize}' mode...")
    selected_tier = next((ARCH_DEFINITIONS[t] for t in sorted(ARCH_DEFINITIONS, key=lambda k: ARCH_DEFINITIONS[k]['vram_threshold'], reverse=True) if vram_gb >= ARCH_DEFINITIONS[t]['vram_threshold']), None)
    if not selected_tier: logger.critical("Could not find a suitable architecture definition."); exit()
    logger.info(f"--- Using '{selected_tier['description']}' ---")
    arch_def = selected_tier['architectures'][mode_to_optimize]
    vram_factor = next((VRAM_SCALING_FACTORS[vt] for vt in sorted(VRAM_SCALING_FACTORS.keys(), reverse=True) if vram_gb >= vt), 1.0)
    logger.info(f"--- Applying VRAM scaling factor of {vram_factor} (based on {vram_gb:.2f} GB) ---")
    num_attention_heads = 8
    hidden_channels_list = sorted([c for c in {((int(c * vram_factor) // num_attention_heads) * num_attention_heads) for c in arch_def['channels']} if c > 0], reverse=True)
    num_layers_list = sorted([l for l in {int(l * vram_factor) for l in arch_def['layers']} if l > 0], reverse=True)
    start_batch_size = int(MODE_PARAMS[mode_to_optimize]['bs'] * vram_factor)
    stress_iterations = MODE_PARAMS[mode_to_optimize]['stress']
    if not num_layers_list or not hidden_channels_list: logger.critical("Search space is empty."); exit()
    logger.info(f"Search space: Layers={num_layers_list}, Channels={hidden_channels_list}")
    logger.info(f"Starting BS: {start_batch_size}, Stress Iterations: {stress_iterations}")
    return start_batch_size, hidden_channels_list, num_layers_list, stress_iterations

def find_max_batch_size_by_stressing(base_config, start_batch_size, stress_iters):
    logger.debug("--- Strategy: Finding max batch size and associated cycle time ---")
    logger.debug("[1/2] Finding OOM ceiling...")
    probe_bs = start_batch_size
    while True:
        config = {**base_config, 'batch_size': probe_bs}
        success, _ = probe_config(config, stress_iterations=5, prefix="	")
        if success:
            logger.debug(f"> OK. Trying batch_size={probe_bs * 2}")
            probe_bs = probe_bs * 2 if probe_bs > 1 else 2
        else:
            oom_ceiling = probe_bs
            logger.debug(f"> OOM at batch_size={oom_ceiling}. Ceiling found."); break
    
    logger.debug(f"[2/2] Binary searching from ceiling ({oom_ceiling}) to find stable edge...")
    bs_candidate, time_at_candidate = 0, float('inf')
    low, high = 1, oom_ceiling - 1
    while low <= high:
        mid = (low + high) // 2
        if mid == 0: break
        config = {**base_config, 'batch_size': mid}
        success, avg_time = probe_config(config, stress_iterations=stress_iters, prefix="	")
        if success:
            bs_candidate, time_at_candidate = mid, avg_time
            low = mid + 1
            logger.debug(f"> Stable at batch_size={mid}. Avg time/iter: {avg_time:.4f}s. Trying higher.")
        else:
            high = mid - 1
            logger.debug(f"> OOM at batch_size={mid}. Trying lower.")

    if bs_candidate == .0:
        logger.warning("Could not find any stable batch size for this architecture.")
        return None, float('inf')
    
    logger.debug(f"> Max stable batch size found: {bs_candidate}.")
    return bs_candidate, time_at_candidate

def save_results(final_profile, project_root_path):
    output_path = os.path.join(project_root_path, 'hardware_profile.json')
    logger.info(f"--- Saving Final Results to '{output_path}' ---")
    with open(output_path, 'w') as f:
        json.dump(final_profile, f, indent=4)
    logger.info("Update complete.")

def get_parameter_cap(current_dataset_size):
    level, cap = "large", PARAMETER_CAPS[float('inf')]
    for size_thresh in sorted(PARAMETER_CAPS.keys()):
        if current_dataset_size < size_thresh:
            cap, level = PARAMETER_CAPS[size_thresh], f"small (<{size_thresh})"
            break
    logger.info(f"--- Dataset size: {current_dataset_size} ({level}). Capping model parameters at ~{cap/1e3:.0f}k to prevent overfitting. ---")
    return cap

def _optimize_for_production(param_combinations, max_params, start_bs, stress_iters):
    logger.info("--- Strategy: Two-stage optimization for 'production' mode. ---")
    logger.info("[1/2] Searching for the highest-quality model architecture...")
    best_architecture = None
    for i, (num_layers, hidden_channels) in enumerate(param_combinations):
        logger.info(f"> Probing Arch {(i+1)}/{len(param_combinations)} (L={num_layers}, C={hidden_channels})...")
        num_params = sum(p.numel() for p in ViSNetPDB(hidden_channels=hidden_channels, num_layers=num_layers).parameters() if p.requires_grad)
        if num_params > max_params:
            logger.info(f"SKIPPED (too complex: {num_params/1e3:.0f}k > {max_params/1e3:.0f}k cap)"); continue
        base_config = {'visnet_hidden_channels': hidden_channels, 'visnet_num_layers': num_layers}
        success, _ = probe_config({**base_config, 'batch_size': 1}, stress_iterations=5)
        if not success:
            logger.info("SKIPPED (too large for VRAM)"); continue
        logger.info(f"FOUND! (Model has {num_params/1e3:.0f}k params and fits VRAM)")
        best_architecture = base_config
        break
    if best_architecture:
        logger.info("[2/2] Architecture selected. Now finding max batch size...")
        optimal_bs, _ = find_max_batch_size_by_stressing(best_architecture, start_bs, stress_iters)
        if optimal_bs: return {**best_architecture, 'batch_size': optimal_bs}
    logger.warning("Could not find any suitable model architecture.")
    return None

def _optimize_for_efficiency_modes(mode, num_layers_list, hidden_channels_list, max_params, start_bs, stress_iters):
    min_time, max_time = TIME_RANGES[mode]
    logger.info(f"--- Strategy: Bayesian Optimization to find most EFFICIENT config in 'sweet spot' ({min_time}-{max_time} mins) for '{mode}' mode. ---")
    
    evaluation_cache = {}

    def evaluate_config(params):
        layer_idx, channel_idx = params
        layers = num_layers_list[layer_idx]
        channels = hidden_channels_list[channel_idx]

        logger.info(f"--- Bayesian Probe: Trying Model Architecture (L={layers}, C={channels}) ---")
        
        num_params = sum(p.numel() for p in ViSNetPDB(hidden_channels=channels, num_layers=layers).parameters() if p.requires_grad)
        if num_params > max_params:
            logger.info(f"> SKIPPED: Model has {num_params/1e3:.1f}k params, exceeding cap."); return 0.0

        base_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers}
        success, _ = probe_config({**base_config, 'batch_size': 1}, stress_iterations=5, prefix="	")
        if not success:
            logger.warning(f"> SKIPPED: Model architecture too large for bs=1."); return 0.0
        
        found_bs, avg_time_per_iter = find_max_batch_size_by_stressing(base_config, start_bs, stress_iters)
        if not found_bs:
            logger.warning("> Could not find a stable batch size for this architecture."); return 0.0

        eval_estimated_time = (avg_time_per_iter * CYCLE_BATCHES) / 60
        eval_efficiency_score = found_bs / eval_estimated_time if eval_estimated_time > 0 else 0.0
        
        # Populate the cache with detailed results
        evaluation_cache[(layers, channels)] = {
            'config': {**base_config, 'batch_size': found_bs},
            'time': eval_estimated_time,
            'efficiency': eval_efficiency_score
        }

        if eval_estimated_time > max_time * 1.2: # Allow some leeway
            logger.info(f"> PRUNING: Estimated time ({eval_estimated_time:.2f}m) too high."); return 0.0

        return -eval_efficiency_score

    search_space = [
        Integer(0, len(num_layers_list) - 1),
        Integer(0, len(hidden_channels_list) - 1),
    ]

    # The optimizer will call evaluate_config and populate the cache
    gp_minimize(
        func=evaluate_config,
        dimensions=search_space,
        n_calls=20,
        n_initial_points=5,
        random_state=1
    )

    best_in_sweet_spot, best_overall = None, None
    best_efficiency_in_sweet_spot, best_overall_efficiency = -1.0, -1.0

    # Iterate through the cache instead of re-running tests
    for eval_data in evaluation_cache.values():
        if eval_data['efficiency'] <= 0: continue

        candidate_config = eval_data['config']
        estimated_time = eval_data['time']
        efficiency_score = eval_data['efficiency']
        
        cfg_layers = candidate_config['visnet_num_layers']
        cfg_channels = candidate_config['visnet_hidden_channels']
        cfg_optimal_bs = candidate_config['batch_size']

        if efficiency_score > best_overall_efficiency:
            best_overall_efficiency, best_overall = efficiency_score, candidate_config

        if min_time <= estimated_time <= max_time:
            if efficiency_score > best_efficiency_in_sweet_spot:
                best_efficiency_in_sweet_spot, best_in_sweet_spot = efficiency_score, candidate_config
                logger.info(f"> NEW BEST in '{mode}' sweet spot: BS={cfg_optimal_bs}, L={cfg_layers}, C={cfg_channels} (Time: {estimated_time:.2f}m, Score: {efficiency_score:.2f})")

    if best_in_sweet_spot:
        logger.info("--- Final Decision: Selecting best configuration from the 'sweet spot'. ---")
        return best_in_sweet_spot
    if best_overall:
        logger.info("--- Final Decision: No config hit the 'sweet spot'. Falling back to best overall efficiency. ---")
        return best_overall
    
    return None

def _optimize_for_smoke_test(param_combinations, max_params):
    logger.info("--- Strategy: Minimal check for 'smoke_test' mode ---")
    num_layers, hidden_channels = param_combinations[-1]
    num_params = sum(p.numel() for p in ViSNetPDB(hidden_channels=hidden_channels, num_layers=num_layers).parameters() if p.requires_grad)
    if num_params > max_params: logger.warning(f"Smallest model ({num_params/1e3:.1f}k params) exceeds data cap ({max_params/1e3:.0f}k).")
    base_config = {'visnet_hidden_channels': hidden_channels, 'visnet_num_layers': num_layers, 'batch_size': 1}
    logger.info(f"--- Trying smallest model (L={num_layers}, C={hidden_channels}, BS=1) for a single iteration ---")
    success, _ = probe_config(base_config, stress_iterations=1, prefix="	")
    if success: return base_config
    return None

def find_optimal_configs(modes_to_optimize, current_dataset_size, project_root_path):
    if not torch.cuda.is_available(): logger.critical("CUDA is not available. Aborting."); return
    logger.warning(f"This optimizer uses a single data sample ('{TEST_SAMPLE_CONFIG['pdb_code']}') for all tests. For best results, ensure this sample is representative of your dataset.")
    prepare_real_data(project_root_path)
    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    logger.info(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")
    max_params = get_parameter_cap(current_dataset_size)
    output_path = os.path.join(project_root_path, 'hardware_profile.json')
    final_profile = {}
    if os.path.exists(output_path):
        try:
            with open(output_path, 'r') as f:
                final_profile = json.load(f)
        except (json.JSONDecodeError, IOError): pass
    modes_to_run = [m for m in OPTIMIZATION_HIERARCHY if m in modes_to_optimize]
    
    optimizer_map = {
        'production': _optimize_for_production,
        'validation': _optimize_for_efficiency_modes,
        'prototyping': _optimize_for_efficiency_modes,
        'smoke_test': _optimize_for_smoke_test
    }

    try:
        for mode in modes_to_run:
            logger.info(f"================ Optimizing for: {mode.upper()} ================")
            start_bs, hidden_channels, num_layers, stress_iters = get_search_space(vram_gb, mode)
            param_combinations = list(product(num_layers, hidden_channels))
            
            optimizer_func = optimizer_map[mode]
            
            if mode == 'production':
                optimizer_args = (param_combinations, max_params, start_bs, stress_iters)
            elif mode in ['validation', 'prototyping']:
                optimizer_args = (mode, num_layers, hidden_channels, max_params, start_bs, stress_iters)
            else: # smoke_test
                optimizer_args = (param_combinations, max_params)

            best_config_for_this_mode = optimizer_func(*optimizer_args)

            if best_config_for_this_mode:
                final_profile[mode] = best_config_for_this_mode
                logger.info(f">>> Found optimal configuration for '{mode}'!")
                pretty_json = json.dumps(best_config_for_this_mode, indent=4)
                for line in pretty_json.split('\n'):
                    logger.info(line)
            else:
                logger.error(f"--- Optimization FAILED for '{mode}'. No stable configuration found. ---")

    except KeyboardInterrupt:
        logger.warning("--- User interrupted optimization. Saving best configuration found so far. ---")
    
    if final_profile: save_results(final_profile, project_root_path)
    else: logger.warning("No valid configurations were found.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Find the optimal hardware configuration for each specified mode.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--modes', nargs='+', default=OPTIMIZATION_HIERARCHY, choices=OPTIMIZATION_HIERARCHY,
        help="A list of development modes to optimize for."
    )
    parser.add_argument(
        '--dataset-size', type=int, default=None,
        help="The total number of samples in the dataset.\nIf not provided, you will be prompted to enter it."
    )
    parser.add_argument(
        '--log-level', type=str, default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level."
    )
    cli_args = parser.parse_args()

    setup_logging(cli_args.log_level)

    # Determine the project root from the script's location within the package structure
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

    dataset_size = cli_args.dataset_size
    if dataset_size is None:
        while True:
            try:
                raw_input = input("Please enter the dataset size (number of samples): ")
                dataset_size = int(raw_input)
                if dataset_size > 0: break
                else: print("Please enter a positive number.")
            except (ValueError, TypeError):
                print("Invalid input. Please enter a valid integer.")
    
    find_optimal_configs(cli_args.modes, dataset_size, project_root)
