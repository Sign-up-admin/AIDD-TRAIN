import logging
from itertools import product
import sys

# --- IMPORTS ---
from ..training.model import ViSNetPDB
from .config import TIME_RANGES, CYCLE_BATCHES
from .core.prober import probe_config, find_max_batch_size_by_stressing

try:
    from skopt import gp_minimize
    from skopt.space import Integer
except ImportError:
    print("ERROR: scikit-optimize library not found. Please install it with: pip install scikit-optimize")
    sys.exit(1)

logger = logging.getLogger("HardwareOptimizer")

def _optimize_for_production(param_combinations, max_params, start_bs, stress_iters, processed_test_data):
    logger.info("--- Strategy: Two-stage optimization for 'production' mode. ---")
    logger.info("[1/2] Searching for the highest-quality model architecture...")
    best_architecture = None
    for i, (num_layers, hidden_channels) in enumerate(param_combinations):
        logger.info(f"> Probing Arch {(i+1)}/{len(param_combinations)} (L={num_layers}, C={hidden_channels})...")
        num_params = sum(p.numel() for p in ViSNetPDB(hidden_channels=hidden_channels, num_layers=num_layers).parameters() if p.requires_grad)
        if num_params > max_params:
            logger.info(f"SKIPPED (too complex: {num_params/1e3:.0f}k > {max_params/1e3:.0f}k cap)"); continue
        base_config = {'visnet_hidden_channels': hidden_channels, 'visnet_num_layers': num_layers}
        success, _ = probe_config({**base_config, 'batch_size': 1}, processed_test_data, stress_iterations=5)
        if not success:
            logger.info("SKIPPED (too large for VRAM)"); continue
        logger.info(f"FOUND! (Model has {num_params/1e3:.0f}k params and fits VRAM)")
        best_architecture = base_config
        break
    if best_architecture:
        logger.info("[2/2] Architecture selected. Now finding max batch size...")
        optimal_bs, _ = find_max_batch_size_by_stressing(best_architecture, processed_test_data, start_bs, stress_iters)
        if optimal_bs: return {**best_architecture, 'batch_size': optimal_bs}
    logger.warning("Could not find any suitable model architecture.")
    return None

def _optimize_for_efficiency_modes(mode, num_layers_list, hidden_channels_list, max_params, start_bs, stress_iters, processed_test_data):
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
        success, _ = probe_config({**base_config, 'batch_size': 1}, processed_test_data, stress_iterations=5, prefix="	")
        if not success:
            logger.warning(f"> SKIPPED: Model architecture too large for bs=1."); return 0.0
        
        found_bs, avg_time_per_iter = find_max_batch_size_by_stressing(base_config, processed_test_data, start_bs, stress_iters)
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

def _optimize_for_smoke_test(param_combinations, max_params, processed_test_data):
    logger.info("--- Strategy: Minimal check for 'smoke_test' mode ---")
    num_layers, hidden_channels = param_combinations[-1]
    num_params = sum(p.numel() for p in ViSNetPDB(hidden_channels=hidden_channels, num_layers=num_layers).parameters() if p.requires_grad)
    if num_params > max_params: logger.warning(f"Smallest model ({num_params/1e3:.1f}k params) exceeds data cap ({max_params/1e3:.0f}k).")
    base_config = {'visnet_hidden_channels': hidden_channels, 'visnet_num_layers': num_layers, 'batch_size': 1}
    logger.info(f"--- Trying smallest model (L={num_layers}, C={hidden_channels}, BS=1) for a single iteration ---")
    success, _ = probe_config(base_config, processed_test_data, stress_iterations=1, prefix="	")
    if success: return base_config
    return None
