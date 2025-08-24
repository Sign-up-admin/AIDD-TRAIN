import torch
import os
import json
import gc
import time
from itertools import product
import argparse

class MockViSNet(torch.nn.Module):
    def __init__(self, hidden_channels, num_layers, lmax, vecnorm_type):
        super().__init__()
        self.layers = torch.nn.ModuleList([
            torch.nn.Linear(hidden_channels, hidden_channels) for _ in range(num_layers)
        ])
        self.embedding = torch.nn.Embedding(100, hidden_channels)

    def forward(self, x, pos, batch):
        h = self.embedding(x)
        for layer in self.layers:
            h = layer(h)
        return torch.sum(h)

def probe_config(config, stress_iterations=20, prefix=""):
    """
    Tests a given configuration for a specific number of iterations to ensure stability.
    """
    base_text = f"{prefix}--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations..."
    bar_length = 20  # Define bar_length here to ensure it's always available
    try:
        model = MockViSNet(
            hidden_channels=config['visnet_hidden_channels'],
            num_layers=config['visnet_num_layers'],
            lmax=2,
            vecnorm_type='max_min'
        ).cuda()
        num_atoms = config['max_atoms']
        batch_size = config['batch_size']
        x = torch.randint(0, 100, (num_atoms * batch_size,)).cuda()
        pos = torch.randn(num_atoms * batch_size, 3).cuda()
        batch = torch.repeat_interleave(torch.arange(batch_size), num_atoms).cuda()
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)

        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(x, pos, batch)
            scaled_output = output * 100
            scaled_output.backward()
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
        # Ensure cleanup happens even if there's an error
        if 'model' in locals():
            del model, x, pos, batch, optimizer
            gc.collect()
            torch.cuda.empty_cache()

def get_search_space(vram_gb, mode_to_optimize):
    """
    Generates a dynamic search space, including stress iterations, based on VRAM and mode.
    """
    print(f"\nGenerating search space for {vram_gb:.2f} GB VRAM and '{mode_to_optimize}' mode...")
    base_hidden_channels = [256, 192, 128, 96, 64]
    base_num_layers = [8, 7, 6, 5, 4, 3]

    if vram_gb > 20: vram_factor, stress_factor = 2.0, 1.5
    elif vram_gb > 10: vram_factor, stress_factor = 1.5, 1.2
    elif vram_gb > 5: vram_factor, stress_factor = 1.2, 1.1
    else: vram_factor, stress_factor = 1.0, 1.0

    hidden_channels_list = sorted([int(c * vram_factor) for c in base_hidden_channels], reverse=True)
    num_layers_list = sorted(list(set(int(l * vram_factor) for l in base_num_layers)), reverse=True)

    if mode_to_optimize == 'validation':
        start_batch_size = int(32 * vram_factor)
        max_atoms_target = int(5000 * vram_factor)
        stress_iterations = int(35 * stress_factor)
    elif mode_to_optimize == 'prototyping':
        start_batch_size = int(64 * vram_factor)
        max_atoms_target = int(1024 * vram_factor)
        stress_iterations = int(25 * stress_factor)
    elif mode_to_optimize == 'smoke_test':
        start_batch_size = int(128 * vram_factor)
        max_atoms_target = 128
        stress_iterations = 15
    else:  # production
        start_batch_size = int(16 * vram_factor)
        max_atoms_target = int(10000 * vram_factor)
        stress_iterations = int(60 * stress_factor)

    return start_batch_size, hidden_channels_list, num_layers_list, max_atoms_target, stress_iterations

def find_max_batch_size_by_stressing(base_config, start_batch_size, stress_iters):
    """
    Finds the maximum stable batch size by first forcing an OOM error to find the ceiling,
    then stepping down to find the highest stable value, and finally confirming stability.
    """
    bs = start_batch_size
    oom_ceiling = -1

    # Phase 1: Find the OOM ceiling by increasing batch size until it fails.
    print("--- Strategy: Forcing OOM to find VRAM ceiling ---")
    while True:
        print(f"[1/3] Finding ceiling...", end='')
        config = {**base_config, 'batch_size': bs}
        # Use a short iteration count just to find the OOM point quickly
        if probe_config(config, stress_iterations=5, prefix="\t"):
            print(f"\t> OK. Trying batch_size={bs*2}")
            if bs == 1: bs = 2
            else: bs *= 2
        else:
            oom_ceiling = bs
            print(f"\t> OOM at batch_size={oom_ceiling}. Ceiling found.")
            break

    # Phase 2: Step down from the ceiling to find the highest stable batch size.
    print(f"\n--- Strategy: Stepping down from ceiling ({oom_ceiling}) to find stable edge ---")
    bs_candidate = 0
    for current_bs in range(oom_ceiling - 1, 0, -1):
        print(f"[2/3] Probing edge...", end='')
        config = {**base_config, 'batch_size': current_bs}
        if probe_config(config, stress_iterations=stress_iters, prefix="\t"):
            bs_candidate = current_bs
            print(f"\t> Stable edge found at batch_size={bs_candidate}.")
            break

    if bs_candidate == 0:
        print("Could not find a stable configuration. This model may be too large for VRAM.")
        return None

    # Phase 3: Confirm the found batch size is stable with a more rigorous test.
    print(f"\n--- Strategy: Final stability check for batch_size={bs_candidate} ---")
    print(f"[3/3] Stability confirmation...", end='')
    config = {**base_config, 'batch_size': bs_candidate}
    if probe_config(config, stress_iterations=stress_iters + 20, prefix="\t"):
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
    Finds the best hardware configuration using a hierarchical, VRAM-aware strategy.
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting.")
        return

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    output_path = 'hardware_profile.json'
    final_profile = {}
    if os.path.exists(output_path):
        try:
            with open(output_path, 'r') as f:
                final_profile = json.load(f)
        except (json.JSONDecodeError, IOError):
            pass

    optimization_hierarchy = ['production', 'validation', 'prototyping', 'smoke_test']
    modes_to_run = [m for m in optimization_hierarchy if m in modes_to_optimize]

    try:
        for mode in modes_to_run:
            print(f"\n================ Optimizing for: {mode.upper()} ================")
            start_bs, hidden_channels, num_layers, max_atoms, stress_iters = get_search_space(vram_gb, mode)
            
            best_config_for_this_mode = None

            seed_config = final_profile.get('production') or final_profile.get('validation')
            if mode in ['validation', 'prototyping', 'smoke_test'] and seed_config:
                print(f"--- Inheriting model architecture from '{ 'production' if 'production' in final_profile else 'validation'}' ---")
                base_config = {
                    'visnet_hidden_channels': seed_config['visnet_hidden_channels'],
                    'visnet_num_layers': seed_config['visnet_num_layers'],
                    'max_atoms': max_atoms
                }
                optimal_bs = find_max_batch_size_by_stressing(base_config, start_bs, stress_iters)
                if optimal_bs:
                    best_config_for_this_mode = {**base_config, 'batch_size': optimal_bs}
            else:
                print(f"--- Strategy: Full search for best model and batch size ---")
                param_combinations = list(product(num_layers, hidden_channels))
                for i, (layers, channels) in enumerate(param_combinations):
                    print(f"\n--- Trying Model Architecture {(i+1)}/{len(param_combinations)} (L={layers}, C={channels}) ---")
                    base_config = {
                        'visnet_hidden_channels': channels, 
                        'visnet_num_layers': layers, 
                        'max_atoms': max_atoms
                    }
                    if not probe_config({**base_config, 'batch_size': 1}, stress_iterations=5, prefix="\t"):
                        print("\t> Model architecture too large for bs=1, skipping...")
                        continue

                    optimal_bs = find_max_batch_size_by_stressing(base_config, start_bs, stress_iters)
                    if optimal_bs:
                        best_config_for_this_mode = {**base_config, 'batch_size': optimal_bs}
                        break
            
            if best_config_for_this_mode:
                final_profile[mode] = best_config_for_this_mode
                print(f"\n>>> Found optimal configuration for '{mode}'!")
                print(json.dumps(best_config_for_this_mode, indent=4))
            else:
                print(f"\n--- Optimization Failed for '{mode}'. No stable configuration found. ---")
                if mode == 'production':
                    print("Aborting further optimization as it depends on this result.")
                    break

    except KeyboardInterrupt:
        print("\n\n--- User interrupted optimization. Saving best configuration found so far. ---")

    if final_profile:
        save_results(final_profile)
    else:
        print("\nNo valid configurations were found.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find the optimal hardware configuration using a hierarchical, VRAM-aware strategy.")
    parser.add_argument(
        '--modes',
        nargs='+',
        default=['production', 'validation', 'prototyping', 'smoke_test'],
        choices=['production', 'validation', 'prototyping', 'smoke_test'],
        help="A list of development modes to optimize for. Runs in hierarchical order."
    )
    args = parser.parse_args()
    find_optimal_configs(args.modes)
