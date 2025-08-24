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

def probe_config(config, stress_iterations=20):
    """
    Tests a given configuration for a specific number of iterations to ensure stability.
    """
    print(f"--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations...", end='')
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

        last_progress_print_time = time.time()
        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(x, pos, batch)
            scaled_output = output * 100
            scaled_output.backward()
            optimizer.step()
            
            current_time = time.time()
            if current_time - last_progress_print_time > 2:
                print(".", end='', flush=True)
                last_progress_print_time = current_time

        print(f" SUCCESS")
        return True
    except torch.cuda.OutOfMemoryError:
        print(" FAILED (OOM)")
        return False
    finally:
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

    if vram_gb > 30: vram_factor, stress_factor = 2.0, 1.5
    elif vram_gb > 16: vram_factor, stress_factor = 1.5, 1.2
    else: vram_factor, stress_factor = 1.0, 1.0

    hidden_channels_list = sorted([int(c * vram_factor) for c in base_hidden_channels], reverse=True)
    num_layers_list = sorted(list(set(int(l * vram_factor) for l in base_num_layers)), reverse=True)

    if mode_to_optimize == 'validation':
        batch_sizes = [int(bs * vram_factor) for bs in [64, 48, 32, 24, 16, 12, 8, 4, 2, 1]]
        max_atoms_target = int(5000 * vram_factor)
        stress_iterations = int(30 * stress_factor)
    elif mode_to_optimize == 'prototyping':
        batch_sizes = [int(bs * vram_factor) for bs in [128, 96, 64, 48, 32, 24, 16, 8]]
        max_atoms_target = int(1024 * vram_factor)
        stress_iterations = int(20 * stress_factor)
    elif mode_to_optimize == 'smoke_test':
        batch_sizes = [int(bs * vram_factor) for bs in [256, 128, 64, 32, 16, 8]]
        max_atoms_target = 128
        stress_iterations = 10
    else:  # production
        batch_sizes = [int(bs * vram_factor) for bs in [32, 24, 16, 12, 8, 6, 4, 2, 1]]
        max_atoms_target = int(10000 * vram_factor)
        stress_iterations = int(50 * stress_factor)

    batch_sizes = sorted(list(set(bs for bs in batch_sizes if bs > 0)), reverse=True)
    return batch_sizes, hidden_channels_list, num_layers_list, max_atoms_target, stress_iterations

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
            batch_sizes, hidden_channels, num_layers, max_atoms, stress_iters = get_search_space(vram_gb, mode)
            
            if mode == 'smoke_test' and 'prototyping' in final_profile:
                print("--- Strategy: Inheriting from 'prototyping' config. ---")
                final_profile[mode] = final_profile['prototyping']
                print(f">>> Directly adopted config for '{mode}'!")
                continue

            seed_config = final_profile.get('production') or final_profile.get('validation')
            best_config_for_this_mode = None

            if mode in ['validation', 'prototyping'] and seed_config:
                print(f"--- Strategy: Seeding with a higher-tier config to find max batch size ---")
                seeded_model = {'visnet_hidden_channels': seed_config['visnet_hidden_channels'], 'visnet_num_layers': seed_config['visnet_num_layers']}
                for i, bs in enumerate(batch_sizes):
                    print(f"[{(i+1):>2}/{len(batch_sizes)}] ", end='')
                    current_config = {**seeded_model, 'batch_size': bs, 'max_atoms': max_atoms}
                    if probe_config(current_config, stress_iterations=stress_iters):
                        best_config_for_this_mode = current_config
                        break
            
            if not best_config_for_this_mode:
                # VRAM-Aware Search Strategy
                search_from_smallest = (vram_gb < 10)
                if search_from_smallest:
                    print(f"--- Strategy: Full search (small-to-large) for low-VRAM GPU ---")
                    param_combinations = list(product(sorted(num_layers), sorted(hidden_channels), sorted(batch_sizes, reverse=True)))
                else:
                    print(f"--- Strategy: Full search (large-to-small) for high-VRAM GPU ---")
                    param_combinations = list(product(num_layers, hidden_channels, batch_sizes))
                
                for i, (layers, channels, bs) in enumerate(param_combinations):
                    print(f"[{(i+1):>3}/{len(param_combinations)}] ", end='')
                    current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': bs, 'max_atoms': max_atoms}
                    if probe_config(current_config, stress_iterations=stress_iters):
                        # For small-to-large search, we save the current best and continue searching for better ones
                        if search_from_smallest:
                            best_config_for_this_mode = current_config
                        else: # For large-to-small, the first success is the best one
                            best_config_for_this_mode = current_config
                            break
            
            if best_config_for_this_mode:
                final_profile[mode] = best_config_for_this_mode
                print(f"\n>>> Found optimal configuration for '{mode}'!")
                print(json.dumps(best_config_for_this_mode, indent=4))
            else:
                print(f"\n--- Optimization Failed for '{mode}'. No stable configuration found. ---")
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
