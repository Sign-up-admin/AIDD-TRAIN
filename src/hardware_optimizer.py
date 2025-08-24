import torch
import os
import json
import gc
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

def probe_config(config, stress_iterations=100):
    """
    Tests a given configuration for OOM errors under sustained load to ensure thermal and power stability.
    Now provides clear, non-optional progress feedback.
    """
    print(f"--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3}...", end='')
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

        # Sustained load simulation with clear progress
        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(x, pos, batch)
            scaled_output = output * 100
            scaled_output.backward()
            optimizer.step()
            if (i + 1) % (max(1, stress_iterations // 5)) == 0:
                print(".", end='', flush=True)

        print(" SUCCESS")
        return True
    except torch.cuda.OutOfMemoryError:
        print(" FAILED (OOM)")
        return False
    finally:
        del model, x, pos, batch, optimizer
        gc.collect()
        torch.cuda.empty_cache()

def find_optimal_config(mode_to_optimize):
    """
    Finds the best hardware configuration by dynamically creating a search space
    based on available VRAM and then using an efficient search strategy.
    Catches KeyboardInterrupt to save progress gracefully.
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting hardware optimization.")
        return

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    # --- 1. Dynamic Search Space Generation ---
    print(f"Generating search space for {vram_gb:.2f} GB VRAM and '{mode_to_optimize}' mode...")
    base_hidden_channels = [256, 192, 128, 96, 64]
    base_num_layers = [8, 7, 6, 5, 4, 3]

    if vram_gb > 30: # For A100/H100
        vram_factor, stress_factor = 2.0, 1.5
    elif vram_gb > 16: # For 3090/4090
        vram_factor, stress_factor = 1.5, 1.2
    else:
        vram_factor, stress_factor = 1.0, 1.0

    hidden_channels_list = sorted([int(c * vram_factor) for c in base_hidden_channels], reverse=True)
    num_layers_list = sorted(list(set(int(l * vram_factor) for l in base_num_layers)), reverse=True)

    if mode_to_optimize == 'validation':
        batch_sizes = [int(bs * vram_factor) for bs in [64, 48, 32, 24, 16, 12, 8, 4, 2, 1]]
        max_atoms_target = int(5000 * vram_factor)
        stress_iterations = int(100 * stress_factor) # Reduced base iterations
    elif mode_to_optimize == 'prototyping':
        batch_sizes = [int(bs * vram_factor) for bs in [128, 96, 64, 48, 32, 24, 16, 8]]
        max_atoms_target = int(1024 * vram_factor)
        stress_iterations = int(50 * stress_factor) # Reduced base iterations
    elif mode_to_optimize == 'smoke_test':
        batch_sizes = [int(bs * vram_factor) for bs in [256, 128, 64, 32, 16, 8]]
        max_atoms_target = 128
        stress_iterations = 20
    else:  # production
        batch_sizes = [int(bs * vram_factor) for bs in [32, 24, 16, 12, 8, 6, 4, 2, 1]]
        max_atoms_target = int(10000 * vram_factor)
        stress_iterations = int(150 * stress_factor) # Reduced base iterations

    batch_sizes = sorted(list(set(bs for bs in batch_sizes if bs > 0)), reverse=True)

    print("\n--- Generated Search Space ---")
    print(f"  - Batch Sizes: {batch_sizes}")
    print(f"  - Hidden Channels: {hidden_channels_list}")
    print(f"  - Num Layers: {num_layers_list}")
    print(f"  - Max Atoms: {max_atoms_target}")
    print(f"  - Stress Iterations: {stress_iterations}\n")

    best_config = None

    # --- 2. Optimized Search Strategy with Graceful Exit ---
    try:
        if mode_to_optimize in ['production', 'smoke_test']:
            strategy = 'Find Most Complex Model' if mode_to_optimize == 'production' else 'Find Quickest Success'
            print(f"--- Optimizing for '{mode_to_optimize}' (Strategy: {strategy}) ---")
            layers_list = num_layers_list if mode_to_optimize == 'production' else sorted(num_layers_list)
            channels_list = hidden_channels_list if mode_to_optimize == 'production' else sorted(hidden_channels_list)
            param_combinations = list(product(layers_list, channels_list, batch_sizes))

            for i, (layers, channels, bs) in enumerate(param_combinations):
                print(f"[{(i+1):>3}/{len(param_combinations)}] ", end='')
                current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': bs, 'max_atoms': max_atoms_target}
                if probe_config(current_config, stress_iterations=stress_iterations):
                    best_config = current_config
                    print(f"\n>>> Found optimal configuration for '{mode_to_optimize}'!")
                    break
        else:  # validation & prototyping with new 2-stage strategy
            print(f"--- Optimizing for '{mode_to_optimize}' (Strategy: 2-Stage Max Batch Size) ---")
            # Stage 1: Find the maximum possible batch size with the simplest model
            print("\n--- Stage 1: Finding Maximum Batch Size ---")
            max_bs = 0
            simplest_model = {'visnet_hidden_channels': min(hidden_channels_list), 'visnet_num_layers': min(num_layers_list)}
            for i, bs in enumerate(batch_sizes):
                print(f"[{(i+1):>2}/{len(batch_sizes)}] ", end='')
                current_config = {**simplest_model, 'batch_size': bs, 'max_atoms': max_atoms_target}
                if probe_config(current_config, stress_iterations=stress_iterations):
                    max_bs = bs
                    print(f"  >>> Max batch size of {max_bs} is stable.")
                    break
            
            if max_bs == 0:
                print("Could not find any stable batch size. Aborting.")
            else:
                # Stage 2: Find the most complex model for that batch size
                print("\n--- Stage 2: Finding Best Model for Max Batch Size ---")
                param_combinations = list(product(num_layers_list, hidden_channels_list))
                for i, (layers, channels) in enumerate(param_combinations):
                    print(f"[{(i+1):>3}/{len(param_combinations)}] ", end='')
                    current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': max_bs, 'max_atoms': max_atoms_target}
                    if probe_config(current_config, stress_iterations=stress_iterations):
                        best_config = current_config
                        print(f"\n>>> Found optimal model for batch size {max_bs}!")
                        break

    except KeyboardInterrupt:
        print("\n\n--- User interrupted optimization. Saving best configuration found so far. ---")

    # --- Display and Save Results ---
    if best_config:
        print(f"\n--- Optimization Complete for '{mode_to_optimize}' ---")
        print("Selected optimal configuration based on mode philosophy:")
        print(json.dumps(best_config, indent=4))
        output_path = 'hardware_profile.json'
        final_profile = {}
        if os.path.exists(output_path):
            try:
                with open(output_path, 'r') as f:
                    final_profile = json.load(f)
            except (json.JSONDecodeError, IOError):
                pass
        final_profile[mode_to_optimize] = best_config
        print(f"\nUpdating settings in '{output_path}'...")
        with open(output_path, 'w') as f:
            json.dump(final_profile, f, indent=4)
        print("Update complete.")
    else:
        print(f"\n--- Optimization Failed for '{mode_to_optimize}' ---")
        print("Could not find any stable configuration.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find the optimal hardware configuration based on the philosophy of the chosen mode.")
    parser.add_argument(
        '--mode',
        type=str,
        required=True,
        choices=['production', 'validation', 'prototyping', 'smoke_test'],
        help="The development mode to optimize for, which determines the optimization strategy."
    )
    args = parser.parse_args()
    find_optimal_config(args.mode)
