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

def probe_config(config, quiet=False, stress_iterations=100):
    """
    Tests a given configuration for OOM errors under sustained load.
    """
    if not quiet:
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

        # Sustained load simulation to heat up the GPU and test stability
        for i in range(stress_iterations):
            optimizer.zero_grad()
            output = model(x, pos, batch)
            scaled_output = output * 100
            scaled_output.backward()
            optimizer.step()
            if not quiet and (i + 1) % 20 == 0:
                print(".", end='', flush=True)

        if not quiet:
            print(" SUCCESS")
        return True
    except torch.cuda.OutOfMemoryError:
        if not quiet:
            print(" FAILED (OOM)")
        return False
    finally:
        del model, x, pos, batch, optimizer
        gc.collect()
        torch.cuda.empty_cache()

def find_optimal_config(mode_to_optimize):
    """
    Finds the best hardware configuration based on the philosophy of the chosen mode.
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting hardware optimization.")
        return

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    # --- Define Search Spaces per Mode ---
    if mode_to_optimize == 'validation':
        batch_sizes = [64, 48, 32, 24, 16, 12, 8, 4, 2, 1]
        max_atoms_target = 5000
    elif mode_to_optimize == 'prototyping':
        batch_sizes = [128, 96, 64, 48, 32, 24, 16, 8]
        max_atoms_target = 1024
    elif mode_to_optimize == 'smoke_test':
        batch_sizes = [256, 128, 64, 32, 16, 8]
        max_atoms_target = 128
    else:  # production
        batch_sizes = [32, 24, 16, 12, 8, 6, 4, 2, 1]
        max_atoms_target = 10000

    hidden_channels_list = [256, 192, 128, 96, 64]
    num_layers_list = [8, 7, 6, 5, 4, 3]

    best_config = None

    # --- Select Optimization Strategy Based on Mode Philosophy ---

    if mode_to_optimize == 'production':
        print(f"\n--- Optimizing for '{mode_to_optimize}' (Strategy: Find Most Complex Model) ---")
        param_combinations = list(product(num_layers_list, hidden_channels_list, batch_sizes))
        for i, (layers, channels, bs) in enumerate(param_combinations):
            print(f"[{(i+1):>3}/{len(param_combinations)}] ", end='')
            current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': bs, 'max_atoms': max_atoms_target}
            if probe_config(current_config, stress_iterations=100):
                best_config = current_config
                print(f"\n>>> Found optimal configuration for '{mode_to_optimize}'!")
                break

    elif mode_to_optimize == 'smoke_test':
        print(f"\n--- Optimizing for '{mode_to_optimize}' (Strategy: Find Quickest Success) ---")
        param_combinations = list(product(sorted(num_layers_list), sorted(hidden_channels_list), batch_sizes))
        for i, (layers, channels, bs) in enumerate(param_combinations):
            print(f"[{(i+1):>3}/{len(param_combinations)}] ", end='')
            current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': bs, 'max_atoms': max_atoms_target}
            if probe_config(current_config, stress_iterations=10, quiet=True):
                best_config = current_config
                print(f"\n>>> Found optimal configuration for '{mode_to_optimize}'!")
                break

    else: # validation & prototyping
        print(f"\n--- Optimizing for '{mode_to_optimize}' (Strategy: Find Max Batch Size) ---")
        successful_configs = []
        model_param_combinations = list(product(num_layers_list, hidden_channels_list))
        for i, (layers, channels) in enumerate(model_param_combinations):
            print(f"\n--- [{(i+1):>2}/{len(model_param_combinations)}] Testing Model: {layers} layers, {channels} channels ---")
            max_bs_for_this_model = 0
            for bs in batch_sizes:
                current_config = {'visnet_hidden_channels': channels, 'visnet_num_layers': layers, 'batch_size': bs, 'max_atoms': max_atoms_target}
                if probe_config(current_config, stress_iterations=100, quiet=True):
                    max_bs_for_this_model = bs
                    successful_configs.append(current_config)
                    print(f"  - Found stable batch size: {max_bs_for_this_model}")
                    break
            if max_bs_for_this_model == 0:
                print("  - This model is too large for any batch size.")
        if successful_configs:
            best_config = sorted(successful_configs, key=lambda c: (c['batch_size'], c['visnet_num_layers'], c['visnet_hidden_channels']), reverse=True)[0]

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
            except json.JSONDecodeError:
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
