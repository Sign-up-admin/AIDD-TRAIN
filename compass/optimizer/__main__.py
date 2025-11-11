import os
import sys

# This must be the first action to ensure all imports are absolute.
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if project_root not in sys.path:
    sys.path.insert(0, project_root)

import torch
import json
import argparse
import logging
from itertools import product

from compass.optimizer.logging_setup import setup_logging
from compass.optimizer.io.data import prepare_real_data
from compass.optimizer.io.results import save_results
from compass.optimizer.search_space.builder import get_parameter_cap, get_search_space
from compass.optimizer.strategies import (
    _optimize_for_production,
    _optimize_for_efficiency_modes,
    _optimize_for_smoke_test,
)
from compass.optimizer.config import OPTIMIZATION_HIERARCHY, TEST_SAMPLE_CONFIG

logger = logging.getLogger("HardwareOptimizer")


def find_optimal_configs(modes_to_optimize, current_dataset_size, project_root_path):
    """
    Finds the optimal hardware configuration for different modes.

    Args:
        modes_to_optimize (list): A list of modes to optimize for.
        current_dataset_size (int): The total number of samples in the dataset.
        project_root_path (str): The root path of the project.
    """
    if not torch.cuda.is_available():
        logger.critical("CUDA is not available. Aborting.")
        return

    logger.warning(
        f"This optimizer uses a single data sample ('{TEST_SAMPLE_CONFIG['pdb_code']}') for all tests. For best results, ensure this sample is representative of your dataset."
    )

    processed_test_data = prepare_real_data(project_root_path)

    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    logger.info(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")

    max_params = get_parameter_cap(current_dataset_size)

    output_path = os.path.join(project_root_path, "hardware_profile.json")
    final_profile = {}
    if os.path.exists(output_path):
        try:
            with open(output_path, "r") as f:
                final_profile = json.load(f)
        except (json.JSONDecodeError, IOError):
            pass

    modes_to_run = [m for m in OPTIMIZATION_HIERARCHY if m in modes_to_optimize]

    optimizer_map = {
        "production": _optimize_for_production,
        "validation": _optimize_for_efficiency_modes,
        "prototyping": _optimize_for_efficiency_modes,
        "smoke_test": _optimize_for_smoke_test,
    }

    try:
        for mode in modes_to_run:
            logger.info(f"================ Optimizing for: {mode.upper()} ================")
            start_bs, hidden_channels, num_layers, stress_iters = get_search_space(vram_gb, mode)
            param_combinations = list(product(num_layers, hidden_channels))

            optimizer_func = optimizer_map[mode]

            if mode == "production":
                optimizer_args = (
                    param_combinations,
                    max_params,
                    start_bs,
                    stress_iters,
                    processed_test_data,
                )
            elif mode in ["validation", "prototyping"]:
                optimizer_args = (
                    mode,
                    num_layers,
                    hidden_channels,
                    max_params,
                    start_bs,
                    stress_iters,
                    processed_test_data,
                )
            else:  # smoke_test
                optimizer_args = (param_combinations, max_params, processed_test_data)

            best_config_for_this_mode = optimizer_func(*optimizer_args)

            if best_config_for_this_mode:
                final_profile[mode] = best_config_for_this_mode
                logger.info(f">>> Found optimal configuration for '{mode}'!")
                pretty_json = json.dumps(best_config_for_this_mode, indent=4)
                for line in pretty_json.split("\n"):
                    logger.info(line)
            else:
                logger.error(
                    f"--- Optimization FAILED for '{mode}'. No stable configuration found. ---"
                )

    except KeyboardInterrupt:
        logger.warning(
            "--- User interrupted optimization. Saving best configuration found so far. ---"
        )

    if final_profile:
        save_results(final_profile, project_root_path)
    else:
        logger.warning("No valid configurations were found.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find the optimal hardware configuration for each specified mode.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--modes",
        nargs="+",
        default=OPTIMIZATION_HIERARCHY,
        choices=OPTIMIZATION_HIERARCHY,
        help="A list of development modes to optimize for.",
    )
    parser.add_argument(
        "--dataset-size",
        type=int,
        default=None,
        help="The total number of samples in the dataset.\nIf not provided, you will be prompted to enter it.",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level.",
    )
    cli_args = parser.parse_args()

    setup_logging(cli_args.log_level)

    dataset_size = cli_args.dataset_size
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

    find_optimal_configs(cli_args.modes, dataset_size, project_root)
