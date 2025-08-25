import os
import json
import logging

# Use a dedicated logger for this utility to avoid interfering with other loggers
# during the config loading phase.
loader_logger = logging.getLogger("ConfigLoader")
if not loader_logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    loader_logger.addHandler(handler)
    loader_logger.setLevel(logging.INFO)
    loader_logger.propagate = False

def load_optimized_settings(config: dict, profile_path: str = 'hardware_profile.json') -> None:
    """
    Checks for a hardware profile and applies optimized settings if they exist for the current mode.

    This function looks for a 'hardware_profile.json' file. If found, and if the
    current execution mode is defined within it, it updates the main
    configuration dictionary with the optimized parameters, logging any overrides.

    Args:
        config (dict): The main configuration dictionary to update.
        profile_path (str): The path to the hardware profile JSON file.
    """
    if not os.path.exists(profile_path):
        loader_logger.info("--- No hardware profile found. Using default configuration. ---")
        loader_logger.info("--- Run the optimizer to generate one: python -m compass.optimizer ---")
        return

    loader_logger.info(f"--- Found hardware profile: '{profile_path}' ---")
    try:
        with open(profile_path, 'r') as f:
            optimized_params = json.load(f)

        current_mode = config.get('execution_mode')
        if current_mode and current_mode in optimized_params:
            mode_settings = optimized_params[current_mode]
            loader_logger.info(f"--- Applying optimized settings for '{current_mode}' mode... ---")
            for key, value in mode_settings.items():
                if key in config and config[key] != value:
                    loader_logger.info(f"    - Overriding '{key}': {config[key]} -> {value}")
                else:
                    loader_logger.info(f"    - Setting '{key}': {value}")
            config.update(mode_settings)
            loader_logger.info("--- Optimized settings applied. ---")
        else:
            loader_logger.info(f"--- No optimized settings found for '{current_mode}' mode in profile. Using default settings. ---")

    except (json.JSONDecodeError, IOError, Exception) as e:
        loader_logger.warning(f"[Warning] Could not load or apply hardware profile: {e}")
