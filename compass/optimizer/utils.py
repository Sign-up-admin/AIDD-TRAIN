import os
import json
import logging

logger = logging.getLogger(__name__)

def load_optimized_settings(config):
    """
    Loads hardware-optimized settings from 'hardware_profile.json' and
    updates the configuration for the current execution mode.
    """
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    profile_path = os.path.join(project_root, 'hardware_profile.json')

    if not os.path.exists(profile_path):
        logger.info("Optimizer profile 'hardware_profile.json' not found. Using default settings.")
        return

    try:
        with open(profile_path, 'r') as f:
            optimized_profiles = json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error reading 'hardware_profile.json': {e}. Using default settings.")
        return

    mode = config.get('execution_mode')
    if not mode:
        logger.warning("'execution_mode' not found in config. Cannot apply optimized settings.")
        return

    mode_settings = optimized_profiles.get(mode)
    if not mode_settings:
        logger.info(f"No optimized settings found for '{mode}' mode in 'hardware_profile.json'. Using default settings.")
        return

    logger.info(f"Applying hardware-optimized settings from 'hardware_profile.json' for '{mode}' mode.")
    
    original_values = {key: config.get(key) for key in mode_settings}

    config.update(mode_settings)

    for key, new_value in mode_settings.items():
        original_value = original_values[key]
        if original_value != new_value:
            logger.info(f"  - Overrode '{key}': from '{original_value}' to '{new_value}'")
        else:
            logger.info(f"  - Set '{key}': '{new_value}' (no change)")
