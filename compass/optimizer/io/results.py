import os
import json
import logging

logger = logging.getLogger("HardwareOptimizer")


def save_results(final_profile, project_root_path):
    """
    Saves the final hardware profile to a JSON file.

    Args:
        final_profile (dict): The final hardware profile to save.
        project_root_path (str): The root path of the project.
    """
    output_path = os.path.join(project_root_path, "hardware_profile.json")
    logger.info(f"--- Saving Final Results to '{output_path}' ---")
    with open(output_path, "w") as f:
        json.dump(final_profile, f, indent=4)
    logger.info("Update complete.")
