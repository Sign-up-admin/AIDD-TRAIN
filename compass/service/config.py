"""
Configuration for COMPASS service.
"""

from typing import Dict
from compass.service.config_manager import get_config_manager


# Backward compatibility: maintain SERVICE_CONFIG dict
def get_service_config() -> Dict:
    """
    Get service configuration (backward compatibility).

    Returns:
        dict: Service configuration
    """
    manager = get_config_manager()
    config = manager.get_dict()

    # Handle upload_max_size conversion (from GB to bytes if needed)
    upload_max_size = config.get("upload_max_size", 10 * 1024 * 1024 * 1024)
    if upload_max_size < 1024 * 1024:  # If less than 1MB, assume it's in GB
        upload_max_size = upload_max_size * 1024 * 1024 * 1024
    config["upload_max_size"] = upload_max_size

    return config


SERVICE_CONFIG = get_service_config()
