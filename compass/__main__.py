"""
This file makes the 'compass' package executable.

You can run the main training script from the project root directory using:
python -m compass
"""

import os
import platform
from multiprocessing import set_start_method

from .main import main
from .config import get_config
from .logger import TrainingLogger

if __name__ == "__main__":
    # --- Config and Logger Setup ---
    config = get_config()
    log_dir = config["log_dir"]
    os.makedirs(log_dir, exist_ok=True)
    logger = TrainingLogger(log_dir=log_dir)
    logger.log("--- Compass Training Process Started ---")

    # Set the multiprocessing start method for non-Linux systems to 'spawn'.
    # This is crucial for preventing deadlocks with CUDA and multiprocessing
    # and must be called within the `if __name__ == '__main__':` block.
    try:
        if platform.system() != "Linux":
            logger.log("Setting multiprocessing start method to 'spawn' for non-Linux OS.")
            set_start_method("spawn")
    except RuntimeError:
        # This can happen if the method has already been set, which is fine.
        logger.log_warning(
            "Multiprocessing start method could not be set, it might have been set already."
        )

    # This entry point is executed when the package is run as a script.
    main(config, logger)
