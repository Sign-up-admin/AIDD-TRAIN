import logging
import os
import sys

# --- GLOBAL LOGGER SETUP ---
logger = logging.getLogger("HardwareOptimizer")
progress_logger = logging.getLogger("ProgressBar")


def setup_logging(log_level="INFO"):
    """
    Sets up the logging for the hardware optimizer.

    Args:
        log_level (str): The logging level to set.
    """
    logger.setLevel(log_level)
    progress_logger.setLevel(log_level)
    logger.propagate = False
    progress_logger.propagate = False
    if not logger.handlers:
        # Ensure the main logs directory exists and place the log file inside it.
        log_dir = 'logs'
        os.makedirs(log_dir, exist_ok=True)
        log_file_path = os.path.join(log_dir, "hardware_optimizer.log")

        file_handler = logging.FileHandler(log_file_path, mode='w')
        file_formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
        console_handler = logging.StreamHandler(sys.stdout)
        console_formatter = logging.Formatter("%(levelname)s: %(message)s")
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)
    if not progress_logger.handlers:
        progress_handler = logging.StreamHandler(sys.stdout)
        progress_handler.setFormatter(logging.Formatter('%(message)s'))
        progress_logger.addHandler(progress_handler)
