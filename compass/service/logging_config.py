"""
Unified logging configuration for COMPASS service.
"""

import logging
import logging.handlers
import os
from pathlib import Path
from typing import Optional


def setup_logging(
    log_dir: str = "logs",
    log_level: str = "INFO",
    service_name: str = "compass",
    enable_file: bool = True,
    enable_console: bool = True,
    max_bytes: int = None,
    backup_count: int = None,
    rotation_strategy: str = "size",
) -> logging.Logger:
    """
    Setup unified logging configuration with rotation support.

    Args:
        log_dir: Directory for log files (can be relative or absolute path)
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        service_name: Name of the service
        enable_file: Enable file logging
        enable_console: Enable console logging
        max_bytes: Maximum log file size before rotation (default: 10MB)
        backup_count: Number of backup files to keep (default: 5)
        rotation_strategy: Rotation strategy - 'size' (RotatingFileHandler) or 'time' (TimedRotatingFileHandler)
            Default: 'size'

    Returns:
        logging.Logger: Configured root logger

    Note:
        - For 'size' rotation: max_bytes and backup_count are used
        - For 'time' rotation: rotation happens daily at midnight, backup_count files are kept
        - Environment variables can override defaults:
            - LOG_MAX_BYTES: Maximum file size in bytes
            - LOG_BACKUP_COUNT: Number of backup files
            - LOG_ROTATION_STRATEGY: 'size' or 'time'
    """
    log_path = Path(log_dir)
    log_path.mkdir(parents=True, exist_ok=True)

    # Get rotation configuration from environment or defaults
    max_bytes = max_bytes or int(os.getenv("LOG_MAX_BYTES", 10 * 1024 * 1024))  # 10MB default
    backup_count = backup_count or int(os.getenv("LOG_BACKUP_COUNT", 5))
    rotation_strategy = rotation_strategy or os.getenv("LOG_ROTATION_STRATEGY", "size").lower()

    # Create formatters
    detailed_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    simple_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%H:%M:%S"
    )

    # Get root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, log_level.upper()))

    # Clear existing handlers
    root_logger.handlers.clear()

    # File handler with rotation
    if enable_file:
        log_file = log_path / f"{service_name}.log"

        if rotation_strategy == "time":
            # Time-based rotation (daily at midnight)
            file_handler = logging.handlers.TimedRotatingFileHandler(
                log_file,
                when="midnight",
                interval=1,
                backupCount=backup_count,
                encoding="utf-8",
                utc=False,
            )
        else:
            # Size-based rotation (default)
            file_handler = logging.handlers.RotatingFileHandler(
                log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
            )

        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(detailed_formatter)
        root_logger.addHandler(file_handler)

        # Error file handler (only errors)
        error_log_file = log_path / f"{service_name}_errors.log"

        if rotation_strategy == "time":
            error_handler = logging.handlers.TimedRotatingFileHandler(
                error_log_file,
                when="midnight",
                interval=1,
                backupCount=backup_count,
                encoding="utf-8",
                utc=False,
            )
        else:
            error_handler = logging.handlers.RotatingFileHandler(
                error_log_file, maxBytes=max_bytes, backupCount=backup_count, encoding="utf-8"
            )

        error_handler.setLevel(logging.ERROR)
        error_handler.setFormatter(detailed_formatter)
        root_logger.addHandler(error_handler)

    # Console handler
    if enable_console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(simple_formatter)
        root_logger.addHandler(console_handler)

    return root_logger
