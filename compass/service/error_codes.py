"""
Standardized error codes for COMPASS service.
"""
from enum import Enum
from typing import Dict, Optional


class ErrorCode(str, Enum):
    """Standardized error codes."""
    
    # General errors (1xxx)
    INTERNAL_ERROR = "ERR_1000"
    INVALID_REQUEST = "ERR_1001"
    SERVICE_UNAVAILABLE = "ERR_1002"
    TIMEOUT = "ERR_1003"
    
    # Validation errors (2xxx)
    VALIDATION_ERROR = "ERR_2000"
    INVALID_PARAMETER = "ERR_2001"
    MISSING_PARAMETER = "ERR_2002"
    INVALID_FILE_TYPE = "ERR_2003"
    FILE_TOO_LARGE = "ERR_2004"
    INVALID_FILE_FORMAT = "ERR_2005"
    
    # Resource errors (3xxx)
    NOT_FOUND = "ERR_3000"
    RESOURCE_NOT_FOUND = "ERR_3001"
    DATASET_NOT_FOUND = "ERR_3002"
    MODEL_NOT_FOUND = "ERR_3003"
    TASK_NOT_FOUND = "ERR_3004"
    
    # Conflict errors (4xxx)
    RESOURCE_CONFLICT = "ERR_4000"
    TASK_ALREADY_RUNNING = "ERR_4001"
    TASK_ALREADY_COMPLETED = "ERR_4002"
    RESOURCE_EXISTS = "ERR_4003"
    
    # Rate limiting (5xxx)
    RATE_LIMIT_EXCEEDED = "ERR_5000"
    TOO_MANY_REQUESTS = "ERR_5001"
    
    # Training errors (6xxx)
    TRAINING_FAILED = "ERR_6000"
    TRAINING_CANCELLED = "ERR_6001"
    INVALID_TRAINING_CONFIG = "ERR_6002"
    DATASET_NOT_READY = "ERR_6003"
    
    # Inference errors (7xxx)
    INFERENCE_FAILED = "ERR_7000"
    MODEL_NOT_LOADED = "ERR_7001"
    INVALID_INPUT_DATA = "ERR_7002"
    BATCH_INFERENCE_FAILED = "ERR_7003"
    
    # Upload errors (8xxx)
    UPLOAD_FAILED = "ERR_8000"
    UPLOAD_QUEUE_FULL = "ERR_8001"
    ZIP_BOMB_DETECTED = "ERR_8002"
    UPLOAD_TIMEOUT = "ERR_8003"
    
    # Service registry errors (9xxx)
    REGISTRY_ERROR = "ERR_9000"
    REGISTRATION_FAILED = "ERR_9001"
    SERVICE_NOT_REGISTERED = "ERR_9002"


# Error code metadata
ERROR_MESSAGES: Dict[ErrorCode, str] = {
    ErrorCode.INTERNAL_ERROR: "Internal server error",
    ErrorCode.INVALID_REQUEST: "Invalid request",
    ErrorCode.SERVICE_UNAVAILABLE: "Service temporarily unavailable",
    ErrorCode.TIMEOUT: "Request timeout",
    
    ErrorCode.VALIDATION_ERROR: "Validation error",
    ErrorCode.INVALID_PARAMETER: "Invalid parameter",
    ErrorCode.MISSING_PARAMETER: "Missing required parameter",
    ErrorCode.INVALID_FILE_TYPE: "Invalid file type",
    ErrorCode.FILE_TOO_LARGE: "File size exceeds limit",
    ErrorCode.INVALID_FILE_FORMAT: "Invalid file format",
    
    ErrorCode.NOT_FOUND: "Resource not found",
    ErrorCode.RESOURCE_NOT_FOUND: "Resource not found",
    ErrorCode.DATASET_NOT_FOUND: "Dataset not found",
    ErrorCode.MODEL_NOT_FOUND: "Model not found",
    ErrorCode.TASK_NOT_FOUND: "Task not found",
    
    ErrorCode.RESOURCE_CONFLICT: "Resource conflict",
    ErrorCode.TASK_ALREADY_RUNNING: "Task is already running",
    ErrorCode.TASK_ALREADY_COMPLETED: "Task is already completed",
    ErrorCode.RESOURCE_EXISTS: "Resource already exists",
    
    ErrorCode.RATE_LIMIT_EXCEEDED: "Rate limit exceeded",
    ErrorCode.TOO_MANY_REQUESTS: "Too many requests",
    
    ErrorCode.TRAINING_FAILED: "Training failed",
    ErrorCode.TRAINING_CANCELLED: "Training cancelled",
    ErrorCode.INVALID_TRAINING_CONFIG: "Invalid training configuration",
    ErrorCode.DATASET_NOT_READY: "Dataset is not ready",
    
    ErrorCode.INFERENCE_FAILED: "Inference failed",
    ErrorCode.MODEL_NOT_LOADED: "Model not loaded",
    ErrorCode.INVALID_INPUT_DATA: "Invalid input data",
    ErrorCode.BATCH_INFERENCE_FAILED: "Batch inference failed",
    
    ErrorCode.UPLOAD_FAILED: "Upload failed",
    ErrorCode.UPLOAD_QUEUE_FULL: "Upload queue is full",
    ErrorCode.ZIP_BOMB_DETECTED: "Zip bomb detected",
    ErrorCode.UPLOAD_TIMEOUT: "Upload timeout",
    
    ErrorCode.REGISTRY_ERROR: "Service registry error",
    ErrorCode.REGISTRATION_FAILED: "Service registration failed",
    ErrorCode.SERVICE_NOT_REGISTERED: "Service not registered",
}


def get_error_message(error_code: ErrorCode, detail: Optional[str] = None) -> str:
    """
    Get error message for error code.
    
    Args:
        error_code: Error code
        detail: Optional additional detail
    
    Returns:
        Error message string
    """
    base_message = ERROR_MESSAGES.get(error_code, "Unknown error")
    if detail:
        return f"{base_message}: {detail}"
    return base_message


def get_http_status_code(error_code: ErrorCode) -> int:
    """
    Map error code to HTTP status code.
    
    Args:
        error_code: Error code
    
    Returns:
        HTTP status code
    """
    code_prefix = error_code.value.split('_')[1][0]
    
    status_map = {
        '1': 500,  # Internal errors
        '2': 400,  # Validation errors
        '3': 404,  # Not found
        '4': 409,  # Conflict
        '5': 429,  # Rate limit
        '6': 400,  # Training errors
        '7': 400,  # Inference errors
        '8': 400,  # Upload errors
        '9': 503,  # Registry errors
    }
    
    return status_map.get(code_prefix, 500)

