"""
Custom exceptions for COMPASS service.
"""

from typing import Optional, Dict
from fastapi import HTTPException, status
from compass.service.error_codes import ErrorCode, get_error_message, get_http_status_code


class ServiceException(Exception):
    """Base exception for service errors."""

    def __init__(
        self,
        message: str,
        status_code: int = status.HTTP_500_INTERNAL_SERVER_ERROR,
        error_code: Optional[ErrorCode] = None,
        detail: Optional[Dict] = None,
    ):
        self.message = message
        self.error_code = error_code
        self.detail = detail if detail is not None else {}

        # Use error code's status code if provided
        if error_code:
            self.status_code = get_http_status_code(error_code)
            # Use standard error message if no custom message provided
            if not message or message == str(Exception()):
                detail_str = str(detail) if detail else None
                self.message = get_error_message(error_code, detail_str)
        else:
            self.status_code = status_code

        super().__init__(self.message)


class ValidationError(ServiceException):
    """Validation error."""

    def __init__(
        self,
        message: str,
        error_code: ErrorCode = ErrorCode.VALIDATION_ERROR,
        detail: Optional[Dict] = None,
    ):
        super().__init__(message=message, error_code=error_code, detail=detail)


class NotFoundError(ServiceException):
    """Resource not found error."""

    def __init__(
        self,
        resource_type: str,
        resource_id: str,
        error_code: ErrorCode = ErrorCode.RESOURCE_NOT_FOUND,
    ):
        # Map to specific error codes
        if resource_type.lower() == "dataset":
            error_code = ErrorCode.DATASET_NOT_FOUND
        elif resource_type.lower() == "model":
            error_code = ErrorCode.MODEL_NOT_FOUND
        elif resource_type.lower() == "task":
            error_code = ErrorCode.TASK_NOT_FOUND

        super().__init__(
            message=f"{resource_type} {resource_id} not found",
            error_code=error_code,
            detail={"resource_type": resource_type, "resource_id": resource_id},
        )


class ConflictError(ServiceException):
    """Resource conflict error."""

    def __init__(
        self,
        message: str,
        error_code: ErrorCode = ErrorCode.RESOURCE_CONFLICT,
        detail: Optional[Dict] = None,
    ):
        super().__init__(message=message, error_code=error_code, detail=detail)


class RateLimitError(ServiceException):
    """Rate limit exceeded error."""

    def __init__(self, message: str = "Rate limit exceeded"):
        super().__init__(message, error_code=ErrorCode.RATE_LIMIT_EXCEEDED)


def sanitize_error_message(error: Exception, include_details: bool = False) -> str:
    """
    Sanitize error message to avoid exposing internal details.

    Args:
        error: The exception
        include_details: Whether to include detailed information (for logging)

    Returns:
        Sanitized error message
    """
    error_type = type(error).__name__
    error_message = str(error)

    # User-friendly messages for common errors
    user_messages = {
        "ValueError": "Invalid input provided",
        "FileNotFoundError": "Required file not found",
        "PermissionError": "Permission denied",
        "ConnectionError": "Failed to connect to service",
        "TimeoutError": "Request timed out",
        "KeyError": "Required parameter missing",
        "TypeError": "Invalid parameter type",
    }

    if include_details:
        # For logging, include full details
        return f"{error_type}: {error_message}"

    # For user response, use friendly message
    user_message = user_messages.get(error_type, "An error occurred")

    return user_message
