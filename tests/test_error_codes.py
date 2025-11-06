"""
Tests for error codes.
"""
import pytest
from compass.service.error_codes import (
    ErrorCode,
    get_error_message,
    get_http_status_code,
    ERROR_MESSAGES
)


def test_error_code_enum():
    """Test ErrorCode enum."""
    assert ErrorCode.INTERNAL_ERROR == "ERR_1000"
    assert ErrorCode.VALIDATION_ERROR == "ERR_2000"
    assert ErrorCode.NOT_FOUND == "ERR_3000"


def test_get_error_message():
    """Test get_error_message function."""
    message = get_error_message(ErrorCode.VALIDATION_ERROR)
    assert message == ERROR_MESSAGES[ErrorCode.VALIDATION_ERROR]
    
    message_with_detail = get_error_message(ErrorCode.VALIDATION_ERROR, "Invalid parameter")
    assert "Invalid parameter" in message_with_detail


def test_get_http_status_code():
    """Test get_http_status_code function."""
    # 1xxx -> 500
    assert get_http_status_code(ErrorCode.INTERNAL_ERROR) == 500
    
    # 2xxx -> 400
    assert get_http_status_code(ErrorCode.VALIDATION_ERROR) == 400
    
    # 3xxx -> 404
    assert get_http_status_code(ErrorCode.NOT_FOUND) == 404
    
    # 4xxx -> 409
    assert get_http_status_code(ErrorCode.RESOURCE_CONFLICT) == 409
    
    # 5xxx -> 429
    assert get_http_status_code(ErrorCode.RATE_LIMIT_EXCEEDED) == 429
    
    # 6xxx -> 400
    assert get_http_status_code(ErrorCode.TRAINING_FAILED) == 400
    
    # 7xxx -> 400
    assert get_http_status_code(ErrorCode.INFERENCE_FAILED) == 400
    
    # 8xxx -> 400
    assert get_http_status_code(ErrorCode.UPLOAD_FAILED) == 400
    
    # 9xxx -> 503
    assert get_http_status_code(ErrorCode.REGISTRY_ERROR) == 503


