"""
Tests for service exceptions.
"""
import pytest
from fastapi import status
from compass.service.exceptions import (
    ServiceException,
    NotFoundError,
    ValidationError,
    ConflictError,
    RateLimitError
)
from compass.service.error_codes import ErrorCode


def test_service_exception():
    """Test ServiceException."""
    exc = ServiceException("Test error", status_code=400)
    assert exc.message == "Test error"
    assert exc.status_code == 400
    assert exc.error_code is None


def test_service_exception_with_error_code():
    """Test ServiceException with error code."""
    exc = ServiceException(
        "Test error",
        error_code=ErrorCode.VALIDATION_ERROR
    )
    assert exc.message == "Test error"
    assert exc.status_code == 400  # Auto-mapped from error code
    assert exc.error_code == ErrorCode.VALIDATION_ERROR


def test_not_found_error():
    """Test NotFoundError."""
    exc = NotFoundError("Task", "test-id")
    assert "Task" in exc.message
    assert "test-id" in exc.message
    assert exc.status_code == 404
    assert exc.error_code == ErrorCode.TASK_NOT_FOUND


def test_not_found_error_dataset():
    """Test NotFoundError for dataset."""
    exc = NotFoundError("Dataset", "test-dataset")
    assert exc.error_code == ErrorCode.DATASET_NOT_FOUND


def test_not_found_error_model():
    """Test NotFoundError for model."""
    exc = NotFoundError("Model", "test-model")
    assert exc.error_code == ErrorCode.MODEL_NOT_FOUND


def test_validation_error():
    """Test ValidationError."""
    exc = ValidationError("Invalid input")
    assert exc.message == "Invalid input"
    assert exc.status_code == 400
    assert exc.error_code == ErrorCode.VALIDATION_ERROR


def test_conflict_error():
    """Test ConflictError."""
    exc = ConflictError("Resource conflict")
    assert exc.message == "Resource conflict"
    assert exc.status_code == 409
    assert exc.error_code == ErrorCode.RESOURCE_CONFLICT


def test_rate_limit_error():
    """Test RateLimitError."""
    exc = RateLimitError()
    assert exc.status_code == 429
    assert exc.error_code == ErrorCode.RATE_LIMIT_EXCEEDED












