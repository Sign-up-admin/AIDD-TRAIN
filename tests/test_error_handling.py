"""
Tests for error handling and exception management.
"""

import pytest
from fastapi.testclient import TestClient
from unittest.mock import Mock, patch, MagicMock
from compass.service.server import app
from compass.service.exceptions import (
    ServiceException,
    NotFoundError,
    ValidationError,
    sanitize_error_message,
)
from compass.service.error_codes import ErrorCode


@pytest.fixture
def client():
    """Create test client."""
    return TestClient(app)


def test_sanitize_error_message():
    """Test error message sanitization."""
    # Test with details (for logging)
    error = ValueError("Invalid input")
    detailed = sanitize_error_message(error, include_details=True)
    assert "ValueError" in detailed
    assert "Invalid input" in detailed
    
    # Test without details (for user response)
    user_friendly = sanitize_error_message(error, include_details=False)
    assert user_friendly == "Invalid input provided"
    
    # Test unknown error type
    unknown_error = Exception("Unknown error")
    user_msg = sanitize_error_message(unknown_error, include_details=False)
    assert user_msg == "An error occurred"


def test_service_exception_handler(client):
    """Test that service exceptions are handled properly."""
    with patch("compass.service.routes.training.training_service") as mock_service:
        mock_service.create_task.side_effect = NotFoundError("Task", "test-id")
        
        response = client.post(
            "/api/v1/training/tasks",
            json={
                "config": {
                    "execution_mode": "prototyping",
                    "epochs": 10,
                    "batch_size": 32,
                    "learning_rate": 0.001,
                    "optimizer": "adam",
                },
            },
        )
        
        # Should return proper error response
        assert response.status_code == 404
        assert "error" in response.json()


def test_validation_exception_handler(client, monkeypatch):
    """Test that validation errors are handled properly."""
    # Ensure authentication is disabled for tests
    monkeypatch.setenv('AUTH_ENABLED', 'false')
    monkeypatch.setenv('FORCE_AUTH_CRITICAL', 'false')
    import importlib
    import compass.service.middleware.auth
    importlib.reload(compass.service.middleware.auth)
    
    # Invalid config should return 422
    response = client.post(
        "/api/v1/training/tasks",
        json={
            "config": {
                "execution_mode": "invalid_mode",
                "epochs": 10,
            },
        },
    )
    
    assert response.status_code == 422
    assert "detail" in response.json()


def test_general_exception_handler(client):
    """Test that unexpected exceptions are handled gracefully."""
    with patch("compass.service.routes.training.training_service") as mock_service:
        mock_service.create_task.side_effect = Exception("Unexpected error")
        
        response = client.post(
            "/api/v1/training/tasks",
            json={
                "config": {
                    "execution_mode": "prototyping",
                    "epochs": 10,
                    "batch_size": 32,
                    "learning_rate": 0.001,
                    "optimizer": "adam",
                },
            },
        )
        
        # Should return 500 but not expose internal details
        assert response.status_code == 500
        data = response.json()
        assert "error" in data
        # Should not contain the actual exception message
        assert "Unexpected error" not in data.get("error", "")


def test_error_codes():
    """Test that error codes are properly mapped."""
    from compass.service.error_codes import get_http_status_code, get_error_message
    
    # Test status code mapping
    assert get_http_status_code(ErrorCode.VALIDATION_ERROR) == 400
    assert get_http_status_code(ErrorCode.NOT_FOUND) == 404
    assert get_http_status_code(ErrorCode.RATE_LIMIT_EXCEEDED) == 429
    assert get_http_status_code(ErrorCode.INTERNAL_ERROR) == 500
    
    # Test error message retrieval
    msg = get_error_message(ErrorCode.VALIDATION_ERROR)
    assert isinstance(msg, str)
    assert len(msg) > 0


def test_exception_with_error_code():
    """Test that exceptions with error codes have proper status codes."""
    exc = ServiceException(
        "Test error",
        error_code=ErrorCode.RATE_LIMIT_EXCEEDED
    )
    assert exc.status_code == 429
    assert exc.error_code == ErrorCode.RATE_LIMIT_EXCEEDED


def test_not_found_error_mapping():
    """Test that NotFoundError maps to correct error codes based on resource type."""
    task_error = NotFoundError("Task", "test-id")
    assert task_error.error_code == ErrorCode.TASK_NOT_FOUND
    
    dataset_error = NotFoundError("Dataset", "test-dataset")
    assert dataset_error.error_code == ErrorCode.DATASET_NOT_FOUND
    
    model_error = NotFoundError("Model", "test-model")
    assert model_error.error_code == ErrorCode.MODEL_NOT_FOUND

