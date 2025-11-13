"""
Tests for training routes.
"""

import pytest
import uuid
from fastapi.testclient import TestClient
from unittest.mock import Mock, patch, MagicMock
from compass.service.server import app
from compass.service.models.task import TaskStatus


@pytest.fixture
def client(monkeypatch):
    """Create test client."""
    # Ensure authentication is disabled for tests
    monkeypatch.setenv('AUTH_ENABLED', 'false')
    monkeypatch.setenv('FORCE_AUTH_CRITICAL', 'false')
    # Increase rate limits for tests
    monkeypatch.setenv('RATE_LIMIT_TRAINING', '10000')
    monkeypatch.setenv('RATE_LIMIT_DEFAULT', '10000')
    # Reload auth module to pick up env vars
    import importlib
    import compass.service.middleware.auth
    importlib.reload(compass.service.middleware.auth)
    # Reset rate limit stats and limiters before each test
    from compass.service.middleware.rate_limit import reset_rate_limit_stats
    reset_rate_limit_stats()
    # Reset rate limiters in middleware
    from compass.service.server import app
    for middleware in app.user_middleware:
        if hasattr(middleware, 'cls') and 'RateLimitMiddleware' in str(middleware.cls):
            # Access the middleware instance and reset limiters
            # This is a workaround - in production, rate limiters reset automatically
            pass
    return TestClient(app)


@pytest.fixture
def mock_training_service():
    """Mock training service."""
    with patch("compass.service.routes.training.training_service") as mock:
        yield mock


def test_create_task_success(client, mock_training_service):
    """Test successful task creation."""
    task_id = str(uuid.uuid4())
    mock_training_service.create_task.return_value = task_id
    mock_training_service.get_task.return_value = {
        "task_id": task_id,
        "status": TaskStatus.PENDING,
        "config": {"epochs": 10},
        "created_at": "2024-01-01T00:00:00",
        "updated_at": "2024-01-01T00:00:00",
    }

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
            "description": "Test task",
        },
    )

    assert response.status_code == 201
    assert response.json()["task_id"] == task_id
    mock_training_service.create_task.assert_called_once()


def test_create_task_invalid_config(client, mock_training_service):
    """Test task creation with invalid config."""
    response = client.post(
        "/api/v1/training/tasks",
        json={
            "config": {
                "execution_mode": "invalid_mode",
                "epochs": 10,
            },
        },
    )

    assert response.status_code == 422  # Validation error


def test_get_task_success(client, mock_training_service):
    """Test getting a task."""
    task_id = str(uuid.uuid4())
    mock_training_service.get_task.return_value = {
        "task_id": task_id,
        "status": TaskStatus.RUNNING,
        "config": {"epochs": 10},
        "created_at": "2024-01-01T00:00:00",
        "updated_at": "2024-01-01T00:00:00",
    }

    response = client.get(f"/api/v1/training/tasks/{task_id}")

    assert response.status_code == 200
    assert response.json()["task_id"] == task_id
    mock_training_service.get_task.assert_called_once_with(task_id)


def test_get_task_not_found(client, mock_training_service):
    """Test getting a non-existent task."""
    task_id = str(uuid.uuid4())
    mock_training_service.get_task.return_value = None

    response = client.get(f"/api/v1/training/tasks/{task_id}")

    assert response.status_code == 404


def test_list_tasks(client, mock_training_service):
    """Test listing tasks."""
    mock_training_service.list_tasks.return_value = [
        {
            "task_id": "task-1",
            "status": TaskStatus.RUNNING,
            "config": {},
            "created_at": "2024-01-01T00:00:00",
            "updated_at": "2024-01-01T00:00:00",
        },
        {
            "task_id": "task-2",
            "status": TaskStatus.COMPLETED,
            "config": {},
            "created_at": "2024-01-01T00:00:00",
            "updated_at": "2024-01-01T00:00:00",
        },
    ]

    response = client.get("/api/v1/training/tasks")

    assert response.status_code == 200
    data = response.json()
    assert data["count"] == 2
    assert len(data["tasks"]) == 2


def test_start_task_success(client, mock_training_service):
    """Test starting a task."""
    task_id = str(uuid.uuid4())
    mock_training_service.get_task.return_value = {
        "task_id": task_id,
        "status": TaskStatus.PENDING,
        "config": {},
        "created_at": "2024-01-01T00:00:00",
        "updated_at": "2024-01-01T00:00:00",
    }
    mock_training_service.start_task.return_value = True

    response = client.post(f"/api/v1/training/tasks/{task_id}/start")

    assert response.status_code == 200
    mock_training_service.start_task.assert_called_once_with(task_id)


def test_stop_task_success(client, mock_training_service):
    """Test stopping a task."""
    task_id = str(uuid.uuid4())
    mock_training_service.get_task.return_value = {
        "task_id": task_id,
        "status": TaskStatus.RUNNING,
        "config": {},
        "created_at": "2024-01-01T00:00:00",
        "updated_at": "2024-01-01T00:00:00",
    }
    mock_training_service.stop_task.return_value = (True, None)

    response = client.post(f"/api/v1/training/tasks/{task_id}/stop")

    assert response.status_code == 200
    mock_training_service.stop_task.assert_called_once_with(task_id)

