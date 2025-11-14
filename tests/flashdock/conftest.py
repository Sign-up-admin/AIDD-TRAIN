"""
FlashDock测试配置和fixtures
"""
import pytest
import sys
import os
from pathlib import Path
from unittest.mock import Mock, MagicMock
from typing import Dict, Any, List

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Add FLASH_DOCK-main to path
flashdock_root = project_root / "FLASH_DOCK-main"
sys.path.insert(0, str(flashdock_root))

# Add FLASH_DOCK-main/services to path
flashdock_services = flashdock_root / "services"
sys.path.insert(0, str(flashdock_services))


@pytest.fixture
def mock_registry_url():
    """Mock registry URL for testing."""
    return "http://localhost:8500"


@pytest.fixture
def mock_compass_url():
    """Mock COMPASS service URL for testing."""
    return "http://localhost:8080"


@pytest.fixture
def mock_flashdock_url():
    """Mock FlashDock URL for testing."""
    return "http://localhost:8501"


@pytest.fixture
def mock_service_info():
    """Create a mock ServiceInfo object."""
    from services.common.service_protocol import ServiceStatus
    
    class MockServiceInfo:
        def __init__(self, service_id="compass-1", host="localhost", port=8080, status="healthy"):
            self.service_id = service_id
            self.host = host
            self.port = port
            self.status = Mock(value=status)
            self.version = "1.0.0"
            self.last_heartbeat = None
    
    return MockServiceInfo


@pytest.fixture
def mock_services_list(mock_service_info):
    """Create a list of mock services."""
    return [
        mock_service_info("compass-1", "localhost", 8080, "healthy"),
        mock_service_info("compass-2", "localhost", 8081, "healthy"),
    ]


@pytest.fixture
def sample_training_config():
    """Sample training configuration for testing."""
    return {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Test training task"
    }


@pytest.fixture
def sample_training_task():
    """Sample training task response."""
    return {
        "task_id": "test-task-123",
        "status": "pending",
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Test training task",
        "created_at": "2024-01-01T00:00:00Z",
        "updated_at": "2024-01-01T00:00:00Z"
    }


@pytest.fixture
def sample_dataset_info():
    """Sample dataset information."""
    return {
        "dataset_id": "test-dataset-123",
        "name": "Test Dataset",
        "description": "Test dataset for testing",
        "file_size": 1024000,
        "created_at": "2024-01-01T00:00:00Z"
    }


@pytest.fixture
def mock_compass_error_response():
    """Mock COMPASS error response."""
    return {
        "error": "Test error message",
        "error_code": "ERR_1000",
        "detail": {
            "field": "test_field",
            "message": "Test detail message"
        }
    }


@pytest.fixture
def temp_test_dir(tmp_path):
    """Create a temporary directory for test files."""
    test_dir = tmp_path / "flashdock_test"
    test_dir.mkdir()
    return test_dir

