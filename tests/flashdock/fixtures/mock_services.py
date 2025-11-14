"""
Mock services and responses for FlashDock testing
"""
from unittest.mock import Mock, MagicMock
from typing import Dict, Any, List, Optional
import json


class MockRegistryClient:
    """Mock registry client for testing."""
    
    def __init__(self, services: Optional[List[Dict[str, Any]]] = None):
        self.services = services or []
        self.registry_url = "http://localhost:8500"
        self.timeout = 5.0
    
    def discover_compass_services(self, healthy_only: bool = True) -> List[Any]:
        """Mock service discovery."""
        if healthy_only:
            return [s for s in self.services if s.get("status") == "healthy"]
        return self.services
    
    def check_registry_health(self) -> bool:
        """Mock health check."""
        return True


class MockCompassService:
    """Mock COMPASS service responses."""
    
    @staticmethod
    def health_check_response():
        """Mock health check response."""
        return {"status": "healthy", "version": "1.0.0"}
    
    @staticmethod
    def create_task_response(task_id: str = "test-task-123"):
        """Mock create task response."""
        return {
            "task_id": task_id,
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
    
    @staticmethod
    def get_task_response(task_id: str, status: str = "pending"):
        """Mock get task response."""
        return {
            "task_id": task_id,
            "status": status,
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
    
    @staticmethod
    def list_tasks_response(tasks: Optional[List[Dict[str, Any]]] = None):
        """Mock list tasks response."""
        if tasks is None:
            tasks = [
                {
                    "task_id": "task-1",
                    "status": "running",
                    "description": "Task 1"
                },
                {
                    "task_id": "task-2",
                    "status": "completed",
                    "description": "Task 2"
                }
            ]
        return {"tasks": tasks}
    
    @staticmethod
    def stop_task_response(task_id: str):
        """Mock stop task response."""
        return {
            "message": f"Task {task_id} stopped successfully",
            "task_id": task_id
        }
    
    @staticmethod
    def error_response(error_code: str = "ERR_1000", message: str = "Test error"):
        """Mock error response."""
        return {
            "error": message,
            "error_code": error_code,
            "detail": {
                "message": message
            }
        }
    
    @staticmethod
    def upload_dataset_response(dataset_id: str = "test-dataset-123"):
        """Mock upload dataset response."""
        return {
            "dataset_id": dataset_id,
            "name": "Test Dataset",
            "description": "Test dataset",
            "file_size": 1024000,
            "created_at": "2024-01-01T00:00:00Z"
        }
    
    @staticmethod
    def list_datasets_response(datasets: Optional[List[Dict[str, Any]]] = None):
        """Mock list datasets response."""
        if datasets is None:
            datasets = [
                {
                    "dataset_id": "dataset-1",
                    "name": "Dataset 1",
                    "file_size": 1024000
                }
            ]
        return {"datasets": datasets}
    
    @staticmethod
    def inference_response(affinity: float = -5.5):
        """Mock inference/affinity prediction response."""
        return {
            "affinity": affinity,
            "model_id": "test-model",
            "protein_path": "test_protein.pdb",
            "ligand_path": "test_ligand.sdf"
        }


def create_mock_response(status_code: int = 200, json_data: Optional[Dict[str, Any]] = None, text: Optional[str] = None):
    """Create a mock requests.Response object."""
    mock_response = Mock()
    mock_response.status_code = status_code
    mock_response.headers = {"content-type": "application/json"}
    
    if json_data:
        mock_response.json.return_value = json_data
        mock_response.text = json.dumps(json_data)
    elif text:
        mock_response.text = text
        mock_response.json.side_effect = ValueError("Not JSON")
    else:
        mock_response.json.return_value = {}
        mock_response.text = ""
    
    mock_response.raise_for_status = Mock()
    if status_code >= 400:
        from requests.exceptions import HTTPError
        mock_response.raise_for_status.side_effect = HTTPError(response=mock_response)
    
    return mock_response

