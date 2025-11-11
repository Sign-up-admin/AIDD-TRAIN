"""
Unit tests for CompassClient.
"""

import unittest
from unittest.mock import Mock, patch, MagicMock
import sys
from pathlib import Path

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from compass_client import CompassClient


class TestCompassClient(unittest.TestCase):
    """Test cases for CompassClient."""

    def setUp(self):
        """Set up test fixtures."""
        self.client = CompassClient()

    @patch("compass_client.requests.request")
    def test_get_service_url(self, mock_request):
        """Test getting service URL."""
        # Mock registry response
        mock_response = Mock()
        mock_response.json.return_value = {
            "services": [
                {
                    "service_id": "compass",
                    "url": "http://localhost:8080",
                    "status": "healthy",
                }
            ]
        }
        mock_request.return_value = mock_response

        url = self.client._get_service_url()
        self.assertIsNotNone(url)
        self.assertIsInstance(url, str)

    @patch("compass_client.requests.request")
    def test_list_training_tasks(self, mock_request):
        """Test listing training tasks."""
        # Mock API response
        mock_response = Mock()
        mock_response.json.return_value = {
            "tasks": [
                {
                    "task_id": "test-task-1",
                    "status": "running",
                    "config": {"execution_mode": "validation"},
                }
            ]
        }
        mock_request.return_value = mock_response

        tasks = self.client.list_training_tasks()
        self.assertIsInstance(tasks, list)
        if tasks:
            self.assertIn("task_id", tasks[0])

    @patch("compass_client.requests.request")
    def test_get_training_task(self, mock_request):
        """Test getting a specific training task."""
        # Mock API response
        mock_response = Mock()
        mock_response.json.return_value = {
            "task_id": "test-task-1",
            "status": "running",
            "config": {"execution_mode": "validation"},
        }
        mock_request.return_value = mock_response

        task = self.client.get_training_task("test-task-1")
        self.assertIsInstance(task, dict)
        self.assertEqual(task["task_id"], "test-task-1")

    def test_task_stream_client_callbacks(self):
        """Test TaskStreamClient callback handling."""
        from compass_client import TaskStreamClient

        # Create mock callbacks
        on_log_called = []
        on_resources_called = []
        on_error_called = []

        def on_log(data):
            on_log_called.append(data)

        def on_resources(data):
            on_resources_called.append(data)

        def on_error(error):
            on_error_called.append(error)

        # Create client with callbacks
        client = TaskStreamClient(
            ws_url="ws://test",
            on_log=on_log,
            on_resources=on_resources,
            on_error=on_error,
        )

        # Test message handling
        client._handle_message('{"type": "log", "data": "test log"}')
        self.assertEqual(len(on_log_called), 1)
        self.assertEqual(on_log_called[0], "test log")

        client._handle_message('{"type": "resources", "data": {"cpu": 50}}')
        self.assertEqual(len(on_resources_called), 1)
        self.assertEqual(on_resources_called[0]["cpu"], 50)

        client._handle_message('{"type": "error", "data": "test error"}')
        self.assertEqual(len(on_error_called), 1)
        self.assertIsInstance(on_error_called[0], Exception)


if __name__ == "__main__":
    unittest.main()











