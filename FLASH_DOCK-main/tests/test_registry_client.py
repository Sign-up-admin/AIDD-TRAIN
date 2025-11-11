"""
Unit tests for FlashDockRegistryClient.
"""

import unittest
import sys
from pathlib import Path
from unittest.mock import Mock, patch

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from registry_client import FlashDockRegistryClient


class TestFlashDockRegistryClient(unittest.TestCase):
    """Test cases for FlashDockRegistryClient."""

    def setUp(self):
        """Set up test fixtures."""
        self.client = FlashDockRegistryClient(registry_url="http://localhost:8500")

    @patch("registry_client.requests.get")
    def test_list_services(self, mock_get):
        """Test listing services."""
        # Mock registry response
        mock_response = Mock()
        mock_response.json.return_value = {
            "services": [
                {
                    "service_id": "compass",
                    "name": "COMPASS Service",
                    "url": "http://localhost:8080",
                    "status": "healthy",
                }
            ]
        }
        mock_get.return_value = mock_response

        services = self.client.list_services()
        self.assertIsInstance(services, list)
        if services:
            self.assertIn("service_id", services[0])

    @patch("registry_client.requests.get")
    def test_get_service(self, mock_get):
        """Test getting a specific service."""
        # Mock registry response
        mock_response = Mock()
        mock_response.json.return_value = {
            "service_id": "compass",
            "name": "COMPASS Service",
            "url": "http://localhost:8080",
            "status": "healthy",
        }
        mock_get.return_value = mock_response

        service = self.client.get_service("compass")
        self.assertIsNotNone(service)
        self.assertEqual(service["service_id"], "compass")

    @patch("registry_client.requests.get")
    def test_service_not_found(self, mock_get):
        """Test behavior when service is not found."""
        # Mock 404 response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response

        service = self.client.get_service("nonexistent")
        self.assertIsNone(service)


if __name__ == "__main__":
    unittest.main()











