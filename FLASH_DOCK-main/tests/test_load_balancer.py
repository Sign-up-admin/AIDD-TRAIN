"""
Unit tests for LoadBalancer.
"""

import unittest
import sys
from pathlib import Path
from unittest.mock import Mock

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from load_balancer import LoadBalancer, LoadBalanceStrategy
from services.common.service_protocol import ServiceInfo


class TestLoadBalancer(unittest.TestCase):
    """Test cases for LoadBalancer."""

    def setUp(self):
        """Set up test fixtures."""
        self.services = [
            ServiceInfo(
                service_id="service1",
                name="Service 1",
                url="http://localhost:8001",
                status="healthy",
            ),
            ServiceInfo(
                service_id="service2",
                name="Service 2",
                url="http://localhost:8002",
                status="healthy",
            ),
            ServiceInfo(
                service_id="service3",
                name="Service 3",
                url="http://localhost:8003",
                status="healthy",
            ),
        ]

    def test_round_robin_selection(self):
        """Test round robin load balancing."""
        balancer = LoadBalancer(strategy=LoadBalanceStrategy.ROUND_ROBIN)
        selected = [balancer.select_service(self.services) for _ in range(6)]
        # Should cycle through services
        self.assertEqual(selected[0].service_id, selected[3].service_id)
        self.assertEqual(selected[1].service_id, selected[4].service_id)
        self.assertEqual(selected[2].service_id, selected[5].service_id)

    def test_random_selection(self):
        """Test random load balancing."""
        balancer = LoadBalancer(strategy=LoadBalanceStrategy.RANDOM)
        selected = balancer.select_service(self.services)
        self.assertIn(selected, self.services)

    def test_least_connections_selection(self):
        """Test least connections load balancing."""
        balancer = LoadBalancer(strategy=LoadBalanceStrategy.LEAST_CONNECTIONS)
        # Set connection counts
        balancer.connection_counts = {
            "service1": 5,
            "service2": 2,
            "service3": 3,
        }
        selected = balancer.select_service(self.services)
        # Should select service with least connections
        self.assertEqual(selected.service_id, "service2")

    def test_empty_services_list(self):
        """Test behavior with empty services list."""
        balancer = LoadBalancer()
        with self.assertRaises(IndexError):
            balancer.select_service([])

    def test_single_service(self):
        """Test behavior with single service."""
        balancer = LoadBalancer()
        selected = balancer.select_service([self.services[0]])
        self.assertEqual(selected.service_id, "service1")


if __name__ == "__main__":
    unittest.main()











