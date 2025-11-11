"""
Load balancer for service selection.
"""

import sys
import random
import logging
from pathlib import Path
from typing import List, Optional, Dict
from enum import Enum

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from services.common.service_protocol import ServiceInfo

logger = logging.getLogger(__name__)


class LoadBalanceStrategy(str, Enum):
    """Load balancing strategies."""

    ROUND_ROBIN = "round_robin"
    RANDOM = "random"
    LEAST_CONNECTIONS = "least_connections"


class LoadBalancer:
    """Load balancer for service selection."""

    def __init__(self, strategy: LoadBalanceStrategy = LoadBalanceStrategy.ROUND_ROBIN):
        """
        Initialize load balancer.

        Args:
            strategy: Load balancing strategy
        """
        self.strategy = strategy
        self.round_robin_index = 0
        self.connection_counts: Dict[str, int] = {}  # Track connections per service

    def select_service(self, services: List[ServiceInfo]) -> Optional[ServiceInfo]:
        """
        Select a service from the list.

        Args:
            services: List of available services

        Returns:
            Optional[ServiceInfo]: Selected service or None if list is empty
        """
        if not services:
            return None

        # Filter to only healthy services
        healthy_services = [s for s in services if s.status.value == "healthy"]
        if not healthy_services:
            logger.warning("No healthy services available, using any available service")
            healthy_services = services

        if self.strategy == LoadBalanceStrategy.ROUND_ROBIN:
            return self._round_robin(healthy_services)
        elif self.strategy == LoadBalanceStrategy.RANDOM:
            return self._random(healthy_services)
        elif self.strategy == LoadBalanceStrategy.LEAST_CONNECTIONS:
            return self._least_connections(healthy_services)
        else:
            # This branch should not be reached as all enum values are covered above
            # However, it serves as a defensive fallback in case of unexpected values
            logger.warning(
                f"Unexpected load balance strategy: {self.strategy}, using first service"
            )
            return healthy_services[0]

    def _round_robin(self, services: List[ServiceInfo]) -> ServiceInfo:
        """Round robin selection."""
        service = services[self.round_robin_index % len(services)]
        self.round_robin_index += 1
        return service

    def _random(self, services: List[ServiceInfo]) -> ServiceInfo:
        """Random selection for load balancing (not security-critical)."""
        # Note: Using random.choice is acceptable for load balancing
        # as it's not used for security/cryptographic purposes
        return random.choice(services)

    def _least_connections(self, services: List[ServiceInfo]) -> ServiceInfo:
        """Least connections selection."""
        # Get connection count for each service (default to 0)
        counts = [self.connection_counts.get(s.service_id, 0) for s in services]
        min_index = counts.index(min(counts))
        return services[min_index]

    def increment_connections(self, service_id: str):
        """Increment connection count for a service."""
        self.connection_counts[service_id] = self.connection_counts.get(service_id, 0) + 1

    def decrement_connections(self, service_id: str):
        """Decrement connection count for a service."""
        if service_id in self.connection_counts:
            self.connection_counts[service_id] = max(0, self.connection_counts[service_id] - 1)
