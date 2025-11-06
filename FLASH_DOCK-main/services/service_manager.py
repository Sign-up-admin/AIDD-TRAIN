"""
Service manager for FLASH-DOCK.
"""
import logging
from typing import Optional, List

import sys
import os
from pathlib import Path

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from services.common.service_protocol import ServiceInfo
# Import from local FLASH_DOCK services directory
sys.path.insert(0, str(Path(__file__).parent))
from registry_client import FlashDockRegistryClient
from load_balancer import LoadBalancer, LoadBalanceStrategy

logger = logging.getLogger(__name__)


class ServiceManager:
    """Manager for service discovery and selection."""
    
    def __init__(
        self,
        registry_url: str = "http://localhost:8500",
        strategy: LoadBalanceStrategy = LoadBalanceStrategy.ROUND_ROBIN
    ):
        """
        Initialize service manager.
        
        Args:
            registry_url: Registry URL
            strategy: Load balancing strategy
        """
        self.registry_client = FlashDockRegistryClient(registry_url)
        self.load_balancer = LoadBalancer(strategy)
        self.cached_services: List[ServiceInfo] = []
        self.cache_ttl = 30  # Cache for 30 seconds
        self.last_cache_update = 0
    
    def get_compass_service(self, force_refresh: bool = False) -> Optional[ServiceInfo]:
        """
        Get a COMPASS service instance.
        
        Args:
            force_refresh: Force refresh of service list
            
        Returns:
            Optional[ServiceInfo]: Selected service or None if unavailable
        """
        import time
        
        # Refresh cache if needed
        if force_refresh or (time.time() - self.last_cache_update) > self.cache_ttl:
            self.cached_services = self.registry_client.discover_compass_services(healthy_only=True)
            self.last_cache_update = time.time()
        
        if not self.cached_services:
            logger.warning("No COMPASS services available")
            return None
        
        # Select service using load balancer
        # Note: Connection count is incremented in compass_client._make_request
        # and decremented in finally block after request completes
        service = self.load_balancer.select_service(self.cached_services)
        if service:
            logger.debug(f"Selected COMPASS service: {service.service_id} at {service.base_url}")
        
        return service
    
    def refresh_services(self):
        """Force refresh service list."""
        import time
        self.cached_services = self.registry_client.discover_compass_services(healthy_only=True)
        self.last_cache_update = time.time()
        logger.info(f"Refreshed service list: {len(self.cached_services)} service(s) available")

