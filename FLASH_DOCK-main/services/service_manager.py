"""
Service manager for FLASH-DOCK.
"""
import sys
import logging
from pathlib import Path
from typing import Optional

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from registry_client import FlashDockRegistryClient

logger = logging.getLogger(__name__)


class ServiceManager:
    """Manager for service discovery and monitoring."""
    
    def __init__(
        self, 
        registry_url: str = "http://localhost:8500",
        timeout: float = 5.0
    ):
        """
        Initialize service manager.
        
        Args:
            registry_url: Registry URL for service discovery
            timeout: Request timeout in seconds (default: 5.0)
        """
        self.registry_client = FlashDockRegistryClient(registry_url=registry_url, timeout=timeout)
        self._services_cache = []
        logger.info(f"ServiceManager initialized with registry at {registry_url}")
    
    def refresh_services(self) -> None:
        """
        Refresh the cached list of COMPASS services.
        
        This method discovers COMPASS services from the registry and updates
        the internal cache. It can be called manually to refresh the service list.
        """
        try:
            self._services_cache = self.registry_client.discover_compass_services(healthy_only=False)
            logger.info(f"Refreshed services: {len(self._services_cache)} COMPASS service(s) found")
        except Exception as e:
            logger.error(f"Failed to refresh services: {e}", exc_info=True)
            # Keep existing cache if refresh fails
            if not self._services_cache:
                self._services_cache = []
