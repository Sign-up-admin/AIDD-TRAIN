"""
Registry client for FLASH-DOCK.
"""
import sys
import logging
from pathlib import Path
from typing import List, Optional

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from services.registry.client import RegistryClient
from services.common.service_protocol import ServiceInfo

logger = logging.getLogger(__name__)


class FlashDockRegistryClient:
    """Registry client for FLASH-DOCK service discovery."""
    
    def __init__(self, registry_url: str = "http://localhost:8500"):
        """
        Initialize registry client.
        
        Args:
            registry_url: Registry URL
        """
        self.registry_url = registry_url
        self.client = RegistryClient(registry_url)
    
    def discover_compass_services(self, healthy_only: bool = True) -> List[ServiceInfo]:
        """
        Discover COMPASS services.
        
        Args:
            healthy_only: Only return healthy services
            
        Returns:
            List[ServiceInfo]: List of discovered services
        """
        try:
            status_filter = "healthy" if healthy_only else None
            services = self.client.discover_services(
                service_name="compass-service",
                status_filter=status_filter
            )
            logger.info(f"Discovered {len(services)} COMPASS service(s)")
            return services
        except Exception as e:
            logger.error(f"Failed to discover services: {e}")
            return []
    
    def check_registry_available(self) -> bool:
        """
        Check if registry is available.
        
        Returns:
            bool: True if registry is available
        """
        return self.client.check_registry_health()

