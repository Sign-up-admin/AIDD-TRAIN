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
import requests

logger = logging.getLogger(__name__)


class FlashDockRegistryClient:
    """Registry client for FLASH-DOCK service discovery."""
    
    def __init__(self, registry_url: str = "http://localhost:8500", timeout: float = 5.0):
        """
        Initialize registry client.
        
        Args:
            registry_url: Registry URL
            timeout: Request timeout in seconds (default: 5.0)
        """
        self.registry_url = registry_url
        self.timeout = timeout
        self.client = RegistryClient(registry_url, timeout=timeout)
    
    def discover_compass_services(self, healthy_only: bool = True) -> List[ServiceInfo]:
        """
        Discover COMPASS services with timeout and error handling.
        
        Args:
            healthy_only: Only return healthy services
            
        Returns:
            List[ServiceInfo]: List of discovered services (empty list on error)
        """
        try:
            status_filter = "healthy" if healthy_only else None
            services = self.client.discover_services(
                service_name="compass-service",
                status_filter=status_filter
            )
            logger.info(f"Discovered {len(services)} COMPASS service(s)")
            return services
        except requests.exceptions.Timeout as e:
            logger.error(f"Timeout while discovering COMPASS services (timeout: {self.timeout}s): {e}")
            logger.warning("Registry may be unavailable or slow to respond")
            return []
        except requests.exceptions.ConnectionError as e:
            logger.error(f"Connection error while discovering COMPASS services: {e}")
            logger.warning(f"Could not connect to registry at {self.registry_url}")
            return []
        except requests.exceptions.RequestException as e:
            logger.error(f"Request error while discovering COMPASS services: {e}")
            return []
        except Exception as e:
            logger.error(f"Unexpected error while discovering services: {e}", exc_info=True)
            return []
    
    def check_registry_available(self) -> bool:
        """
        Check if registry is available.
        
        Returns:
            bool: True if registry is available
        """
        try:
            logger.debug(f"Checking registry health at {self.registry_url} (timeout: {self.timeout}s)")
            result = self.client.check_registry_health()
            logger.debug(f"Registry health check result: {result}")
            return result
        except Exception as e:
            logger.warning(f"Error checking registry health: {e}", exc_info=True)
            return False

