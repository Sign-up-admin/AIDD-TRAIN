"""
Debug and monitoring utilities for FLASH-DOCK services.
"""
import logging
import traceback
import sys
from typing import Optional, Dict, Any
from datetime import datetime
import requests

# Configure detailed logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler('flashdock_debug.log', encoding='utf-8')
    ]
)

logger = logging.getLogger(__name__)


class ServiceMonitor:
    """Monitor service health and connection status."""
    
    def __init__(self, registry_url: str = "http://localhost:8500"):
        """Initialize monitor."""
        self.registry_url = registry_url
        self.session = requests.Session()
        self.session.timeout = 5
    
    def check_registry(self) -> Dict[str, Any]:
        """Check registry health."""
        result = {
            "available": False,
            "url": self.registry_url,
            "error": None,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            response = self.session.get(f"{self.registry_url}/health", timeout=5)
            if response.status_code == 200:
                result["available"] = True
                result["data"] = response.json()
            else:
                result["error"] = f"HTTP {response.status_code}"
        except requests.exceptions.ConnectionError as e:
            result["error"] = f"Connection failed: {str(e)}"
        except requests.exceptions.Timeout:
            result["error"] = "Request timeout"
        except Exception as e:
            result["error"] = f"Unexpected error: {str(e)}"
            logger.error(f"Registry check error: {traceback.format_exc()}")
        
        return result
    
    def check_compass_services(self) -> Dict[str, Any]:
        """Check COMPASS services in registry."""
        result = {
            "available": False,
            "services": [],
            "count": 0,
            "error": None,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            response = self.session.get(
                f"{self.registry_url}/api/v1/services",
                params={"service_name": "compass-service"},
                timeout=5
            )
            if response.status_code == 200:
                data = response.json()
                result["services"] = data.get("services", [])
                result["count"] = data.get("count", 0)
                result["available"] = result["count"] > 0
            else:
                result["error"] = f"HTTP {response.status_code}"
        except Exception as e:
            result["error"] = f"Error: {str(e)}"
            logger.error(f"COMPASS service check error: {traceback.format_exc()}")
        
        return result
    
    def check_compass_health(self, base_url: str) -> Dict[str, Any]:
        """Check COMPASS service health directly."""
        result = {
            "available": False,
            "url": base_url,
            "error": None,
            "timestamp": datetime.now().isoformat()
        }
        
        try:
            response = self.session.get(f"{base_url}/health", timeout=5)
            if response.status_code == 200:
                result["available"] = True
                result["data"] = response.json()
            else:
                result["error"] = f"HTTP {response.status_code}"
        except Exception as e:
            result["error"] = f"Error: {str(e)}"
            logger.debug(f"COMPASS health check error: {str(e)}")
        
        return result
    
    def full_diagnostic(self) -> Dict[str, Any]:
        """Run full diagnostic check."""
        logger.info("Running full service diagnostic...")
        
        diagnostic = {
            "timestamp": datetime.now().isoformat(),
            "registry": self.check_registry(),
            "compass_services": self.check_compass_services(),
            "compass_health": []
        }
        
        # Check health of each COMPASS service
        for service in diagnostic["compass_services"]["services"]:
            health = self.check_compass_health(service.get("base_url", ""))
            diagnostic["compass_health"].append(health)
        
        return diagnostic


def format_diagnostic(diagnostic: Dict[str, Any]) -> str:
    """Format diagnostic for display."""
    lines = []
    lines.append("=" * 60)
    lines.append("Service Diagnostic Report")
    lines.append("=" * 60)
    lines.append(f"Timestamp: {diagnostic['timestamp']}")
    lines.append("")
    
    # Registry status
    registry = diagnostic["registry"]
    lines.append("Registry Status:")
    lines.append(f"  URL: {registry['url']}")
    lines.append(f"  Available: {'YES' if registry['available'] else 'NO'}")
    if registry['error']:
        lines.append(f"  Error: {registry['error']}")
    lines.append("")
    
    # COMPASS services
    compass = diagnostic["compass_services"]
    lines.append("COMPASS Services in Registry:")
    lines.append(f"  Count: {compass['count']}")
    if compass['error']:
        lines.append(f"  Error: {compass['error']}")
    else:
        for service in compass['services']:
            lines.append(f"  - {service.get('service_id', 'unknown')}")
            lines.append(f"    URL: {service.get('base_url', 'unknown')}")
            lines.append(f"    Status: {service.get('status', 'unknown')}")
    lines.append("")
    
    # Direct health checks
    lines.append("Direct Health Checks:")
    for health in diagnostic["compass_health"]:
        lines.append(f"  {health['url']}: {'OK' if health['available'] else 'FAIL'}")
        if health['error']:
            lines.append(f"    Error: {health['error']}")
    lines.append("=" * 60)
    
    return "\n".join(lines)

