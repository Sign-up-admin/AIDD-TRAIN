"""
Service Registry Client.
Provides client interface for service registration and discovery.
"""
import logging
import requests
from typing import List, Optional, Dict
from datetime import datetime

from services.common.service_protocol import ServiceInfo, ServiceStatus
from services.registry.models import ServiceRegistrationRequest

logger = logging.getLogger(__name__)


class RegistryClient:
    """Client for interacting with service registry."""
    
    def __init__(self, registry_url: str = "http://localhost:8500"):
        """
        Initialize registry client.
        
        Args:
            registry_url: URL of the service registry
        """
        self.registry_url = registry_url.rstrip('/')
        self.session = requests.Session()
        self.session.headers.update({'Content-Type': 'application/json'})
    
    def register_service(
        self,
        service_name: str,
        host: str,
        port: int,
        metadata: Optional[Dict] = None,
        version: str = "1.0.0"
    ) -> str:
        """
        Register a service with the registry.
        
        Args:
            service_name: Name of the service
            host: Host address
            port: Port number
            metadata: Optional metadata dictionary
            version: Service version
            
        Returns:
            str: Service ID assigned by registry
        """
        request = ServiceRegistrationRequest(
            service_name=service_name,
            host=host,
            port=port,
            metadata=metadata or {},
            version=version
        )
        
        try:
            response = self.session.post(
                f"{self.registry_url}/api/v1/services/register",
                json=request.dict()
            )
            response.raise_for_status()
            data = response.json()
            service_id = data['service_id']
            logger.info(f"Registered service {service_name} with ID: {service_id}")
            return service_id
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to register service: {e}")
            raise
    
    def send_heartbeat(self, service_id: str) -> bool:
        """
        Send heartbeat to registry.
        
        Args:
            service_id: Service ID
            
        Returns:
            bool: True if successful
        """
        try:
            response = self.session.post(
                f"{self.registry_url}/api/v1/services/{service_id}/heartbeat"
            )
            response.raise_for_status()
            return True
        except requests.exceptions.RequestException as e:
            logger.warning(f"Failed to send heartbeat: {e}")
            return False
    
    def deregister_service(self, service_id: str) -> bool:
        """
        Deregister a service from the registry.
        
        Args:
            service_id: Service ID
            
        Returns:
            bool: True if successful
        """
        try:
            response = self.session.delete(
                f"{self.registry_url}/api/v1/services/{service_id}"
            )
            response.raise_for_status()
            logger.info(f"Deregistered service {service_id}")
            return True
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to deregister service: {e}")
            return False
    
    def discover_services(
        self,
        service_name: Optional[str] = None,
        status_filter: Optional[str] = None
    ) -> List[ServiceInfo]:
        """
        Discover services from the registry.
        
        Args:
            service_name: Filter by service name
            status_filter: Filter by status (healthy, unhealthy, unknown)
            
        Returns:
            List[ServiceInfo]: List of discovered services
        """
        try:
            params = {}
            if service_name:
                params['service_name'] = service_name
            if status_filter:
                params['status_filter'] = status_filter
            
            response = self.session.get(
                f"{self.registry_url}/api/v1/services",
                params=params
            )
            response.raise_for_status()
            data = response.json()
            
            services = []
            for service_dict in data['services']:
                try:
                    service = ServiceInfo.from_dict(service_dict)
                    services.append(service)
                except Exception as e:
                    logger.warning(f"Failed to parse service info: {e}")
                    continue
            
            logger.debug(f"Discovered {len(services)} services")
            return services
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to discover services: {e}")
            return []
    
    def get_service(self, service_id: str) -> Optional[ServiceInfo]:
        """
        Get service information by ID.
        
        Args:
            service_id: Service ID
            
        Returns:
            Optional[ServiceInfo]: Service information or None if not found
        """
        try:
            response = self.session.get(
                f"{self.registry_url}/api/v1/services/{service_id}"
            )
            response.raise_for_status()
            data = response.json()
            return ServiceInfo.from_dict(data)
        except requests.exceptions.RequestException as e:
            logger.error(f"Failed to get service: {e}")
            return None
    
    def check_registry_health(self) -> bool:
        """
        Check if registry is healthy.
        
        Returns:
            bool: True if registry is healthy
        """
        try:
            response = self.session.get(
                f"{self.registry_url}/health",
                timeout=5
            )
            return response.status_code == 200
        except requests.exceptions.RequestException:
            return False

