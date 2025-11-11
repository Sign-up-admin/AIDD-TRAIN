"""
Health checker for services.
"""

import os
import time
import requests
import logging
from datetime import datetime, timedelta
from typing import Dict, List, Optional
from threading import Thread

from services.common.service_protocol import ServiceInfo, ServiceStatus

logger = logging.getLogger(__name__)


class HealthChecker:
    """Health checker for registered services."""

    def __init__(self, check_interval: int = 10, timeout: int = 5):
        """
        Initialize health checker.

        Args:
            check_interval: Interval between health checks in seconds
            timeout: Request timeout in seconds
        """
        self.check_interval = check_interval
        self.timeout = timeout
        self.running = False
        self.thread: Optional[Thread] = None

    def check_service_health(self, service: ServiceInfo) -> ServiceStatus:
        """
        Check health of a single service.

        Args:
            service: Service information

        Returns:
            ServiceStatus: Health status
        """
        try:
            health_url = f"{service.base_url}/health"
            response = requests.get(health_url, timeout=self.timeout)
            if response.status_code == 200:
                return ServiceStatus.HEALTHY
            else:
                return ServiceStatus.UNHEALTHY
        except requests.exceptions.Timeout:
            logger.warning(f"Health check timeout for {service.service_id}")
            return ServiceStatus.UNHEALTHY
        except requests.exceptions.ConnectionError:
            logger.warning(f"Health check connection error for {service.service_id}")
            return ServiceStatus.UNHEALTHY
        except Exception as e:
            logger.error(f"Health check error for {service.service_id}: {e}")
            return ServiceStatus.UNHEALTHY

    def check_all_services(self, services: Dict[str, ServiceInfo]) -> None:
        """
        Check health of all services.

        Args:
            services: Dictionary of service_id -> ServiceInfo
        """
        from services.registry.storage import ServiceRegistryStorage

        # Get storage instance (will be injected or use default)
        storage = getattr(self, "_storage", None)
        if storage is None:
            # Use default storage path
            db_path = os.getenv("REGISTRY_DB_PATH", "registry.db")
            storage = ServiceRegistryStorage(db_path=db_path)

        for service_id, service in list(services.items()):
            status = self.check_service_health(service)
            service.status = status
            now = datetime.now()
            service.last_heartbeat = now

            # Update in database
            try:
                storage.update_service_status(service_id, status)
                storage.update_heartbeat(service_id, now)
            except Exception as e:
                logger.warning(f"Failed to update service {service_id} in database: {e}")

            logger.debug(f"Service {service_id} health: {status.value}")

    def start_background_checking(self, services: Dict[str, ServiceInfo]):
        """
        Start background health checking in a separate thread.

        Args:
            services: Dictionary of service_id -> ServiceInfo
        """
        if self.running:
            return

        self.running = True

        def check_loop():
            while self.running:
                try:
                    self.check_all_services(services)
                except Exception as e:
                    logger.error(f"Error in health check loop: {e}")
                time.sleep(self.check_interval)

        self.thread = Thread(target=check_loop, daemon=True)
        self.thread.start()
        logger.info(f"Started background health checking (interval: {self.check_interval}s)")

    def stop(self):
        """Stop background health checking."""
        self.running = False
        if self.thread:
            self.thread.join(timeout=5)
        logger.info("Stopped background health checking")
