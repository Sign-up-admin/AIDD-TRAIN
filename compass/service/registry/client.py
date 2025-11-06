"""
Registry client for COMPASS service registration.
"""

import logging
import threading
import time
from typing import Optional

from services.registry.client import RegistryClient
from services.common.utils import get_local_ip
from compass.service.config import SERVICE_CONFIG

logger = logging.getLogger(__name__)


class CompassRegistryClient:
    """Registry client for COMPASS service."""

    def __init__(self, registry_url: Optional[str] = None):
        """
        Initialize registry client.

        Args:
            registry_url: Registry URL (defaults to config)
        """
        self.registry_url = registry_url or SERVICE_CONFIG["registry_url"]
        self.client = RegistryClient(self.registry_url)
        self.service_id: Optional[str] = None
        self.heartbeat_thread: Optional[threading.Thread] = None
        self.running = False

    def register(
        self,
        host: str,
        port: int,
        metadata: Optional[dict] = None,
        max_retries: int = 5,
        initial_delay: float = 1.0,
        max_delay: float = 60.0,
        backoff_factor: float = 2.0,
    ) -> str:
        """
        Register COMPASS service with registry with retry mechanism.

        Args:
            host: Service host
            port: Service port
            metadata: Optional metadata
            max_retries: Maximum number of retry attempts
            initial_delay: Initial delay between retries (seconds)
            max_delay: Maximum delay between retries (seconds)
            backoff_factor: Exponential backoff factor

        Returns:
            str: Service ID

        Raises:
            Exception: If registration fails after all retries
        """
        delay = initial_delay
        last_exception = None

        for attempt in range(max_retries + 1):
            try:
                service_id = self.client.register_service(
                    service_name=SERVICE_CONFIG["service_name"],
                    host=host,
                    port=port,
                    metadata=metadata or {},
                    version="1.0.0",
                )
                self.service_id = service_id
                logger.info(f"Registered COMPASS service with ID: {service_id}")

                # Start heartbeat thread
                self.start_heartbeat()

                return service_id
            except Exception as e:
                last_exception = e
                if attempt < max_retries:
                    logger.warning(
                        f"Failed to register service (attempt {attempt + 1}/{max_retries + 1}): {e}. "
                        f"Retrying in {delay:.1f}s..."
                    )
                    time.sleep(delay)
                    delay = min(delay * backoff_factor, max_delay)
                else:
                    logger.error(
                        f"Failed to register service after {max_retries + 1} attempts: {e}"
                    )

        # All retries failed
        raise Exception(
            f"Failed to register service after {max_retries + 1} attempts"
        ) from last_exception

    def start_heartbeat(self):
        """Start heartbeat thread."""
        if self.heartbeat_thread and self.heartbeat_thread.is_alive():
            logger.warning("Heartbeat thread already running")
            return

        self.running = True

        def heartbeat_loop():
            logger.debug("Heartbeat thread started")
            while self.running and self.service_id:
                try:
                    self.client.send_heartbeat(self.service_id)
                    # Sleep with periodic checks to allow quick shutdown
                    sleep_interval = SERVICE_CONFIG["heartbeat_interval"]
                    slept = 0
                    while slept < sleep_interval and self.running and self.service_id:
                        time.sleep(min(1.0, sleep_interval - slept))
                        slept += 1.0
                except Exception as e:
                    logger.warning(f"Heartbeat error: {e}")
                    if self.running:
                        # Only retry if still running
                        time.sleep(5)  # Retry after 5 seconds

            logger.debug("Heartbeat thread stopped")

        self.heartbeat_thread = threading.Thread(
            target=heartbeat_loop, daemon=True, name="HeartbeatThread"
        )
        self.heartbeat_thread.start()
        logger.info("Started heartbeat thread")

    def stop_heartbeat(self):
        """Stop heartbeat thread."""
        if not self.heartbeat_thread:
            return

        if self.running:
            logger.info("Stopping heartbeat thread...")
            self.running = False

            # Wait for thread to finish with timeout
            if self.heartbeat_thread.is_alive():
                self.heartbeat_thread.join(timeout=10)
                if self.heartbeat_thread.is_alive():
                    logger.warning("Heartbeat thread did not stop within timeout")
                else:
                    logger.info("Heartbeat thread stopped successfully")

        # Clear thread reference
        self.heartbeat_thread = None

    def deregister(self):
        """Deregister service from registry."""
        if self.service_id:
            # Stop heartbeat first
            self.stop_heartbeat()

            # Deregister service
            try:
                self.client.deregister_service(self.service_id)
                logger.info(f"Deregistered COMPASS service: {self.service_id}")
            except Exception as e:
                logger.error(f"Error during deregistration: {e}")

            self.service_id = None
