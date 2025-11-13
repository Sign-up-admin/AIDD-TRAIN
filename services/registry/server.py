"""
Service Registry Server.
Provides service registration, discovery, and health monitoring.
"""

import os
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).resolve().parent.parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

import logging
import argparse
from contextlib import asynccontextmanager
from typing import Dict, List, Optional
from datetime import datetime
from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware
import uvicorn

from services.registry.models import (
    ServiceRegistrationRequest,
    ServiceRegistrationResponse,
    ServiceQueryResponse,
    HealthCheckResponse,
)
from services.registry.health_checker import HealthChecker
from services.registry.storage import ServiceRegistryStorage
from services.common.service_protocol import ServiceInfo, ServiceStatus
from services.common.utils import generate_service_id

# Configure unified logging
try:
    from compass.service.logging_config import setup_logging

    log_dir = os.getenv("REGISTRY_LOG_DIR", "logs")
    log_level = os.getenv("LOG_LEVEL", "INFO")
    setup_logging(
        log_dir=log_dir,
        log_level=log_level,
        service_name="registry",
        enable_file=True,
        enable_console=True,
    )
except ImportError:
    # Fallback to basic logging if unified config not available
    log_dir = os.getenv("REGISTRY_LOG_DIR", "logs")
    log_level = os.getenv("LOG_LEVEL", "INFO")
    os.makedirs(log_dir, exist_ok=True)
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s",
        handlers=[
            logging.FileHandler(os.path.join(log_dir, "registry.log"), encoding="utf-8"),
            logging.StreamHandler(),
        ],
    )

logger = logging.getLogger(__name__)

# Persistent storage for services
db_path = os.getenv("REGISTRY_DB_PATH", "registry.db")
storage = ServiceRegistryStorage(db_path=db_path)

# In-memory service cache (loaded from database)
services: Dict[str, ServiceInfo] = {}

# Health checker
health_checker: Optional[HealthChecker] = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Lifespan context manager for startup and shutdown events."""
    # Startup
    global health_checker

    # Load services from database
    services.update(storage.load_all_services())
    logger.info(f"Loaded {len(services)} services from database")

    # Cleanup stale services
    stale_count = storage.cleanup_stale_services(max_age_seconds=300)
    if stale_count > 0:
        # Reload after cleanup
        services.clear()
        services.update(storage.load_all_services())

    check_interval = int(os.getenv("REGISTRY_CHECK_INTERVAL", "10"))
    timeout = int(os.getenv("REGISTRY_TIMEOUT", "5"))
    health_checker = HealthChecker(check_interval=check_interval, timeout=timeout)
    health_checker.start_background_checking(services)
    logger.info("Service Registry started")

    yield

    # Shutdown - persist current state
    for service_info in services.values():
        try:
            storage.register_service(service_info)
        except Exception as e:
            logger.warning(f"Failed to persist service {service_info.service_id} on shutdown: {e}")

    if health_checker:
        health_checker.stop()
    logger.info("Service Registry stopped")


# Create FastAPI app with lifespan
app = FastAPI(
    title="Service Registry",
    description="Service registration and discovery service",
    version="1.0.0",
    lifespan=lifespan,
)

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def root():
    """Root endpoint."""
    return {"service": "Service Registry", "version": "1.0.0", "registered_services": len(services)}


@app.get("/health", response_model=HealthCheckResponse)
async def health():
    """Health check endpoint."""
    return HealthCheckResponse(status="healthy", message="Service registry is running")


@app.post("/api/v1/services/register", response_model=ServiceRegistrationResponse)
async def register_service(request: ServiceRegistrationRequest):
    """
    Register a new service.

    Args:
        request: Service registration request

    Returns:
        ServiceRegistrationResponse: Registration response with service ID
    """
    # Generate service ID
    service_id = generate_service_id(request.service_name)

    # Construct base URL
    base_url = f"http://{request.host}:{request.port}"

    # Create service info
    service_info = ServiceInfo(
        service_name=request.service_name,
        service_id=service_id,
        host=request.host,
        port=request.port,
        base_url=base_url,
        status=ServiceStatus.UNKNOWN,
        metadata=request.metadata,
        version=request.version,
        registered_at=datetime.now(),
        last_heartbeat=datetime.now(),
    )

    # Register service in memory and persist to database
    services[service_id] = service_info
    storage.register_service(service_info)

    logger.info(f"Registered service: {service_id} at {base_url}")

    # Ensure registered_at is not None (it's set in __post_init__ but MyPy doesn't know)
    registered_at = service_info.registered_at
    if registered_at is None:
        registered_at = datetime.now()

    return ServiceRegistrationResponse(
        service_id=service_id,
        message=f"Service {request.service_name} registered successfully",
        registered_at=registered_at,
    )


@app.post("/api/v1/services/{service_id}/heartbeat")
async def heartbeat(service_id: str):
    """
    Update service heartbeat.

    Args:
        service_id: Service ID
    """
    if service_id not in services:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Service {service_id} not found"
        )

    now = datetime.now()
    services[service_id].last_heartbeat = now
    storage.update_heartbeat(service_id, now)
    logger.debug(f"Heartbeat updated for service {service_id}")

    return {"message": "Heartbeat updated", "timestamp": now.isoformat()}


@app.delete("/api/v1/services/{service_id}")
async def deregister_service(service_id: str):
    """
    Deregister a service.

    Args:
        service_id: Service ID
    """
    if service_id not in services:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Service {service_id} not found"
        )

    services.pop(service_id)
    storage.delete_service(service_id)
    logger.info(f"Deregistered service: {service_id}")

    return {"message": f"Service {service_id} deregistered successfully"}


@app.get("/api/v1/services", response_model=ServiceQueryResponse)
async def list_services(service_name: Optional[str] = None, status_filter: Optional[str] = None):
    """
    List all registered services.

    Args:
        service_name: Filter by service name
        status_filter: Filter by status (healthy, unhealthy, unknown)

    Returns:
        ServiceQueryResponse: List of services
    """
    filtered_services = list(services.values())

    # Filter by service name
    if service_name:
        filtered_services = [s for s in filtered_services if s.service_name == service_name]

    # Filter by status
    if status_filter:
        try:
            status_enum = ServiceStatus(status_filter.lower())
            filtered_services = [s for s in filtered_services if s.status == status_enum]
        except ValueError:
            pass  # Invalid status filter, ignore

    # Convert to dict for JSON serialization
    service_list = [s.to_dict() for s in filtered_services]

    return ServiceQueryResponse(services=service_list, count=len(service_list))


@app.get("/api/v1/services/{service_id}")
async def get_service(service_id: str):
    """
    Get service information by ID.

    Args:
        service_id: Service ID

    Returns:
        ServiceInfo: Service information
    """
    if service_id not in services:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Service {service_id} not found"
        )

    return services[service_id].to_dict()


def main():
    """Main entry point for registry server."""
    parser = argparse.ArgumentParser(description="Service Registry Server")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8500, help="Port to bind to")
    parser.add_argument(
        "--check-interval", type=int, default=10, help="Health check interval (seconds)"
    )
    parser.add_argument("--timeout", type=int, default=5, help="Health check timeout (seconds)")
    args = parser.parse_args()

    # Set environment variables for health checker
    os.environ["REGISTRY_CHECK_INTERVAL"] = str(args.check_interval)
    os.environ["REGISTRY_TIMEOUT"] = str(args.timeout)

    logger.info(f"Starting Service Registry on {args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
