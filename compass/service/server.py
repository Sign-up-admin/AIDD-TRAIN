"""
COMPASS Service main server.
"""

import os
import logging
import argparse
import atexit
from fastapi import FastAPI, Request, status, WebSocket
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
import uvicorn

from compass.service.config import SERVICE_CONFIG
from compass.service.registry.client import CompassRegistryClient
from compass.service.routes import health, training, data, models, inference
from compass.service.exceptions import ServiceException, sanitize_error_message
from compass.service.logging_config import setup_logging
from services.common.utils import get_local_ip

# Configure unified logging
setup_logging(
    log_dir=SERVICE_CONFIG["log_dir"],
    log_level=os.getenv("LOG_LEVEL", "INFO"),
    service_name="compass-service",
)
logger = logging.getLogger(__name__)

# Create FastAPI app with enhanced documentation
app = FastAPI(
    title="COMPASS Service",
    description="""
    COMPASS (Computational Protein-Ligand Affinity Service) - Training and Inference Service
    
    ## Features
    
    * **Training Management**: Create, start, pause, stop, and monitor training tasks
    * **Model Inference**: Perform single and batch predictions for protein-ligand binding affinity
    * **Dataset Management**: Upload and manage training datasets
    * **Model Management**: List and manage trained models
    
    ## API Versions
    
    Current API version: **v1**
    
    ## Authentication
    
    Currently, the API does not require authentication. In production, consider implementing:
    - API key authentication
    - OAuth2/JWT tokens
    - Role-based access control
    
    ## Rate Limiting
    
    Default rate limits (configurable):
    - General endpoints: 100 requests/minute
    - Upload endpoints: 10 requests/minute
    - Training endpoints: 5 requests/minute
    
    ## Error Codes
    
    * `400`: Bad Request - Invalid input parameters
    * `404`: Not Found - Resource not found
    * `409`: Conflict - Resource conflict
    * `413`: Payload Too Large - File exceeds size limit
    * `429`: Too Many Requests - Rate limit exceeded
    * `500`: Internal Server Error - Server error
    * `503`: Service Unavailable - Service temporarily unavailable
    
    ## Support
    
    For issues and questions, please contact the development team.
    """,
    version="1.0.0",
    contact={
        "name": "COMPASS Development Team",
        "email": "support@compass.example.com",
    },
    license_info={
        "name": "MIT",
    },
    tags_metadata=[
        {
            "name": "health",
            "description": "Health check endpoints for service monitoring",
        },
        {
            "name": "training",
            "description": "Training task management - create, start, pause, stop, and monitor training tasks",
        },
        {
            "name": "data",
            "description": "Dataset management - upload, list, and manage training datasets",
        },
        {
            "name": "models",
            "description": "Model management - list trained models and their metadata",
        },
        {
            "name": "inference",
            "description": "Inference operations - perform single and batch predictions",
        },
    ],
)

# Enable CORS with restricted origins
# In production, replace with specific allowed origins
# For development, allow localhost origins
cors_origins = os.getenv(
    "CORS_ORIGINS",
    "http://localhost:8501,http://127.0.0.1:8501,http://localhost:3000,http://127.0.0.1:3000"
).split(",")
cors_origins = [origin.strip() for origin in cors_origins if origin.strip()]

# Allow all origins only in development mode (not recommended for production)
allow_all_origins = os.getenv("CORS_ALLOW_ALL", "false").lower() == "true"
if allow_all_origins:
    logger.warning("CORS is configured to allow all origins. This is not recommended for production!")
    cors_origins = ["*"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=cors_origins,
    allow_credentials=True,
    allow_methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
    allow_headers=["*"],
    expose_headers=["*"],
)

# Add rate limiting middleware
from compass.service.middleware.rate_limit import RateLimitMiddleware
import os

# Configure rate limits
default_limit = int(os.getenv("RATE_LIMIT_DEFAULT", "100"))
default_window = int(os.getenv("RATE_LIMIT_WINDOW", "60"))

# Per-endpoint rate limits
per_endpoint_limits = {
    "/api/v1/training/tasks": {
        "limit": int(os.getenv("RATE_LIMIT_TRAINING", "30")),  # 从5增加到30，允许更频繁的请求
        "window": int(os.getenv("RATE_LIMIT_TRAINING_WINDOW", "60")),
    },
    "/api/v1/data/upload": {
        "limit": int(os.getenv("RATE_LIMIT_UPLOAD", "10")),
        "window": int(os.getenv("RATE_LIMIT_UPLOAD_WINDOW", "60")),
    },
    "/api/v1/inference": {
        "limit": int(os.getenv("RATE_LIMIT_INFERENCE", "50")),
        "window": int(os.getenv("RATE_LIMIT_INFERENCE_WINDOW", "60")),
    },
}

app.add_middleware(
    RateLimitMiddleware,
    default_limit=default_limit,
    default_window=default_window,
    per_endpoint_limits=per_endpoint_limits,
)

# Add security headers middleware
from compass.service.middleware.security_headers import SecurityHeadersMiddleware

app.add_middleware(SecurityHeadersMiddleware)

# Add authentication middleware (optional, enabled via AUTH_ENABLED env var)
from compass.service.middleware.auth import AuthMiddleware

app.add_middleware(AuthMiddleware)

# Add metrics middleware
from compass.service.middleware.metrics import MetricsMiddleware, get_metrics_collector

app.add_middleware(MetricsMiddleware)


# Exception handlers
@app.exception_handler(ServiceException)
async def service_exception_handler(request: Request, exc: ServiceException):
    """Handle service exceptions."""
    logger.error(
        f"Service error: {exc.message}",
        extra={
            "status_code": exc.status_code,
            "error_code": exc.error_code.value if exc.error_code else None,
            "path": request.url.path,
        },
    )

    response_content = {"error": exc.message, "status_code": exc.status_code}

    # Add error code if available
    if exc.error_code:
        response_content["error_code"] = exc.error_code.value

    # Add detail if available
    if exc.detail:
        response_content["detail"] = exc.detail

    return JSONResponse(status_code=exc.status_code, content=response_content)


@app.exception_handler(RequestValidationError)
async def validation_exception_handler(request: Request, exc: RequestValidationError):
    """Handle validation errors."""
    logger.warning(f"Validation error: {exc.errors()}", extra={"path": str(request.url)})
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content={
            "error": "Validation error",
            "detail": exc.errors(),
            "path": str(request.url.path),
        },
    )


@app.exception_handler(Exception)
async def general_exception_handler(request: Request, exc: Exception):
    """Handle unexpected exceptions."""
    error_message = sanitize_error_message(exc, include_details=True)
    logger.error(
        f"Unexpected error: {error_message}", extra={"path": str(request.url)}, exc_info=True
    )
    # Don't expose internal error details to users
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content={"error": "Internal server error", "path": str(request.url.path)},
    )


# Include routers
app.include_router(health.router)
app.include_router(training.router)
app.include_router(data.router)
app.include_router(models.router)
app.include_router(inference.router)

# Debug: Test if we can access the WebSocket route from training router
# This is just for debugging - the route should be registered via include_router
try:
    # Check if WebSocket route exists in training router
    training_ws_routes = [r for r in training.router.routes if hasattr(r, 'path') and 'stream' in r.path]
    if training_ws_routes:
        logger.info(f"Found {len(training_ws_routes)} WebSocket route(s) in training router")
        for route in training_ws_routes:
            logger.info(f"  WebSocket route: {route.path}, type: {type(route).__name__}")
    else:
        logger.warning("No WebSocket routes found in training router")
        # List all routes in training router for debugging
        all_training_routes = [r for r in training.router.routes if hasattr(r, 'path')]
        logger.info(f"Training router has {len(all_training_routes)} total routes")
        for route in all_training_routes[:10]:  # Show first 10
            logger.info(f"  Route: {route.path}, type: {type(route).__name__}")
except Exception as e:
    logger.error(f"Error checking WebSocket routes: {e}", exc_info=True)

# Debug: Log all registered routes including WebSocket
logger.info(f"Registered routes in app (total: {len(app.routes)}):")
websocket_found = False
for route in app.routes:
    if hasattr(route, 'path'):
        # Check if it's a WebSocket route
        route_type = "HTTP"
        route_class = type(route).__name__
        if 'WebSocket' in route_class or 'websocket' in route_class.lower():
            route_type = "WebSocket"
            websocket_found = True
        if 'stream' in route.path.lower():
            logger.info(f"  {route_type} ({route_class}): {route.path}")
            websocket_found = True
logger.info(f"WebSocket routes found in app: {websocket_found}")

# Registry client
registry_client: CompassRegistryClient = None


@app.on_event("startup")
async def startup_event():
    """Initialize service on startup."""
    global registry_client

    # Create backup on startup if configured
    backup_on_startup = os.getenv("BACKUP_ON_STARTUP", "false").lower() == "true"
    if backup_on_startup:
        try:
            from compass.service.utils.backup import get_backup_manager

            backup_manager = get_backup_manager()
            if backup_manager:
                backup_path = backup_manager.create_backup(
                    datasets_dir=SERVICE_CONFIG.get("data_dir"),
                    models_dir=SERVICE_CONFIG.get("checkpoint_dir"),
                    database_file=os.getenv("REGISTRY_DB_PATH", "registry.db"),
                )
                if backup_path:
                    logger.info(f"Startup backup created: {backup_path}")
        except Exception as e:
            logger.warning(f"Failed to create startup backup: {e}")

    # Register with service registry (with retry)
    try:
        host = SERVICE_CONFIG["host"]
        if host == "0.0.0.0":
            host = get_local_ip()
        port = SERVICE_CONFIG["port"]

        registry_client = CompassRegistryClient()

        # Register with retry mechanism
        max_retries = int(os.getenv("REGISTRY_RETRY_MAX", "5"))
        service_id = registry_client.register(
            host=host,
            port=port,
            metadata={
                "max_workers": SERVICE_CONFIG["max_workers"],
                "device": "cuda" if os.getenv("CUDA_VISIBLE_DEVICES") else "cpu",
            },
            max_retries=max_retries,
            initial_delay=1.0,
            max_delay=60.0,
            backoff_factor=2.0,
        )
        SERVICE_CONFIG["service_id"] = service_id

        # Register cleanup on exit
        atexit.register(cleanup)

        logger.info(f"COMPASS Service started on {host}:{port}")
        logger.info(f"Registered with service registry: {service_id}")
    except Exception as e:
        logger.error(f"Failed to register with service registry after retries: {e}")
        logger.warning("Service will continue without registry registration")
        # Still set registry_client for potential later registration
        registry_client = CompassRegistryClient()


@app.on_event("shutdown")
async def shutdown_event():
    """Cleanup on shutdown."""
    cleanup()
    logger.info("COMPASS Service stopped")


def cleanup():
    """Cleanup function for service deregistration."""
    global registry_client
    if registry_client:
        try:
            # Ensure heartbeat thread is stopped
            registry_client.stop_heartbeat()
            registry_client.deregister()
        except Exception as e:
            logger.error(f"Failed to deregister from service registry: {e}")
        finally:
            # Ensure running flag is cleared
            if registry_client:
                registry_client.running = False


@app.get("/")
async def root():
    """Root endpoint."""
    return {
        "service": "COMPASS Service",
        "version": "1.0.0",
        "service_id": SERVICE_CONFIG.get("service_id"),
        "status": "running",
    }


def main():
    """Main entry point for COMPASS service."""
    parser = argparse.ArgumentParser(description="COMPASS Service")
    parser.add_argument("--host", default=SERVICE_CONFIG["host"], help="Host to bind to")
    parser.add_argument("--port", type=int, default=SERVICE_CONFIG["port"], help="Port to bind to")
    parser.add_argument(
        "--registry-url", default=SERVICE_CONFIG["registry_url"], help="Service registry URL"
    )
    args = parser.parse_args()

    # Update config
    SERVICE_CONFIG["host"] = args.host
    SERVICE_CONFIG["port"] = args.port
    SERVICE_CONFIG["registry_url"] = args.registry_url

    logger.info(f"Starting COMPASS Service on {args.host}:{args.port}")
    logger.info(f"Registry URL: {args.registry_url}")

    uvicorn.run(app, host=args.host, port=args.port)


if __name__ == "__main__":
    main()
