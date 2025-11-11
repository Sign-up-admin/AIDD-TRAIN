"""
Simple API key authentication middleware.
"""

import os
import logging
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse
from fastapi import status

logger = logging.getLogger(__name__)

# API key from environment variable
API_KEY = os.getenv("API_KEY", "")
# Enable/disable authentication
AUTH_ENABLED = os.getenv("AUTH_ENABLED", "false").lower() == "true"

# Public endpoints that don't require authentication
PUBLIC_ENDPOINTS = [
    "/health",
    "/health/ready",
    "/docs",
    "/openapi.json",
    "/redoc",
    "/metrics",
]


class AuthMiddleware(BaseHTTPMiddleware):
    """Simple API key authentication middleware."""

    async def dispatch(self, request: Request, call_next):
        """Check API key for protected endpoints."""
        # Skip authentication if disabled
        if not AUTH_ENABLED:
            return await call_next(request)

        # Skip authentication for public endpoints
        if request.url.path in PUBLIC_ENDPOINTS:
            return await call_next(request)

        # Skip authentication for OPTIONS requests (CORS preflight)
        if request.method == "OPTIONS":
            return await call_next(request)

        # Get API key from header
        api_key = request.headers.get("X-API-Key") or request.headers.get("Authorization")

        # Extract key from "Bearer <key>" format
        if api_key and api_key.startswith("Bearer "):
            api_key = api_key[7:]

        # Check API key
        if not api_key or api_key != API_KEY:
            logger.warning(f"Authentication failed for {request.url.path} from {request.client.host}")
            return JSONResponse(
                status_code=status.HTTP_401_UNAUTHORIZED,
                content={
                    "error": "Unauthorized",
                    "message": "Invalid or missing API key",
                    "status_code": 401,
                },
            )

        return await call_next(request)





