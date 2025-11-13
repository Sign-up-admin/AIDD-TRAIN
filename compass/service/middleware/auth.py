"""
API key authentication middleware with support for multiple keys and key rotation.
"""

import os
import logging
import secrets
import time
from typing import List, Set, Optional
from collections import defaultdict, deque
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse
from fastapi import status

logger = logging.getLogger(__name__)

# API key from environment variable (single key for backward compatibility)
API_KEY = os.getenv("API_KEY", "")
# Multiple API keys (comma-separated)
API_KEYS_STR = os.getenv("API_KEYS", "")
# Enable/disable authentication
AUTH_ENABLED = os.getenv("AUTH_ENABLED", "false").lower() == "true"
# Check if we're in production
IS_PRODUCTION = os.getenv("ENVIRONMENT", "development").lower() == "production"
# Force authentication for critical endpoints in production
FORCE_AUTH_CRITICAL = os.getenv("FORCE_AUTH_CRITICAL", "true").lower() == "true"

# Critical endpoints that require authentication in production
CRITICAL_ENDPOINTS = [
    "/api/v1/training/tasks",
    "/api/v1/data/upload",
    "/api/v1/data/datasets",
    "/api/v1/inference",
    "/api/v1/models",
]

# Public endpoints that don't require authentication
PUBLIC_ENDPOINTS = [
    "/health",
    "/health/ready",
    "/docs",
    "/openapi.json",
    "/redoc",
    "/metrics",
    "/",  # Root endpoint
]

# Authentication failure tracking for monitoring
_auth_failures: deque = deque(maxlen=1000)  # Track last 1000 failures
_failure_counts: defaultdict = defaultdict(int)  # Per-IP failure counts
_failure_window_seconds = 300  # 5 minutes window


def load_api_keys() -> Set[str]:
    """
    Load API keys from environment variables.
    Supports both single API_KEY and multiple API_KEYS.

    Returns:
        Set of valid API keys
    """
    keys = set()

    # Add single API key if set
    if API_KEY:
        keys.add(API_KEY)

    # Add multiple API keys if set
    if API_KEYS_STR:
        for key in API_KEYS_STR.split(","):
            key = key.strip()
            if key:
                keys.add(key)

    return keys


# Load API keys
VALID_API_KEYS = load_api_keys()

# Warn if authentication is disabled in production
if IS_PRODUCTION and not AUTH_ENABLED:
    logger.warning(
        "SECURITY WARNING: Authentication is disabled in production environment! "
        "This is a critical security risk. Set AUTH_ENABLED=true and provide API_KEY."
    )

# Warn if API key is not set when auth is enabled
if AUTH_ENABLED and not VALID_API_KEYS:
    logger.error(
        "SECURITY ERROR: AUTH_ENABLED is true but no API keys are set! "
        "Authentication will fail for all requests. Set API_KEY or API_KEYS environment variable."
    )

# In production, force authentication for critical endpoints
if IS_PRODUCTION and FORCE_AUTH_CRITICAL and not AUTH_ENABLED:
    logger.error(
        "SECURITY ERROR: Production environment requires authentication for critical endpoints, "
        "but AUTH_ENABLED is false. Critical endpoints will be blocked."
    )


def is_critical_endpoint(path: str) -> bool:
    """Check if endpoint is critical and requires authentication."""
    return any(path.startswith(endpoint) for endpoint in CRITICAL_ENDPOINTS)


def track_auth_failure(ip: str, path: str):
    """Track authentication failure for monitoring."""
    now = time.time()
    _auth_failures.append((now, ip, path))
    _failure_counts[ip] += 1

    # Log warning if too many failures from same IP
    if _failure_counts[ip] >= 10:
        logger.warning(
            f"Multiple authentication failures from {ip} (count: {_failure_counts[ip]})",
            extra={"ip": ip, "failure_count": _failure_counts[ip], "path": path},
        )


def get_recent_failures(window_seconds: int = 300) -> List[tuple]:
    """Get recent authentication failures within time window."""
    now = time.time()
    cutoff = now - window_seconds
    return [(t, ip, path) for t, ip, path in _auth_failures if t >= cutoff]


class AuthMiddleware(BaseHTTPMiddleware):
    """
    API key authentication middleware with support for multiple keys.

    Features:
    - Supports single API_KEY or multiple API_KEYS (comma-separated)
    - Constant-time key comparison to prevent timing attacks
    - Authentication failure tracking and monitoring
    - Force authentication for critical endpoints in production
    """

    def _should_skip_auth(self, path: str, method: str) -> bool:
        """Check if authentication should be skipped for this request."""
        return path in PUBLIC_ENDPOINTS or method == "OPTIONS"

    def _check_critical_endpoint_auth(self, path: str) -> Optional[JSONResponse]:
        """Check if critical endpoint requires authentication in production."""
        if IS_PRODUCTION and FORCE_AUTH_CRITICAL and is_critical_endpoint(path):
            if not AUTH_ENABLED or not VALID_API_KEYS:
                logger.error(
                    f"SECURITY ERROR: Critical endpoint {path} requires authentication "
                    f"but authentication is disabled or no keys configured. Blocking request."
                )
                return JSONResponse(
                    status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
                    content={
                        "error": "Service Unavailable",
                        "message": "Authentication required but not configured",
                        "status_code": 503,
                    },
                )
        return None

    def _extract_api_key(self, request: Request) -> Optional[str]:
        """Extract API key from request headers."""
        api_key = request.headers.get("X-API-Key") or request.headers.get("Authorization")
        if api_key and api_key.startswith("Bearer "):
            api_key = api_key[7:]
        return api_key

    def _validate_api_key(self, api_key: str) -> bool:
        """Validate API key against valid keys using constant-time comparison."""
        for valid_key in VALID_API_KEYS:
            if secrets.compare_digest(api_key, valid_key):
                return True
        return False

    def _create_unauthorized_response(self, message: str) -> JSONResponse:
        """Create unauthorized response."""
        return JSONResponse(
            status_code=status.HTTP_401_UNAUTHORIZED,
            content={
                "error": "Unauthorized",
                "message": message,
                "status_code": 401,
            },
            headers={"WWW-Authenticate": "Bearer"},
        )

    async def dispatch(self, request: Request, call_next):
        """Check API key for protected endpoints."""
        client_ip = request.client.host if request.client else "unknown"
        path = request.url.path

        # Skip authentication for public endpoints
        if self._should_skip_auth(path, request.method):
            return await call_next(request)

        # Check critical endpoint authentication in production
        critical_response = self._check_critical_endpoint_auth(path)
        if critical_response:
            return critical_response

        # Skip authentication if disabled (but still check critical endpoints above)
        if not AUTH_ENABLED:
            return await call_next(request)

        # Extract and validate API key
        api_key = self._extract_api_key(request)
        if not api_key:
            track_auth_failure(client_ip, path)
            logger.warning(
                f"Authentication failed: missing API key for {path} from {client_ip}",
                extra={"path": path, "ip": client_ip},
            )
            return self._create_unauthorized_response("Missing API key")

        if not self._validate_api_key(api_key):
            track_auth_failure(client_ip, path)
            logger.warning(
                f"Authentication failed: invalid API key for {path} from {client_ip}",
                extra={"path": path, "ip": client_ip},
            )
            return self._create_unauthorized_response("Invalid API key")

        return await call_next(request)
