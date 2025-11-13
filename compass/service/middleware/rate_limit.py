"""
Rate limiting middleware for COMPASS service.
"""

import time
import logging
from typing import Dict, Tuple, Optional
from collections import defaultdict, deque
from fastapi import Request, HTTPException, status
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from compass.service.error_codes import ErrorCode

logger = logging.getLogger(__name__)

# Rate limit statistics for monitoring
_rate_limit_stats = {
    "total_requests": 0,
    "rate_limited_requests": 0,
    "by_endpoint": defaultdict(int),
    "by_ip": defaultdict(int),
}


class RateLimiter:
    """Simple in-memory rate limiter using sliding window."""

    def __init__(self, max_requests: int = 100, window_seconds: int = 60):
        """
        Initialize rate limiter.

        Args:
            max_requests: Maximum number of requests per window
            window_seconds: Time window in seconds
        """
        self.max_requests = max_requests
        self.window_seconds = window_seconds
        self.requests: Dict[str, deque] = defaultdict(lambda: deque())

    def is_allowed(self, identifier: str) -> Tuple[bool, Optional[int]]:
        """
        Check if request is allowed.

        Args:
            identifier: Request identifier (IP address, user ID, etc.)

        Returns:
            Tuple of (is_allowed, remaining_requests)
        """
        now = time.time()
        request_times = self.requests[identifier]

        # Remove old requests outside the window
        while request_times and request_times[0] < now - self.window_seconds:
            request_times.popleft()

        # Check if limit exceeded
        if len(request_times) >= self.max_requests:
            return False, 0

        # Add current request
        request_times.append(now)

        # Calculate remaining requests
        remaining = max(0, self.max_requests - len(request_times))
        return True, remaining

    def get_remaining(self, identifier: str) -> int:
        """Get remaining requests for identifier."""
        now = time.time()
        request_times = self.requests[identifier]

        # Remove old requests
        while request_times and request_times[0] < now - self.window_seconds:
            request_times.popleft()

        return max(0, self.max_requests - len(request_times))

    def reset(self, identifier: Optional[str] = None):
        """Reset rate limit for identifier or all."""
        if identifier:
            self.requests.pop(identifier, None)
        else:
            self.requests.clear()


class RateLimitMiddleware(BaseHTTPMiddleware):
    """Rate limiting middleware for FastAPI."""

    def __init__(
        self,
        app,
        default_limit: int = 100,
        default_window: int = 60,
        per_endpoint_limits: Optional[Dict[str, Dict[str, int]]] = None,
        identifier_func: Optional[callable] = None,
    ):
        global _middleware_instance
        _middleware_instance = self
        """
        Initialize rate limit middleware.

        Args:
            app: FastAPI application
            default_limit: Default requests per window
            default_window: Default window in seconds
            per_endpoint_limits: Per-endpoint limits in format:
                {
                    "/api/v1/training/tasks": {"limit": 5, "window": 60},
                    "/api/v1/data/upload": {"limit": 10, "window": 60}
                }
            identifier_func: Function to extract identifier from request
                Default: uses client IP address
        """
        super().__init__(app)
        self.default_limiter = RateLimiter(default_limit, default_window)
        self.endpoint_limiters: Dict[str, RateLimiter] = {}

        # Setup per-endpoint limiters
        if per_endpoint_limits:
            for endpoint, config in per_endpoint_limits.items():
                limit = config.get("limit", default_limit)
                window = config.get("window", default_window)
                self.endpoint_limiters[endpoint] = RateLimiter(limit, window)

        self.identifier_func = identifier_func or self._get_client_ip

    def _get_client_ip(self, request: Request) -> str:
        """Extract client IP address from request."""
        # Check X-Forwarded-For header (for proxies)
        forwarded = request.headers.get("X-Forwarded-For")
        if forwarded:
            # Take the first IP in the chain
            return forwarded.split(",")[0].strip()

        # Check X-Real-IP header
        real_ip = request.headers.get("X-Real-IP")
        if real_ip:
            return real_ip.strip()

        # Fall back to client host
        if request.client:
            return request.client.host

        return "unknown"

    async def dispatch(self, request: Request, call_next):
        """Process request with rate limiting."""
        # Skip rate limiting for health checks and public endpoints
        skip_paths = ["/health", "/health/ready", "/docs", "/openapi.json", "/redoc", "/"]
        if request.url.path in skip_paths:
            return await call_next(request)

        # Get identifier
        identifier = self.identifier_func(request)

        # Get appropriate limiter
        limiter = self.default_limiter
        for endpoint, endpoint_limiter in self.endpoint_limiters.items():
            if request.url.path.startswith(endpoint):
                limiter = endpoint_limiter
                break

        # Update statistics
        _rate_limit_stats["total_requests"] += 1
        _rate_limit_stats["by_endpoint"][request.url.path] += 1
        _rate_limit_stats["by_ip"][identifier] += 1

        # Check rate limit
        is_allowed, remaining = limiter.is_allowed(identifier)

        if not is_allowed:
            # Update rate limit statistics
            _rate_limit_stats["rate_limited_requests"] += 1

            logger.warning(
                f"Rate limit exceeded for {identifier} on {request.url.path} "
                f"(limit: {limiter.max_requests}/{limiter.window_seconds}s)",
                extra={
                    "ip": identifier,
                    "path": request.url.path,
                    "method": request.method,
                    "limit": limiter.max_requests,
                    "window": limiter.window_seconds,
                    "user_agent": request.headers.get("User-Agent", "unknown"),
                },
            )
            reset_time = int(time.time()) + limiter.window_seconds
            return JSONResponse(
                status_code=status.HTTP_429_TOO_MANY_REQUESTS,
                content={
                    "error": "Rate limit exceeded",
                    "error_code": ErrorCode.RATE_LIMIT_EXCEEDED.value,
                    "status_code": 429,
                    "detail": f"Too many requests. Limit: {limiter.max_requests} per {limiter.window_seconds}s. Please try again later.",
                },
                headers={
                    "X-RateLimit-Limit": str(limiter.max_requests),
                    "X-RateLimit-Remaining": "0",
                    "X-RateLimit-Reset": str(reset_time),
                    "Retry-After": str(limiter.window_seconds),
                },
            )

        # Add rate limit headers
        response = await call_next(request)
        reset_time = int(time.time()) + limiter.window_seconds
        response.headers["X-RateLimit-Limit"] = str(limiter.max_requests)
        response.headers["X-RateLimit-Remaining"] = str(remaining)
        response.headers["X-RateLimit-Reset"] = str(reset_time)

        return response


def get_rate_limit_stats() -> Dict:
    """Get rate limit statistics for monitoring."""
    return {
        "total_requests": _rate_limit_stats["total_requests"],
        "rate_limited_requests": _rate_limit_stats["rate_limited_requests"],
        "rate_limit_percentage": (
            (_rate_limit_stats["rate_limited_requests"] / _rate_limit_stats["total_requests"] * 100)
            if _rate_limit_stats["total_requests"] > 0
            else 0.0
        ),
        "top_endpoints": dict(
            sorted(_rate_limit_stats["by_endpoint"].items(), key=lambda x: x[1], reverse=True)[:10]
        ),
        "top_ips": dict(
            sorted(_rate_limit_stats["by_ip"].items(), key=lambda x: x[1], reverse=True)[:10]
        ),
    }


# Global reference to middleware instance for testing
_middleware_instance = None


def reset_rate_limit_stats():
    """Reset rate limit statistics."""
    global _middleware_instance  # noqa: F824
    _rate_limit_stats["total_requests"] = 0
    _rate_limit_stats["rate_limited_requests"] = 0
    _rate_limit_stats["by_endpoint"].clear()
    _rate_limit_stats["by_ip"].clear()
    # Also reset rate limiters if middleware instance is available
    if _middleware_instance:
        _middleware_instance.default_limiter.reset()
        for limiter in _middleware_instance.endpoint_limiters.values():
            limiter.reset()
