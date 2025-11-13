"""
Security headers middleware for COMPASS service.
"""

from starlette.middleware.base import BaseHTTPMiddleware
from starlette.requests import Request
from starlette.responses import Response
import logging

logger = logging.getLogger(__name__)


class SecurityHeadersMiddleware(BaseHTTPMiddleware):
    """Middleware to add security headers to responses."""

    async def dispatch(self, request: Request, call_next):
        """Add security headers to response."""
        response: Response = await call_next(request)

        # Content-Security-Policy - strict policy to prevent XSS
        # For API endpoints, we use a restrictive policy
        # Frontend endpoints may need to adjust this based on their requirements
        # Check if this is an API endpoint
        is_api_endpoint = request.url.path.startswith("/api/")

        if is_api_endpoint:
            # Stricter CSP for API endpoints (JSON responses, no scripts)
            csp_policy = (
                "default-src 'none'; "
                "frame-ancestors 'none'; "
                "base-uri 'none'; "
                "form-action 'none';"
            )
        else:
            # More permissive CSP for frontend/documentation endpoints
            csp_policy = (
                "default-src 'self'; "
                "script-src 'self'; "  # Removed unsafe-inline and unsafe-eval for better security
                "style-src 'self' 'unsafe-inline'; "  # Keep unsafe-inline for inline styles (common in web apps)
                "img-src 'self' data: https:; "
                "font-src 'self' data:; "
                "connect-src 'self' ws: wss:; "
                "frame-ancestors 'none'; "
                "base-uri 'self'; "
                "form-action 'self';"
            )
        response.headers["Content-Security-Policy"] = csp_policy

        # X-Frame-Options
        response.headers["X-Frame-Options"] = "DENY"

        # X-Content-Type-Options
        response.headers["X-Content-Type-Options"] = "nosniff"

        # X-XSS-Protection
        response.headers["X-XSS-Protection"] = "1; mode=block"

        # Referrer-Policy
        response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"

        # Permissions-Policy
        response.headers["Permissions-Policy"] = (
            "geolocation=(), " "microphone=(), " "camera=(), " "payment=(), " "usb=()"
        )

        # Strict-Transport-Security (only if HTTPS)
        # Uncomment in production when HTTPS is enabled
        # response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"

        return response
