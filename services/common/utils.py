"""
Common utility functions for services.
"""
import socket
import uuid
from typing import Optional


def get_local_ip() -> str:
    """Get local IP address."""
    try:
        # Connect to a remote server to determine local IP
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"


def generate_service_id(service_name: str) -> str:
    """Generate a unique service ID."""
    return f"{service_name}-{uuid.uuid4().hex[:8]}"


def validate_url(url: str) -> bool:
    """Validate URL format."""
    try:
        from urllib.parse import urlparse
        result = urlparse(url)
        return all([result.scheme, result.netloc])
    except Exception:
        return False


