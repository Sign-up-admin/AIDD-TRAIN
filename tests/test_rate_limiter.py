"""
Tests for rate limiter.
"""
import pytest
import time
from compass.service.middleware.rate_limit import RateLimiter


def test_rate_limiter_init():
    """Test RateLimiter initialization."""
    limiter = RateLimiter(max_requests=10, window_seconds=60)
    assert limiter.max_requests == 10
    assert limiter.window_seconds == 60


def test_rate_limiter_allows_requests():
    """Test rate limiter allows requests within limit."""
    limiter = RateLimiter(max_requests=5, window_seconds=60)
    
    for i in range(5):
        allowed, remaining = limiter.is_allowed("test-ip")
        assert allowed is True
        assert remaining == 5 - i - 1


def test_rate_limiter_blocks_excess_requests():
    """Test rate limiter blocks requests exceeding limit."""
    limiter = RateLimiter(max_requests=3, window_seconds=60)
    
    # Use up all requests
    for _ in range(3):
        limiter.is_allowed("test-ip")
    
    # Next request should be blocked
    allowed, remaining = limiter.is_allowed("test-ip")
    assert allowed is False
    assert remaining == 0


def test_rate_limiter_get_remaining():
    """Test getting remaining requests."""
    limiter = RateLimiter(max_requests=5, window_seconds=60)
    
    limiter.is_allowed("test-ip")
    assert limiter.get_remaining("test-ip") == 4
    
    limiter.is_allowed("test-ip")
    assert limiter.get_remaining("test-ip") == 3


def test_rate_limiter_reset():
    """Test resetting rate limiter."""
    limiter = RateLimiter(max_requests=3, window_seconds=60)
    
    # Use up all requests
    for _ in range(3):
        limiter.is_allowed("test-ip")
    
    # Reset
    limiter.reset("test-ip")
    
    # Should be able to make requests again
    allowed, remaining = limiter.is_allowed("test-ip")
    assert allowed is True
    assert remaining == 2


