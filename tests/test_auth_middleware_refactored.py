"""
Tests for refactored AuthMiddleware.
"""

import pytest
from unittest.mock import Mock, patch
from starlette.requests import Request
from starlette.responses import JSONResponse

from compass.service.middleware.auth import AuthMiddleware


@pytest.fixture
def auth_middleware():
    """Create AuthMiddleware instance."""
    return AuthMiddleware(Mock())


@pytest.fixture
def mock_request():
    """Create a mock request."""
    request = Mock(spec=Request)
    request.url.path = "/api/v1/training/tasks"
    request.method = "GET"
    request.client.host = "127.0.0.1"
    request.headers = {}
    return request


class TestAuthMiddlewareHelpers:
    """Test AuthMiddleware helper methods."""

    def test_should_skip_auth_public_endpoint(self, auth_middleware):
        """Test skipping auth for public endpoints."""
        result = auth_middleware._should_skip_auth("/health", "GET")
        assert result is True

    def test_should_skip_auth_options(self, auth_middleware):
        """Test skipping auth for OPTIONS requests."""
        result = auth_middleware._should_skip_auth("/api/v1/training/tasks", "OPTIONS")
        assert result is True

    def test_should_not_skip_auth_protected(self, auth_middleware):
        """Test not skipping auth for protected endpoints."""
        result = auth_middleware._should_skip_auth("/api/v1/training/tasks", "GET")
        assert result is False

    def test_extract_api_key_from_header(self, auth_middleware):
        """Test extracting API key from header."""
        request = Mock()
        request.headers = {"X-API-Key": "test-key"}
        key = auth_middleware._extract_api_key(request)
        assert key == "test-key"

    def test_extract_api_key_from_bearer(self, auth_middleware):
        """Test extracting API key from Bearer token."""
        request = Mock()
        request.headers = {"Authorization": "Bearer test-key"}
        key = auth_middleware._extract_api_key(request)
        assert key == "test-key"

    def test_extract_api_key_missing(self, auth_middleware):
        """Test extracting API key when missing."""
        request = Mock()
        request.headers = {}
        key = auth_middleware._extract_api_key(request)
        assert key is None

    @patch("compass.service.middleware.auth.VALID_API_KEYS", {"valid-key"})
    def test_validate_api_key_valid(self, auth_middleware):
        """Test validating valid API key."""
        result = auth_middleware._validate_api_key("valid-key")
        assert result is True

    @patch("compass.service.middleware.auth.VALID_API_KEYS", {"valid-key"})
    def test_validate_api_key_invalid(self, auth_middleware):
        """Test validating invalid API key."""
        result = auth_middleware._validate_api_key("invalid-key")
        assert result is False

    def test_create_unauthorized_response(self, auth_middleware):
        """Test creating unauthorized response."""
        response = auth_middleware._create_unauthorized_response("Missing API key")
        assert isinstance(response, JSONResponse)
        assert response.status_code == 401


