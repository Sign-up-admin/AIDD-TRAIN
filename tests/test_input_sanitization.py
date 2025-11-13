"""
Tests for input sanitization.
"""

import pytest
from compass.service.utils.input_sanitizer import (
    sanitize_string,
    sanitize_description,
    sanitize_task_id,
    sanitize_filename,
)


def test_sanitize_string_basic():
    """Test basic string sanitization."""
    result = sanitize_string("test string")
    assert result == "test string"


def test_sanitize_string_xss_prevention():
    """Test XSS prevention in string sanitization."""
    malicious = "<script>alert('XSS')</script>"
    result = sanitize_string(malicious)
    assert "<script>" not in result
    assert "&lt;script&gt;" in result


def test_sanitize_string_max_length():
    """Test string length limiting."""
    long_string = "a" * 200
    result = sanitize_string(long_string, max_length=100)
    assert len(result) == 100


def test_sanitize_description_none():
    """Test sanitizing None description."""
    result = sanitize_description(None)
    assert result is None


def test_sanitize_description_xss():
    """Test XSS prevention in description."""
    malicious = "<img src=x onerror=alert('XSS')>"
    result = sanitize_description(malicious)
    assert "<img" not in result
    assert "&lt;img" in result


def test_sanitize_task_id_valid():
    """Test sanitizing valid task ID."""
    task_id = "550e8400-e29b-41d4-a716-446655440000"
    result = sanitize_task_id(task_id)
    assert result == task_id.lower()


def test_sanitize_task_id_invalid():
    """Test sanitizing invalid task ID."""
    with pytest.raises(ValueError):
        sanitize_task_id("invalid-task-id")


def test_sanitize_filename_path_traversal():
    """Test path traversal prevention in filename."""
    malicious = "../../etc/passwd"
    result = sanitize_filename(malicious)
    assert ".." not in result
    assert "/" not in result
    assert "\\" not in result


def test_sanitize_filename_control_chars():
    """Test control character removal in filename."""
    filename = "test\x00file.txt"
    result = sanitize_filename(filename)
    assert "\x00" not in result

