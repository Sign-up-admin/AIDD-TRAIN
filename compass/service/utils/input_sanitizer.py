"""
Input sanitization utilities for XSS prevention and security.
"""

import re
import html
from typing import Optional


def sanitize_string(value: str, max_length: Optional[int] = None) -> str:
    """
    Sanitize a string input to prevent XSS and other security issues.

    Args:
        value: Input string to sanitize
        max_length: Maximum allowed length (None for no limit)

    Returns:
        Sanitized string
    """
    if not isinstance(value, str):
        raise ValueError("Input must be a string")

    # Remove null bytes and control characters (except newline and tab)
    value = re.sub(r"[\x00-\x08\x0b-\x0c\x0e-\x1f]", "", value)

    # Limit length if specified
    if max_length and len(value) > max_length:
        value = value[:max_length]

    # HTML escape to prevent XSS (for display purposes)
    # Note: This is defensive - JSON responses don't execute HTML, but this adds extra protection
    value = html.escape(value, quote=True)

    return value


def sanitize_description(description: Optional[str]) -> Optional[str]:
    """
    Sanitize a description field.

    Args:
        description: Description string to sanitize

    Returns:
        Sanitized description or None
    """
    if description is None:
        return None

    if not isinstance(description, str):
        raise ValueError("Description must be a string")

    # Remove null bytes and control characters
    description = re.sub(r"[\x00-\x08\x0b-\x0c\x0e-\x1f]", "", description)

    # Remove javascript: protocol and other dangerous protocols
    dangerous_protocols = ["javascript:", "data:", "vbscript:", "file:", "about:"]
    for protocol in dangerous_protocols:
        # Case-insensitive replacement
        description = re.sub(re.escape(protocol), "", description, flags=re.IGNORECASE)

    # Limit length (10000 characters)
    if len(description) > 10000:
        description = description[:10000]

    # HTML escape to prevent XSS
    description = html.escape(description, quote=True)

    return description


def sanitize_task_id(task_id: str) -> str:
    """
    Sanitize and validate a task ID (UUID format).

    Args:
        task_id: Task ID to sanitize

    Returns:
        Sanitized task ID

    Raises:
        ValueError: If task_id is not a valid UUID format
    """
    if not isinstance(task_id, str):
        raise ValueError("Task ID must be a string")

    # Remove any whitespace
    task_id = task_id.strip()

    # Validate UUID format (basic check)
    uuid_pattern = re.compile(
        r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$", re.IGNORECASE
    )
    if not uuid_pattern.match(task_id):
        raise ValueError("Invalid task ID format")

    return task_id.lower()


def sanitize_filename(filename: str) -> str:
    """
    Sanitize a filename to prevent path traversal and other issues.

    Args:
        filename: Filename to sanitize

    Returns:
        Sanitized filename
    """
    if not isinstance(filename, str):
        raise ValueError("Filename must be a string")

    # Remove path components
    filename = filename.replace("..", "").replace("/", "").replace("\\", "")

    # Remove null bytes and control characters
    filename = re.sub(r"[\x00-\x1f]", "", filename)

    # Limit length
    if len(filename) > 255:
        filename = filename[:255]

    return filename
