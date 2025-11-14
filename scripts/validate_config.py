#!/usr/bin/env python3
"""
Configuration validation script.

Validates that all required configuration values are set and have valid formats.
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


def validate_env_var(
    name: str, required: bool = False, validator: callable = None
) -> Tuple[bool, str]:
    """
    Validate an environment variable.

    Args:
        name: Environment variable name
        required: Whether the variable is required
        validator: Optional validation function

    Returns:
        Tuple of (is_valid, error_message)
    """
    value = os.getenv(name)

    if required and not value:
        return False, f"Required environment variable {name} is not set"

    if value and validator:
        try:
            validator(value)
        except Exception as e:
            return False, f"Invalid value for {name}: {str(e)}"

    return True, ""


def validate_port(value: str) -> None:
    """Validate port number."""
    port = int(value)
    if port < 1 or port > 65535:
        raise ValueError(f"Port must be between 1 and 65535, got {port}")


def validate_boolean(value: str) -> None:
    """Validate boolean value."""
    if value.lower() not in ["true", "false"]:
        raise ValueError(f"Boolean must be 'true' or 'false', got {value}")


def validate_positive_int(value: str) -> None:
    """Validate positive integer."""
    num = int(value)
    if num <= 0:
        raise ValueError(f"Must be a positive integer, got {num}")


def validate_positive_float(value: str) -> None:
    """Validate positive float."""
    num = float(value)
    if num <= 0:
        raise ValueError(f"Must be a positive number, got {num}")


def validate_cors_origins(value: str) -> None:
    """Validate CORS origins."""
    origins = [o.strip() for o in value.split(",")]
    for origin in origins:
        if origin == "*":
            raise ValueError("Wildcard origin '*' is not allowed for security")
        if not (origin.startswith("http://") or origin.startswith("https://")):
            raise ValueError(f"Origin must start with http:// or https://, got {origin}")


def validate_config() -> Tuple[bool, List[str]]:
    """
    Validate all configuration values.

    Returns:
        Tuple of (is_valid, list_of_errors)
    """
    errors = []

    # Required configuration
    required_vars = [
        ("COMPASS_HOST", True),
        ("COMPASS_PORT", True, validate_port),
        ("ENVIRONMENT", False),
    ]

    # Optional configuration with validators
    optional_vars = [
        ("AUTH_ENABLED", False, validate_boolean),
        ("FORCE_AUTH_CRITICAL", False, validate_boolean),
        ("CORS_ORIGINS", False, validate_cors_origins),
        ("COMPASS_MAX_WORKERS", False, validate_positive_int),
        ("MAX_CONCURRENT_UPLOADS", False, validate_positive_int),
        ("MAX_CONCURRENT_TASKS", False, validate_positive_int),
        ("DB_CONNECTION_TIMEOUT", False, validate_positive_float),
        ("DB_BUSY_TIMEOUT", False, validate_positive_int),
        ("RATE_LIMIT_DEFAULT", False, validate_positive_int),
    ]

    # Check required variables
    for var_info in required_vars:
        name = var_info[0]
        required = var_info[1] if len(var_info) > 1 else False
        validator = var_info[2] if len(var_info) > 2 else None

        is_valid, error = validate_env_var(name, required, validator)
        if not is_valid:
            errors.append(error)

    # Check optional variables
    for var_info in optional_vars:
        name = var_info[0]
        validator = var_info[1] if len(var_info) > 1 else None

        is_valid, error = validate_env_var(name, False, validator)
        if not is_valid:
            errors.append(error)

    # Production-specific checks
    if os.getenv("ENVIRONMENT", "").lower() == "production":
        if os.getenv("AUTH_ENABLED", "").lower() != "true":
            errors.append("AUTH_ENABLED must be 'true' in production environment")

        if not os.getenv("API_KEY") and not os.getenv("API_KEYS"):
            errors.append("API_KEY or API_KEYS must be set in production environment")

        if os.getenv("CORS_ORIGINS", "").strip() == "":
            errors.append("CORS_ORIGINS must be set in production environment")

    return len(errors) == 0, errors


def main():
    """Main entry point."""
    is_valid, errors = validate_config()

    if is_valid:
        print("✓ Configuration is valid")
        sys.exit(0)
    else:
        print("✗ Configuration validation failed:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)


if __name__ == "__main__":
    main()



