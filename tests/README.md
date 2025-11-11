# COMPASS Service Tests

This directory contains unit and integration tests for the COMPASS service.

## Test Structure

```
tests/
├── __init__.py
├── conftest.py              # Pytest configuration and fixtures
├── test_service_exceptions.py
├── test_error_codes.py
├── test_config_manager.py
├── test_upload_queue.py
├── test_rate_limiter.py
├── test_progress_tracker.py
├── test_models_validation.py
├── test_backup_manager.py
└── test_benchmark.py
```

## Running Tests

### Run all tests
```bash
pytest
```

### Run with coverage
```bash
pytest --cov=compass --cov-report=html
```

### Run specific test file
```bash
pytest tests/test_service_exceptions.py
```

### Run specific test
```bash
pytest tests/test_service_exceptions.py::test_service_exception
```

### Run with markers
```bash
pytest -m unit
pytest -m integration
pytest -m slow
```

## Test Coverage

Current test coverage includes:
- Service exceptions and error handling
- Error codes and status mapping
- Configuration management
- Upload queue management
- Rate limiting
- Progress tracking
- Model validation
- Backup management
- Benchmark utilities

## Adding New Tests

When adding new tests:
1. Follow the naming convention: `test_*.py` for test files
2. Use descriptive test function names: `test_<feature>_<scenario>`
3. Add appropriate fixtures in `conftest.py` if needed
4. Use markers for test categorization (`@pytest.mark.unit`, `@pytest.mark.integration`, etc.)
5. Ensure tests are independent and can run in any order












