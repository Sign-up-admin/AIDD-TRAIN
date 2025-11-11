# Code Quality Guide

This document describes the code quality standards and tools used in the AIDD-TRAIN project.

## Overview

The project uses a comprehensive set of code quality tools to ensure consistency, maintainability, and security:

- **Black**: Code formatting
- **Flake8**: Code style and complexity checking
- **Pylint**: Code quality analysis
- **MyPy**: Static type checking
- **Bandit**: Security vulnerability scanning
- **Pytest**: Unit testing and coverage

## Quick Start

### Installation

Install all development dependencies:

```bash
pip install -r requirements-dev.txt
```

### Running Checks

#### Run All Checks

```bash
# Using Python script
python scripts/run_all_checks.py

# Using Windows batch script
scripts\run_lint.bat

# CI mode (non-interactive)
python scripts/run_all_checks.py --ci
```

#### Run Specific Tools

```bash
# Run only Black
python scripts/run_all_checks.py --tool black

# Run only Flake8
python scripts/run_all_checks.py --tool flake8

# Run only tests (skip other checks)
python scripts/run_all_checks.py --no-tests
python -m pytest
```

#### Auto-fix Formatting

```bash
# Auto-fix formatting issues
python scripts/run_all_checks.py --format

# Or use the batch script
scripts\format_code.bat

# Or run Black directly
python -m black compass services FLASH_DOCK-main/services FLASH_DOCK-main/pages
```

#### Run Tests

```bash
# Run all tests with coverage
python scripts/run_tests.bat

# Or directly
python -m pytest tests/ --cov=compass --cov-report=html --cov-report=term
```

## Tool Configuration

### Black (Code Formatting)

**Configuration:** `pyproject.toml`

- Line length: 100 characters
- Target Python versions: 3.8-3.11
- Auto-formats code to PEP 8 style

**Usage:**
```bash
# Check formatting
black --check .

# Format code
black .
```

**Excluded directories:**
- `.eggs`, `.git`, `.hg`
- `.mypy_cache`, `.tox`, `.venv`, `venv`
- `_build`, `buck-out`, `build`, `dist`
- `checkpoints`, `logs`, `data`, `processed_data`
- `FLASH_DOCK-main/others`

### Flake8 (Code Style)

**Configuration:** `.flake8`

- Max line length: 100
- Max complexity: 10
- Ignored errors: E203, E501, W503, E402, F401 (for __init__.py)

**Common Issues:**
- `W293`: Blank line contains whitespace (auto-fixable by Black)
- `E241`: Multiple spaces after punctuation (auto-fixable by Black)
- `E261`: At least two spaces before inline comment (auto-fixable by Black)
- `C901`: Function too complex (requires manual refactoring)

**Usage:**
```bash
flake8 compass services FLASH_DOCK-main/services FLASH_DOCK-main/pages
```

### Pylint (Code Quality)

**Configuration:** `.pylintrc`

- Max line length: 100
- Disabled checks: missing docstrings, too-few-public-methods, etc.
- Jobs: 4 (parallel processing)

**Usage:**
```bash
pylint compass services FLASH_DOCK-main/services FLASH_DOCK-main/pages
```

### MyPy (Type Checking)

**Configuration:** `pyproject.toml`

- Python version: 3.10
- Strict type checking disabled (gradual typing)
- Ignores missing imports for third-party libraries (rdkit, torch, etc.)

**Usage:**
```bash
mypy compass services
```

**Common Issues:**
- Missing type annotations
- Type mismatches
- Missing type stubs (install with `pip install types-requests`)

### Bandit (Security)

**Configuration:** Default settings

- Scans for common security issues
- Reports high and medium severity findings

**Usage:**
```bash
bandit -r compass services FLASH_DOCK-main
```

### Pytest (Testing)

**Configuration:** `pytest.ini` and `pyproject.toml`

- Test paths: `tests/`
- Coverage: compass package
- Reports: HTML, XML, terminal

**Usage:**
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=compass --cov-report=html

# Run specific test
pytest tests/test_config_manager.py

# Run with markers
pytest -m unit
pytest -m integration
```

## Code Quality Standards

### Success Criteria

- ✅ All code passes Black formatting check
- ✅ Flake8 errors < 50 (target: < 20)
- ✅ Pylint score > 7.0/10
- ✅ MyPy errors < 20 (target: < 5)
- ✅ Bandit: No high-severity security issues
- ✅ Test coverage > 60%
- ✅ All tests passing

### Priority Levels

- **P0**: Security issues, critical bugs, failing tests
- **P1**: Code quality issues, type errors, high complexity
- **P2**: Style issues, formatting problems
- **P3**: Optimization suggestions, documentation improvements

## Workflow

### Before Committing

1. Run auto-formatter:
   ```bash
   python scripts/run_all_checks.py --format
   ```

2. Run all checks:
   ```bash
   python scripts/run_all_checks.py
   ```

3. Fix any P0 and P1 issues

4. Run tests:
   ```bash
   python -m pytest
   ```

### CI/CD Integration

The scripts support CI mode for automated checks:

```bash
# Non-interactive mode, exits with error code on failure
python scripts/run_all_checks.py --ci
```

Exit codes:
- `0`: All checks passed
- `1-5`: Number of failed tools

## Reports

All reports are saved to `lint_reports/` directory:

- `{timestamp}_black_report.txt`
- `{timestamp}_flake8_report.txt`
- `{timestamp}_pylint_report.txt`
- `{timestamp}_mypy_report.txt`
- `{timestamp}_bandit_report.txt`
- `{timestamp}_pytest_report.txt`
- `{timestamp}_summary.md`

Generate unified HTML report:
```bash
python scripts/generate_reports.py
```

## Common Issues and Fixes

### Formatting Issues

**Problem:** Black reports files need reformatting

**Solution:**
```bash
python scripts/run_all_checks.py --format
# or
python -m black .
```

### Style Issues (Flake8)

**Problem:** W293 (blank line contains whitespace), E241 (multiple spaces)

**Solution:** Run Black to auto-fix most issues:
```bash
python -m black .
```

### Complexity Issues (C901)

**Problem:** Function too complex (complexity > 10)

**Solution:** Refactor into smaller functions:
```python
# Before
def complex_function():
    # 50 lines of complex logic
    pass

# After
def complex_function():
    step1()
    step2()
    step3()

def step1():
    # Specific logic
    pass
```

### Type Errors (MyPy)

**Problem:** Missing type annotations or type mismatches

**Solution:**
```python
# Before
def get_user(id):
    return database.get(id)

# After
def get_user(id: int) -> User | None:
    return database.get(id)
```

### Missing Type Stubs

**Problem:** MyPy reports "Library stubs not installed"

**Solution:**
```bash
pip install types-requests types-pyyaml
# or
mypy --install-types
```

## Best Practices

1. **Run checks frequently**: Don't let issues accumulate
2. **Fix P0 issues immediately**: Security and critical bugs
3. **Auto-fix when possible**: Use Black for formatting
4. **Incremental improvements**: Fix a few issues at a time
5. **Test before committing**: Ensure tests pass
6. **Review reports**: Understand what needs fixing

## Troubleshooting

### Tools not found

If a tool is not found, install it:
```bash
pip install -r requirements-dev.txt
```

### Encoding errors on Windows

The scripts handle Windows console encoding automatically. If you encounter issues:
```bash
# Set console encoding
chcp 65001
```

### Tests failing due to missing dependencies

Install project dependencies:
```bash
pip install -r requirements.txt
pip install -r requirements-service.txt
```

## Resources

- [Black Documentation](https://black.readthedocs.io/)
- [Flake8 Documentation](https://flake8.pycqa.org/)
- [Pylint Documentation](https://pylint.pycqa.org/)
- [MyPy Documentation](https://mypy.readthedocs.io/)
- [Bandit Documentation](https://bandit.readthedocs.io/)
- [Pytest Documentation](https://docs.pytest.org/)

## Getting Help

If you encounter issues:

1. Check the tool's documentation
2. Review the generated reports in `lint_reports/`
3. Check the `issues_analysis.md` for known issues
4. Ask the team for guidance











