@echo off
REM Run tests with coverage
REM Usage: run_tests.bat [additional pytest arguments]

cd /d "%~dp0\.."
set PYTHONPATH=%CD%;%PYTHONPATH%
python -m pytest tests\ --cov=compass --cov-report=html --cov-report=xml --cov-report=term %*

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Tests completed with failures. See htmlcov\ for coverage report.
    exit /b %ERRORLEVEL%
) else (
    echo.
    echo All tests passed!
    exit /b 0
)





