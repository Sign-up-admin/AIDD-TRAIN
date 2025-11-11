@echo off
REM Run all code quality checks
REM Usage: run_lint.bat [--format] [--tool TOOL_NAME] [--no-tests] [--ci]

cd /d "%~dp0\.."
set PYTHONPATH=%CD%;%PYTHONPATH%
python scripts\run_all_checks.py %*

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Checks completed with errors. See lint_reports\ for details.
    exit /b %ERRORLEVEL%
) else (
    echo.
    echo All checks passed!
    exit /b 0
)





