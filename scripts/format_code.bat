@echo off
REM Format code with Black
REM Usage: format_code.bat

cd /d "%~dp0\.."
set PYTHONPATH=%CD%;%PYTHONPATH%
echo Formatting code with Black...
echo.

python -m black compass services FLASH_DOCK-main\services FLASH_DOCK-main\pages

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Formatting completed with errors.
    exit /b %ERRORLEVEL%
) else (
    echo.
    echo Code formatting completed successfully!
    exit /b 0
)





