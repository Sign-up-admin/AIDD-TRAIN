@echo off
chcp 65001 >nul 2>&1
cd /d "%~dp0"

echo Testing start script variables...
echo.

set "PROJECT_ROOT=%CD%"
echo Project Root: %PROJECT_ROOT%

set "PYTHON_CMD=python"
if exist "C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe" (
    set "PYTHON_CMD=C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe"
    echo Found AIDDTRAIN Python: %PYTHON_CMD%
) else (
    if exist "%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe" (
        set "PYTHON_CMD=%USERPROFILE%\anaconda3\envs\AIDDTRAIN\python.exe"
        echo Found AIDDTRAIN Python: %PYTHON_CMD%
    ) else (
        echo Using default python
    )
)

set "FLASHDOCK_PYTHON=python"
if exist "C:\ProgramData\Anaconda3\envs\flash_dock\python.exe" (
    set "FLASHDOCK_PYTHON=C:\ProgramData\Anaconda3\envs\flash_dock\python.exe"
    echo Found flash_dock Python: %FLASHDOCK_PYTHON%
) else (
    if exist "%USERPROFILE%\anaconda3\envs\flash_dock\python.exe" (
        set "FLASHDOCK_PYTHON=%USERPROFILE%\anaconda3\envs\flash_dock\python.exe"
        echo Found flash_dock Python: %FLASHDOCK_PYTHON%
    ) else (
        echo Using default python for FLASH-DOCK
    )
)

echo.
echo All variables set correctly!
pause

