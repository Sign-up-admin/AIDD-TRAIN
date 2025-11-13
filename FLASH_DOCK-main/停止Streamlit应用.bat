@echo off
echo 正在停止Streamlit应用...
echo.

REM 查找并停止所有Streamlit进程
taskkill /F /IM streamlit.exe 2>nul
if %errorlevel% equ 0 (
    echo Streamlit进程已停止
) else (
    echo 未找到streamlit.exe进程
)

REM 查找并停止Python进程（可能运行Streamlit）
for /f "tokens=2" %%a in ('netstat -ano ^| findstr :8501 ^| findstr LISTENING') do (
    echo 找到占用8501端口的进程: %%a
    taskkill /F /PID %%a 2>nul
    if %errorlevel% equ 0 (
        echo 进程 %%a 已停止
    )
)

echo.
echo 完成！
pause








