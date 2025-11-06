@echo off
echo ============================================
echo FlashDock 项目文件清理脚本
echo ============================================
echo.
echo 将要删除以下文件：
echo.
echo [1] 重复的启动脚本
echo     - 启动FlashDock.bat
echo.
echo [2] SSL修复相关临时文件
echo     - fix_ssl_and_test.py
echo     - fix_ssl_cert.bat
echo     - fix_ssl_permanent.md
echo     - fix_windows_sh.py
echo.
echo [3] 多余的FlashDock版本（旧版本）
echo     - FlashDock_web.py
echo     - FlashDock_web_v1.py
echo.
echo [4] 测试和临时文档
echo     - Python3.12环境测试报告.md
echo     - Python3.12环境运行指南.md
echo     - Python版本兼容性分析.md
echo     - SSL修复总结.md
echo.
echo ============================================
echo.
set /p confirm="确认删除这些文件吗？(Y/N): "
if /i not "%confirm%"=="Y" (
    echo 已取消删除操作
    pause
    exit /b 0
)

echo.
echo 开始清理...
echo.

REM 删除重复的启动脚本
if exist "启动FlashDock.bat" (
    del /f /q "启动FlashDock.bat"
    echo [删除] 启动FlashDock.bat
)

REM 删除SSL修复相关文件
if exist "fix_ssl_and_test.py" (
    del /f /q "fix_ssl_and_test.py"
    echo [删除] fix_ssl_and_test.py
)
if exist "fix_ssl_cert.bat" (
    del /f /q "fix_ssl_cert.bat"
    echo [删除] fix_ssl_cert.bat
)
if exist "fix_ssl_permanent.md" (
    del /f /q "fix_ssl_permanent.md"
    echo [删除] fix_ssl_permanent.md
)
if exist "fix_windows_sh.py" (
    del /f /q "fix_windows_sh.py"
    echo [删除] fix_windows_sh.py
)

REM 删除多余的FlashDock版本
if exist "FlashDock_web.py" (
    del /f /q "FlashDock_web.py"
    echo [删除] FlashDock_web.py
)
if exist "FlashDock_web_v1.py" (
    del /f /q "FlashDock_web_v1.py"
    echo [删除] FlashDock_web_v1.py
)

REM 删除测试和临时文档
if exist "Python3.12环境测试报告.md" (
    del /f /q "Python3.12环境测试报告.md"
    echo [删除] Python3.12环境测试报告.md
)
if exist "Python3.12环境运行指南.md" (
    del /f /q "Python3.12环境运行指南.md"
    echo [删除] Python3.12环境运行指南.md
)
if exist "Python版本兼容性分析.md" (
    del /f /q "Python版本兼容性分析.md"
    echo [删除] Python版本兼容性分析.md
)
if exist "SSL修复总结.md" (
    del /f /q "SSL修复总结.md"
    echo [删除] SSL修复总结.md
)

echo.
echo ============================================
echo 清理完成！
echo ============================================
echo.
echo 保留的重要文件：
echo   - FlashDock.py (主程序)
echo   - start_flashdock_fixed.bat (启动脚本)
echo   - README.md (项目说明)
echo   - requirements.txt (依赖列表)
echo.
pause


