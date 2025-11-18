"""检查FLASH-DOCK日志和错误"""
import subprocess
import sys

def check_logs():
    print("=" * 60)
    print("检查 FLASH-DOCK 日志和错误")
    print("=" * 60)
    print()
    
    # 检查Streamlit日志
    print("1. 检查Streamlit日志...")
    log_paths = [
        "~/.streamlit/logs/",
        "/root/.streamlit/logs/",
    ]
    
    for log_path in log_paths:
        cmd = f"ls -la {log_path} 2>/dev/null | head -10 || echo 'Log directory not found'"
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0 and "not found" not in result.stdout:
            print(f"   [INFO] 日志目录: {log_path}")
            print(f"   {result.stdout}")
        else:
            print(f"   [INFO] {log_path} 不存在")
    
    # 检查Python错误
    print("\n2. 检查Python进程状态...")
    cmd = "ps aux | grep 'streamlit.*FlashDock' | grep -v grep | head -1"
    result = subprocess.run(
        ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", cmd],
        capture_output=True,
        text=True,
        encoding='utf-8',
        errors='ignore',
        timeout=10
    )
    if result.returncode == 0 and result.stdout.strip():
        print("   [OK] Streamlit进程信息:")
        print(f"   {result.stdout.strip()}")
    else:
        print("   [WARNING] 未找到Streamlit进程")
    
    # 尝试获取stderr输出
    print("\n3. 检查可能的错误...")
    print("   提示: 请查看启动FLASH-DOCK的WSL窗口")
    print("   如果窗口已关闭，可以重新运行启动脚本查看错误")
    
    print("\n" + "=" * 60)
    print("建议操作:")
    print("=" * 60)
    print("1. 查看启动FLASH-DOCK的WSL窗口中的错误信息")
    print("2. 如果窗口已关闭，运行: start_flashdock_debug.bat")
    print("3. 在浏览器中按F12查看控制台错误")
    print("4. 检查浏览器Network标签中的请求状态")

if __name__ == "__main__":
    check_logs()

