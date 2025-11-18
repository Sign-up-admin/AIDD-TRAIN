"""
测试FLASH-DOCK启动
"""

import sys
import subprocess
import time
import requests

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def check_service(url, timeout=3):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=timeout)
        return r.status_code < 500
    except:
        return False

def main():
    print("=" * 60)
    print("FLASH-DOCK 启动测试")
    print("=" * 60)
    print()
    
    # 检查当前状态
    print("1. 检查当前服务状态...")
    if check_service("http://localhost:8501"):
        print("   [OK] FLASH-DOCK 已在运行")
        print("   访问地址: http://localhost:8501")
        return 0
    else:
        print("   [INFO] FLASH-DOCK 未运行")
    
    print()
    print("2. 建议的启动方法:")
    print("   - 方法1: start_flashdock_fixed.bat (推荐)")
    print("   - 方法2: start_flashdock_wsl.bat (已修复环境名称)")
    print("   - 方法3: python fix_and_start_services.py")
    print()
    print("3. 环境信息:")
    print("   - WSL发行版: Ubuntu-24.04")
    print("   - Conda环境: flash_dock")
    print("   - 端口: 8501")
    print()
    
    # 询问是否启动
    print("是否现在启动 FLASH-DOCK? (y/n): ", end="")
    try:
        choice = input().strip().lower()
        if choice == 'y':
            print()
            print("正在启动 FLASH-DOCK...")
            print("提示: 服务将在新窗口中启动")
            print()
            
            # 使用修复后的启动脚本
            try:
                subprocess.Popen(
                    ["start_flashdock_fixed.bat"],
                    shell=True,
                    cwd="."
                )
                print("[OK] 启动命令已执行")
                print()
                print("等待服务启动（约15秒）...")
                
                # 等待服务启动
                for i in range(15):
                    time.sleep(1)
                    if check_service("http://localhost:8501"):
                        print()
                        print("[SUCCESS] FLASH-DOCK 已成功启动!")
                        print("访问地址: http://localhost:8501")
                        return 0
                    if (i + 1) % 5 == 0:
                        print(f"  等待中... ({i+1}/15)")
                
                print()
                print("[WARNING] 服务可能还在启动中")
                print("请检查新打开的窗口中的状态")
                print("或手动访问: http://localhost:8501")
                return 1
            except Exception as e:
                print(f"[ERROR] 启动失败: {e}")
                return 1
        else:
            print("已取消启动")
            return 0
    except KeyboardInterrupt:
        print("\n已取消")
        return 0

if __name__ == "__main__":
    sys.exit(main())

