"""检查并启动所有服务"""
import subprocess
import sys
import time
import requests
from pathlib import Path

def check_service(url, name, timeout=2):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=timeout)
        return r.status_code == 200 or r.status_code < 500
    except:
        return False

def start_service(name, cmd, cwd):
    """启动服务"""
    print(f"\n启动 {name}...")
    try:
        if sys.platform == 'win32':
            subprocess.Popen(
                cmd,
                cwd=cwd,
                creationflags=subprocess.CREATE_NEW_CONSOLE
            )
        else:
            subprocess.Popen(cmd, cwd=cwd)
        print(f"[OK] {name} 已在新窗口启动")
        return True
    except Exception as e:
        print(f"[ERROR] 启动 {name} 失败: {e}")
        return False

def main():
    project_root = Path(__file__).parent.resolve()
    python_cmd = sys.executable
    
    print("=" * 60)
    print("检查并启动所有服务")
    print("=" * 60)
    
    # 检查服务状态
    print("\n检查服务状态...")
    registry_ok = check_service("http://localhost:8500/health", "Registry")
    compass_ok = check_service("http://localhost:8080/health", "COMPASS")
    flashdock_ok = check_service("http://localhost:8501", "FLASH-DOCK")
    
    print(f"Registry (8500): {'[运行中]' if registry_ok else '[未运行]'}")
    print(f"COMPASS (8080): {'[运行中]' if compass_ok else '[未运行]'}")
    print(f"FLASH-DOCK (8501): {'[运行中]' if flashdock_ok else '[未运行]'}")
    
    # 启动服务注册中心
    if not registry_ok:
        env = os.environ.copy()
        env['PYTHONPATH'] = str(project_root)
        cmd = [python_cmd, str(project_root / "services" / "registry" / "server.py"),
               "--host", "0.0.0.0", "--port", "8500"]
        start_service("服务注册中心", cmd, str(project_root))
        time.sleep(3)
    
    # 启动 COMPASS 服务
    if not compass_ok:
        env = os.environ.copy()
        env['PYTHONPATH'] = str(project_root)
        cmd = [python_cmd, str(project_root / "compass" / "service_main.py"),
               "--host", "0.0.0.0", "--port", "8080",
               "--registry-url", "http://localhost:8500"]
        start_service("COMPASS 服务", cmd, str(project_root))
        time.sleep(3)
    
    # 启动 FLASH-DOCK
    if not flashdock_ok:
        wsl_cmd = (
            'wsl -d Ubuntu-24.04 bash -c "'
            'source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || '
            'source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; '
            'conda activate flash_dock && '
            'export PYTHONPATH=/mnt/e/Qinchaojun/AIDD-TRAIN && '
            'cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main && '
            'streamlit run FlashDock.py --server.port 8501 --server.address 0.0.0.0"'
        )
        if sys.platform == 'win32':
            subprocess.Popen(
                wsl_cmd,
                shell=True,
                creationflags=subprocess.CREATE_NEW_CONSOLE
            )
        else:
            subprocess.Popen(wsl_cmd, shell=True)
        print("[OK] FLASH-DOCK 已在新窗口启动")
        time.sleep(5)
    
    # 最终检查
    print("\n" + "=" * 60)
    print("最终服务状态")
    print("=" * 60)
    time.sleep(5)
    
    registry_ok = check_service("http://localhost:8500/health", "Registry")
    compass_ok = check_service("http://localhost:8080/health", "COMPASS")
    flashdock_ok = check_service("http://localhost:8501", "FLASH-DOCK")
    
    print(f"Registry (8500): {'[OK]' if registry_ok else '[FAIL]'}")
    print(f"COMPASS (8080): {'[OK]' if compass_ok else '[FAIL]'}")
    print(f"FLASH-DOCK (8501): {'[OK]' if flashdock_ok else '[FAIL]'}")
    
    if registry_ok and compass_ok and flashdock_ok:
        print("\n[SUCCESS] 所有服务已启动！")
        print("\n服务地址:")
        print("  - Registry: http://localhost:8500")
        print("  - COMPASS: http://localhost:8080")
        print("  - FLASH-DOCK: http://localhost:8501")
    else:
        print("\n[WARNING] 部分服务可能还在启动中，请检查新窗口中的状态")

if __name__ == "__main__":
    import os
    sys.exit(main())

