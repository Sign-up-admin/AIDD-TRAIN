"""直接启动所有服务（用于调试）"""
import os
import sys
import subprocess
import time
from pathlib import Path

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

def find_python_aiddtrain():
    """查找 AIDDTRAIN 环境的 Python"""
    possible_paths = [
        r"D:\conda_envs\AIDDTRAIN\python.exe",
        r"C:\ProgramData\Anaconda3\envs\AIDDTRAIN\python.exe",
        os.path.expanduser(r"~\anaconda3\envs\AIDDTRAIN\python.exe"),
        os.path.expanduser(r"~\miniconda3\envs\AIDDTRAIN\python.exe"),
    ]
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    
    return None

def main():
    project_root = Path(__file__).parent.resolve()
    
    # 尝试使用 AIDDTRAIN 环境的 Python
    python_aiddtrain = find_python_aiddtrain()
    if python_aiddtrain:
        python_cmd = python_aiddtrain
        print(f"[INFO] 使用 AIDDTRAIN 环境: {python_cmd}")
    else:
        python_cmd = sys.executable
        print(f"[WARNING] 未找到 AIDDTRAIN 环境，使用当前 Python: {python_cmd}")
    
    print("=" * 60)
    print("启动所有服务")
    print("=" * 60)
    print(f"项目根目录: {project_root}")
    print(f"Python: {python_cmd}")
    print()
    
    # 设置环境变量
    env = os.environ.copy()
    env['PYTHONPATH'] = str(project_root)
    
    # 启动服务注册中心
    print("启动服务注册中心 (端口 8500)...")
    registry_cmd = [
        python_cmd,
        str(project_root / "services" / "registry" / "server.py"),
        "--host", "0.0.0.0",
        "--port", "8500"
    ]
    print(f"命令: {' '.join(registry_cmd)}")
    
    registry_process = subprocess.Popen(
        registry_cmd,
        cwd=str(project_root),
        env=env,
        creationflags=subprocess.CREATE_NEW_CONSOLE if sys.platform == 'win32' else 0
    )
    print(f"[OK] 服务注册中心已启动 (PID: {registry_process.pid})")
    print("等待服务注册中心启动...")
    time.sleep(5)
    
    # 启动 COMPASS 服务
    print("\n启动 COMPASS 服务 (端口 8080)...")
    compass_cmd = [
        python_cmd,
        str(project_root / "compass" / "service_main.py"),
        "--host", "0.0.0.0",
        "--port", "8080",
        "--registry-url", "http://localhost:8500"
    ]
    print(f"命令: {' '.join(compass_cmd)}")
    
    compass_process = subprocess.Popen(
        compass_cmd,
        cwd=str(project_root),
        env=env,
        creationflags=subprocess.CREATE_NEW_CONSOLE if sys.platform == 'win32' else 0
    )
    print(f"[OK] COMPASS 服务已启动 (PID: {compass_process.pid})")
    print("等待 COMPASS 服务启动...")
    time.sleep(5)
    
    # 启动 FLASH-DOCK (在 WSL 中运行)
    print("\n启动 FLASH-DOCK (端口 8501)...")
    print("检查 FLASH-DOCK 状态...")
    try:
        import requests
        flashdock_check = requests.get("http://localhost:8501", timeout=2)
        if flashdock_check.status_code == 200:
            print("[INFO] FLASH-DOCK 已在运行")
        else:
            raise Exception("FLASH-DOCK 未正常运行")
    except:
        print("FLASH-DOCK 未运行，正在启动...")
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
            flashdock_process = subprocess.Popen(
                wsl_cmd,
                shell=True,
                creationflags=subprocess.CREATE_NEW_CONSOLE
            )
            print(f"[OK] FLASH-DOCK 已在新窗口启动 (PID: {flashdock_process.pid})")
        else:
            flashdock_process = subprocess.Popen(wsl_cmd, shell=True)
            print(f"[OK] FLASH-DOCK 已启动 (PID: {flashdock_process.pid})")
        print("等待 FLASH-DOCK 启动...")
        time.sleep(5)
    
    print("\n" + "=" * 60)
    print("服务启动完成")
    print("=" * 60)
    print("服务注册中心: http://localhost:8500")
    print("COMPASS 服务: http://localhost:8080")
    print("FLASH-DOCK: http://localhost:8501")
    print()
    print("提示: 服务在新窗口中运行，请检查新窗口中的状态")
    print("如果服务启动失败，请查看新窗口中的错误信息")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n中断启动")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] 启动失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

