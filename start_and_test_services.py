"""
启动所有服务并验证服务运行状态
"""
import sys
import os
import time
import requests
import subprocess
from pathlib import Path
from typing import Tuple, Optional

# 服务配置
REGISTRY_PORT = 8500
COMPASS_PORT = 8080
FLASHDOCK_PORT = 8501

REGISTRY_URL = f"http://localhost:{REGISTRY_PORT}"
COMPASS_URL = f"http://localhost:{COMPASS_PORT}"
FLASHDOCK_URL = f"http://localhost:{FLASHDOCK_PORT}"

def check_service(url: str, name: str, timeout: int = 5) -> Tuple[bool, Optional[str]]:
    """检查服务是否运行"""
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code == 200:
            return True, None
        else:
            return False, f"Status code: {response.status_code}"
    except requests.exceptions.ConnectionError:
        return False, "Connection refused"
    except requests.exceptions.Timeout:
        return False, "Timeout"
    except Exception as e:
        return False, str(e)

def find_python_executable(env_name: str) -> Optional[str]:
    """查找Conda环境的Python可执行文件"""
    # 检查常见的Conda安装位置
    conda_locations = [
        os.path.join(os.environ.get("USERPROFILE", ""), "anaconda3", "envs", env_name, "python.exe"),
        os.path.join("C:\\ProgramData\\Anaconda3", "envs", env_name, "python.exe"),
        os.path.join(os.environ.get("USERPROFILE", ""), "miniconda3", "envs", env_name, "python.exe"),
        os.path.join("C:\\ProgramData\\Miniconda3", "envs", env_name, "python.exe"),
    ]
    
    for path in conda_locations:
        if os.path.exists(path):
            return path
    return None

def start_registry_service(project_root: Path) -> Optional[subprocess.Popen]:
    """启动服务注册中心"""
    python_cmd = find_python_executable("AIDDTRAIN") or sys.executable
    
    env = os.environ.copy()
    env['PYTHONPATH'] = str(project_root)
    
    cmd = [
        python_cmd,
        str(project_root / "services" / "registry" / "server.py"),
        "--host", "0.0.0.0",
        "--port", str(REGISTRY_PORT)
    ]
    
    try:
        if sys.platform == 'win32':
            # 在新窗口中启动，不使用 PIPE（这样可以看到错误信息）
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                creationflags=subprocess.CREATE_NEW_CONSOLE
            )
        else:
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
        print(f"[INFO] 服务注册中心进程已启动 (PID: {process.pid})")
        return process
    except Exception as e:
        print(f"Failed to start registry service: {e}")
        return None

def start_compass_service(project_root: Path) -> Optional[subprocess.Popen]:
    """启动COMPASS服务"""
    python_cmd = find_python_executable("AIDDTRAIN") or sys.executable
    
    env = os.environ.copy()
    env['PYTHONPATH'] = str(project_root)
    
    cmd = [
        python_cmd,
        str(project_root / "compass" / "service_main.py"),
        "--host", "0.0.0.0",
        "--port", str(COMPASS_PORT),
        "--registry-url", REGISTRY_URL
    ]
    
    try:
        if sys.platform == 'win32':
            # 在新窗口中启动，不使用 PIPE（这样可以看到错误信息）
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                creationflags=subprocess.CREATE_NEW_CONSOLE
            )
        else:
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            )
        print(f"[INFO] COMPASS 服务进程已启动 (PID: {process.pid})")
        print(f"[INFO] 服务将在新窗口中运行，请查看新窗口中的日志")
        return process
    except Exception as e:
        print(f"Failed to start COMPASS service: {e}")
        return None

def check_wsl_flashdock_ready(wsl_distro: str = "Ubuntu-24.04", env_name: str = "flash_dock") -> bool:
    """检查WSL FlashDock环境是否可用"""
    if sys.platform != "win32":
        return False
    
    try:
        # 检查WSL是否可用
        result = subprocess.run(
            ["wsl", "--status"],
            capture_output=True,
            timeout=5
        )
        if result.returncode != 0:
            return False
        
        # 检查WSL发行版是否存在（即使停止也会自动启动）
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "--", "echo", "test"],
            capture_output=True,
            timeout=10  # 增加超时时间，因为可能需要启动WSL
        )
        if result.returncode != 0:
            return False
        
        # 检查conda环境是否存在（需要先初始化conda）
        check_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda env list | grep -q '^{env_name} '"
        )
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", check_cmd],
            capture_output=True,
            timeout=10
        )
        return result.returncode == 0
    except Exception as e:
        print(f"[DEBUG] WSL检测失败: {e}")
        return False


def start_flashdock_service(project_root: Path) -> Optional[subprocess.Popen]:
    """启动FLASH-DOCK前端服务（优先使用WSL环境）"""
    flashdock_dir = project_root / "FLASH_DOCK-main"
    
    # WSL配置
    wsl_distro = "Ubuntu-24.04"
    wsl_env_name = "flash_dock"
    wsl_project_root = "/mnt/e/Qinchaojun/AIDD-TRAIN"
    wsl_flashdock_dir = f"{wsl_project_root}/FLASH_DOCK-main"
    wsl_port = str(FLASHDOCK_PORT)
    
    # 优先尝试使用WSL环境
    if sys.platform == "win32" and check_wsl_flashdock_ready(wsl_distro, wsl_env_name):
        print(f"[INFO] Using WSL environment to start FLASH-DOCK")
        print(f"  WSL Distribution: {wsl_distro}")
        print(f"  Environment: {wsl_env_name}")
        print(f"  Project Path: {wsl_project_root}")
        
        # 构建WSL命令
        wsl_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {wsl_env_name} && "
            f"export PYTHONPATH={wsl_project_root} && "
            f"cd {wsl_flashdock_dir} && "
            f"streamlit run FlashDock.py --server.port {wsl_port} --server.address 0.0.0.0"
        )
        
        cmd = ["wsl", "-d", wsl_distro, "bash", "-c", wsl_cmd]
        
        try:
            if sys.platform == 'win32':
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    creationflags=subprocess.CREATE_NEW_CONSOLE
                )
            else:
                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
            return process
        except Exception as e:
            print(f"Failed to start FLASH-DOCK service (WSL): {e}")
            return None
    else:
        # WSL环境不可用，直接失败
        print("=" * 60)
        print("[ERROR] WSL FlashDock environment not available")
        print("=" * 60)
        print("FlashDock must run in WSL environment")
        print()
        print("Please follow these steps to set up WSL environment:")
        print("1. Ensure WSL is installed and enabled")
        print(f"2. Ensure WSL distribution '{wsl_distro}' exists")
        print("3. Manually create conda environment and install dependencies in WSL")
        print()
        return None

def wait_for_service(url: str, name: str, max_wait: int = 30, interval: int = 2) -> bool:
    """等待服务启动"""
    print(f"Waiting for {name} to start...")
    start_time = time.time()
    
    while time.time() - start_time < max_wait:
        is_running, error = check_service(url, name, timeout=2)
        if is_running:
            print(f"[OK] {name} is running")
            return True
        elapsed = int(time.time() - start_time)
        time.sleep(interval)
        # 每10秒显示一次进度，避免输出过多
        if elapsed % 10 == 0 or elapsed < 5:
            print(f"  Still waiting... ({elapsed}/{max_wait}s)")
    
    print(f"[WARNING] {name} may still be starting (waited {max_wait} seconds)")
    return False

def main():
    """主函数"""
    print("=" * 70)
    print("启动所有服务并验证")
    print("=" * 70)
    print()
    
    project_root = Path(__file__).parent.resolve()
    
    # 检查服务是否已经运行
    print("检查现有服务状态...")
    registry_running, _ = check_service(f"{REGISTRY_URL}/health", "Registry")
    compass_running, _ = check_service(f"{COMPASS_URL}/health", "COMPASS")
    flashdock_running, _ = check_service(FLASHDOCK_URL, "FLASH-DOCK")
    
    processes = []
    
    # 启动服务注册中心
    if not registry_running:
        print("\n启动服务注册中心...")
        process = start_registry_service(project_root)
        if process:
            processes.append(("Registry", process))
            if not wait_for_service(f"{REGISTRY_URL}/health", "Registry", max_wait=30):
                print("[ERROR] Failed to start Registry service")
                return 1
        else:
            print("[ERROR] Failed to start Registry service")
            return 1
    else:
        print("[OK] Registry service is already running")
    
    # 启动COMPASS服务
    if not compass_running:
        print("\n启动COMPASS服务...")
        print("提示: COMPASS 服务启动可能需要更长时间（30-60秒），请耐心等待")
        process = start_compass_service(project_root)
        if process:
            processes.append(("COMPASS", process))
            # 增加等待时间到60秒，因为 COMPASS 服务启动较慢
            if not wait_for_service(f"{COMPASS_URL}/health", "COMPASS", max_wait=60):
                print("[WARNING] COMPASS 服务可能仍在启动中")
                print("提示: 请检查新窗口中的服务状态和错误信息")
                print("如果服务正常启动，可以稍后手动检查 http://localhost:8080/health")
                # 不直接返回错误，让服务继续尝试启动
        else:
            print("[ERROR] Failed to start COMPASS service")
            return 1
    else:
        print("[OK] COMPASS service is already running")
    
    # 启动FLASH-DOCK前端
    if not flashdock_running:
        print("\n启动FLASH-DOCK前端...")
        process = start_flashdock_service(project_root)
        if process:
            processes.append(("FLASH-DOCK", process))
            # FLASH-DOCK没有health端点，检查主页
            if not wait_for_service(FLASHDOCK_URL, "FLASH-DOCK", max_wait=60):
                print("[WARNING] FLASH-DOCK may still be starting...")
        else:
            print("[ERROR] Failed to start FLASH-DOCK service")
            return 1
    else:
        print("[OK] FLASH-DOCK service is already running")
    
    # 最终验证
    print("\n" + "=" * 70)
    print("服务状态验证")
    print("=" * 70)
    
    all_ok = True
    registry_ok, registry_error = check_service(f"{REGISTRY_URL}/health", "Registry")
    compass_ok, compass_error = check_service(f"{COMPASS_URL}/health", "COMPASS")
    flashdock_ok, flashdock_error = check_service(FLASHDOCK_URL, "FLASH-DOCK")
    
    print(f"Registry (8500): {'[OK]' if registry_ok else f'[FAIL] {registry_error}'}")
    print(f"COMPASS (8080): {'[OK]' if compass_ok else f'[FAIL] {compass_error}'}")
    print(f"FLASH-DOCK (8501): {'[OK]' if flashdock_ok else f'[FAIL] {flashdock_error}'}")
    
    if not (registry_ok and compass_ok and flashdock_ok):
        all_ok = False
    
    if all_ok:
        print("\n[SUCCESS] All services are running!")
        print("\n服务地址:")
        print(f"  - Registry: {REGISTRY_URL}")
        print(f"  - COMPASS: {COMPASS_URL}")
        print(f"  - COMPASS API Docs: {COMPASS_URL}/docs")
        print(f"  - FLASH-DOCK: {FLASHDOCK_URL}")
        print("\n可以开始测试了！")
        return 0
    else:
        print("\n[ERROR] Some services failed to start")
        return 1

if __name__ == "__main__":
    sys.exit(main())

