"""
修复并启动所有服务
"""

import sys
import os
import time
import requests
import subprocess
from pathlib import Path

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def check_service(url, name, timeout=3):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=timeout)
        return r.status_code == 200
    except:
        return False


def check_port(port):
    """检查端口占用的进程"""
    try:
        result = subprocess.run(
            ["netstat", "-ano"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=5
        )
        processes = []
        for line in result.stdout.splitlines():
            if f":{port}" in line and "LISTENING" in line:
                parts = line.split()
                if len(parts) >= 5:
                    pid = parts[-1]
                    processes.append(pid)
        return processes
    except Exception as e:
        print(f"[ERROR] 检查端口 {port} 时出错: {e}")
        return []


def stop_process(pid, service_name):
    """停止进程"""
    try:
        result = subprocess.run(
            ["taskkill", "/F", "/PID", pid],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=5
        )
        if result.returncode == 0:
            print(f"  [OK] 已停止进程 {pid} ({service_name})")
            return True
        else:
            print(f"  [WARNING] 进程 {pid} 已停止或不存在")
            return False
    except Exception as e:
        print(f"  [ERROR] 停止进程 {pid} 时出错: {e}")
        return False


def stop_service_port(port, service_name):
    """停止指定端口的服务"""
    print(f"\n[步骤] 停止 {service_name} (端口 {port})...")
    processes = check_port(port)
    if not processes:
        print(f"  [INFO] 端口 {port} 上没有运行的服务")
        return
    
    for pid in processes:
        stop_process(pid, service_name)
    time.sleep(1)


def stop_wsl_flashdock(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
    """停止 WSL 中的 FlashDock streamlit 进程"""
    if sys.platform != "win32":
        return False
    
    print(f"\n[步骤] 停止 WSL 中的 FlashDock 服务...")
    print(f"  WSL 发行版: {wsl_distro}")
    
    try:
        # 检查 WSL 是否可用
        result = subprocess.run(
            ["wsl", "--status"],
            capture_output=True,
            timeout=5
        )
        if result.returncode != 0:
            print("  [INFO] WSL 不可用，跳过 WSL 进程停止")
            return False
        
        # 查找 WSL 中运行 streamlit 的进程
        find_cmd = (
            f"ps aux | grep -E 'streamlit.*FlashDock|streamlit.*8501' | "
            f"grep -v grep | awk '{{print $2}}'"
        )
        
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", find_cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        
        if result.returncode == 0 and result.stdout.strip():
            pids = [pid.strip() for pid in result.stdout.strip().split('\n') if pid.strip()]
            if pids:
                print(f"  找到 {len(pids)} 个 streamlit 进程: {', '.join(pids)}")
                
                # 停止所有找到的进程
                for pid in pids:
                    kill_cmd = f"kill -9 {pid} 2>/dev/null || true"
                    kill_result = subprocess.run(
                        ["wsl", "-d", wsl_distro, "bash", "-c", kill_cmd],
                        capture_output=True,
                        text=True,
                        encoding='utf-8',
                        errors='ignore',
                        timeout=5
                    )
                    if kill_result.returncode == 0:
                        print(f"  [OK] 已停止 WSL 进程 {pid}")
                
                time.sleep(2)
                print("  [OK] WSL 中的 FlashDock 进程已停止")
                return True
            else:
                print("  [INFO] WSL 中没有找到运行中的 streamlit 进程")
                return True
        else:
            print("  [INFO] WSL 中没有找到运行中的 streamlit 进程")
            return True
            
    except subprocess.TimeoutExpired:
        print("  [ERROR] WSL 命令执行超时")
        return False
    except Exception as e:
        print(f"  [ERROR] 停止 WSL FlashDock 进程时出错: {e}")
        return False


def stop_all_services():
    """停止所有服务（包括 WSL 中的）"""
    print("=" * 60)
    print("停止所有现有服务")
    print("=" * 60)
    
    # 停止服务注册中心 (端口 8500)
    stop_service_port(8500, "服务注册中心")
    
    # 停止 COMPASS 服务 (端口 8080)
    stop_service_port(8080, "COMPASS 服务")
    
    # 停止 FLASH-DOCK 前端 (端口 8501)
    stop_service_port(8501, "FLASH-DOCK 前端")
    
    # 停止 WSL 中的 FlashDock 进程
    stop_wsl_flashdock()
    
    # 等待所有进程完全停止
    print("\n等待所有进程完全停止...")
    time.sleep(3)
    
    print("\n" + "=" * 60)
    print("所有服务已成功停止")
    print("=" * 60)
    print()


def find_python_executable(env_name="AIDDTRAIN"):
    """查找conda环境中的Python可执行文件"""
    # 支持多个环境名称（如 flash_dock 和 flashdocker）
    env_names = [env_name]
    if env_name == "flash_dock":
        env_names = ["flash_dock", "flashdocker"]
    elif env_name == "flashdocker":
        env_names = ["flashdocker", "flash_dock"]
    
    possible_paths = []
    for name in env_names:
        possible_paths.extend([
            f"D:\\conda_envs\\{name}\\python.exe",
            f"C:\\ProgramData\\Anaconda3\\envs\\{name}\\python.exe",
            f"{os.environ.get('USERPROFILE', '')}\\anaconda3\\envs\\{name}\\python.exe",
            f"{os.environ.get('USERPROFILE', '')}\\miniconda3\\envs\\{name}\\python.exe",
        ])
    
    for path in possible_paths:
        if os.path.exists(path):
            return path
    return None


def start_registry_service(project_root, python_exe=None):
    """启动服务注册中心（使用Windows环境）"""
    if python_exe is None:
        python_exe = sys.executable
    
    env = os.environ.copy()
    env["PYTHONPATH"] = str(project_root)
    
    cmd = [
        python_exe,
        str(project_root / "services" / "registry" / "server.py"),
        "--host", "0.0.0.0",
        "--port", "8500",
    ]
    
    print(f"启动服务注册中心: {' '.join(cmd)}")
    print(f"工作目录: {project_root}")
    print(f"PYTHONPATH: {env['PYTHONPATH']}")
    
    try:
        if sys.platform == "win32":
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                creationflags=subprocess.CREATE_NEW_CONSOLE,
            )
        else:
            process = subprocess.Popen(cmd, env=env, cwd=str(project_root))
        
        print(f"[INFO] 服务进程已启动 (PID: {process.pid})")
    except Exception as e:
        print(f"[ERROR] 启动服务注册中心失败: {e}")
        return False
    
    # 等待服务启动（增加等待时间和重试次数）
    print("等待服务注册中心启动...")
    for i in range(15):  # 增加到15秒
        time.sleep(1)
        if check_service("http://localhost:8500/health", "服务注册中心", timeout=2):
            print("[OK] 服务注册中心已启动")
            return True
        if (i + 1) % 3 == 0:  # 每3秒显示一次进度
            print(f"  等待中... ({i+1}/15)")
    
    print("[WARNING] 服务注册中心可能未完全启动，请检查新窗口中的错误信息")
    print("提示: 服务可能仍在启动中，请稍后手动检查 http://localhost:8500/health")
    return False


def start_compass_service(project_root, python_exe=None):
    """启动COMPASS服务（使用Windows环境）"""
    if python_exe is None:
        python_exe = sys.executable
    
    env = os.environ.copy()
    env["PYTHONPATH"] = str(project_root)
    
    cmd = [
        python_exe,
        str(project_root / "compass" / "service_main.py"),
        "--host", "0.0.0.0",
        "--port", "8080",
        "--registry-url", "http://localhost:8500",
    ]
    
    print(f"启动COMPASS服务: {' '.join(cmd)}")
    print(f"工作目录: {project_root}")
    print(f"PYTHONPATH: {env['PYTHONPATH']}")
    
    try:
        if sys.platform == "win32":
            process = subprocess.Popen(
                cmd,
                env=env,
                cwd=str(project_root),
                creationflags=subprocess.CREATE_NEW_CONSOLE,
            )
        else:
            process = subprocess.Popen(cmd, env=env, cwd=str(project_root))
        
        print(f"[INFO] 服务进程已启动 (PID: {process.pid})")
    except Exception as e:
        print(f"[ERROR] 启动COMPASS服务失败: {e}")
        return False
    
    # 等待服务启动（增加等待时间）
    print("等待COMPASS服务启动...")
    print("提示: COMPASS服务启动可能需要更长时间，请耐心等待")
    for i in range(20):  # 增加到20秒
        time.sleep(1)
        if check_service("http://localhost:8080/health", "COMPASS服务", timeout=2):
            print("[OK] COMPASS服务已启动")
            return True
        if (i + 1) % 5 == 0:  # 每5秒显示一次进度
            print(f"  等待中... ({i+1}/20)")
    
    print("[WARNING] COMPASS服务可能未完全启动，请检查新窗口中的错误信息")
    print("提示: 服务可能仍在启动中，请稍后手动检查 http://localhost:8080/health")
    return False


def check_wsl_flashdock_ready(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
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
        
        # 检查WSL发行版是否存在（即使停止也会返回成功）
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "--", "echo", "test"],
            capture_output=True,
            timeout=10  # 增加超时时间，因为可能需要启动WSL
        )
        if result.returncode != 0:
            return False
        
        # 检查conda环境是否存在
        check_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda env list | grep -q '^{env_name} '"
        )
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", check_cmd],
            capture_output=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        return result.returncode == 0
    except Exception as e:
        print(f"[DEBUG] WSL检测失败: {e}")
        return False


def verify_wsl_dependencies(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
    """验证 WSL 中的关键依赖"""
    if sys.platform != "win32":
        return False, "非 Windows 系统"
    
    print(f"\n[步骤] 验证 WSL 中的依赖...")
    print(f"  WSL 发行版: {wsl_distro}")
    print(f"  环境名称: {env_name}")
    
    try:
        # 验证关键依赖
        verify_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {env_name} && "
            f"python -c \"import sh; import streamlit; import streamlit_molstar; print('OK')\""
        )
        
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", verify_cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=15
        )
        
        if result.returncode == 0 and "OK" in result.stdout:
            print("  [OK] 关键依赖验证通过 (sh, streamlit, streamlit_molstar)")
            return True, "依赖验证通过"
        else:
            error_msg = result.stderr.strip() if result.stderr else result.stdout.strip()
            print(f"  [WARNING] 依赖验证失败: {error_msg}")
            
            # 尝试安装缺失的 sh 模块
            if "No module named 'sh'" in error_msg or "ModuleNotFoundError" in error_msg:
                print("  [INFO] 尝试安装缺失的 sh 模块...")
                install_cmd = (
                    f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
                    f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
                    f"conda activate {env_name} && "
                    f"pip install sh>=1.14.0 -q"
                )
                
                install_result = subprocess.run(
                    ["wsl", "-d", wsl_distro, "bash", "-c", install_cmd],
                    capture_output=True,
                    text=True,
                    encoding='utf-8',
                    errors='ignore',
                    timeout=60
                )
                
                if install_result.returncode == 0:
                    print("  [OK] sh 模块安装成功")
                    # 再次验证
                    verify_result = subprocess.run(
                        ["wsl", "-d", wsl_distro, "bash", "-c", verify_cmd],
                        capture_output=True,
                        text=True,
                        encoding='utf-8',
                        errors='ignore',
                        timeout=15
                    )
                    if verify_result.returncode == 0 and "OK" in verify_result.stdout:
                        print("  [OK] 依赖验证通过")
                        return True, "依赖验证通过（已安装 sh 模块）"
                    else:
                        return False, "依赖验证失败（即使安装了 sh 模块）"
                else:
                    return False, f"无法安装 sh 模块: {install_result.stderr.strip()}"
            
            return False, f"依赖验证失败: {error_msg}"
            
    except subprocess.TimeoutExpired:
        return False, "依赖验证超时"
    except Exception as e:
        return False, f"依赖验证出错: {e}"


def start_flashdock_service(project_root, python_exe=None):
    """启动FLASH-DOCK服务（仅此服务使用WSL环境，WSL不可用时直接失败）"""
    flashdock_dir = project_root / "FLASH_DOCK-main"
    if not flashdock_dir.exists():
        print(f"[ERROR] FLASH-DOCK目录不存在: {flashdock_dir}")
        return False
    
    # WSL配置
    wsl_distro = "Ubuntu-24.04"
    wsl_env_name = "flash_dock"
    wsl_project_root = "/mnt/e/Qinchaojun/AIDD-TRAIN"
    wsl_flashdock_dir = f"{wsl_project_root}/FLASH_DOCK-main"
    wsl_port = "8501"
    
    # 检查WSL环境是否可用
    if sys.platform != "win32":
        print("[ERROR] FlashDock 只能在 Windows 上通过 WSL 运行")
        return False
    
    if not check_wsl_flashdock_ready(wsl_distro, wsl_env_name):
        print("=" * 60)
        print("[ERROR] WSL FlashDock 环境不可用")
        print("=" * 60)
        print("FlashDock 必须使用 WSL 环境运行")
        print()
        print("请执行以下步骤设置 WSL 环境：")
        print("1. 确保 WSL 已安装并启用")
        print(f"2. 确保 WSL 发行版 '{wsl_distro}' 存在")
        print("3. 在 WSL 中手动创建 conda 环境并安装依赖")
        print()
        return False
    
    # 验证依赖
    deps_ok, deps_msg = verify_wsl_dependencies(wsl_distro, wsl_env_name)
    if not deps_ok:
        print("=" * 60)
        print("[ERROR] WSL 依赖验证失败")
        print("=" * 60)
        print(f"错误信息: {deps_msg}")
        print()
        print("请手动在 WSL 中安装缺失的依赖：")
        print(f"  wsl -d {wsl_distro}")
        print(f"  conda activate {wsl_env_name}")
        print(f"  pip install sh>=1.14.0")
        print()
        return False
    
    # 使用WSL环境启动
    print(f"[INFO] 使用WSL环境启动FLASH-DOCK")
    print(f"  WSL发行版: {wsl_distro}")
    print(f"  环境名称: {wsl_env_name}")
    print(f"  项目路径: {wsl_project_root}")
    
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
    
    print(f"启动FLASH-DOCK服务 (WSL): wsl -d {wsl_distro} bash -c ...")
    print(f"工作目录: {wsl_flashdock_dir}")
    print()
    
    try:
        if sys.platform == "win32":
            # 使用UTF-8编码避免编码错误
            subprocess.Popen(
                cmd,
                creationflags=subprocess.CREATE_NEW_CONSOLE,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        else:
            subprocess.Popen(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
    except Exception as e:
        print(f"[ERROR] 启动 FlashDock 服务失败: {e}")
        return False
    
    # 等待服务启动（增加等待时间，因为WSL启动可能较慢）
    print("等待FLASH-DOCK服务启动...")
    print("提示: WSL环境启动可能需要更长时间，请耐心等待")
    for i in range(30):  # 增加到30秒
        time.sleep(1)
        # 尝试多个健康检查端点
        if (check_service("http://localhost:8501/_stcore/health", "FLASH-DOCK服务", timeout=2) or
            check_service("http://localhost:8501", "FLASH-DOCK服务", timeout=2)):
            print("[OK] FLASH-DOCK服务已启动")
            return True
        if (i + 1) % 5 == 0:  # 每5秒显示一次进度
            print(f"  等待中... ({i+1}/30)")
    
    print("[WARNING] FLASH-DOCK服务可能未完全启动")
    print("提示: 请检查新打开的WSL窗口中的错误信息")
    print("或者手动启动: start_flashdock_wsl.bat")
    return False


def main():
    print("=" * 60)
    print("修复并启动所有服务")
    print("=" * 60)
    print()
    
    # 获取项目根目录
    project_root = Path(__file__).parent.resolve()
    print(f"项目根目录: {project_root}")
    print()
    
    # 查找Python可执行文件
    python_exe = find_python_executable("AIDDTRAIN")
    if python_exe:
        print(f"找到AIDDTRAIN环境: {python_exe}")
    else:
        print("未找到AIDDTRAIN环境，使用默认Python")
        python_exe = sys.executable
    print()
    
    # 检查FLASH-DOCK环境（必须使用WSL）
    wsl_ready = False
    if sys.platform == "win32":
        wsl_ready = check_wsl_flashdock_ready()
        if wsl_ready:
            print("[OK] WSL FlashDock环境已就绪")
        else:
            print("[WARNING] WSL FlashDock环境不可用，FlashDock服务将无法启动")
    else:
        print("[WARNING] 非Windows系统，FlashDock服务将无法启动（需要WSL）")
    print()
    
    # 先停止所有现有服务（包括 WSL 中的）
    stop_all_services()
    
    # 检查当前服务状态
    print("检查当前服务状态...")
    registry_ok = check_service("http://localhost:8500/health", "服务注册中心", timeout=2)
    compass_ok = check_service("http://localhost:8080/health", "COMPASS服务", timeout=2)
    flashdock_ok = check_service("http://localhost:8501/_stcore/health", "FLASH-DOCK服务", timeout=2)
    
    print(f"服务注册中心 (8500): {'[OK] 运行中' if registry_ok else '[FAIL] 未运行'}")
    print(f"COMPASS服务 (8080): {'[OK] 运行中' if compass_ok else '[FAIL] 未运行'}")
    print(f"FLASH-DOCK服务 (8501): {'[OK] 运行中' if flashdock_ok else '[FAIL] 未运行'}")
    print()
    
    # 启动服务注册中心
    if not registry_ok:
        print("=" * 60)
        print("启动服务注册中心...")
        print("=" * 60)
        if not start_registry_service(project_root, python_exe):
            print("[ERROR] 服务注册中心启动失败")
            return 1
        print()
    else:
        print("[INFO] 服务注册中心已在运行")
        print()
    
    # 等待一下确保注册中心完全启动
    if not registry_ok:
        time.sleep(3)
    
    # 启动COMPASS服务
    if not compass_ok:
        print("=" * 60)
        print("启动COMPASS服务...")
        print("=" * 60)
        if not start_compass_service(project_root, python_exe):
            print("[ERROR] COMPASS服务启动失败")
            return 1
        print()
    else:
        print("[INFO] COMPASS服务已在运行")
        print()
    
    # 启动FLASH-DOCK服务（必须使用WSL）
    if not flashdock_ok and wsl_ready:
        print("=" * 60)
        print("启动FLASH-DOCK服务...")
        print("=" * 60)
        if not start_flashdock_service(project_root):
            print("[ERROR] FLASH-DOCK服务启动失败")
            print("请检查WSL环境配置")
        print()
    elif flashdock_ok:
        print("[INFO] FLASH-DOCK服务已在运行")
        print()
    elif not wsl_ready:
        print("[WARNING] 跳过FLASH-DOCK服务（WSL环境不可用）")
        print("请手动设置WSL环境：创建conda环境并安装依赖")
        print()
    
    # 最终检查
    print("=" * 60)
    print("最终服务状态检查")
    print("=" * 60)
    time.sleep(2)
    
    registry_ok = check_service("http://localhost:8500/health", "服务注册中心", timeout=3)
    compass_ok = check_service("http://localhost:8080/health", "COMPASS服务", timeout=3)
    flashdock_ok = check_service("http://localhost:8501/_stcore/health", "FLASH-DOCK服务", timeout=3)
    
    print(f"服务注册中心 (8500): {'[OK] 运行中' if registry_ok else '[FAIL] 未运行'}")
    print(f"COMPASS服务 (8080): {'[OK] 运行中' if compass_ok else '[FAIL] 未运行'}")
    if wsl_ready:
        print(f"FLASH-DOCK服务 (8501): {'[OK] 运行中' if flashdock_ok else '[FAIL] 未运行'}")
    else:
        print(f"FLASH-DOCK服务 (8501): [SKIP] WSL环境不可用")
    print()
    
    if registry_ok and compass_ok:
        print("[SUCCESS] 核心服务都已成功启动!")
        print()
        print("服务地址:")
        print("  - 服务注册中心: http://localhost:8500")
        print("  - COMPASS服务: http://localhost:8080")
        print("  - API文档: http://localhost:8080/docs")
        if flashdock_ok:
            print("  - FLASH-DOCK前端: http://localhost:8501")
        return 0
    else:
        print("[WARNING] 部分服务可能未完全启动")
        print("请检查新打开的窗口中的错误信息")
        return 1


if __name__ == "__main__":
    sys.exit(main())

