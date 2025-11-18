"""
停止所有服务
"""
import subprocess
import sys
import time
import requests
import os

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


def check_port(port: int) -> list:
    """检查端口占用的进程"""
    try:
        # 使用netstat查找占用端口的进程
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


def stop_process(pid: str, service_name: str) -> bool:
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


def stop_service(port: int, service_name: str):
    """停止指定端口的服务"""
    print(f"\n[步骤] 停止 {service_name} (端口 {port})...")
    
    processes = check_port(port)
    if not processes:
        print(f"  [INFO] 端口 {port} 上没有运行的服务")
        return
    
    for pid in processes:
        stop_process(pid, service_name)
    
    # 等待进程完全停止
    time.sleep(1)


def verify_service_stopped(url: str, service_name: str) -> bool:
    """验证服务是否已停止"""
    try:
        response = requests.get(url, timeout=2)
        return False  # 服务仍在运行
    except:
        return True  # 服务已停止


def stop_wsl_flashdock(wsl_distro: str = "Ubuntu-24.04", env_name: str = "flash_dock") -> bool:
    """停止 WSL 中的 FlashDock streamlit 进程"""
    if sys.platform != "win32":
        return False
    
    print(f"\n[步骤] 停止 WSL 中的 FlashDock 服务...")
    print(f"  WSL 发行版: {wsl_distro}")
    print(f"  环境名称: {env_name}")
    
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
                        timeout=5
                    )
                    if kill_result.returncode == 0:
                        print(f"  [OK] 已停止 WSL 进程 {pid}")
                    else:
                        print(f"  [WARNING] 停止 WSL 进程 {pid} 失败")
                
                # 等待进程完全停止
                time.sleep(2)
                
                # 再次检查是否还有进程
                verify_result = subprocess.run(
                    ["wsl", "-d", wsl_distro, "bash", "-c", find_cmd],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                
                if verify_result.returncode == 0 and verify_result.stdout.strip():
                    print("  [WARNING] 仍有 streamlit 进程在运行")
                    return False
                else:
                    print("  [OK] WSL 中的 FlashDock 进程已全部停止")
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


def main():
    """主函数"""
    print("=" * 60)
    print("停止 COMPASS 项目所有服务")
    print("=" * 60)
    
    # 停止服务注册中心 (端口 8500)
    stop_service(8500, "服务注册中心")
    
    # 停止 COMPASS 服务 (端口 8080)
    stop_service(8080, "COMPASS 服务")
    
    # 停止 FLASH-DOCK 前端 (端口 8501)
    # 先停止 Windows 上的进程
    stop_service(8501, "FLASH-DOCK 前端")
    
    # 停止 WSL 中的 FlashDock 进程
    stop_wsl_flashdock()
    
    # 等待所有进程完全停止
    print("\n等待所有进程完全停止...")
    time.sleep(2)
    
    # 验证服务已停止
    print("\n验证服务状态...")
    services = [
        ("http://localhost:8500/health", "服务注册中心"),
        ("http://localhost:8080/health", "COMPASS 服务"),
        ("http://localhost:8501/_stcore/health", "FLASH-DOCK 前端"),
    ]
    
    all_stopped = True
    for url, name in services:
        if verify_service_stopped(url, name):
            print(f"  [OK] {name} 已停止")
        else:
            print(f"  [WARNING] {name} 可能仍在运行")
            all_stopped = False
    
    # 检查端口占用
    print("\n检查端口占用...")
    ports = [8500, 8080, 8501]
    for port in ports:
        processes = check_port(port)
        if processes:
            print(f"  [WARNING] 端口 {port} 仍被占用，进程: {', '.join(processes)}")
            all_stopped = False
        else:
            print(f"  [OK] 端口 {port} 已释放")
    
    if all_stopped:
        print("\n" + "=" * 60)
        print("所有服务已成功停止！")
        print("=" * 60)
        return 0
    else:
        print("\n" + "=" * 60)
        print("警告：部分服务可能仍在运行")
        print("=" * 60)
        return 1


if __name__ == "__main__":
    sys.exit(main())




