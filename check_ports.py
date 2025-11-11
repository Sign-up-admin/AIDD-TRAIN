"""
Port checker utility for COMPASS service system.
Checks if required ports are available and provides information about processes using them.
"""
import socket
import sys
import subprocess
import platform
from typing import List, Tuple, Optional


def check_port(host: str = "localhost", port: int = 8500) -> Tuple[bool, Optional[str]]:
    """
    Check if a port is available.
    
    Args:
        host: Host address
        port: Port number
        
    Returns:
        Tuple of (is_available, process_info)
    """
    try:
        # Try to bind to the port
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(1)
        result = sock.connect_ex((host, port))
        sock.close()
        
        if result == 0:
            # Port is in use, try to get process info
            process_info = get_process_using_port(port)
            return False, process_info
        else:
            return True, None
    except Exception as e:
        # Return False for availability when check fails (port status unknown)
        return False, f"Error checking port: {e}"


def get_process_using_port(port: int) -> Optional[str]:
    """
    Get information about the process using a port.
    
    Args:
        port: Port number
        
    Returns:
        Process information string or None
    """
    try:
        system = platform.system()
        if system == "Windows":
            # Windows: use netstat
            result = subprocess.run(
                ["netstat", "-ano"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if f":{port}" in line and "LISTENING" in line:
                        parts = line.split()
                        if len(parts) >= 5:
                            pid = parts[-1]
                            # Try to get process name
                            try:
                                tasklist_result = subprocess.run(
                                    ["tasklist", "/FI", f"PID eq {pid}", "/FO", "CSV"],
                                    capture_output=True,
                                    text=True,
                                    timeout=5
                                )
                                if tasklist_result.returncode == 0:
                                    lines = tasklist_result.stdout.split('\n')
                                    if len(lines) > 1:
                                        process_info = lines[1].split(',')[0].strip('"')
                                        return f"PID {pid}: {process_info}"
                            except Exception:
                                pass
                            return f"PID {pid}"
        else:
            # Linux/Mac: use lsof or netstat
            try:
                result = subprocess.run(
                    ["lsof", "-i", f":{port}"],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                if result.returncode == 0 and result.stdout:
                    lines = result.stdout.split('\n')
                    if len(lines) > 1:
                        parts = lines[1].split()
                        if len(parts) >= 2:
                            return f"PID {parts[1]}: {parts[0]}"
            except FileNotFoundError:
                # Try netstat as fallback
                try:
                    result = subprocess.run(
                        ["netstat", "-tulpn"],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    if result.returncode == 0:
                        for line in result.stdout.split('\n'):
                            if f":{port}" in line:
                                return line.strip()
                except Exception:
                    pass
    except Exception as e:
        return f"Error getting process info: {e}"
    
    return "Unknown process"


def check_required_ports() -> List[Tuple[int, bool, Optional[str], str]]:
    """
    Check all required ports for the COMPASS service system.
    
    Returns:
        List of tuples: (port, is_available, process_info, description)
    """
    ports_to_check = [
        (8500, "服务注册中心 (Service Registry)"),
        (8080, "COMPASS服务 (COMPASS Service)"),
    ]
    
    results = []
    for port, description in ports_to_check:
        is_available, process_info = check_port("localhost", port)
        results.append((port, is_available, process_info, description))
    
    return results


def main():
    """Main entry point."""
    import io
    import sys
    
    # Set UTF-8 encoding for Windows
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    
    print("=" * 60)
    print("COMPASS服务系统端口检查工具")
    print("Port Checker for COMPASS Service System")
    print("=" * 60)
    print()
    
    results = check_required_ports()
    
    all_available = True
    for port, is_available, process_info, description in results:
        status = "[OK] 可用 (Available)" if is_available else "[IN USE] 被占用 (In Use)"
        print(f"端口 {port} ({description}): {status}")
        
        if not is_available:
            all_available = False
            if process_info:
                print(f"  占用进程: {process_info}")
            print(f"  建议: 停止占用该端口的进程，或使用其他端口")
            print()
    
    print()
    print("=" * 60)
    
    if all_available:
        print("[OK] 所有端口都可用，可以启动服务")
        return 0
    else:
        print("[ERROR] 部分端口被占用，需要解决端口冲突")
        print()
        print("解决方案:")
        print("1. 停止占用端口的进程")
        print("2. 或修改服务配置使用其他端口")
        print("3. Windows: 使用 'netstat -ano | findstr :PORT' 查看端口占用")
        print("4. Linux/Mac: 使用 'lsof -i :PORT' 查看端口占用")
        return 1


if __name__ == "__main__":
    sys.exit(main())

