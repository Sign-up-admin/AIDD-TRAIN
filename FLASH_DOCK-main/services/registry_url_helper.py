"""
注册中心 URL 辅助函数
在 WSL 中自动检测 Windows 主机 IP 地址
"""
import os
import subprocess
import requests


def get_registry_url():
    """
    获取注册中心 URL，在 WSL 中自动检测 Windows 主机 IP
    
    Returns:
        str: 注册中心 URL
    """
    # 首先检查环境变量
    registry_url = os.getenv("REGISTRY_URL", None)
    if registry_url:
        return registry_url
    
    # 如果在 WSL 中，尝试获取 Windows 主机 IP
    if os.path.exists("/etc/resolv.conf"):
        # 方法1: 从路由表获取默认网关（最可靠）
        try:
            result = subprocess.run(
                ["ip", "route", "show", "default"],
                capture_output=True,
                text=True,
                timeout=2
            )
            if result.returncode == 0:
                for line in result.stdout.split('\n'):
                    if 'default via' in line:
                        parts = line.split()
                        if len(parts) >= 3:
                            gateway_ip = parts[2]
                            # 测试网关 IP 是否可访问注册中心
                            test_url = f"http://{gateway_ip}:8500/health"
                            try:
                                response = requests.get(test_url, timeout=2)
                                if response.status_code == 200:
                                    return f"http://{gateway_ip}:8500"
                            except:
                                pass
        except:
            pass
        
        # 方法2: 从 /etc/resolv.conf 获取 nameserver
        try:
            with open("/etc/resolv.conf", "r") as f:
                for line in f:
                    if line.startswith("nameserver"):
                        windows_ip = line.split()[1]
                        # 测试 Windows IP 是否可访问
                        test_url = f"http://{windows_ip}:8500/health"
                        try:
                            response = requests.get(test_url, timeout=2)
                            if response.status_code == 200:
                                return f"http://{windows_ip}:8500"
                        except:
                            pass
        except:
            pass
    
    # 默认使用 localhost
    return "http://localhost:8500"

