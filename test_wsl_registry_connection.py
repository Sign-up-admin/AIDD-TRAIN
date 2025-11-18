#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""测试 WSL 中访问 Windows 注册中心"""

import subprocess
import requests
import sys

def get_windows_host_ip():
    """获取 Windows 主机在 WSL 中的 IP 地址"""
    try:
        # 方法1: 从 /etc/resolv.conf 获取
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", "cat /etc/resolv.conf | grep nameserver | head -1 | awk '{print $2}'"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except:
        pass
    
    try:
        # 方法2: 从 ip route 获取
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", "ip route show | grep default | head -1 | awk '{print $3}'"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0 and result.stdout.strip():
            return result.stdout.strip()
    except:
        pass
    
    return None

def test_registry_connection(host):
    """测试注册中心连接"""
    urls = [
        f"http://{host}:8500/health",
        "http://localhost:8500/health",
    ]
    
    for url in urls:
        try:
            response = requests.get(url, timeout=3)
            if response.status_code == 200:
                print(f"[OK] 可以访问: {url}")
                return url.replace("/health", "")
            else:
                print(f"[FAIL] {url}: HTTP {response.status_code}")
        except Exception as e:
            print(f"[FAIL] {url}: {e}")
    
    return None

if __name__ == "__main__":
    print("=" * 60)
    print("测试 WSL 中访问 Windows 注册中心")
    print("=" * 60)
    print()
    
    # 获取 Windows 主机 IP
    print("1. 获取 Windows 主机 IP...")
    windows_ip = get_windows_host_ip()
    if windows_ip:
        print(f"   Windows 主机 IP: {windows_ip}")
    else:
        print("   无法获取 Windows 主机 IP")
        sys.exit(1)
    
    print()
    print("2. 测试注册中心连接...")
    registry_url = test_registry_connection(windows_ip)
    
    print()
    print("=" * 60)
    if registry_url:
        print(f"[SUCCESS] 找到可用的注册中心 URL: {registry_url}")
        print()
        print("解决方案:")
        print(f"  在 FlashDock 中配置 registry_url = '{registry_url}'")
    else:
        print("[ERROR] 无法连接到注册中心")
        print()
        print("可能的原因:")
        print("  1. 注册中心未启动")
        print("  2. Windows 防火墙阻止了连接")
        print("  3. 网络配置问题")

