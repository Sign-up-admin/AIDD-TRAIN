#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""测试 Windows 主机 IP 地址"""

import requests
import os

def test_ips():
    """测试多个可能的 Windows IP 地址"""
    # 可能的 Windows 主机 IP 地址
    test_ips = [
        "172.30.224.1",  # WSL 默认网关（最常见）
        "192.168.1.150",  # Windows 实际 IP
        "10.255.255.254",  # WSL DNS 服务器
    ]
    
    # 从路由表获取默认网关
    try:
        import subprocess
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
                        if gateway_ip not in test_ips:
                            test_ips.insert(0, gateway_ip)  # 优先测试网关 IP
    except:
        pass
    
    print("=" * 60)
    print("测试 Windows 主机 IP 地址")
    print("=" * 60)
    print()
    
    working_ip = None
    
    for ip in test_ips:
        print(f"测试 {ip}...")
        try:
            url = f"http://{ip}:8500/health"
            response = requests.get(url, timeout=3)
            if response.status_code == 200:
                print(f"  [OK] 可以访问注册中心!")
                print(f"  URL: {url}")
                working_ip = ip
                break
            else:
                print(f"  [FAIL] HTTP {response.status_code}")
        except Exception as e:
            print(f"  [FAIL] {e}")
        print()
    
    if working_ip:
        print("=" * 60)
        print(f"[SUCCESS] 找到可用的 Windows 主机 IP: {working_ip}")
        print(f"注册中心 URL: http://{working_ip}:8500")
        print("=" * 60)
        return working_ip
    else:
        print("=" * 60)
        print("[ERROR] 无法找到可用的 Windows 主机 IP")
        print("=" * 60)
        return None

if __name__ == "__main__":
    test_ips()

