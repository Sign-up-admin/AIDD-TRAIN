#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""在 WSL 中测试访问 Windows 注册中心"""

import requests
import os

def get_windows_host_ip():
    """获取 Windows 主机 IP"""
    try:
        if os.path.exists("/etc/resolv.conf"):
            with open("/etc/resolv.conf", "r") as f:
                for line in f:
                    if line.startswith("nameserver"):
                        return line.split()[1]
    except:
        pass
    return None

def test_registry_connection():
    """测试注册中心连接"""
    print("=" * 60)
    print("测试 WSL 中访问 Windows 注册中心")
    print("=" * 60)
    print()
    
    # 获取 Windows IP
    windows_ip = get_windows_host_ip()
    print(f"1. Windows 主机 IP: {windows_ip}")
    print()
    
    # 测试不同的 URL
    test_urls = []
    if windows_ip:
        test_urls.append(("Windows IP", f"http://{windows_ip}:8500"))
    test_urls.append(("localhost", "http://localhost:8500"))
    
    print("2. 测试注册中心连接:")
    working_url = None
    
    for name, base_url in test_urls:
        print(f"   测试 {name} ({base_url})...")
        try:
            # 测试健康检查
            health_url = f"{base_url}/health"
            response = requests.get(health_url, timeout=3)
            if response.status_code == 200:
                print(f"   [OK] 健康检查成功")
                
                # 测试服务查询
                services_url = f"{base_url}/api/v1/services"
                params = {"service_name": "compass-service"}
                response = requests.get(services_url, params=params, timeout=3)
                if response.status_code == 200:
                    data = response.json()
                    count = data.get("count", 0)
                    print(f"   [OK] 找到 {count} 个 COMPASS 服务")
                    if count > 0:
                        print(f"   [SUCCESS] 可用的注册中心 URL: {base_url}")
                        working_url = base_url
                        break
                else:
                    print(f"   [FAIL] 服务查询失败: HTTP {response.status_code}")
            else:
                print(f"   [FAIL] 健康检查失败: HTTP {response.status_code}")
        except Exception as e:
            print(f"   [FAIL] 连接失败: {e}")
        print()
    
    print("=" * 60)
    if working_url:
        print(f"[SUCCESS] 找到可用的注册中心: {working_url}")
        print()
        print("建议:")
        print(f"  在 FlashDock 中使用: {working_url}")
        return working_url
    else:
        print("[ERROR] 无法连接到注册中心")
        print()
        print("可能的原因:")
        print("  1. 注册中心未启动")
        print("  2. Windows 防火墙阻止了连接")
        print("  3. 网络配置问题")
        return None

if __name__ == "__main__":
    test_registry_connection()

