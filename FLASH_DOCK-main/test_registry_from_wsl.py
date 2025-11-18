#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""在 WSL 中测试访问 Windows 注册中心"""

import requests
import os

def get_windows_host_ip():
    """获取 Windows 主机 IP"""
    try:
        with open('/etc/resolv.conf', 'r') as f:
            for line in f:
                if line.startswith('nameserver'):
                    return line.split()[1]
    except:
        pass
    return None

def test_registry():
    """测试注册中心连接"""
    windows_ip = get_windows_host_ip()
    print(f"Windows 主机 IP: {windows_ip}")
    print()
    
    # 测试不同的 URL
    test_urls = []
    if windows_ip:
        test_urls.append(f"http://{windows_ip}:8500")
    test_urls.append("http://localhost:8500")
    
    for base_url in test_urls:
        print(f"测试: {base_url}")
        try:
            # 测试健康检查
            health_url = f"{base_url}/health"
            response = requests.get(health_url, timeout=3)
            if response.status_code == 200:
                print(f"  [OK] 健康检查成功")
                
                # 测试服务查询
                services_url = f"{base_url}/api/v1/services"
                params = {"service_name": "compass-service"}
                response = requests.get(services_url, params=params, timeout=3)
                if response.status_code == 200:
                    data = response.json()
                    count = data.get("count", 0)
                    print(f"  [OK] 找到 {count} 个 COMPASS 服务")
                    if count > 0:
                        print(f"  [SUCCESS] 可用的注册中心 URL: {base_url}")
                        return base_url
                else:
                    print(f"  [FAIL] 服务查询失败: {response.status_code}")
            else:
                print(f"  [FAIL] 健康检查失败: {response.status_code}")
        except Exception as e:
            print(f"  [FAIL] 连接失败: {e}")
        print()
    
    return None

if __name__ == "__main__":
    print("=" * 60)
    print("在 WSL 中测试访问 Windows 注册中心")
    print("=" * 60)
    print()
    
    registry_url = test_registry()
    
    if registry_url:
        print("=" * 60)
        print(f"建议的注册中心 URL: {registry_url}")
        print("=" * 60)
        print()
        print("如果 FlashDock 使用 localhost:8500 无法连接，")
        print(f"请尝试使用: {registry_url}")

