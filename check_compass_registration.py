#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""检查 COMPASS 服务注册状态"""

import requests
import sys

def check_registry_services():
    """检查注册中心中的服务"""
    try:
        url = "http://localhost:8500/api/v1/services"
        response = requests.get(url, params={"service_name": "compass-service"}, timeout=5)
        if response.status_code == 200:
            data = response.json()
            services = data.get("services", [])
            count = data.get("count", 0)
            print(f"注册中心中的 COMPASS 服务数量: {count}")
            if services:
                for svc in services:
                    print(f"  - 服务ID: {svc.get('service_id')}")
                    print(f"    地址: {svc.get('host')}:{svc.get('port')}")
                    print(f"    状态: {svc.get('status', 'unknown')}")
            else:
                print("  未找到 COMPASS 服务")
            return count > 0
        else:
            print(f"注册中心响应错误: {response.status_code}")
            return False
    except Exception as e:
        print(f"检查注册中心失败: {e}")
        return False

def check_compass_health():
    """检查 COMPASS 服务健康状态"""
    try:
        url = "http://localhost:8080/health"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            print(f"COMPASS 服务健康检查: OK")
            print(f"  状态: {data}")
            return True
        else:
            print(f"COMPASS 服务健康检查失败: {response.status_code}")
            return False
    except Exception as e:
        print(f"COMPASS 服务不可访问: {e}")
        return False

def check_registry_health():
    """检查注册中心健康状态"""
    try:
        url = "http://localhost:8500/health"
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            print(f"服务注册中心: OK")
            return True
        else:
            print(f"服务注册中心响应错误: {response.status_code}")
            return False
    except Exception as e:
        print(f"服务注册中心不可访问: {e}")
        return False

if __name__ == "__main__":
    print("=" * 60)
    print("COMPASS 服务注册状态检查")
    print("=" * 60)
    print()
    
    # 检查注册中心
    print("1. 检查服务注册中心...")
    registry_ok = check_registry_health()
    print()
    
    if not registry_ok:
        print("错误: 服务注册中心未运行")
        sys.exit(1)
    
    # 检查 COMPASS 服务健康
    print("2. 检查 COMPASS 服务健康状态...")
    compass_health_ok = check_compass_health()
    print()
    
    # 检查注册中心中的服务
    print("3. 检查注册中心中的 COMPASS 服务...")
    registered = check_registry_services()
    print()
    
    # 总结
    print("=" * 60)
    print("诊断结果:")
    print("=" * 60)
    if compass_health_ok and registered:
        print("[OK] COMPASS 服务正常运行并已注册")
    elif compass_health_ok and not registered:
        print("[WARNING] COMPASS 服务运行正常但未注册到注册中心")
        print("  可能原因: 注册失败或服务刚启动")
    elif not compass_health_ok and registered:
        print("[WARNING] COMPASS 服务已注册但健康检查失败")
        print("  可能原因: 服务启动失败或崩溃")
    else:
        print("[ERROR] COMPASS 服务未运行且未注册")
        print("  需要重新启动 COMPASS 服务")

