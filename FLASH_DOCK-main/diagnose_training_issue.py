#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
诊断训练管理页面的问题
检查注册中心连接、服务发现、数据集列表加载等
"""

import sys
import os
from pathlib import Path

# Add parent directory to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))
# Add FLASH_DOCK-main/services to path
flashdock_services = Path(__file__).parent / "services"
sys.path.insert(0, str(flashdock_services))

from registry_url_helper import get_registry_url
from compass_client import CompassClient
from service_manager import ServiceManager
import requests

def test_registry_connection(registry_url):
    """测试注册中心连接"""
    print(f"\n1. 测试注册中心连接: {registry_url}")
    try:
        response = requests.get(f"{registry_url}/health", timeout=5)
        if response.status_code == 200:
            print(f"   [OK] 注册中心可访问")
            return True
        else:
            print(f"   [FAIL] HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"   [FAIL] {e}")
        return False

def test_service_discovery(registry_url):
    """测试服务发现"""
    print(f"\n2. 测试服务发现: {registry_url}")
    try:
        response = requests.get(
            f"{registry_url}/api/v1/services",
            params={"service_name": "compass-service"},
            timeout=5
        )
        if response.status_code == 200:
            data = response.json()
            count = data.get("count", 0)
            services = data.get("services", [])
            print(f"   [OK] 找到 {count} 个 COMPASS 服务")
            if services:
                for svc in services:
                    print(f"      - {svc.get('service_id')}: {svc.get('host')}:{svc.get('port')} (状态: {svc.get('status')})")
            return count > 0, services
        else:
            print(f"   [FAIL] HTTP {response.status_code}")
            return False, []
    except Exception as e:
        print(f"   [FAIL] {e}")
        return False, []

def test_compass_client(registry_url):
    """测试 CompassClient"""
    print(f"\n3. 测试 CompassClient 初始化: {registry_url}")
    try:
        client = CompassClient(registry_url=registry_url)
        print(f"   [OK] CompassClient 初始化成功")
        
        # 测试获取服务 URL
        try:
            service_url = client._get_service_url()
            print(f"   [OK] 获取到服务 URL: {service_url}")
            return True, client
        except Exception as e:
            print(f"   [FAIL] 无法获取服务 URL: {e}")
            return False, None
    except Exception as e:
        print(f"   [FAIL] CompassClient 初始化失败: {e}")
        return False, None

def test_list_datasets(client):
    """测试数据集列表加载"""
    print(f"\n4. 测试数据集列表加载")
    try:
        datasets = client.list_datasets()
        print(f"   [OK] 成功加载数据集列表，共 {len(datasets)} 个数据集")
        if datasets:
            for ds in datasets[:3]:  # 只显示前3个
                print(f"      - {ds.get('dataset_id')}: {ds.get('name')}")
        return True
    except Exception as e:
        print(f"   [FAIL] 加载数据集列表失败: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    print("=" * 70)
    print("训练管理页面问题诊断")
    print("=" * 70)
    
    # 获取注册中心 URL
    registry_url = get_registry_url()
    print(f"\n检测到的注册中心 URL: {registry_url}")
    
    # 测试1: 注册中心连接
    if not test_registry_connection(registry_url):
        print("\n[ERROR] 注册中心不可访问，无法继续测试")
        return 1
    
    # 测试2: 服务发现
    has_services, services = test_service_discovery(registry_url)
    if not has_services:
        print("\n[ERROR] 未找到 COMPASS 服务，无法继续测试")
        print("\n可能的原因:")
        print("  1. COMPASS 服务未启动")
        print("  2. COMPASS 服务未注册到注册中心")
        print("  3. 服务状态不是 'healthy'")
        return 1
    
    # 检查服务状态
    healthy_services = [s for s in services if s.get('status') == 'healthy']
    if not healthy_services:
        print("\n[WARNING] 没有健康状态的 COMPASS 服务")
        print("  这可能是因为 CompassClient 使用 healthy_only=True")
        print("  建议: 等待服务完全启动，或检查服务健康状态")
    
    # 测试3: CompassClient
    client_ok, client = test_compass_client(registry_url)
    if not client_ok or not client:
        print("\n[ERROR] CompassClient 无法获取服务 URL")
        return 1
    
    # 测试4: 数据集列表
    if not test_list_datasets(client):
        print("\n[ERROR] 数据集列表加载失败")
        return 1
    
    print("\n" + "=" * 70)
    print("[SUCCESS] 所有测试通过！")
    print("=" * 70)
    return 0

if __name__ == "__main__":
    sys.exit(main())

