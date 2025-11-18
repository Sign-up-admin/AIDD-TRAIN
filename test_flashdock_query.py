#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""测试 FlashDock 的查询方式"""

import requests
import json

def test_flashdock_query():
    """模拟 FlashDock 的查询"""
    registry_url = "http://localhost:8500"
    
    # FlashDock 使用的查询方式
    url = f"{registry_url}/api/v1/services"
    params = {"service_name": "compass-service"}
    
    print("=" * 60)
    print("测试 FlashDock 查询方式")
    print("=" * 60)
    print(f"URL: {url}")
    print(f"参数: {params}")
    print()
    
    try:
        response = requests.get(url, params=params, timeout=5)
        print(f"状态码: {response.status_code}")
        print()
        
        if response.status_code == 200:
            data = response.json()
            print("响应数据:")
            print(json.dumps(data, indent=2, ensure_ascii=False))
            print()
            
            services = data.get("services", [])
            count = data.get("count", 0)
            
            print(f"服务数量: {count}")
            print(f"服务列表长度: {len(services)}")
            
            if count > 0:
                print("\n[OK] 找到 COMPASS 服务")
                for i, svc in enumerate(services, 1):
                    print(f"\n服务 {i}:")
                    print(f"  service_id: {svc.get('service_id')}")
                    print(f"  service_name: {svc.get('service_name')}")
                    print(f"  host: {svc.get('host')}")
                    print(f"  port: {svc.get('port')}")
                    print(f"  status: {svc.get('status')}")
            else:
                print("\n[ERROR] 未找到 COMPASS 服务")
        else:
            print(f"[ERROR] 请求失败: {response.status_code}")
            print(f"响应内容: {response.text}")
            
    except Exception as e:
        print(f"[ERROR] 查询失败: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_flashdock_query()

