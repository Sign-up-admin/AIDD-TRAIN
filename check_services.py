#!/usr/bin/env python3
"""
检查服务状态
"""

import requests
import time

def check_service(url, name, timeout=3):
    try:
        r = requests.get(url, timeout=timeout)
        return r.status_code == 200, r.status_code
    except requests.exceptions.RequestException as e:
        return False, str(e)

if __name__ == "__main__":
    print('=' * 50)
    print('检查当前服务状态')
    print('=' * 50)

    registry_ok, registry_status = check_service('http://localhost:8500/health', '服务注册中心')
    compass_ok, compass_status = check_service('http://localhost:8080/health', 'COMPASS服务')
    flashdock_ok, flashdock_status = check_service('http://localhost:8501/_stcore/health', 'FLASH-DOCK服务')

    print(f'服务注册中心 (8500): {"[OK] 运行中" if registry_ok else "[FAIL] " + str(registry_status)}')
    print(f'COMPASS服务 (8080): {"[OK] 运行中" if compass_ok else "[FAIL] " + str(compass_status)}')
    print(f'FLASH-DOCK服务 (8501): {"[OK] 运行中" if flashdock_ok else "[FAIL] " + str(flashdock_status)}')

    print()
    if not compass_ok:
        print('COMPASS服务未运行，可能启动失败')
        print('需要重新启动COMPASS服务')
    else:
        print('所有服务运行正常')
