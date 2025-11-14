"""
验证所有服务是否正常运行
"""

import sys
import requests

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def check_service(url, name):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=3)
        if r.status_code == 200:
            try:
                data = r.json()
                return True, data
            except:
                return True, {"status": "ok"}
        else:
            return False, {"error": f"Status code: {r.status_code}"}
    except Exception as e:
        return False, {"error": str(e)}


print("=" * 60)
print("验证服务状态")
print("=" * 60)
print()

# 检查服务注册中心
registry_ok, registry_data = check_service("http://localhost:8500/health", "服务注册中心")
print(f"服务注册中心 (8500): {'[OK]' if registry_ok else '[FAIL]'}")
if registry_ok:
    print(f"  响应: {registry_data}")
else:
    print(f"  错误: {registry_data}")

print()

# 检查COMPASS服务
compass_ok, compass_data = check_service("http://localhost:8080/health", "COMPASS服务")
print(f"COMPASS服务 (8080): {'[OK]' if compass_ok else '[FAIL]'}")
if compass_ok:
    print(f"  响应: {compass_data}")
else:
    print(f"  错误: {compass_data}")

print()

# 检查服务注册中心的服务列表
if registry_ok:
    try:
        r = requests.get("http://localhost:8500/api/v1/services", timeout=3)
        if r.status_code == 200:
            services = r.json()
            print(f"已注册的服务数量: {services.get('count', 0)}")
            if services.get('services'):
                print("已注册的服务:")
                for svc in services['services']:
                    print(f"  - {svc.get('service_name')} ({svc.get('service_id')}) - {svc.get('status')}")
    except Exception as e:
        print(f"无法获取服务列表: {e}")

print()

if registry_ok and compass_ok:
    print("[SUCCESS] 所有服务都正常运行!")
    sys.exit(0)
else:
    print("[WARNING] 部分服务未正常运行")
    sys.exit(1)



