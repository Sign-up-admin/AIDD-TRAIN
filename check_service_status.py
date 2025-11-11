"""检查服务状态"""
import sys
import io
import requests

# 设置UTF-8编码
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

def check_service(url, name):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=3)
        return r.status_code == 200
    except:
        return False

print("=" * 60)
print("服务状态检查")
print("=" * 60)
print()

registry_ok = check_service("http://localhost:8500/health", "服务注册中心")
compass_ok = check_service("http://localhost:8080/health", "COMPASS服务")

print(f"服务注册中心 (8500): {'[OK] 运行中' if registry_ok else '[FAIL] 未运行'}")
print(f"COMPASS服务 (8080): {'[OK] 运行中' if compass_ok else '[FAIL] 未运行'}")
print()

if registry_ok and compass_ok:
    print("[OK] 所有服务都在运行!")
    print()
    print("服务地址:")
    print("  - 服务注册中心: http://localhost:8500")
    print("  - COMPASS服务: http://localhost:8080")
    print("  - API文档: http://localhost:8080/docs")
    sys.exit(0)
else:
    print("[WARNING] 部分服务未运行")
    if not registry_ok:
        print("  - 服务注册中心未运行")
    if not compass_ok:
        print("  - COMPASS服务未运行")
    sys.exit(1)












