"""检查所有服务状态"""
import requests
import sys
import io

# 设置UTF-8编码
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def check_service(url, name, timeout=3):
    """检查服务是否运行"""
    try:
        r = requests.get(url, timeout=timeout)
        return r.status_code == 200 or r.status_code < 500
    except:
        return False

print("=" * 60)
print("所有服务状态检查")
print("=" * 60)
print()

services = [
    ("服务注册中心", "http://localhost:8500/health", 8500),
    ("COMPASS服务", "http://localhost:8080/health", 8080),
    ("FLASH-DOCK", "http://localhost:8501", 8501),
]

all_ok = True
for name, url, port in services:
    status = check_service(url, name)
    status_text = "✅ 运行中" if status else "❌ 未运行"
    print(f"{name} (端口 {port}): {status_text}")
    if not status:
        all_ok = False

print()
print("=" * 60)
if all_ok:
    print("✅ 所有服务都在运行！")
    print()
    print("服务地址:")
    print("  - 服务注册中心: http://localhost:8500")
    print("  - COMPASS服务: http://localhost:8080")
    print("  - FLASH-DOCK: http://localhost:8501")
else:
    print("⚠️  部分服务未运行，请检查服务状态")
print("=" * 60)

