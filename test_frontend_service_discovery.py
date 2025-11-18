"""测试前端服务发现功能"""
import sys
from pathlib import Path

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

# 添加项目路径
project_root = Path(__file__).parent.resolve()
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main"))

from FLASH_DOCK_main.services.registry_client import FlashDockRegistryClient
from FLASH_DOCK_main.services.compass_client import CompassClient

def main():
    print("=" * 60)
    print("测试前端服务发现功能")
    print("=" * 60)
    
    # 测试注册中心客户端
    print("\n1. 测试注册中心客户端...")
    registry_client = FlashDockRegistryClient(registry_url="http://localhost:8500", timeout=5.0)
    
    # 检查注册中心是否可用
    if not registry_client.check_registry_available():
        print("[ERROR] 服务注册中心不可用")
        return 1
    print("[OK] 服务注册中心可用")
    
    # 发现COMPASS服务
    print("\n2. 发现COMPASS服务...")
    services = registry_client.discover_compass_services(healthy_only=True)
    print(f"找到 {len(services)} 个健康的COMPASS服务:")
    for service in services:
        print(f"  - {service.service_id}: {service.host}:{service.port} (状态: {service.status.value})")
    
    if not services:
        print("[ERROR] 没有找到COMPASS服务")
        return 1
    
    # 测试COMPASS客户端
    print("\n3. 测试COMPASS客户端...")
    compass_client = CompassClient(registry_url="http://localhost:8500", timeout=5.0)
    
    # 尝试获取服务URL
    try:
        service_url = compass_client._get_service_url()
        print(f"[OK] 成功获取服务URL: {service_url}")
    except Exception as e:
        print(f"[ERROR] 获取服务URL失败: {e}")
        return 1
    
    # 测试API调用
    print("\n4. 测试API调用...")
    try:
        datasets = compass_client.list_datasets()
        print(f"[OK] 成功获取数据集列表: {len(datasets)} 个数据集")
    except Exception as e:
        print(f"[ERROR] 获取数据集列表失败: {e}")
        return 1
    
    print("\n" + "=" * 60)
    print("所有测试通过！")
    print("=" * 60)
    return 0

if __name__ == "__main__":
    sys.exit(main())

