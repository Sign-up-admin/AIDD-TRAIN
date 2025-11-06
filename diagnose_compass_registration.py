"""
诊断脚本来检查COMPASS服务注册问题
"""
import sys
import requests
from pathlib import Path
import io

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# Add parent directory to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

def check_registry():
    """检查服务注册中心"""
    try:
        response = requests.get("http://localhost:8500/health", timeout=5)
        if response.status_code == 200:
            print("[OK] 服务注册中心运行正常")
            return True
        else:
            print(f"[FAIL] 服务注册中心响应异常: {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] 无法连接到服务注册中心: {e}")
        return False

def check_compass_service():
    """检查COMPASS服务"""
    try:
        response = requests.get("http://localhost:8080/health", timeout=5)
        if response.status_code == 200:
            print("[OK] COMPASS服务运行正常")
            data = response.json()
            print(f"  状态: {data.get('status', 'N/A')}")
            return True
        else:
            print(f"[FAIL] COMPASS服务响应异常: {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] 无法连接到COMPASS服务: {e}")
        return False

def check_registered_services():
    """检查已注册的服务"""
    try:
        response = requests.get("http://localhost:8500/api/v1/services?service_name=compass-service&status_filter=healthy", timeout=5)
        if response.status_code == 200:
            data = response.json()
            services = data.get('services', [])
            if services:
                print(f"[OK] 找到 {len(services)} 个已注册的COMPASS服务:")
                for svc in services:
                    print(f"  - Service ID: {svc.get('service_id')}")
                    print(f"    地址: {svc.get('base_url')}")
                    print(f"    状态: {svc.get('status')}")
                return True
            else:
                print("[FAIL] 未找到已注册的COMPASS服务")
                return False
        else:
            print(f"[FAIL] 查询服务列表失败: {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] 查询服务列表时出错: {e}")
        return False

def check_all_services():
    """检查所有已注册的服务"""
    try:
        response = requests.get("http://localhost:8500/api/v1/services", timeout=5)
        if response.status_code == 200:
            data = response.json()
            services = data.get('services', [])
            print(f"\n所有已注册的服务 ({len(services)} 个):")
            for svc in services:
                print(f"  - {svc.get('service_name')} ({svc.get('service_id')})")
                print(f"    地址: {svc.get('base_url')}")
                print(f"    状态: {svc.get('status')}")
            return services
        else:
            print(f"[FAIL] 查询所有服务失败: {response.status_code}")
            return []
    except Exception as e:
        print(f"[FAIL] 查询所有服务时出错: {e}")
        return []

def main():
    print("=" * 60)
    print("COMPASS服务注册诊断工具")
    print("=" * 60)
    print()
    
    # 检查服务注册中心
    print("1. 检查服务注册中心...")
    registry_ok = check_registry()
    print()
    
    # 检查COMPASS服务
    print("2. 检查COMPASS服务...")
    compass_ok = check_compass_service()
    print()
    
    # 检查已注册的服务
    print("3. 检查已注册的COMPASS服务...")
    registered = check_registered_services()
    print()
    
    # 检查所有服务
    print("4. 检查所有已注册的服务...")
    all_services = check_all_services()
    print()
    
    # 诊断结果
    print("=" * 60)
    print("诊断结果:")
    print("=" * 60)
    if registry_ok and compass_ok and not registered:
        print("[ERROR] 问题: COMPASS服务正在运行，但未注册到服务注册中心")
        print("\n可能的原因:")
        print("  1. COMPASS服务启动时注册失败")
        print("  2. 注册后心跳失败，服务被标记为不健康")
        print("  3. 服务名称不匹配")
        print("\n解决方案:")
        print("  1. 查看COMPASS服务的启动日志，查找注册错误")
        print("  2. 重启COMPASS服务")
        print("  3. 检查COMPASS服务的注册配置")
    elif not registry_ok:
        print("[ERROR] 问题: 服务注册中心不可用")
    elif not compass_ok:
        print("[ERROR] 问题: COMPASS服务不可用")
    elif registered:
        print("[OK] 一切正常: COMPASS服务已正确注册")
    print()

if __name__ == "__main__":
    main()

