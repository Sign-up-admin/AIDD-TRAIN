"""检查服务注册中心中的服务列表"""
import requests
import json
import sys

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def main():
    registry_url = "http://localhost:8500"
    
    print("=" * 60)
    print("检查服务注册中心")
    print("=" * 60)
    
    # 检查注册中心健康状态
    try:
        health_response = requests.get(f"{registry_url}/health", timeout=5)
        if health_response.status_code == 200:
            print("[OK] 服务注册中心运行正常")
        else:
            print(f"[WARNING] 服务注册中心健康检查失败: {health_response.status_code}")
    except Exception as e:
        print(f"[ERROR] 无法连接到服务注册中心: {e}")
        return 1
    
    # 获取所有服务列表
    try:
        services_response = requests.get(f"{registry_url}/api/v1/services", timeout=5)
        if services_response.status_code == 200:
            services_data = services_response.json()
            print(f"\n服务列表 (总数: {services_data.get('total', 0)}):")
            print("-" * 60)
            
            services = services_data.get("services", [])
            if not services:
                print("  [WARNING] 没有注册的服务")
            else:
                for service in services:
                    print(f"\n服务名称: {service.get('name', 'N/A')}")
                    print(f"  服务ID: {service.get('service_id', 'N/A')}")
                    print(f"  地址: {service.get('host', 'N/A')}:{service.get('port', 'N/A')}")
                    print(f"  状态: {service.get('status', 'N/A')}")
                    print(f"  版本: {service.get('version', 'N/A')}")
                    if service.get('metadata'):
                        print(f"  元数据: {json.dumps(service['metadata'], indent=4, ensure_ascii=False)}")
            
            # 特别检查COMPASS服务（使用service_name参数查询）
            print("\n" + "=" * 60)
            print("检查COMPASS服务（使用service_name参数）")
            print("=" * 60)
            
            try:
                compass_response = requests.get(
                    f"{registry_url}/api/v1/services",
                    params={"service_name": "compass-service"},
                    timeout=5
                )
                if compass_response.status_code == 200:
                    compass_data = compass_response.json()
                    compass_services = compass_data.get("services", [])
                    if compass_services:
                        print(f"[OK] 找到 {len(compass_services)} 个COMPASS服务:")
                        for service in compass_services:
                            print(f"  - {service.get('service_id')}: {service.get('host')}:{service.get('port')} (状态: {service.get('status')})")
                    else:
                        print("[ERROR] 没有找到COMPASS服务（使用service_name参数查询）")
                        return 1
                else:
                    print(f"[ERROR] 查询COMPASS服务失败: {compass_response.status_code}")
                    return 1
            except Exception as e:
                print(f"[ERROR] 查询COMPASS服务时出错: {e}")
                return 1
            
            # 检查健康状态过滤
            print("\n检查健康的COMPASS服务:")
            try:
                healthy_response = requests.get(
                    f"{registry_url}/api/v1/services",
                    params={"service_name": "compass-service", "status_filter": "healthy"},
                    timeout=5
                )
                if healthy_response.status_code == 200:
                    healthy_data = healthy_response.json()
                    healthy_services = healthy_data.get("services", [])
                    if healthy_services:
                        print(f"[OK] 找到 {len(healthy_services)} 个健康的COMPASS服务")
                    else:
                        print("[WARNING] 没有找到健康的COMPASS服务")
                        return 1
            except Exception as e:
                print(f"[ERROR] 查询健康服务时出错: {e}")
                return 1
        else:
            print(f"[ERROR] 获取服务列表失败: {services_response.status_code}")
            print(f"响应: {services_response.text}")
            return 1
    except Exception as e:
        print(f"[ERROR] 获取服务列表时出错: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print("\n" + "=" * 60)
    print("检查完成")
    print("=" * 60)
    return 0

if __name__ == "__main__":
    sys.exit(main())

