"""清理重复注册的服务"""
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
    print("清理重复注册的COMPASS服务")
    print("=" * 60)
    
    # 获取所有COMPASS服务
    try:
        response = requests.get(
            f"{registry_url}/api/v1/services",
            params={"service_name": "compass-service"},
            timeout=5
        )
        if response.status_code == 200:
            data = response.json()
            services = data.get("services", [])
            print(f"\n找到 {len(services)} 个COMPASS服务")
            
            if len(services) > 1:
                print("\n检测到多个COMPASS服务，保留最新的一个，删除其他的...")
                
                # 按注册时间排序，保留最新的
                services.sort(key=lambda x: x.get("registered_at", ""), reverse=True)
                keep_service = services[0]
                services_to_delete = services[1:]
                
                print(f"\n保留服务: {keep_service.get('service_id')} (注册时间: {keep_service.get('registered_at')})")
                
                for service in services_to_delete:
                    service_id = service.get("service_id")
                    print(f"\n删除服务: {service_id} (注册时间: {service.get('registered_at')})")
                    
                    try:
                        delete_response = requests.delete(
                            f"{registry_url}/api/v1/services/{service_id}",
                            timeout=5
                        )
                        if delete_response.status_code == 200:
                            print(f"  [OK] 已删除服务 {service_id}")
                        else:
                            print(f"  [ERROR] 删除失败: {delete_response.status_code}")
                    except Exception as e:
                        print(f"  [ERROR] 删除服务时出错: {e}")
                
                print("\n" + "=" * 60)
                print("清理完成")
                print("=" * 60)
            else:
                print("\n只有一个COMPASS服务，无需清理")
        else:
            print(f"[ERROR] 获取服务列表失败: {response.status_code}")
            return 1
    except Exception as e:
        print(f"[ERROR] 清理服务时出错: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

