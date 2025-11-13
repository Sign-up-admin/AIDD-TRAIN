"""检查服务状态"""
import requests
import time
import sys

def check_service(url, name):
    """检查单个服务状态"""
    try:
        response = requests.get(url, timeout=2)
        if response.status_code == 200:
            return True, "运行中"
        else:
            return False, f"状态码: {response.status_code}"
    except requests.exceptions.ConnectionError:
        return False, "未运行"
    except Exception as e:
        return False, f"错误: {str(e)}"

def main():
    """主函数"""
    print("检查服务状态...")
    print("=" * 60)
    
    services = [
        ("http://localhost:8500/health", "服务注册中心"),
        ("http://localhost:8080/health", "COMPASS服务"),
        ("http://localhost:8501/_stcore/health", "FLASH-DOCK前端"),
    ]
    
    all_running = True
    for url, name in services:
        is_running, status = check_service(url, name)
        status_symbol = "[OK]" if is_running else "[FAIL]"
        print(f"{status_symbol} {name}: {status}")
        if not is_running:
            all_running = False
    
    print("=" * 60)
    if all_running:
        print("所有服务运行正常！")
        return 0
    else:
        print("部分服务未运行，请检查服务窗口")
        return 1

if __name__ == "__main__":
    sys.exit(main())
