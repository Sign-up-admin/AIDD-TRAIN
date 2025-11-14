"""检查COMPASS服务状态"""
import requests
import sys

try:
    response = requests.get('http://localhost:8080/health', timeout=5)
    if response.status_code == 200:
        print("COMPASS服务状态: 运行中 (HTTP 200)")
        print(f"服务地址: http://localhost:8080")
        sys.exit(0)
    else:
        print(f"COMPASS服务状态: 异常 (HTTP {response.status_code})")
        sys.exit(1)
except requests.exceptions.ConnectionError:
    print("COMPASS服务状态: 未运行 (连接被拒绝)")
    print("请检查服务是否已启动")
    sys.exit(1)
except Exception as e:
    print(f"检查服务时出错: {e}")
    sys.exit(1)




