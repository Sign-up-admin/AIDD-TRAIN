"""检查FLASH-DOCK服务状态"""
import requests
import sys

try:
    r = requests.get('http://localhost:8501', timeout=2)
    if r.status_code < 500:
        print("FLASH-DOCK运行中")
        print(f"状态码: {r.status_code}")
        sys.exit(0)
    else:
        print("FLASH-DOCK响应异常")
        print(f"状态码: {r.status_code}")
        sys.exit(1)
except requests.exceptions.ConnectionError:
    print("FLASH-DOCK未运行或还在启动中")
    sys.exit(1)
except Exception as e:
    print(f"检查失败: {e}")
    sys.exit(1)

