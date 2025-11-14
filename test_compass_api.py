"""测试COMPASS API"""
import requests

try:
    response = requests.get('http://localhost:8080/api/v1/training/tasks', timeout=5)
    print(f"API测试: HTTP {response.status_code}")
    if response.status_code == 200:
        data = response.json()
        tasks = data.get("tasks", [])
        print(f"任务数量: {len(tasks)}")
        print("COMPASS服务API正常工作！")
    else:
        print(f"API返回异常状态码: {response.status_code}")
except Exception as e:
    print(f"API测试失败: {e}")




