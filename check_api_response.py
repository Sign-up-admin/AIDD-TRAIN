"""检查API返回格式"""
import requests
import json

response = requests.get("http://localhost:8080/api/v1/training/tasks", timeout=10)
print(f"Status Code: {response.status_code}")
print(f"Content-Type: {response.headers.get('content-type')}")
print(f"\nResponse Type: {type(response.json())}")
print(f"\nResponse Content:")
print(json.dumps(response.json(), indent=2, ensure_ascii=False))










