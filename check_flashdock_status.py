"""检查FLASH-DOCK服务状态"""
import sys
import io
import requests

# 设置UTF-8编码
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

try:
    r = requests.get("http://localhost:8501", timeout=3)
    if r.status_code == 200:
        print("FLASH-DOCK (8501): [OK] 运行中")
        sys.exit(0)
    else:
        print(f"FLASH-DOCK (8501): [FAIL] 状态码: {r.status_code}")
        sys.exit(1)
except Exception as e:
    print(f"FLASH-DOCK (8501): [FAIL] 未运行 - {e}")
    sys.exit(1)
