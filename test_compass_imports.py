#!/usr/bin/env python3
"""
测试COMPASS服务的关键模块导入
"""

import sys
import os
from pathlib import Path

# 设置项目根目录
project_root = Path(__file__).parent.resolve()
sys.path.insert(0, str(project_root))
os.environ["PYTHONPATH"] = str(project_root)

print("=" * 60)
print("测试COMPASS服务关键模块导入")
print("=" * 60)
print(f"项目根目录: {project_root}")
print(f"Python路径: {sys.path[:3]}")
print()

errors = []

# 测试基础依赖
print("1. 测试基础依赖...")
try:
    import fastapi
    print(f"   [OK] fastapi {fastapi.__version__}")
except ImportError as e:
    print(f"   [FAIL] fastapi: {e}")
    errors.append(f"fastapi: {e}")

try:
    import uvicorn
    print(f"   [OK] uvicorn {uvicorn.__version__}")
except ImportError as e:
    print(f"   [FAIL] uvicorn: {e}")
    errors.append(f"uvicorn: {e}")

try:
    import pydantic
    print(f"   [OK] pydantic {pydantic.__version__}")
except ImportError as e:
    print(f"   [FAIL] pydantic: {e}")
    errors.append(f"pydantic: {e}")

# 测试services.common.utils
print("\n2. 测试services.common.utils...")
try:
    from services.common.utils import get_local_ip
    ip = get_local_ip()
    print(f"   [OK] services.common.utils.get_local_ip -> {ip}")
except ImportError as e:
    print(f"   [FAIL] services.common.utils: {e}")
    errors.append(f"services.common.utils: {e}")

# 测试compass.service.config
print("\n3. 测试compass.service.config...")
try:
    from compass.service.config import SERVICE_CONFIG
    print(f"   [OK] compass.service.config.SERVICE_CONFIG")
    print(f"        host: {SERVICE_CONFIG.get('host')}")
    print(f"        port: {SERVICE_CONFIG.get('port')}")
    print(f"        registry_url: {SERVICE_CONFIG.get('registry_url')}")
except ImportError as e:
    print(f"   [FAIL] compass.service.config: {e}")
    errors.append(f"compass.service.config: {e}")
except Exception as e:
    print(f"   [FAIL] compass.service.config加载错误: {e}")
    errors.append(f"compass.service.config加载: {e}")

# 测试compass.service.server
print("\n4. 测试compass.service.server...")
try:
    from compass.service.server import app, main
    print(f"   [OK] compass.service.server (app, main)")
except ImportError as e:
    print(f"   [FAIL] compass.service.server: {e}")
    errors.append(f"compass.service.server: {e}")
    import traceback
    traceback.print_exc()

# 测试compass.service_main
print("\n5. 测试compass.service_main...")
try:
    import compass.service_main
    print(f"   [OK] compass.service_main")
except ImportError as e:
    print(f"   [FAIL] compass.service_main: {e}")
    errors.append(f"compass.service_main: {e}")
    import traceback
    traceback.print_exc()

# 总结
print("\n" + "=" * 60)
if errors:
    print(f"[FAIL] 发现 {len(errors)} 个导入错误:")
    for error in errors:
        print(f"  - {error}")
    sys.exit(1)
else:
    print("[OK] 所有关键模块导入成功!")
    sys.exit(0)

