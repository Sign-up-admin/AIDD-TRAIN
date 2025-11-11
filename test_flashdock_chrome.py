"""
FlashDock前端Chrome测试辅助脚本
用于自动化测试和验证FlashDock前端功能
"""

import requests
import time
import json
from pathlib import Path
from typing import Dict, List, Optional

# 服务地址
REGISTRY_URL = "http://localhost:8500"
COMPASS_URL = "http://localhost:8080"
FLASHDOCK_URL = "http://localhost:8501"


def check_service(url: str, name: str) -> bool:
    """检查服务是否运行"""
    try:
        response = requests.get(url, timeout=5)
        return response.status_code == 200
    except:
        return False


def test_homepage():
    """测试主页功能"""
    print("\n=== 测试主页功能 ===")
    try:
        response = requests.get(FLASHDOCK_URL, timeout=10)
        if response.status_code == 200:
            print("[OK] 主页加载成功")
            # 检查页面内容
            if "FlashDock" in response.text or "streamlit" in response.text.lower():
                print("[OK] 页面内容正常")
            return True
        else:
            print(f"[FAIL] 主页加载失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] 主页测试失败: {e}")
        return False


def test_compass_api():
    """测试COMPASS API连接"""
    print("\n=== 测试COMPASS API ===")
    try:
        # 测试健康检查
        health_url = f"{COMPASS_URL}/health"
        response = requests.get(health_url, timeout=5)
        if response.status_code == 200:
            print("[OK] COMPASS服务健康检查通过")
        else:
            print(f"[FAIL] COMPASS健康检查失败: HTTP {response.status_code}")
            return False

        # 测试任务列表API
        tasks_url = f"{COMPASS_URL}/api/v1/training/tasks"
        response = requests.get(tasks_url, timeout=5)
        if response.status_code == 200:
            print("[OK] 任务列表API正常")
            tasks = response.json().get("tasks", [])
            print(f"  当前任务数: {len(tasks)}")
            return True
        else:
            print(f"[FAIL] 任务列表API失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] COMPASS API测试失败: {e}")
        return False


def test_registry_api():
    """测试服务注册中心API"""
    print("\n=== 测试服务注册中心 ===")
    try:
        health_url = f"{REGISTRY_URL}/health"
        response = requests.get(health_url, timeout=5)
        if response.status_code == 200:
            print("[OK] 服务注册中心健康检查通过")
            return True
        else:
            print(f"[FAIL] 服务注册中心健康检查失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"[FAIL] 服务注册中心测试失败: {e}")
        return False


def create_test_task() -> Optional[str]:
    """创建测试任务"""
    print("\n=== 创建测试任务 ===")
    try:
        task_url = f"{COMPASS_URL}/api/v1/training/tasks"
        task_config = {
            "config": {
                "execution_mode": "smoke_test",
                "epochs": 10,
                "batch_size": 2,
                "learning_rate": 0.001,
                "optimizer": "adam"
            },
            "description": "Chrome测试任务"
        }
        response = requests.post(task_url, json=task_config, timeout=10)
        if response.status_code == 201:
            task = response.json()
            task_id = task.get("task_id")
            print(f"[OK] 测试任务创建成功: {task_id}")
            return task_id
        else:
            print(f"[FAIL] 创建任务失败: HTTP {response.status_code}")
            print(f"  响应: {response.text}")
            return None
    except Exception as e:
        print(f"[FAIL] 创建任务失败: {e}")
        return None


def test_task_operations(task_id: str):
    """测试任务操作（启动、停止）"""
    print(f"\n=== 测试任务操作 (任务ID: {task_id}) ===")
    
    # 启动任务
    print("\n1. 测试启动任务...")
    try:
        start_url = f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start"
        response = requests.post(start_url, timeout=10)
        if response.status_code == 200:
            print("[OK] 任务启动请求成功")
            # 等待任务状态更新
            time.sleep(2)
            # 检查任务状态
            task_url = f"{COMPASS_URL}/api/v1/training/tasks/{task_id}"
            response = requests.get(task_url, timeout=5)
            if response.status_code == 200:
                task = response.json()
                status = task.get("status")
                print(f"  任务状态: {status}")
                if status in ["running", "initializing"]:
                    print("[OK] 任务已启动")
                    
                    # 测试停止任务
                    print("\n2. 测试停止任务...")
                    stop_url = f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop"
                    response = requests.post(stop_url, timeout=10)
                    if response.status_code == 200:
                        print("[OK] 任务停止请求成功")
                        time.sleep(2)
                        # 检查任务状态
                        response = requests.get(task_url, timeout=5)
                        if response.status_code == 200:
                            task = response.json()
                            status = task.get("status")
                            print(f"  任务状态: {status}")
                            if status in ["cancelled", "completed", "failed"]:
                                print("[OK] 任务已停止")
                                return True
                            else:
                                print(f"[WARN] 任务状态: {status} (可能仍在停止中)")
                                return True
                        else:
                            print("[FAIL] 无法获取任务状态")
                            return False
                    else:
                        print(f"[FAIL] 停止任务失败: HTTP {response.status_code}")
                        print(f"  响应: {response.text}")
                        return False
                else:
                    print(f"[WARN] 任务状态: {status} (可能仍在启动中)")
                    return True
            else:
                print("[FAIL] 无法获取任务状态")
                return False
        else:
            print(f"[FAIL] 启动任务失败: HTTP {response.status_code}")
            print(f"  响应: {response.text}")
            return False
    except Exception as e:
        print(f"[FAIL] 任务操作测试失败: {e}")
        return False


def main():
    """主测试函数"""
    print("=" * 70)
    print("FlashDock前端Chrome测试辅助脚本")
    print("=" * 70)
    
    # 检查服务状态
    print("\n检查服务状态...")
    registry_ok = check_service(REGISTRY_URL, "服务注册中心")
    compass_ok = check_service(f"{COMPASS_URL}/health", "COMPASS服务")
    flashdock_ok = check_service(FLASHDOCK_URL, "FLASH-DOCK前端")
    
    print(f"服务注册中心 (8500): {'[OK]' if registry_ok else '[FAIL]'}")
    print(f"COMPASS服务 (8080): {'[OK]' if compass_ok else '[FAIL]'}")
    print(f"FLASH-DOCK前端 (8501): {'[OK]' if flashdock_ok else '[FAIL]'}")
    
    if not (registry_ok and compass_ok and flashdock_ok):
        print("\n[WARN] 警告: 部分服务未运行，请先启动所有服务")
        print("可以使用以下命令启动服务:")
        print("  python check_and_start_services.py")
        print("  或")
        print("  start_all_services.bat")
        return
    
    # 运行测试
    results = []
    
    # 1. 测试主页
    results.append(("主页功能", test_homepage()))
    
    # 2. 测试API
    results.append(("服务注册中心", test_registry_api()))
    results.append(("COMPASS API", test_compass_api()))
    
    # 3. 测试任务操作
    task_id = create_test_task()
    if task_id:
        results.append(("任务操作", test_task_operations(task_id)))
    
    # 输出测试结果
    print("\n" + "=" * 70)
    print("测试结果汇总")
    print("=" * 70)
    for name, result in results:
        status = "[PASS] 通过" if result else "[FAIL] 失败"
        print(f"{name}: {status}")
    
    print("\n" + "=" * 70)
    print("Chrome浏览器测试指南")
    print("=" * 70)
    print("1. 打开Chrome浏览器，访问:", FLASHDOCK_URL)
    print("2. 按F12打开开发者工具")
    print("3. Network标签页: 勾选'Preserve log'，清空请求记录")
    print("4. Console标签页: 清空日志，准备观察错误")
    print("5. 测试以下功能:")
    print("   - 页面导航切换")
    print("   - 创建训练任务")
    print("   - 启动任务")
    print("   - 停止任务（重点测试）")
    print("   - WebSocket终端显示")
    print("   - 任务状态更新")
    print("6. 观察Network和Console中的错误信息")
    print("\n详细测试指南请参考: CHROME_DEBUG_GUIDE.md")


if __name__ == "__main__":
    main()

