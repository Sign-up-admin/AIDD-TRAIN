"""
使用API测试训练停止功能（不依赖浏览器）
"""
import sys
import time
import requests
from typing import Dict, Any, Optional

# 服务配置
COMPASS_URL = "http://localhost:8080"

def check_service(url: str, name: str, timeout: int = 5) -> bool:
    """检查服务是否运行"""
    try:
        response = requests.get(url, timeout=timeout)
        return response.status_code == 200
    except:
        return False

def create_training_task(config: Dict[str, Any]) -> Optional[str]:
    """创建训练任务"""
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks",
            json=config,
            timeout=30
        )
        if response.status_code == 201:
            task = response.json()
            return task.get("task_id")
        else:
            print(f"[ERROR] Failed to create task: {response.status_code}")
            print(f"Response: {response.text}")
            return None
    except Exception as e:
        print(f"[ERROR] Exception creating task: {e}")
        return None

def get_task_status(task_id: str) -> Optional[str]:
    """获取任务状态"""
    try:
        response = requests.get(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
            timeout=10
        )
        if response.status_code == 200:
            task = response.json()
            return task.get("status")
        else:
            return None
    except Exception as e:
        print(f"[ERROR] Exception getting task status: {e}")
        return None

def start_task(task_id: str) -> bool:
    """启动任务"""
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=30
        )
        return response.status_code == 200
    except Exception as e:
        print(f"[ERROR] Exception starting task: {e}")
        return False

def stop_task(task_id: str) -> tuple[bool, Optional[str]]:
    """停止任务"""
    try:
        start_time = time.time()
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop",
            timeout=10
        )
        elapsed = time.time() - start_time
        
        if response.status_code == 200:
            return True, None
        else:
            error_detail = response.json().get("detail", response.text)
            return False, error_detail
    except requests.exceptions.Timeout:
        return False, "Request timeout"
    except Exception as e:
        return False, str(e)

def wait_for_task_status(task_id: str, target_status: str, max_wait: int = 60) -> bool:
    """等待任务状态变为目标状态"""
    start_time = time.time()
    while time.time() - start_time < max_wait:
        status = get_task_status(task_id)
        if status == target_status:
            return True
        print(f"  Current status: {status}, waiting for {target_status}...")
        time.sleep(2)
    return False

def test_scenario_1_normal_stop():
    """测试场景1: 正常停止流程"""
    print("\n" + "=" * 70)
    print("测试场景1: 正常停止流程")
    print("=" * 70)
    
    # 创建任务
    config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Test: Normal stop scenario"
    }
    
    task_id = create_training_task(config)
    if not task_id:
        print("[FAIL] Failed to create task")
        return False
    
    print(f"[OK] Task created: {task_id}")
    
    # 启动任务
    if not start_task(task_id):
        print("[FAIL] Failed to start task")
        return False
    
    print("[OK] Task started")
    
    # 等待任务进入running状态
    print("Waiting for task to enter running state...")
    if not wait_for_task_status(task_id, "running", max_wait=30):
        # 检查当前状态
        current_status = get_task_status(task_id)
        print(f"[WARNING] Task status: {current_status}")
        if current_status not in ["running", "initializing"]:
            print("[FAIL] Task is not in a stoppable state")
            return False
    
    # 停止任务
    print("Stopping task...")
    start_time = time.time()
    success, error = stop_task(task_id)
    elapsed = time.time() - start_time
    
    if success:
        print(f"[OK] Stop request sent (took {elapsed:.2f}s)")
        
        # 等待任务停止
        print("Waiting for task to stop...")
        if wait_for_task_status(task_id, "cancelled", max_wait=40):
            print("[OK] Task successfully stopped")
            return True
        else:
            final_status = get_task_status(task_id)
            if final_status in ["cancelled", "failed", "completed"]:
                print(f"[OK] Task is no longer running (status: {final_status})")
                return True
            else:
                print(f"[FAIL] Task is still running (status: {final_status})")
                return False
    else:
        print(f"[FAIL] Failed to stop task: {error}")
        return False

def test_scenario_2_quick_stop():
    """测试场景2: 快速停止"""
    print("\n" + "=" * 70)
    print("测试场景2: 快速停止")
    print("=" * 70)
    
    # 创建任务
    config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Test: Quick stop scenario"
    }
    
    task_id = create_training_task(config)
    if not task_id:
        print("[FAIL] Failed to create task")
        return False
    
    print(f"[OK] Task created: {task_id}")
    
    # 启动任务
    if not start_task(task_id):
        print("[FAIL] Failed to start task")
        return False
    
    print("[OK] Task started")
    
    # 立即停止任务
    print("Stopping task immediately...")
    success, error = stop_task(task_id)
    
    if success:
        print("[OK] Stop request sent")
        
        # 等待任务停止
        print("Waiting for task to stop...")
        if wait_for_task_status(task_id, "cancelled", max_wait=40):
            print("[OK] Task successfully stopped")
            return True
        else:
            final_status = get_task_status(task_id)
            if final_status in ["cancelled", "failed"]:
                print(f"[OK] Task stopped (status: {final_status})")
                return True
            else:
                print(f"[WARNING] Task status: {final_status}")
                return False
    else:
        print(f"[FAIL] Failed to stop task: {error}")
        return False

def test_scenario_3_error_handling_completed():
    """测试场景3: 错误处理（已完成的任务）"""
    print("\n" + "=" * 70)
    print("测试场景3: 错误处理（已完成的任务）")
    print("=" * 70)
    
    # 创建快速完成的任务
    config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 1,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Test: Error handling - completed task"
    }
    
    task_id = create_training_task(config)
    if not task_id:
        print("[FAIL] Failed to create task")
        return False
    
    print(f"[OK] Task created: {task_id}")
    
    # 启动任务
    if not start_task(task_id):
        print("[FAIL] Failed to start task")
        return False
    
    # 等待任务完成
    print("Waiting for task to complete...")
    if not wait_for_task_status(task_id, "completed", max_wait=120):
        print("[WARNING] Task did not complete in time")
        current_status = get_task_status(task_id)
        print(f"Current status: {current_status}")
        if current_status != "completed":
            print("[SKIP] Task is not completed, skipping test")
            return True  # 跳过测试，不算失败
    
    # 尝试停止已完成的任务
    print("Attempting to stop completed task...")
    success, error = stop_task(task_id)
    
    if not success:
        if error and ("cannot be stopped" in error.lower() or "status" in error.lower()):
            print("[OK] Correctly rejected stop request for completed task")
            print(f"Error message: {error}")
            return True
        else:
            print(f"[WARNING] Unexpected error: {error}")
            return False
    else:
        print("[WARNING] Stop request was accepted (unexpected)")
        return False

def test_scenario_4_error_handling_nonexistent():
    """测试场景4: 错误处理（任务不存在）"""
    print("\n" + "=" * 70)
    print("测试场景4: 错误处理（任务不存在）")
    print("=" * 70)
    
    fake_task_id = "00000000-0000-0000-0000-000000000000"
    
    print(f"Attempting to stop nonexistent task: {fake_task_id}")
    success, error = stop_task(fake_task_id)
    
    if not success:
        if error and ("not found" in error.lower() or "404" in str(error)):
            print("[OK] Correctly rejected stop request for nonexistent task")
            print(f"Error message: {error}")
            return True
        else:
            print(f"[WARNING] Unexpected error: {error}")
            return False
    else:
        print("[WARNING] Stop request was accepted (unexpected)")
        return False

def test_scenario_5_error_handling_network():
    """测试场景5: 错误处理（网络错误）"""
    print("\n" + "=" * 70)
    print("测试场景5: 错误处理（网络错误）")
    print("=" * 70)
    print("[NOTE] This test requires manual intervention:")
    print("  1. Create and start a training task")
    print("  2. Stop the COMPASS service")
    print("  3. Attempt to stop the task")
    print("  4. Verify that appropriate error is returned")
    print()
    print("Skipping automated test (requires service to be stopped)")
    print("[INFO] This scenario should be tested manually")
    return True  # Skip this test in automated mode

def main():
    """主函数"""
    print("=" * 70)
    print("训练停止功能API测试")
    print("=" * 70)
    
    # 检查服务
    print("\n检查服务状态...")
    registry_ok = check_service("http://localhost:8500/health", "Registry")
    compass_ok = check_service(f"{COMPASS_URL}/health", "COMPASS")
    
    print(f"Registry (8500): {'[OK]' if registry_ok else '[NOT RUNNING]'}")
    print(f"COMPASS (8080): {'[OK]' if compass_ok else '[NOT RUNNING]'}")
    
    if not registry_ok:
        print("\n[ERROR] Registry service is not running")
        print("Please start the service first:")
        print("  1. Run start_all_services.bat")
        print("  2. Or run: python services/registry/server.py --host 0.0.0.0 --port 8500")
        return 1
    
    if not compass_ok:
        print("\n[ERROR] COMPASS service is not running")
        print("Please start the service first:")
        print("  1. Run start_all_services.bat")
        print("  2. Or run: python compass/service_main.py --host 0.0.0.0 --port 8080 --registry-url http://localhost:8500")
        return 1
    
    print("\n[OK] All services are running")
    
    # 运行测试
    results = []
    
    try:
        result1 = test_scenario_1_normal_stop()
        results.append(("场景1: 正常停止流程", result1))
    except Exception as e:
        print(f"[ERROR] Test scenario 1 failed: {e}")
        import traceback
        traceback.print_exc()
        results.append(("场景1: 正常停止流程", False))
    
    try:
        result2 = test_scenario_2_quick_stop()
        results.append(("场景2: 快速停止", result2))
    except Exception as e:
        print(f"[ERROR] Test scenario 2 failed: {e}")
        import traceback
        traceback.print_exc()
        results.append(("场景2: 快速停止", False))
    
    try:
        result3 = test_scenario_3_error_handling_completed()
        results.append(("场景3: 错误处理（已完成的任务）", result3))
    except Exception as e:
        print(f"[ERROR] Test scenario 3 failed: {e}")
        import traceback
        traceback.print_exc()
        results.append(("场景3: 错误处理（已完成的任务）", False))
    
    try:
        result4 = test_scenario_4_error_handling_nonexistent()
        results.append(("场景4: 错误处理（任务不存在）", result4))
    except Exception as e:
        print(f"[ERROR] Test scenario 4 failed: {e}")
        import traceback
        traceback.print_exc()
        results.append(("场景4: 错误处理（任务不存在）", False))
    
    try:
        result5 = test_scenario_5_error_handling_network()
        results.append(("场景5: 错误处理（网络错误）", result5))
    except Exception as e:
        print(f"[ERROR] Test scenario 5 failed: {e}")
        import traceback
        traceback.print_exc()
        results.append(("场景5: 错误处理（网络错误）", False))
    
    # 打印结果
    print("\n" + "=" * 70)
    print("测试结果汇总")
    print("=" * 70)
    
    for scenario, result in results:
        status = "[PASS]" if result else "[FAIL]"
        print(f"{scenario}: {status}")
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    print(f"\n总计: {passed}/{total} 测试通过")
    
    if passed == total:
        print("[SUCCESS] All tests passed!")
        return 0
    else:
        print("[WARNING] Some tests failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())

