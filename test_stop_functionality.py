"""
测试前端停止训练功能
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
        print(f"[OK] {name} is running (status: {response.status_code})")
        return response.status_code == 200
    except Exception as e:
        print(f"[FAIL] {name} is not running: {e}")
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
            task_id = task.get("task_id")
            print(f"[OK] Task created: {task_id}")
            return task_id
        else:
            print(f"[FAIL] Failed to create task: {response.status_code}")
            print(f"Response: {response.text}")
            return None
    except Exception as e:
        print(f"[FAIL] Exception creating task: {e}")
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
        if response.status_code == 200:
            print(f"[OK] Task started: {task_id}")
            return True
        else:
            print(f"[FAIL] Failed to start task: {response.status_code}")
            print(f"Response: {response.text}")
            return False
    except Exception as e:
        print(f"[FAIL] Exception starting task: {e}")
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
            print(f"[OK] Stop request sent successfully (took {elapsed:.2f}s)")
            return True, None
        else:
            error_detail = response.json().get("detail", response.text) if response.headers.get("content-type", "").startswith("application/json") else response.text
            print(f"[FAIL] Failed to stop task: {response.status_code}")
            print(f"Error: {error_detail}")
            return False, error_detail
    except requests.exceptions.Timeout:
        print(f"[FAIL] Request timeout")
        return False, "Request timeout"
    except Exception as e:
        print(f"[FAIL] Exception: {e}")
        return False, str(e)

def wait_for_task_status(task_id: str, target_status: str, max_wait: int = 60) -> bool:
    """等待任务状态变为目标状态"""
    start_time = time.time()
    while time.time() - start_time < max_wait:
        status = get_task_status(task_id)
        if status == target_status:
            print(f"[OK] Task status changed to {target_status}")
            return True
        print(f"  Current status: {status}, waiting for {target_status}...")
        time.sleep(2)
    print(f"[FAIL] Timeout waiting for status {target_status}")
    return False

def main():
    print("=" * 70)
    print("Testing Training Stop Functionality")
    print("=" * 70)
    
    # 1. 检查服务状态
    print("\n1. Checking services...")
    compass_ok = check_service(f"{COMPASS_URL}/health", "COMPASS Service")
    registry_ok = check_service("http://localhost:8500/health", "Registry Service")
    
    if not compass_ok:
        print("\n[ERROR] COMPASS service is not running. Please start it first.")
        return
    
    # 2. 创建训练任务
    print("\n2. Creating training task...")
    config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "Frontend stop functionality test"
    }
    
    task_id = create_training_task(config)
    if not task_id:
        print("\n[ERROR] Failed to create task")
        return
    
    # 3. 启动任务
    print("\n3. Starting task...")
    if not start_task(task_id):
        print("\n[ERROR] Failed to start task")
        return
    
    # 4. 等待任务进入running状态
    print("\n4. Waiting for task to start running...")
    if not wait_for_task_status(task_id, "running", max_wait=30):
        current_status = get_task_status(task_id)
        print(f"\n[WARNING] Task did not reach running status. Current status: {current_status}")
        if current_status not in ["running", "initializing"]:
            print("[ERROR] Task is not in a stoppable state")
            return
    
    # 5. 停止任务
    print("\n5. Stopping task...")
    success, error = stop_task(task_id)
    if not success:
        print(f"\n[ERROR] Failed to send stop request: {error}")
        return
    
    # 6. 等待任务状态变为cancelled
    print("\n6. Waiting for task to be cancelled...")
    if wait_for_task_status(task_id, "cancelled", max_wait=30):
        print("\n" + "=" * 70)
        print("[SUCCESS] Test passed! Task was successfully stopped.")
        print("=" * 70)
    else:
        current_status = get_task_status(task_id)
        print(f"\n[WARNING] Task status is {current_status}, expected cancelled")
        if current_status in ["completed", "failed"]:
            print("[INFO] Task ended before cancellation could take effect")
        else:
            print("[ERROR] Task stop may have failed")

if __name__ == "__main__":
    main()


