"""测试停止现有运行中的任务"""
import requests
import time
import json
from datetime import datetime

COMPASS_URL = "http://localhost:8080"

def get_task_status(task_id):
    """获取任务状态"""
    try:
        response = requests.get(f"{COMPASS_URL}/api/v1/training/tasks/{task_id}", timeout=10)
        if response.status_code == 200:
            return response.json()
        return None
    except Exception as e:
        print(f"Error getting task status: {e}")
        return None

def stop_task(task_id):
    """停止任务"""
    try:
        start_time = time.time()
        response = requests.post(f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop", timeout=10)
        elapsed = time.time() - start_time
        
        if response.status_code == 200:
            return True, None, elapsed, response.json()
        else:
            error_detail = response.json().get("detail", response.text) if response.headers.get("content-type", "").startswith("application/json") else response.text
            return False, error_detail, elapsed, None
    except Exception as e:
        return False, str(e), 0.0, None

def monitor_status(task_id, duration=30, interval=0.5):
    """监控任务状态"""
    print(f"\n监控任务状态变化（持续 {duration} 秒）...")
    start_time = time.time()
    last_status = None
    status_changes = []
    
    while time.time() - start_time < duration:
        task = get_task_status(task_id)
        if task:
            current_status = task.get("status")
            progress = task.get("progress", {})
            cancelled = progress.get("cancelled", False) if isinstance(progress, dict) else False
            
            if current_status != last_status:
                timestamp = datetime.now().strftime('%H:%M:%S')
                status_changes.append((timestamp, last_status, current_status))
                print(f"\n[{timestamp}] 状态变化: {last_status} -> {current_status}")
                if isinstance(progress, dict):
                    stage = progress.get("stage", "unknown")
                    print(f"  阶段: {stage}, 已取消标志: {cancelled}")
                last_status = current_status
            
            if current_status in ['cancelled', 'completed', 'failed']:
                print(f"\n任务已进入最终状态: {current_status}")
                break
        
        time.sleep(interval)
    
    return status_changes

# 查找运行中的任务
print("=" * 70)
print("测试停止现有任务功能")
print("=" * 70)

print("\n查找运行中的任务...")
try:
    response = requests.get(f"{COMPASS_URL}/api/v1/training/tasks", timeout=10)
    if response.status_code == 200:
        data = response.json()
        running_tasks = []
        
        # Handle different response formats
        if isinstance(data, dict) and 'tasks' in data:
            tasks = data['tasks']
        elif isinstance(data, list):
            tasks = data
        elif isinstance(data, dict):
            tasks = [data]
        else:
            tasks = []
        
        for task in tasks:
            if isinstance(task, dict):
                task_id = task.get('task_id', '')
                status = task.get('status', '')
                if status in ['running', 'initializing']:
                    running_tasks.append((task_id, status))
        
        if not running_tasks:
            print("没有找到运行中的任务")
            print("\n注意: CPU使用率过高，无法创建新任务")
            print("请等待CPU使用率降低后再测试")
            exit(1)
        
        # 使用第一个运行中的任务
        task_id, status = running_tasks[0]
        print(f"\n找到运行中的任务: {task_id[:8]}... (状态: {status})")
        
        # 获取当前状态
        print("\n获取任务当前状态...")
        task = get_task_status(task_id)
        if task:
            print(f"当前状态: {task.get('status')}")
            progress = task.get('progress', {})
            if isinstance(progress, dict):
                print(f"阶段: {progress.get('stage', 'unknown')}")
                print(f"已取消标志: {progress.get('cancelled', False)}")
        
        # 发送停止请求
        print(f"\n发送停止请求...")
        success, error, elapsed, result = stop_task(task_id)
        
        if success:
            print(f"[OK] 停止请求已发送（响应时间: {elapsed:.2f}秒）")
            if result:
                print(f"响应: {json.dumps(result, indent=2, ensure_ascii=False)}")
            
            # 监控状态变化
            status_changes = monitor_status(task_id, duration=30, interval=0.5)
            
            # 最终状态检查
            print("\n最终状态检查...")
            final_task = get_task_status(task_id)
            if final_task:
                final_status = final_task.get("status")
                final_progress = final_task.get("progress", {})
                print(f"最终状态: {final_status}")
                
                if isinstance(final_progress, dict):
                    print(f"最终进度信息:")
                    print(f"  阶段: {final_progress.get('stage', 'unknown')}")
                    print(f"  已取消标志: {final_progress.get('cancelled', False)}")
                
                print("\n" + "=" * 70)
                if final_status == 'cancelled':
                    print("[SUCCESS] 任务已成功停止！")
                    print("\n状态变化记录:")
                    for timestamp, old_status, new_status in status_changes:
                        print(f"  [{timestamp}] {old_status} -> {new_status}")
                    exit(0)
                elif final_status in ['completed', 'failed']:
                    print(f"[INFO] 任务已进入最终状态: {final_status}")
                    exit(0)
                else:
                    print(f"[WARNING] 任务状态仍为: {final_status}")
                    print("任务可能仍在停止中，请稍后检查")
                    exit(1)
            else:
                print("[ERROR] 无法获取最终状态")
                exit(1)
        else:
            print(f"[FAIL] 停止请求失败: {error}")
            print(f"响应时间: {elapsed:.2f}秒")
            exit(1)
    else:
        print(f"获取任务列表失败: {response.status_code}")
        exit(1)
except Exception as e:
    print(f"错误: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

