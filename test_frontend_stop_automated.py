"""
自动化测试前端停止功能 - 模拟完整的前端操作流程
"""
import sys
import os
import time
import requests
import json
from datetime import datetime

# 添加路径以便导入compass_client
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'FLASH_DOCK-main'))

COMPASS_URL = "http://localhost:8080"

def check_service():
    """检查服务是否运行"""
    try:
        response = requests.get(f"{COMPASS_URL}/health", timeout=2)
        return response.status_code == 200
    except:
        return False

def create_test_task():
    """创建一个测试任务"""
    task_config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 8,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "自动化测试-前端停止功能"
    }
    
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks",
            json=task_config,
            timeout=10
        )
        if response.status_code == 201:
            task = response.json()
            return task.get("task_id")
        else:
            print(f"[ERROR] 创建任务失败: {response.status_code}")
            print(f"响应: {response.text}")
            return None
    except Exception as e:
        print(f"[ERROR] 创建任务时出错: {e}")
        return None

def start_task(task_id):
    """启动任务"""
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=10
        )
        return response.status_code == 200
    except Exception as e:
        print(f"[ERROR] 启动任务时出错: {e}")
        return False

def get_task_status(task_id):
    """获取任务状态"""
    try:
        response = requests.get(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
            timeout=10
        )
        if response.status_code == 200:
            return response.json()
        return None
    except:
        return None

def stop_task(task_id):
    """停止任务（模拟前端调用）"""
    try:
        start_time = time.time()
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop",
            timeout=10
        )
        elapsed = time.time() - start_time
        
        if response.status_code == 200:
            result = response.json()
            return True, result, elapsed
        else:
            error_detail = response.json().get("detail", response.text) if response.headers.get("content-type", "").startswith("application/json") else response.text
            return False, error_detail, elapsed
    except requests.exceptions.Timeout:
        return False, "Request timeout", 10.0
    except Exception as e:
        return False, str(e), 0.0

def simulate_frontend_stop_flow(task_id):
    """模拟前端停止流程"""
    print("\n" + "="*70)
    print("模拟前端停止操作流程")
    print("="*70)
    
    # 步骤1: 重新获取最新状态（模拟前端状态检查）
    print("\n[步骤1] 重新获取最新状态（模拟前端状态检查）...")
    try:
        current_task = get_task_status(task_id)
        if not current_task:
            print("  [FAIL] 无法获取任务状态")
            return False
        
        current_status = current_task.get("status")
        print(f"  [OK] 当前状态: {current_status}")
        
        if current_status != "running" and current_status != "initializing":
            print(f"  [WARN] 任务状态为 '{current_status}'，无法停止")
            return False
    except Exception as e:
        print(f"  [ERROR] 状态检查失败: {e}")
        return False
    
    # 步骤2: 发送停止请求（模拟前端点击停止按钮）
    print("\n[步骤2] 发送停止请求（模拟前端点击停止按钮）...")
    stop_success, stop_result, stop_elapsed = stop_task(task_id)
    
    if not stop_success:
        print(f"  [FAIL] 停止请求失败: {stop_result}")
        print(f"  响应时间: {stop_elapsed:.2f}秒")
        return False
    
    print(f"  [OK] 停止请求已发送（响应时间: {stop_elapsed:.2f}秒）")
    if isinstance(stop_result, dict):
        print(f"  响应内容: {json.dumps(stop_result, indent=2, ensure_ascii=False)}")
    
    # 步骤3: 轮询检查任务状态（模拟前端轮询）
    print("\n[步骤3] 轮询检查任务状态（模拟前端轮询，最多30秒）...")
    max_poll_time = 30.0
    poll_interval = 1.0
    poll_elapsed = 0.0
    final_status = None
    status_updates = []
    
    while poll_elapsed < max_poll_time:
        try:
            updated_task = get_task_status(task_id)
            if updated_task:
                final_status = updated_task.get('status')
                progress = updated_task.get('progress', {})
                cancelled = progress.get('cancelled', False)
                
                status_info = {
                    'time': poll_elapsed,
                    'status': final_status,
                    'cancelled': cancelled,
                    'stage': progress.get('stage', 'unknown')
                }
                status_updates.append(status_info)
                
                # 显示状态更新
                print(f"  [{poll_elapsed:.1f}s] 状态: {final_status}, 取消标志: {cancelled}, 阶段: {progress.get('stage', 'unknown')}")
                
                if final_status in ['cancelled', 'completed', 'failed']:
                    print(f"  [SUCCESS] 任务状态已更新为: {final_status}")
                    break
                elif final_status not in ['running', 'initializing']:
                    print(f"  [WARN] 任务状态: {final_status} (非预期状态)")
                    break
        except Exception as poll_e:
            print(f"  [WARN] 轮询时出错: {str(poll_e)[:50]}... (继续轮询)")
        
        time.sleep(poll_interval)
        poll_elapsed += poll_interval
    
    # 步骤4: 验证结果
    print("\n[步骤4] 验证停止结果...")
    if final_status == 'cancelled':
        print("  [SUCCESS] 任务已成功停止！")
        
        # 显示状态变化时间线
        print("\n  状态变化时间线:")
        for update in status_updates:
            print(f"    {update['time']:.1f}s: {update['status']} (cancelled={update['cancelled']}, stage={update['stage']})")
        
        # 计算响应时间
        if status_updates:
            first_update = status_updates[0]
            last_update = status_updates[-1]
            total_time = last_update['time']
            print(f"\n  总响应时间: {total_time:.1f}秒")
            print(f"  API响应时间: {stop_elapsed:.2f}秒")
            print(f"  状态更新延迟: {total_time - stop_elapsed:.1f}秒")
        
        return True
    elif final_status in ['completed', 'failed']:
        print(f"  [INFO] 任务已进入最终状态: {final_status}")
        return True
    elif final_status in ['running', 'initializing']:
        print(f"  [FAIL] 任务状态仍为: {final_status}")
        print("  停止功能可能存在问题")
        return False
    else:
        print(f"  [WARN] 未知状态: {final_status}")
        return False

def main():
    """主测试函数"""
    print("="*70)
    print("前端停止功能自动化测试")
    print("="*70)
    
    # 检查服务
    print("\n[准备] 检查COMPASS服务...")
    if not check_service():
        print("[ERROR] COMPASS服务未运行，请先启动服务")
        print("运行: python check_and_start_services.py")
        return 1
    print("[OK] COMPASS服务运行正常")
    
    # 创建任务
    print("\n[准备] 创建测试任务...")
    task_id = create_test_task()
    if not task_id:
        print("[ERROR] 无法创建测试任务")
        return 1
    print(f"[OK] 任务创建成功")
    print(f"任务ID: {task_id}")
    
    # 启动任务
    print("\n[准备] 启动测试任务...")
    if not start_task(task_id):
        print("[ERROR] 无法启动测试任务")
        return 1
    print("[OK] 任务启动请求已发送")
    
    # 等待任务进入running状态
    print("\n[准备] 等待任务进入running状态...")
    max_wait = 60
    waited = 0
    while waited < max_wait:
        task = get_task_status(task_id)
        if task:
            status = task.get("status")
            if waited % 5 == 0:  # 每5秒打印一次
                print(f"  [{waited}s] 状态: {status}")
            if status == "running":
                print(f"  [OK] 任务已进入running状态（等待了{waited}秒）")
                break
            elif status in ["failed", "cancelled", "completed"]:
                print(f"  [WARN] 任务已进入最终状态: {status}")
                return 1
        time.sleep(2)
        waited += 2
    
    if waited >= max_wait:
        print("[ERROR] 超时：任务未在60秒内进入running状态")
        return 1
    
    # 等待几秒让训练真正开始
    print("\n[准备] 等待5秒让训练真正开始...")
    time.sleep(5)
    
    # 执行前端停止流程测试
    success = simulate_frontend_stop_flow(task_id)
    
    # 最终结果
    print("\n" + "="*70)
    if success:
        print("测试结果: [SUCCESS] 前端停止功能测试通过！")
        print("="*70)
        return 0
    else:
        print("测试结果: [FAIL] 前端停止功能测试失败")
        print("="*70)
        return 1

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\n用户中断测试")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] 未预期的错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)










