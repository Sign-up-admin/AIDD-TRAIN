"""
训练停止功能实践调试脚本
用于创建训练任务、启动训练、尝试停止并观察问题
"""
import sys
import time
import requests
import json
from typing import Dict, Optional, Any
from datetime import datetime

# 服务配置
COMPASS_URL = "http://localhost:8080"
REGISTRY_URL = "http://localhost:8500"


def check_service(url: str, name: str, timeout: int = 5) -> bool:
    """检查服务是否运行"""
    try:
        response = requests.get(url, timeout=timeout)
        return response.status_code == 200
    except:
        return False


def create_training_task() -> Optional[str]:
    """创建一个简单的训练任务（smoke_test模式）"""
    print("\n" + "="*70)
    print("创建训练任务")
    print("="*70)
    
    task_config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 5,  # 少量epochs用于快速测试
            "batch_size": 8,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "调试停止功能测试任务"
    }
    
    try:
        print(f"发送创建任务请求到 {COMPASS_URL}/api/v1/training/tasks")
        print(f"配置: {json.dumps(task_config, indent=2, ensure_ascii=False)}")
        
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks",
            json=task_config,
            timeout=10
        )
        
        if response.status_code == 201:
            task = response.json()
            task_id = task.get("task_id")
            print(f"\n[OK] 任务创建成功")
            print(f"任务ID: {task_id}")
            print(f"初始状态: {task.get('status')}")
            return task_id
        else:
            print(f"\n[ERROR] 创建任务失败: {response.status_code}")
            print(f"响应: {response.text}")
            return None
    except Exception as e:
        print(f"\n[ERROR] 创建任务时出错: {e}")
        import traceback
        traceback.print_exc()
        return None


def start_task(task_id: str) -> bool:
    """启动训练任务"""
    print(f"\n启动任务 {task_id}...")
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=10
        )
        if response.status_code == 200:
            print(f"[OK] 任务启动请求已发送")
            return True
        else:
            print(f"[ERROR] 启动失败: {response.status_code}")
            print(f"响应: {response.text}")
            return False
    except Exception as e:
        print(f"[ERROR] 启动任务时出错: {e}")
        return False


def get_task_status(task_id: str) -> Optional[Dict[str, Any]]:
    """获取任务状态"""
    try:
        response = requests.get(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
            timeout=10
        )
        if response.status_code == 200:
            return response.json()
        else:
            print(f"[ERROR] 获取任务状态失败: {response.status_code}")
            return None
    except Exception as e:
        print(f"[ERROR] 获取任务状态时出错: {e}")
        return None


def wait_for_running(task_id: str, max_wait: int = 60) -> bool:
    """等待任务进入running状态"""
    print(f"\n等待任务进入running状态（最多等待{max_wait}秒）...")
    start_time = time.time()
    
    while time.time() - start_time < max_wait:
        task = get_task_status(task_id)
        if task:
            status = task.get("status")
            progress = task.get("progress", {})
            stage = progress.get("stage", "unknown")
            
            print(f"[{datetime.now().strftime('%H:%M:%S')}] 状态: {status}, 阶段: {stage}")
            
            if status == "running":
                print(f"[OK] 任务已进入running状态")
                return True
            elif status in ["failed", "cancelled", "completed"]:
                print(f"[WARNING] 任务已进入最终状态: {status}")
                return False
        
        time.sleep(2)
    
    print(f"[WARNING] 超时：任务未在{max_wait}秒内进入running状态")
    return False


def stop_task(task_id: str) -> tuple[bool, Optional[str], float]:
    """停止任务并记录响应时间"""
    print(f"\n" + "="*70)
    print(f"发送停止请求")
    print("="*70)
    
    try:
        start_time = time.time()
        print(f"[{datetime.now().strftime('%H:%M:%S')}] 发送停止请求...")
        
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop",
            timeout=10
        )
        
        elapsed = time.time() - start_time
        
        print(f"[{datetime.now().strftime('%H:%M:%S')}] API响应时间: {elapsed:.2f}秒")
        print(f"状态码: {response.status_code}")
        
        if response.status_code == 200:
            result = response.json()
            print(f"响应内容: {json.dumps(result, indent=2, ensure_ascii=False)}")
            return True, None, elapsed
        else:
            error_detail = response.json().get("detail", response.text) if response.headers.get("content-type", "").startswith("application/json") else response.text
            print(f"错误详情: {error_detail}")
            return False, error_detail, elapsed
    except requests.exceptions.Timeout:
        print(f"[ERROR] 请求超时")
        return False, "Request timeout", 10.0
    except Exception as e:
        print(f"[ERROR] 停止请求异常: {e}")
        import traceback
        traceback.print_exc()
        return False, str(e), 0.0


def monitor_task_status(task_id: str, duration: float = 60.0, interval: float = 1.0):
    """监控任务状态变化"""
    print(f"\n" + "="*70)
    print(f"监控任务状态变化（持续 {duration} 秒，间隔 {interval} 秒）")
    print("="*70)
    
    start_time = time.time()
    last_status = None
    last_stage = None
    last_cancelled = None
    
    while time.time() - start_time < duration:
        task = get_task_status(task_id)
        if task:
            current_status = task.get("status")
            progress = task.get("progress", {})
            current_stage = progress.get("stage", "unknown")
            current_cancelled = progress.get("cancelled", False)
            
            # 检测状态变化
            status_changed = current_status != last_status
            stage_changed = current_stage != last_stage
            cancelled_changed = current_cancelled != last_cancelled
            
            if status_changed or stage_changed or cancelled_changed:
                timestamp = datetime.now().strftime('%H:%M:%S')
                print(f"\n[{timestamp}] 状态变化检测:")
                if status_changed:
                    print(f"  状态: {last_status} -> {current_status}")
                if stage_changed:
                    print(f"  阶段: {last_stage} -> {current_stage}")
                if cancelled_changed:
                    print(f"  取消标志: {last_cancelled} -> {current_cancelled}")
                
                # 显示完整进度信息
                print(f"  完整进度信息:")
                print(f"    {json.dumps(progress, indent=4, ensure_ascii=False)}")
                
                last_status = current_status
                last_stage = current_stage
                last_cancelled = current_cancelled
            
            # 如果任务已经停止，提前结束监控
            if current_status in ['cancelled', 'completed', 'failed']:
                print(f"\n[{datetime.now().strftime('%H:%M:%S')}] 任务已进入最终状态: {current_status}")
                break
        
        time.sleep(interval)
    
    print(f"\n监控结束")


def main():
    """主函数"""
    print("="*70)
    print("训练停止功能实践调试")
    print("="*70)
    
    # 检查服务
    print("\n检查服务状态...")
    registry_ok = check_service(f"{REGISTRY_URL}/health", "Registry")
    compass_ok = check_service(f"{COMPASS_URL}/health", "COMPASS")
    
    print(f"Registry (8500): {'[OK]' if registry_ok else '[NOT RUNNING]'}")
    print(f"COMPASS (8080): {'[OK]' if compass_ok else '[NOT RUNNING]'}")
    
    if not compass_ok:
        print("\n[ERROR] COMPASS服务未运行，请先启动服务")
        return 1
    
    # 创建训练任务
    task_id = create_training_task()
    if not task_id:
        print("\n[ERROR] 无法创建训练任务")
        return 1
    
    # 启动任务
    if not start_task(task_id):
        print("\n[ERROR] 无法启动训练任务")
        return 1
    
    # 等待任务进入running状态
    if not wait_for_running(task_id, max_wait=120):
        print("\n[ERROR] 任务未进入running状态")
        return 1
    
    # 等待几秒让训练真正开始
    print("\n等待5秒让训练真正开始...")
    time.sleep(5)
    
    # 发送停止请求
    success, error, elapsed = stop_task(task_id)
    
    if success:
        print(f"\n[OK] 停止请求已发送（响应时间: {elapsed:.2f}秒）")
    else:
        print(f"\n[FAIL] 停止请求失败: {error}")
        print(f"响应时间: {elapsed:.2f}秒")
    
    # 监控状态变化
    monitor_task_status(task_id, duration=60.0, interval=1.0)
    
    # 最终状态检查
    print("\n" + "="*70)
    print("最终状态检查")
    print("="*70)
    final_task = get_task_status(task_id)
    if final_task:
        final_status = final_task.get("status")
        final_progress = final_task.get("progress", {})
        print(f"最终状态: {final_status}")
        print(f"最终进度信息:")
        print(json.dumps(final_progress, indent=2, ensure_ascii=False))
        
        if final_status == 'cancelled':
            print("\n[SUCCESS] 任务已成功停止")
            return 0
        elif final_status in ['completed', 'failed']:
            print(f"\n[INFO] 任务已进入最终状态: {final_status}")
            return 0
        else:
            print(f"\n[WARNING] 任务状态仍为: {final_status}")
            print("任务可能仍在停止中，请检查日志")
            return 1
    else:
        print("[ERROR] 无法获取最终状态")
        return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\n用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] 未预期的错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)










