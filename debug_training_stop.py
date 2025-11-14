"""
训练停止功能调试工具

用于诊断和监控训练任务的停止功能。
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
            print(f"[ERROR] Failed to get task status: {response.status_code}")
            print(f"Response: {response.text}")
            return None
    except Exception as e:
        print(f"[ERROR] Exception getting task status: {e}")
        return None


def stop_task(task_id: str) -> tuple[bool, Optional[str], float]:
    """停止任务并记录响应时间"""
    try:
        start_time = time.time()
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/stop",
            timeout=10
        )
        elapsed = time.time() - start_time
        
        if response.status_code == 200:
            return True, None, elapsed
        else:
            error_detail = response.json().get("detail", response.text) if response.headers.get("content-type", "").startswith("application/json") else response.text
            return False, error_detail, elapsed
    except requests.exceptions.Timeout:
        return False, "Request timeout", 10.0
    except Exception as e:
        return False, str(e), 0.0


def monitor_task_status(task_id: str, duration: float = 30.0, interval: float = 1.0):
    """监控任务状态变化"""
    print(f"\n开始监控任务 {task_id} 的状态变化（持续 {duration} 秒）...")
    start_time = time.time()
    last_status = None
    
    while time.time() - start_time < duration:
        task = get_task_status(task_id)
        if task:
            current_status = task.get("status")
            progress = task.get("progress", {})
            
            if current_status != last_status:
                print(f"\n[{datetime.now().strftime('%H:%M:%S')}] 状态变化: {last_status} -> {current_status}")
                if progress:
                    stage = progress.get("stage", "unknown")
                    cancelled = progress.get("cancelled", False)
                    print(f"  阶段: {stage}, 已取消标志: {cancelled}")
                last_status = current_status
            
            # 如果任务已经停止，提前结束监控
            if current_status in ['cancelled', 'completed', 'failed']:
                print(f"\n任务已进入最终状态: {current_status}")
                break
        
        time.sleep(interval)
    
    print("\n监控结束")


def test_stop_functionality():
    """测试停止功能"""
    print("=" * 70)
    print("训练停止功能调试工具")
    print("=" * 70)
    
    # 检查服务
    print("\n检查服务状态...")
    registry_ok = check_service(f"{REGISTRY_URL}/health", "Registry")
    compass_ok = check_service(f"{COMPASS_URL}/health", "COMPASS")
    
    print(f"Registry (8500): {'[OK]' if registry_ok else '[NOT RUNNING]'}")
    print(f"COMPASS (8080): {'[OK]' if compass_ok else '[NOT RUNNING]'}")
    
    if not compass_ok:
        print("\n[ERROR] COMPASS服务未运行，请先启动服务")
        return 1
    
    # 获取任务ID
    print("\n请输入要停止的任务ID（或按Enter跳过，使用交互模式）:")
    task_id_input = input().strip()
    
    if not task_id_input:
        # 交互模式：列出所有任务
        print("\n获取任务列表...")
        try:
            response = requests.get(
                f"{COMPASS_URL}/api/v1/training/tasks",
                timeout=10
            )
            if response.status_code == 200:
                tasks = response.json()
                if isinstance(tasks, list) and len(tasks) > 0:
                    print("\n可用任务:")
                    running_tasks = []
                    for i, task in enumerate(tasks, 1):
                        task_id = task.get("task_id", "unknown")
                        status = task.get("status", "unknown")
                        description = task.get("description", "")
                        print(f"  {i}. {task_id[:8]}... - {status} - {description}")
                        if status in ['running', 'initializing']:
                            running_tasks.append((i, task_id, status))
                    
                    if running_tasks:
                        print("\n运行中的任务:")
                        for idx, task_id, status in running_tasks:
                            print(f"  {idx}. {task_id} ({status})")
                        
                        print("\n请选择要停止的任务编号（或输入任务ID）:")
                        choice = input().strip()
                        
                        # 尝试解析为数字
                        try:
                            choice_num = int(choice)
                            for idx, tid, _ in running_tasks:
                                if idx == choice_num:
                                    task_id_input = tid
                                    break
                        except ValueError:
                            # 不是数字，当作任务ID
                            task_id_input = choice
                    else:
                        print("\n没有运行中的任务")
                        return 0
                else:
                    print("\n没有可用任务")
                    return 0
            else:
                print(f"[ERROR] 获取任务列表失败: {response.status_code}")
                return 1
        except Exception as e:
            print(f"[ERROR] 获取任务列表时出错: {e}")
            return 1
    
    if not task_id_input:
        print("\n未指定任务ID，退出")
        return 0
    
    # 获取任务当前状态
    print(f"\n获取任务 {task_id_input} 的当前状态...")
    task = get_task_status(task_id_input)
    if not task:
        print("[ERROR] 无法获取任务状态")
        return 1
    
    current_status = task.get("status")
    progress = task.get("progress", {})
    
    print(f"当前状态: {current_status}")
    print(f"进度信息: {json.dumps(progress, indent=2, ensure_ascii=False)}")
    
    if current_status not in ['running', 'initializing']:
        print(f"\n[WARNING] 任务状态为 '{current_status}'，无法停止")
        print("只有运行中或初始化中的任务可以停止")
        return 0
    
    # 确认停止
    print(f"\n确认要停止任务 {task_id_input} 吗？(y/n): ", end="")
    confirm = input().strip().lower()
    if confirm != 'y':
        print("已取消")
        return 0
    
    # 发送停止请求
    print(f"\n发送停止请求...")
    success, error, elapsed = stop_task(task_id_input)
    
    if success:
        print(f"[OK] 停止请求已发送（响应时间: {elapsed:.2f}秒）")
        
        # 监控状态变化
        monitor_task_status(task_id_input, duration=30.0, interval=0.5)
        
        # 最终状态检查
        print("\n最终状态检查...")
        final_task = get_task_status(task_id_input)
        if final_task:
            final_status = final_task.get("status")
            final_progress = final_task.get("progress", {})
            print(f"最终状态: {final_status}")
            print(f"最终进度: {json.dumps(final_progress, indent=2, ensure_ascii=False)}")
            
            if final_status == 'cancelled':
                print("\n[SUCCESS] 任务已成功停止")
                return 0
            elif final_status in ['completed', 'failed']:
                print(f"\n[INFO] 任务已进入最终状态: {final_status}")
                return 0
            else:
                print(f"\n[WARNING] 任务状态仍为: {final_status}")
                print("任务可能仍在停止中，请稍后检查")
                return 1
        else:
            print("[ERROR] 无法获取最终状态")
            return 1
    else:
        print(f"[FAIL] 停止请求失败: {error}")
        print(f"响应时间: {elapsed:.2f}秒")
        return 1


if __name__ == "__main__":
    try:
        sys.exit(test_stop_functionality())
    except KeyboardInterrupt:
        print("\n\n用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] 未预期的错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)










