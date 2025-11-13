"""
Chrome浏览器调试自动化脚本
自动创建任务、启动任务，然后等待用户在Chrome中点击停止按钮
"""
import requests
import time
import json
from datetime import datetime
from typing import Optional

# 服务配置
COMPASS_URL = "http://localhost:8080"

def create_and_start_task():
    """创建并启动训练任务"""
    print("=" * 70)
    print("Chrome浏览器调试自动化脚本")
    print("=" * 70)
    print()
    
    # 创建任务
    print("1. 创建训练任务...")
    config = {
        "config": {
            "execution_mode": "smoke_test",
            "epochs": 10,
            "batch_size": 2,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": f"Chrome调试测试 - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    }
    
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks",
            json=config,
            timeout=30
        )
        
        if response.status_code != 201:
            print(f"[FAIL] 创建任务失败: {response.status_code}")
            print(f"响应: {response.text}")
            return None
        
        task = response.json()
        task_id = task.get("task_id")
        print(f"[OK] 任务创建成功: {task_id}")
        print()
        
        # 启动任务
        print("2. 启动训练任务...")
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=30
        )
        
        if response.status_code != 200:
            print(f"[FAIL] 启动任务失败: {response.status_code}")
            print(f"响应: {response.text}")
            return None
        
        print("[OK] 任务启动请求已发送")
        print()
        
        # 等待running状态
        print("3. 等待任务进入running状态...")
        max_wait = 60
        start_time = time.time()
        
        while time.time() - start_time < max_wait:
            try:
                response = requests.get(
                    f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
                    timeout=10
                )
                if response.status_code == 200:
                    task = response.json()
                    status = task.get("status")
                    print(f"  当前状态: {status} (等待了 {int(time.time() - start_time)} 秒)")
                    
                    if status == "running":
                        print("[OK] 任务已进入running状态")
                        print()
                        return task_id
                    elif status in ["cancelled", "completed", "failed"]:
                        print(f"[WARNING] 任务已结束，状态: {status}")
                        return None
            except Exception as e:
                print(f"[WARNING] 获取状态失败: {e}")
            
            time.sleep(2)
        
        print("[WARNING] 等待超时，但任务可能仍在运行")
        print()
        return task_id
        
    except Exception as e:
        print(f"[FAIL] 错误: {e}")
        import traceback
        traceback.print_exc()
        return None

def monitor_task(task_id: str):
    """监控任务状态"""
    print("=" * 70)
    print("任务监控中...")
    print("=" * 70)
    print()
    print(f"任务ID: {task_id}")
    print()
    print("请在Chrome浏览器中:")
    print("1. 打开 http://localhost:8501")
    print("2. 按F12打开开发者工具")
    print("3. 切换到 Network 标签页，勾选 'Preserve log'")
    print("4. 切换到 Console 标签页")
    print("5. 在训练管理页面找到任务，点击 '停止任务' 按钮")
    print("6. 观察 Network 和 Console 标签页")
    print()
    print("监控任务状态变化（按Ctrl+C停止监控）...")
    print()
    
    last_status = None
    start_time = time.time()
    
    try:
        while True:
            try:
                response = requests.get(
                    f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
                    timeout=10
                )
                if response.status_code == 200:
                    task = response.json()
                    status = task.get("status")
                    
                    if status != last_status:
                        elapsed = int(time.time() - start_time)
                        print(f"[{elapsed}s] 任务状态: {status}")
                        last_status = status
                        
                        if status in ["cancelled", "completed", "failed"]:
                            print()
                            print(f"[OK] 任务已结束，最终状态: {status}")
                            break
            except Exception as e:
                print(f"[WARNING] 获取状态失败: {e}")
            
            time.sleep(2)
    except KeyboardInterrupt:
        print()
        print("[INFO] 监控被用户中断")
        if last_status:
            print(f"最后状态: {last_status}")
    
    print()
    print("=" * 70)
    print("调试信息记录")
    print("=" * 70)
    print()
    print("请记录以下信息:")
    print("1. Network标签页中的停止请求:")
    print("   - 请求URL")
    print("   - 请求状态码 (200/400/404/500/超时)")
    print("   - 响应时间")
    print("   - 响应内容")
    print()
    print("2. Console标签页中的错误:")
    print("   - JavaScript错误")
    print("   - 网络错误")
    print()
    print("3. 页面显示:")
    print("   - 是否显示加载状态")
    print("   - 是否显示成功/错误消息")
    print("   - 错误消息是否可见")
    print("   - 任务状态是否更新")
    print()
    print("4. 后端服务窗口中的日志:")
    print("   - 查找 [STOP] 标记的日志")
    print("   - 检查错误或异常")
    print()

if __name__ == "__main__":
    task_id = create_and_start_task()
    if task_id:
        monitor_task(task_id)
    else:
        print("[FAIL] 无法创建或启动任务")








