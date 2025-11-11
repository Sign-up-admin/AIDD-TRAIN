"""
Chrome浏览器调试辅助脚本
帮助创建训练任务并监控停止请求
"""
import requests
import time
import json
from datetime import datetime
from typing import Optional, Dict, Any

# 服务配置
COMPASS_URL = "http://localhost:8080"
REGISTRY_URL = "http://localhost:8500"
FLASHDOCK_URL = "http://localhost:8501"

def check_services():
    """检查所有服务是否运行"""
    print("=" * 70)
    print("检查服务状态")
    print("=" * 70)
    
    services = {
        "服务注册中心": REGISTRY_URL + "/health",
        "COMPASS服务": COMPASS_URL + "/health",
        "FLASH-DOCK前端": FLASHDOCK_URL,
    }
    
    all_ok = True
    for name, url in services.items():
        try:
            r = requests.get(url, timeout=5)
            if r.status_code == 200 or (name == "FLASH-DOCK前端" and r.status_code < 500):
                print(f"[OK] {name}: 运行中")
            else:
                print(f"[FAIL] {name}: 状态码 {r.status_code}")
                all_ok = False
        except Exception as e:
            print(f"[FAIL] {name}: 无法连接 - {e}")
            all_ok = False
    
    print()
    return all_ok

def create_training_task(config: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """创建训练任务"""
    if config is None:
        config = {
            "config": {
                "execution_mode": "smoke_test",
                "epochs": 10,
                "batch_size": 2,
                "learning_rate": 0.001,
                "optimizer": "adam"
            },
            "description": f"Chrome调试测试任务 - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        }
    
    print("=" * 70)
    print("创建训练任务")
    print("=" * 70)
    print(f"配置: {json.dumps(config, indent=2, ensure_ascii=False)}")
    print()
    
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks",
            json=config,
            timeout=30
        )
        
        if response.status_code == 201:
            task = response.json()
            task_id = task.get("task_id")
            print(f"[OK] 任务创建成功!")
            print(f"  任务ID: {task_id}")
            print(f"  状态: {task.get('status')}")
            print()
            return task_id
        else:
            print(f"[FAIL] 任务创建失败: {response.status_code}")
            print(f"  响应: {response.text}")
            return None
    except Exception as e:
        print(f"[FAIL] 创建任务时出错: {e}")
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
        return None
    except Exception as e:
        print(f"获取任务状态失败: {e}")
        return None

def start_task(task_id: str) -> bool:
    """启动训练任务"""
    print("=" * 70)
    print("启动训练任务")
    print("=" * 70)
    print(f"任务ID: {task_id}")
    print()
    
    try:
        response = requests.post(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=30
        )
        
        if response.status_code == 200:
            print("[OK] 任务启动请求已发送")
            print()
            return True
        else:
            print(f"[FAIL] 任务启动失败: {response.status_code}")
            print(f"  响应: {response.text}")
            return False
    except Exception as e:
        print(f"[FAIL] 启动任务时出错: {e}")
        return False

def wait_for_task_status(task_id: str, target_status: str, max_wait: int = 60, check_interval: int = 2):
    """等待任务状态变为目标状态"""
    print(f"等待任务状态变为 '{target_status}'...")
    print("(按Ctrl+C可以中断等待)")
    print()
    
    start_time = time.time()
    last_status = None
    
    try:
        while time.time() - start_time < max_wait:
            status = get_task_status(task_id)
            if status:
                if status != last_status:
                    print(f"  当前状态: {status} (等待了 {int(time.time() - start_time)} 秒)")
                    last_status = status
                
                if status == target_status:
                    print(f"[OK] 任务状态已变为 '{target_status}'")
                    print()
                    return True
            
            time.sleep(check_interval)
        
        print(f"[WARNING] 等待超时 ({max_wait}秒)，当前状态: {last_status}")
        print()
        return False
    except KeyboardInterrupt:
        print("\n[WARNING] 等待被用户中断")
        print(f"  当前状态: {last_status}")
        print()
        return False

def monitor_stop_request(task_id: str, duration: int = 30):
    """监控停止请求（模拟）"""
    print("=" * 70)
    print("监控停止请求")
    print("=" * 70)
    print("请在Chrome浏览器中:")
    print("1. 打开开发者工具 (F12)")
    print("2. 切换到 Network 标签页")
    print("3. 勾选 'Preserve log'")
    print("4. 点击 '停止任务' 按钮")
    print("5. 观察 Network 标签页中的请求")
    print()
    print(f"任务ID: {task_id}")
    print(f"监控时长: {duration} 秒")
    print()
    
    # 定期检查任务状态
    start_time = time.time()
    last_status = None
    
    try:
        while time.time() - start_time < duration:
            status = get_task_status(task_id)
            if status and status != last_status:
                print(f"[{int(time.time() - start_time)}s] 任务状态: {status}")
                last_status = status
                
                if status in ["cancelled", "completed", "failed"]:
                    print(f"[OK] 任务已结束，状态: {status}")
                    break
            
            time.sleep(2)
    except KeyboardInterrupt:
        print("\n[INFO] 监控被用户中断")
    
    # 获取最终状态
    final_status = get_task_status(task_id)
    print()
    print(f"最终任务状态: {final_status}")
    print()

def print_chrome_debug_instructions():
    """打印Chrome调试说明"""
    print("=" * 70)
    print("Chrome浏览器调试步骤")
    print("=" * 70)
    print()
    print("1. 打开Chrome浏览器，访问: http://localhost:8501")
    print("2. 按F12打开开发者工具")
    print("3. 切换到以下标签页:")
    print("   - Network (网络请求监控)")
    print("   - Console (控制台日志)")
    print()
    print("4. 在训练管理页面:")
    print("   - 找到刚创建的任务")
    print("   - 点击 '停止任务' 按钮")
    print()
    print("5. 在Network标签页中观察:")
    print("   - 查找 POST 请求到 /api/v1/training/tasks/{task_id}/stop")
    print("   - 检查请求状态码 (200/400/404/500/超时)")
    print("   - 检查响应时间")
    print("   - 检查响应内容")
    print()
    print("6. 在Console标签页中观察:")
    print("   - 是否有JavaScript错误")
    print("   - 是否有网络错误")
    print()
    print("7. 观察页面显示:")
    print("   - 是否显示加载状态")
    print("   - 是否显示成功/错误消息")
    print("   - 错误消息是否可见")
    print("   - 任务状态是否更新")
    print()
    print("=" * 70)
    print()

def main():
    """主函数"""
    print()
    print("=" * 70)
    print("Chrome浏览器调试辅助工具")
    print("=" * 70)
    print()
    
    # 检查服务
    if not check_services():
        print("[WARNING] 警告: 部分服务未运行，请先启动服务")
        print()
        return
    
    print()
    
    # 创建任务
    task_id = create_training_task()
    if not task_id:
        print("无法创建任务，退出")
        return
    
    print()
    print("=" * 70)
    print("下一步操作")
    print("=" * 70)
    print()
    print("选择操作:")
    print("1. 自动启动任务并等待running状态")
    print("2. 手动启动任务（在浏览器中）")
    print()
    
    choice = input("请选择 (1/2，默认1): ").strip() or "1"
    
    if choice == "1":
        # 启动任务
        if start_task(task_id):
            # 等待running状态
            if wait_for_task_status(task_id, "running", max_wait=60):
                print()
                print_chrome_debug_instructions()
                print()
                print("现在请在Chrome浏览器中点击 '停止任务' 按钮")
                print("同时观察Network和Console标签页")
                print()
                # 监控停止请求
                monitor_stop_request(task_id, duration=60)
            else:
                print("任务未能进入running状态")
                current_status = get_task_status(task_id)
                print(f"当前状态: {current_status}")
                if current_status in ["running", "initializing"]:
                    print()
                    print_chrome_debug_instructions()
                    monitor_stop_request(task_id, duration=60)
    else:
        print()
        print(f"任务ID: {task_id}")
        print("请在浏览器中手动启动任务，然后点击停止按钮")
        print()
        print_chrome_debug_instructions()
        monitor_stop_request(task_id, duration=120)
    
    print()
    print("=" * 70)
    print("调试完成")
    print("=" * 70)
    print()
    print("请记录以下信息:")
    print("1. Network标签页中的请求详情")
    print("2. Console标签页中的错误信息")
    print("3. 页面显示的行为")
    print("4. 后端服务窗口中的日志")
    print()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n程序被用户中断")
    except Exception as e:
        print(f"\n\n错误: {e}")
        import traceback
        traceback.print_exc()

