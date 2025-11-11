"""
准备前端测试环境 - 创建并启动一个测试任务
"""
import sys
import time
import requests
import json

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
            "epochs": 10,  # 足够的epochs以便有时间测试停止
            "batch_size": 8,
            "learning_rate": 0.001,
            "optimizer": "adam"
        },
        "description": "前端停止功能测试任务"
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
            print(f"创建任务失败: {response.status_code}")
            print(f"响应: {response.text}")
            return None
    except Exception as e:
        print(f"创建任务时出错: {e}")
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
        print(f"启动任务时出错: {e}")
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

def main():
    print("="*70)
    print("准备前端测试环境")
    print("="*70)
    
    # 检查服务
    print("\n1. 检查COMPASS服务...")
    if not check_service():
        print("   ❌ COMPASS服务未运行，请先启动服务")
        print("   运行: python check_and_start_services.py")
        return 1
    print("   [OK] COMPASS服务运行正常")
    
    # 创建任务
    print("\n2. 创建测试任务...")
    task_id = create_test_task()
    if not task_id:
        print("   ❌ 无法创建测试任务")
        return 1
    print(f"   [OK] 任务创建成功")
    print(f"   任务ID: {task_id}")
    
    # 启动任务
    print("\n3. 启动测试任务...")
    if not start_task(task_id):
        print("   ❌ 无法启动测试任务")
        return 1
    print("   [OK] 任务启动请求已发送")
    
    # 等待任务进入running状态
    print("\n4. 等待任务进入running状态...")
    max_wait = 60
    waited = 0
    while waited < max_wait:
        task = get_task_status(task_id)
        if task:
            status = task.get("status")
            print(f"   [{waited}s] 状态: {status}")
            if status == "running":
                print("   [OK] 任务已进入running状态")
                break
            elif status in ["failed", "cancelled", "completed"]:
                print(f"   ⚠️ 任务已进入最终状态: {status}")
                return 1
        time.sleep(2)
        waited += 2
    
    if waited >= max_wait:
        print("   ⚠️ 超时：任务未在60秒内进入running状态")
        return 1
    
    # 显示测试说明
    print("\n" + "="*70)
    print("测试任务已准备就绪！")
    print("="*70)
    print(f"\n任务ID: {task_id}")
    print(f"任务状态: running")
    print("\n现在可以在前端测试停止功能：")
    print("\n1. 打开浏览器访问前端页面（通常是 http://localhost:8501 或 Streamlit 配置的端口）")
    print("2. 进入'训练管理'页面")
    print("3. 在'任务列表'标签页中找到这个任务")
    print("4. 点击'停止任务'按钮")
    print("5. 观察：")
    print("   - 应该看到'任务停止请求已发送'的成功消息")
    print("   - 应该看到进度条显示轮询进度")
    print("   - 应该看到状态实时更新（从running变为cancelled）")
    print("   - 整个过程应该在30秒内完成")
    print("\n如果遇到问题，请检查：")
    print("  - 前端页面是否正确加载")
    print("  - 网络连接是否正常")
    print("  - 后端服务是否仍在运行")
    print("\n" + "="*70)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

