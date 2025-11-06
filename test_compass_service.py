"""
独立测试COMPASS服务功能的脚本
用于诊断COMPASS服务是否正常工作
"""
import sys
import requests
import time
from pathlib import Path
import io
from typing import Dict, Optional

# Fix Windows console encoding
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

COMPASS_BASE_URL = "http://localhost:8080"
TIMEOUT = 30


def print_header(text: str):
    """打印标题"""
    print("\n" + "=" * 60)
    print(text)
    print("=" * 60)


def print_result(success: bool, message: str):
    """打印测试结果"""
    status = "[OK]" if success else "[FAIL]"
    print(f"{status} {message}")


def test_health_check() -> bool:
    """测试健康检查"""
    print_header("1. 测试健康检查 (/health)")
    try:
        response = requests.get(f"{COMPASS_BASE_URL}/health", timeout=5)
        if response.status_code == 200:
            data = response.json()
            print_result(True, f"健康检查成功: {data.get('status')}")
            print(f"  时间戳: {data.get('timestamp')}")
            return True
        else:
            print_result(False, f"健康检查失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print_result(False, f"健康检查失败: {e}")
        return False


def test_readiness_check() -> bool:
    """测试就绪检查"""
    print_header("2. 测试就绪检查 (/health/ready)")
    try:
        response = requests.get(f"{COMPASS_BASE_URL}/health/ready", timeout=5)
        if response.status_code == 200:
            data = response.json()
            print_result(True, f"就绪检查成功: {data.get('status')}")
            return True
        else:
            print_result(False, f"就绪检查失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print_result(False, f"就绪检查失败: {e}")
        return False


def test_create_task() -> Optional[str]:
    """测试创建任务"""
    print_header("3. 测试创建任务 (POST /api/v1/training/tasks)")
    
    config = {
        "execution_mode": "smoke_test",  # 使用快速测试模式
        "epochs": 2,
        "batch_size": 1,
        "learning_rate": 0.0001
    }
    
    data = {
        "config": config,
        "description": "测试任务 - COMPASS服务诊断"
    }
    
    try:
        print(f"发送创建任务请求...")
        print(f"  配置: {config}")
        response = requests.post(
            f"{COMPASS_BASE_URL}/api/v1/training/tasks",
            json=data,
            timeout=TIMEOUT
        )
        
        if response.status_code == 201:
            result = response.json()
            task_id = result.get('task_id')
            status = result.get('status')
            print_result(True, f"任务创建成功")
            print(f"  任务ID: {task_id}")
            print(f"  状态: {status}")
            return task_id
        else:
            print_result(False, f"任务创建失败: HTTP {response.status_code}")
            try:
                error_data = response.json()
                print(f"  错误详情: {error_data}")
            except:
                print(f"  响应内容: {response.text}")
            return None
    except Exception as e:
        print_result(False, f"任务创建失败: {e}")
        return None


def test_get_task_status(task_id: str) -> Optional[Dict]:
    """测试获取任务状态"""
    try:
        response = requests.get(
            f"{COMPASS_BASE_URL}/api/v1/training/tasks/{task_id}",
            timeout=10
        )
        if response.status_code == 200:
            return response.json()
        else:
            print(f"  获取状态失败: HTTP {response.status_code}")
            return None
    except Exception as e:
        print(f"  获取状态失败: {e}")
        return None


def test_start_task(task_id: str) -> bool:
    """测试启动任务"""
    print_header("4. 测试启动任务 (POST /api/v1/training/tasks/{task_id}/start)")
    
    try:
        # 先获取初始状态
        initial_status = None
        task_status = test_get_task_status(task_id)
        if task_status:
            initial_status = task_status.get('status')
            print(f"初始状态: {initial_status}")
        
        # 发送启动请求
        print(f"发送启动任务请求...")
        start_time = time.time()
        response = requests.post(
            f"{COMPASS_BASE_URL}/api/v1/training/tasks/{task_id}/start",
            timeout=TIMEOUT
        )
        
        if response.status_code == 200:
            result = response.json()
            print_result(True, f"启动命令发送成功")
            print(f"  响应: {result.get('message')}")
            
            # 监控状态变化
            print(f"\n监控任务状态变化...")
            max_checks = 10
            check_interval = 0.5
            status_changed = False
            
            for i in range(max_checks):
                time.sleep(check_interval)
                task_status = test_get_task_status(task_id)
                if task_status:
                    current_status = task_status.get('status')
                    elapsed = time.time() - start_time
                    print(f"  检查 {i+1}/{max_checks}: 状态={current_status} (耗时 {elapsed:.1f}s)")
                    
                    # 检查状态是否变化
                    if initial_status:
                        if current_status != initial_status:
                            status_changed = True
                            print_result(True, f"状态已变化: {initial_status} -> {current_status}")
                            if current_status in ['initializing', 'running']:
                                print_result(True, "任务已成功启动")
                                return True
                            elif current_status == 'failed':
                                error = task_status.get('error', 'Unknown error')
                                print_result(False, f"任务启动失败: {error}")
                                return False
                    
                    # 如果已经是运行状态，直接返回成功
                    if current_status in ['initializing', 'running']:
                        print_result(True, "任务正在运行")
                        return True
                    elif current_status == 'failed':
                        error = task_status.get('error', 'Unknown error')
                        print_result(False, f"任务失败: {error}")
                        return False
                else:
                    print(f"  检查 {i+1}/{max_checks}: 无法获取状态")
            
            if not status_changed:
                print_result(False, f"状态未在 {max_checks * check_interval} 秒内变化")
                print(f"  初始状态: {initial_status}")
                print(f"  当前状态: {task_status.get('status') if task_status else 'Unknown'}")
                return False
            
            return True
        else:
            print_result(False, f"启动任务失败: HTTP {response.status_code}")
            try:
                error_data = response.json()
                print(f"  错误详情: {error_data}")
            except:
                print(f"  响应内容: {response.text}")
            return False
    except Exception as e:
        print_result(False, f"启动任务失败: {e}")
        return False


def test_task_logs(task_id: str):
    """测试获取任务日志"""
    print_header("5. 测试获取任务日志 (GET /api/v1/training/tasks/{task_id}/logs)")
    
    try:
        response = requests.get(
            f"{COMPASS_BASE_URL}/api/v1/training/tasks/{task_id}/logs",
            params={"limit": 20},
            timeout=10
        )
        if response.status_code == 200:
            result = response.json()
            logs = result.get('logs', [])
            print_result(True, f"获取日志成功，共 {len(logs)} 条")
            if logs:
                print("\n最近日志:")
                for log in logs[-5:]:  # 只显示最后5条
                    print(f"  {log}")
            return True
        else:
            print_result(False, f"获取日志失败: HTTP {response.status_code}")
            return False
    except Exception as e:
        print_result(False, f"获取日志失败: {e}")
        return False


def test_task_progress(task_id: str):
    """测试获取任务进度"""
    print_header("6. 测试获取任务进度")
    
    task_status = test_get_task_status(task_id)
    if task_status:
        progress = task_status.get('progress', {})
        if progress:
            stage = progress.get('stage', 'unknown')
            progress_pct = progress.get('progress', 0.0)
            message = progress.get('message', '')
            print_result(True, f"获取进度成功")
            print(f"  阶段: {stage}")
            print(f"  进度: {progress_pct * 100:.1f}%")
            print(f"  消息: {message}")
            
            # 显示详细进度信息
            if stage == 'data_processing':
                dp_info = progress.get('data_processing', {})
                if dp_info:
                    completed = dp_info.get('completed', 0)
                    total = dp_info.get('total', 0)
                    percentage = dp_info.get('percentage', 0)
                    print(f"  数据处理: {completed}/{total} ({percentage:.1f}%)")
            elif stage == 'training':
                train_info = progress.get('training', {})
                if train_info:
                    epoch = train_info.get('current_epoch', 0)
                    total_epochs = train_info.get('total_epochs', 0)
                    print(f"  训练: Epoch {epoch}/{total_epochs}")
            
            return True
        else:
            print("  暂无进度信息")
            return False
    else:
        print_result(False, "无法获取任务状态")
        return False


def monitor_task_for_processing(task_id: str, max_wait: int = 30):
    """监控任务，等待数据处理开始"""
    print_header("7. 监控任务 - 等待数据处理开始")
    
    print(f"监控任务 {max_wait} 秒，等待数据处理开始...")
    start_time = time.time()
    check_interval = 2
    
    while (time.time() - start_time) < max_wait:
        task_status = test_get_task_status(task_id)
        if task_status:
            status = task_status.get('status')
            progress = task_status.get('progress', {})
            stage = progress.get('stage', 'unknown')
            elapsed = time.time() - start_time
            
            print(f"  时间: {elapsed:.1f}s | 状态: {status} | 阶段: {stage}")
            
            # 检查是否开始数据处理
            if stage == 'data_processing':
                print_result(True, "数据处理已开始！")
                dp_info = progress.get('data_processing', {})
                if dp_info:
                    completed = dp_info.get('completed', 0)
                    total = dp_info.get('total', 0)
                    print(f"  已处理: {completed}/{total}")
                return True
            
            # 检查是否失败
            if status == 'failed':
                error = task_status.get('error', 'Unknown error')
                print_result(False, f"任务失败: {error}")
                return False
            
            # 检查是否完成（可能很快完成）
            if status == 'completed':
                print_result(True, "任务已完成")
                return True
        
        time.sleep(check_interval)
    
    print_result(False, f"在 {max_wait} 秒内未检测到数据处理开始")
    return False


def main():
    """主测试函数"""
    print_header("COMPASS服务独立测试")
    print(f"测试目标: {COMPASS_BASE_URL}")
    print(f"超时设置: {TIMEOUT}秒")
    
    results = {
        "health_check": False,
        "readiness_check": False,
        "create_task": False,
        "start_task": False,
        "task_logs": False,
        "task_progress": False,
        "data_processing_started": False
    }
    
    task_id = None
    
    # 1. 健康检查
    results["health_check"] = test_health_check()
    if not results["health_check"]:
        print("\n[错误] 健康检查失败，无法继续测试")
        return 1
    
    # 2. 就绪检查
    results["readiness_check"] = test_readiness_check()
    if not results["readiness_check"]:
        print("\n[警告] 就绪检查失败，但继续测试...")
    
    # 3. 创建任务
    task_id = test_create_task()
    if task_id:
        results["create_task"] = True
    else:
        print("\n[错误] 任务创建失败，无法继续测试")
        return 1
    
    # 4. 启动任务
    if task_id:
        results["start_task"] = test_start_task(task_id)
    
    # 5. 获取日志
    if task_id:
        results["task_logs"] = test_task_logs(task_id)
    
    # 6. 获取进度
    if task_id:
        results["task_progress"] = test_task_progress(task_id)
    
    # 7. 监控数据处理
    if task_id and results["start_task"]:
        results["data_processing_started"] = monitor_task_for_processing(task_id, max_wait=30)
    
    # 总结
    print_header("测试总结")
    print("\n测试结果:")
    for test_name, result in results.items():
        status = "✓" if result else "✗"
        print(f"  {status} {test_name}")
    
    passed = sum(1 for r in results.values() if r)
    total = len(results)
    print(f"\n通过: {passed}/{total}")
    
    if passed == total:
        print("\n[成功] 所有测试通过！COMPASS服务工作正常。")
        return 0
    else:
        print("\n[警告] 部分测试失败，请检查COMPASS服务。")
        return 1


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n[中断] 测试被用户中断")
        sys.exit(130)
    except Exception as e:
        print(f"\n\n[错误] 测试过程中发生异常: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

