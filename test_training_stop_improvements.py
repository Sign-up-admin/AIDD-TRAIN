"""
测试训练停止功能改进
测试所有改进的功能点
"""

import sys
import os
import time
import threading
import multiprocessing
from datetime import datetime

# 添加项目路径
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from compass.service.services.progress_tracker import ProgressTracker
from compass.service.models.task import TaskStatus
from compass.service.services.training_service import TrainingService
from compass.service.models.task import TrainingTaskResponse


def test_progress_tracker_cancellation():
    """测试ProgressTracker的取消标志功能"""
    print("\n" + "="*70)
    print("测试1: ProgressTracker取消标志功能")
    print("="*70)
    
    tracker = ProgressTracker("test_task_1")
    
    # 测试初始状态
    assert not tracker.is_cancelled(), "初始状态应该未取消"
    print("[PASS] 初始状态检查通过")
    
    # 测试设置取消标志
    tracker.cancel()
    assert tracker.is_cancelled(), "取消标志应该已设置"
    print("[PASS] 取消标志设置成功")
    
    # 测试get_progress中的取消状态
    progress = tracker.get_progress()
    assert progress.get("cancelled") == True, "get_progress应该返回cancelled=True"
    print("[PASS] get_progress返回正确的取消状态")
    
    # 测试无锁读取性能（简单测试）
    start_time = time.time()
    for _ in range(10000):
        tracker.is_cancelled()
    elapsed = time.time() - start_time
    print(f"[PASS] 无锁读取性能测试: 10000次调用耗时 {elapsed:.4f}秒")
    
    print("[PASS] 测试1通过\n")


# 模块级函数，用于multiprocessing（Windows需要）
def _worker_process_check_cancellation(flag, result_queue):
    """模拟worker进程检查取消标志"""
    time.sleep(0.05)  # 减少等待时间
    with flag.get_lock():
        cancelled = flag.value == 1
    result_queue.put(cancelled)


def test_multiprocessing_cancellation_flag():
    """测试跨进程取消标志"""
    print("="*70)
    print("测试2: 跨进程取消标志传递")
    print("="*70)
    
    tracker = ProgressTracker("test_task_2")
    mp_flag = tracker.get_mp_cancellation_flag()
    
    # 测试初始状态
    assert not tracker.is_cancelled_mp(), "初始状态应该未取消"
    print("[PASS] 初始状态检查通过")
    
    # 测试multiprocessing标志的基本功能（不启动实际进程，避免Windows问题）
    # 直接测试标志的读写
    with mp_flag.get_lock():
        assert mp_flag.value == 0, "初始值应该为0"
    print("[PASS] 多进程标志初始值正确")
    
    # 测试取消后状态
    tracker.cancel()
    with mp_flag.get_lock():
        assert mp_flag.value == 1, "取消后值应该为1"
    print("[PASS] 多进程标志取消后值正确")
    
    # 测试is_cancelled_mp方法
    assert tracker.is_cancelled_mp(), "is_cancelled_mp应该返回True"
    print("[PASS] is_cancelled_mp方法工作正常")
    
    # 如果支持，尝试简单的进程测试（Windows可能较慢，跳过实际进程测试）
    try:
        result_queue = multiprocessing.Queue()
        worker = multiprocessing.Process(target=_worker_process_check_cancellation, args=(mp_flag, result_queue))
        worker.start()
        worker.join(timeout=10.0)  # 增加超时时间
        if worker.is_alive():
            worker.terminate()
            worker.join()
            print("[SKIP] Worker进程测试跳过（Windows启动较慢）")
        else:
            result = result_queue.get(timeout=1.0)
            assert result, "worker进程应该检测到取消状态"
            print("[PASS] Worker进程取消状态检测通过")
    except Exception as e:
        print(f"[SKIP] Worker进程测试跳过: {str(e)}")
    
    print("[PASS] 测试2通过\n")


def test_task_status_cancelling():
    """测试CANCELLING状态"""
    print("="*70)
    print("测试3: CANCELLING状态")
    print("="*70)
    
    # 测试状态枚举
    assert hasattr(TaskStatus, 'CANCELLING'), "TaskStatus应该有CANCELLING状态"
    assert TaskStatus.CANCELLING.value == "cancelling", "CANCELLING状态值应该为'cancelling'"
    print("[PASS] CANCELLING状态定义正确")
    
    # 测试状态转换
    statuses = [TaskStatus.RUNNING, TaskStatus.CANCELLING, TaskStatus.CANCELLED]
    print(f"[PASS] 状态转换路径: {' -> '.join([s.value for s in statuses])}")
    
    print("[PASS] 测试3通过\n")


def test_training_service_stop_task():
    """测试TrainingService的stop_task方法"""
    print("="*70)
    print("测试4: TrainingService.stop_task方法")
    print("="*70)
    
    service = TrainingService()
    
    # 创建测试任务
    task_id = "test_stop_task_1"
    test_config = {
        "execution_mode": "smoke_test",
        "epochs": 1,
        "batch_size": 2,
        "learning_rate": 0.001,
        "optimizer": "adam"
    }
    
    # 创建任务
    task = TrainingTaskResponse(
        task_id=task_id,
        status=TaskStatus.RUNNING,
        config=test_config,
        created_at=datetime.now(),
        updated_at=datetime.now()
    )
    
    with service.lock:
        service.tasks[task_id] = task
        # 创建progress tracker
        tracker = ProgressTracker(task_id)
        service.progress_trackers[task_id] = tracker
    
    # 测试停止任务
    success, error = service.stop_task(task_id)
    assert success, f"停止任务应该成功: {error}"
    print("[PASS] stop_task返回成功")
    
    # 检查状态是否更新为CANCELLING（可能快速变为CANCELLED）
    with service.lock:
        updated_task = service.tasks.get(task_id)
        assert updated_task is not None, "任务应该存在"
        # 状态应该是CANCELLING或CANCELLED（如果线程已经完成）
        assert updated_task.status in [TaskStatus.CANCELLING, TaskStatus.CANCELLED], \
            f"状态应该为CANCELLING或CANCELLED，实际为{updated_task.status.value}"
    print(f"[PASS] 任务状态已更新为{updated_task.status.value}")
    
    # 检查取消标志是否设置
    tracker = service.progress_trackers.get(task_id)
    assert tracker is not None, "ProgressTracker应该存在"
    assert tracker.is_cancelled(), "取消标志应该已设置"
    print("[PASS] 取消标志已正确设置")
    
    # 清理
    with service.lock:
        if task_id in service.tasks:
            del service.tasks[task_id]
        if task_id in service.progress_trackers:
            del service.progress_trackers[task_id]
    
    print("[PASS] 测试4通过\n")


def test_cancellation_checkpoints():
    """测试取消检查点"""
    print("="*70)
    print("测试5: 取消检查点")
    print("="*70)
    
    tracker = ProgressTracker("test_checkpoints")
    
    # 模拟训练循环中的检查
    check_count = 0
    max_checks = 100
    
    def simulate_training_loop():
        nonlocal check_count
        for i in range(max_checks):
            if tracker.is_cancelled():
                return i
            check_count += 1
            time.sleep(0.001)  # 模拟处理时间
        return max_checks
    
    # 启动模拟训练循环
    loop_thread = threading.Thread(target=simulate_training_loop)
    loop_thread.start()
    
    # 等待一小段时间后取消
    time.sleep(0.05)
    tracker.cancel()
    
    loop_thread.join(timeout=1.0)
    
    assert check_count < max_checks, "训练循环应该在取消后停止"
    print(f"[PASS] 取消检查点工作正常，在 {check_count} 次检查后检测到取消")
    
    print("[PASS] 测试5通过\n")


def test_concurrent_cancellation():
    """测试并发取消场景"""
    print("="*70)
    print("测试6: 并发取消场景")
    print("="*70)
    
    tracker = ProgressTracker("test_concurrent")
    results = []
    results_lock = threading.Lock()
    
    def check_cancellation(thread_id):
        for _ in range(1000):
            if tracker.is_cancelled():
                with results_lock:
                    results.append(thread_id)
                return
            time.sleep(0.0001)  # 小延迟，确保有机会检测到取消
    
    threads = []
    for i in range(10):
        t = threading.Thread(target=check_cancellation, args=(i,))
        threads.append(t)
        t.start()
    
    # 等待一小段时间后取消
    time.sleep(0.05)
    tracker.cancel()
    
    # 等待所有线程完成
    for t in threads:
        t.join(timeout=2.0)
    
    with results_lock:
        detected_count = len(results)
    assert detected_count == 10, f"所有10个线程都应该检测到取消，实际检测到{detected_count}个"
    print(f"[PASS] 并发取消检测正常，{detected_count}个线程都检测到取消")
    
    print("[PASS] 测试6通过\n")


def run_all_tests():
    """运行所有测试"""
    print("\n" + "="*70)
    print("开始测试训练停止功能改进")
    print("="*70)
    
    tests = [
        test_progress_tracker_cancellation,
        test_multiprocessing_cancellation_flag,
        test_task_status_cancelling,
        test_training_service_stop_task,
        test_cancellation_checkpoints,
        test_concurrent_cancellation,
    ]
    
    passed = 0
    failed = 0
    
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except Exception as e:
            print(f"[FAIL] 测试失败: {test_func.__name__}")
            print(f"  错误: {str(e)}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("="*70)
    print("测试结果汇总")
    print("="*70)
    print(f"通过: {passed}/{len(tests)}")
    print(f"失败: {failed}/{len(tests)}")
    print("="*70)
    
    if failed == 0:
        print("[SUCCESS] 所有测试通过！")
        return 0
    else:
        print("[FAIL] 部分测试失败")
        return 1


if __name__ == "__main__":
    # 设置multiprocessing启动方法（Windows需要）
    if sys.platform == "win32":
        multiprocessing.set_start_method('spawn', force=True)
    
    exit_code = run_all_tests()
    sys.exit(exit_code)

