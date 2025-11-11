"""
详细测试训练停止功能 - 检查数据加载阶段的阻塞问题
"""
import sys
import time
import requests
import json
from typing import Dict, Optional, Any
from datetime import datetime

COMPASS_URL = "http://localhost:8080"


def get_task_logs(task_id: str) -> Optional[list]:
    """获取任务日志"""
    try:
        response = requests.get(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}/logs",
            timeout=10
        )
        if response.status_code == 200:
            logs = response.json()
            return logs.get("logs", [])
        return None
    except Exception as e:
        print(f"[ERROR] 获取日志失败: {e}")
        return None


def get_task_status(task_id: str) -> Optional[Dict[str, Any]]:
    """获取任务状态"""
    try:
        response = requests.get(
            f"{COMPASS_URL}/api/v1/training/tasks/{task_id}",
            timeout=10
        )
        if response.status_code == 200:
            return response.json()
        return None
    except Exception as e:
        print(f"[ERROR] 获取状态失败: {e}")
        return None


def analyze_stop_timing(task_id: str):
    """分析停止请求的时序"""
    print("\n" + "="*70)
    print("分析停止请求时序")
    print("="*70)
    
    # 获取日志
    logs = get_task_logs(task_id)
    if logs:
        print(f"\n获取到 {len(logs)} 条日志")
        print("\n关键日志（包含取消相关）:")
        for log in logs[-50:]:  # 只看最后50条
            log_text = log.get("message", "") if isinstance(log, dict) else str(log)
            if any(keyword in log_text.lower() for keyword in ["cancel", "stop", "cancellation", "停止", "取消"]):
                timestamp = log.get("timestamp", "") if isinstance(log, dict) else ""
                print(f"  [{timestamp}] {log_text}")
    
    # 获取当前状态
    task = get_task_status(task_id)
    if task:
        print(f"\n当前状态: {task.get('status')}")
        progress = task.get("progress", {})
        print(f"阶段: {progress.get('stage')}")
        print(f"取消标志: {progress.get('cancelled')}")
        print(f"最后更新: {progress.get('last_update')}")


def main():
    """主函数"""
    if len(sys.argv) < 2:
        print("用法: python test_stop_detailed.py <task_id>")
        print("\n示例:")
        print("  python test_stop_detailed.py 705be14a-0e5a-4ee7-a287-9be0eb8601e3")
        return 1
    
    task_id = sys.argv[1]
    print(f"分析任务: {task_id}")
    analyze_stop_timing(task_id)
    return 0


if __name__ == "__main__":
    sys.exit(main())







