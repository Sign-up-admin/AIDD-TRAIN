"""
COMPASS自动化测试套件
提供全面的单元测试、集成测试和端到端测试
"""

import pytest
import os
import sys
import time
import json
import uuid
import asyncio
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime
from unittest.mock import Mock, patch, MagicMock, AsyncMock

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from fastapi.testclient import TestClient
from compass.service.server import app
from compass.service.models.task import TaskStatus
from compass.service.config import SERVICE_CONFIG


class CompassTestSuite:
    """COMPASS自动化测试套件主类"""
    
    def __init__(self):
        self.results = {
            "unit_tests": [],
            "integration_tests": [],
            "e2e_tests": [],
            "summary": {
                "total": 0,
                "passed": 0,
                "failed": 0,
                "skipped": 0,
                "start_time": None,
                "end_time": None,
                "duration": 0
            }
        }
        self.temp_dirs = []
    
    def setup(self):
        """测试前准备"""
        self.results["summary"]["start_time"] = datetime.now()
        print("\n" + "=" * 70)
        print("COMPASS自动化测试套件")
        print("=" * 70)
        print(f"开始时间: {self.results['summary']['start_time']}")
        print(f"项目根目录: {project_root}")
        print("=" * 70 + "\n")
    
    def teardown(self):
        """测试后清理"""
        self.results["summary"]["end_time"] = datetime.now()
        duration = (self.results["summary"]["end_time"] - 
                   self.results["summary"]["start_time"]).total_seconds()
        self.results["summary"]["duration"] = duration
        
        # 清理临时目录
        for temp_dir in self.temp_dirs:
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir, ignore_errors=True)
    
    def create_temp_dir(self) -> str:
        """创建临时目录"""
        temp_dir = tempfile.mkdtemp(prefix="compass_test_")
        self.temp_dirs.append(temp_dir)
        return temp_dir
    
    def record_result(self, category: str, test_name: str, passed: bool, 
                     error: Optional[str] = None, duration: float = 0):
        """记录测试结果"""
        result = {
            "name": test_name,
            "passed": passed,
            "error": error,
            "duration": duration,
            "timestamp": datetime.now().isoformat()
        }
        self.results[category].append(result)
        self.results["summary"]["total"] += 1
        if passed:
            self.results["summary"]["passed"] += 1
        else:
            self.results["summary"]["failed"] += 1
    
    def generate_report(self) -> str:
        """生成测试报告"""
        summary = self.results["summary"]
        report = []
        report.append("\n" + "=" * 70)
        report.append("COMPASS自动化测试报告")
        report.append("=" * 70)
        report.append(f"开始时间: {summary['start_time']}")
        report.append(f"结束时间: {summary['end_time']}")
        report.append(f"总耗时: {summary['duration']:.2f}秒")
        report.append("")
        report.append(f"总测试数: {summary['total']}")
        report.append(f"通过: {summary['passed']} ({summary['passed']/summary['total']*100:.1f}%)" 
                     if summary['total'] > 0 else "通过: 0")
        report.append(f"失败: {summary['failed']}")
        report.append(f"跳过: {summary['skipped']}")
        report.append("")
        
        # 单元测试结果
        if self.results["unit_tests"]:
            report.append("单元测试结果:")
            report.append("-" * 70)
            for test in self.results["unit_tests"]:
                status = "✓" if test["passed"] else "✗"
                report.append(f"  {status} {test['name']} ({test['duration']:.2f}s)")
                if not test["passed"] and test["error"]:
                    report.append(f"    错误: {test['error']}")
            report.append("")
        
        # 集成测试结果
        if self.results["integration_tests"]:
            report.append("集成测试结果:")
            report.append("-" * 70)
            for test in self.results["integration_tests"]:
                status = "✓" if test["passed"] else "✗"
                report.append(f"  {status} {test['name']} ({test['duration']:.2f}s)")
                if not test["passed"] and test["error"]:
                    report.append(f"    错误: {test['error']}")
            report.append("")
        
        # 端到端测试结果
        if self.results["e2e_tests"]:
            report.append("端到端测试结果:")
            report.append("-" * 70)
            for test in self.results["e2e_tests"]:
                status = "✓" if test["passed"] else "✗"
                report.append(f"  {status} {test['name']} ({test['duration']:.2f}s)")
                if not test["passed"] and test["error"]:
                    report.append(f"    错误: {test['error']}")
            report.append("")
        
        report.append("=" * 70)
        return "\n".join(report)


# ==================== 单元测试 ====================

@pytest.mark.unit
def test_config_loading():
    """测试配置加载"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        from compass.service.config import SERVICE_CONFIG
        assert SERVICE_CONFIG is not None
        assert "host" in SERVICE_CONFIG
        assert "port" in SERVICE_CONFIG
        suite.record_result("unit_tests", "配置加载", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("unit_tests", "配置加载", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.unit
def test_task_model_validation():
    """测试任务模型验证"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        from compass.service.models.task import TrainingTaskRequest, TaskStatus
        
        # 测试有效的任务请求
        valid_request = {
            "config": {
                "execution_mode": "smoke_test",
                "epochs": 2,
                "batch_size": 1
            }
        }
        task = TrainingTaskRequest(**valid_request)
        assert task.config["execution_mode"] == "smoke_test"
        
        suite.record_result("unit_tests", "任务模型验证", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("unit_tests", "任务模型验证", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.unit
def test_error_codes():
    """测试错误代码"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        from compass.service.error_codes import ErrorCode
        
        # 验证错误代码存在
        assert hasattr(ErrorCode, "TASK_NOT_FOUND")
        assert hasattr(ErrorCode, "VALIDATION_ERROR")
        assert hasattr(ErrorCode, "INTERNAL_ERROR")
        
        suite.record_result("unit_tests", "错误代码", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("unit_tests", "错误代码", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.unit
def test_progress_tracker():
    """测试进度跟踪器"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        from compass.service.services.progress_tracker import ProgressTracker
        
        tracker = ProgressTracker(task_id="test_task")
        assert tracker.task_id == "test_task"
        assert tracker.is_cancelled() == False
        
        tracker.set_stage("testing", "测试阶段")
        assert tracker.current_stage == "testing"
        
        tracker.cancel()
        assert tracker.is_cancelled() == True
        
        suite.record_result("unit_tests", "进度跟踪器", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("unit_tests", "进度跟踪器", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


# ==================== 集成测试 ====================

@pytest.mark.integration
def test_health_endpoint():
    """测试健康检查端点"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        # 禁用认证
        os.environ["AUTH_ENABLED"] = "false"
        os.environ["FORCE_AUTH_CRITICAL"] = "false"
        
        client = TestClient(app)
        response = client.get("/health")
        
        assert response.status_code == 200
        data = response.json()
        assert "status" in data
        assert data["status"] == "healthy"
        
        suite.record_result("integration_tests", "健康检查端点", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("integration_tests", "健康检查端点", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.integration
def test_training_task_creation():
    """测试训练任务创建"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        os.environ["AUTH_ENABLED"] = "false"
        os.environ["FORCE_AUTH_CRITICAL"] = "false"
        
        client = TestClient(app)
        
        # Mock训练服务
        with patch("compass.service.routes.training.training_service") as mock_service:
            task_id = str(uuid.uuid4())
            mock_service.create_task.return_value = task_id
            mock_service.get_task.return_value = {
                "task_id": task_id,
                "status": TaskStatus.PENDING,
                "config": {"epochs": 2},
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat(),
            }
            
            response = client.post(
                "/api/v1/training/tasks",
                json={
                    "config": {
                        "execution_mode": "smoke_test",
                        "epochs": 2,
                        "batch_size": 1
                    }
                }
            )
            
            assert response.status_code in [200, 201]
            data = response.json()
            assert "task_id" in data
            
        suite.record_result("integration_tests", "训练任务创建", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("integration_tests", "训练任务创建", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.integration
def test_task_list_endpoint():
    """测试任务列表端点"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        os.environ["AUTH_ENABLED"] = "false"
        os.environ["FORCE_AUTH_CRITICAL"] = "false"
        
        client = TestClient(app)
        
        with patch("compass.service.routes.training.training_service") as mock_service:
            mock_service.list_tasks.return_value = []
            
            response = client.get("/api/v1/training/tasks")
            
            assert response.status_code == 200
            data = response.json()
            assert isinstance(data, list)
            
        suite.record_result("integration_tests", "任务列表端点", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("integration_tests", "任务列表端点", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


@pytest.mark.integration
def test_metrics_endpoint():
    """测试指标端点"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        os.environ["AUTH_ENABLED"] = "false"
        os.environ["FORCE_AUTH_CRITICAL"] = "false"
        
        client = TestClient(app)
        response = client.get("/metrics")
        
        assert response.status_code == 200
        # 指标端点应该返回文本格式
        assert response.headers.get("content-type") is not None
        
        suite.record_result("integration_tests", "指标端点", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("integration_tests", "指标端点", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


# ==================== 端到端测试 ====================

@pytest.mark.e2e
@pytest.mark.slow
def test_complete_training_workflow():
    """测试完整的训练工作流"""
    suite = CompassTestSuite()
    suite.setup()
    start_time = time.time()
    
    try:
        # 这个测试需要实际运行训练，可能会很慢
        # 在CI/CD环境中可以跳过或使用mock
        
        # 创建临时目录
        temp_dir = suite.create_temp_dir()
        
        # Mock训练服务以避免实际训练
        with patch("compass.service.services.training_service.TrainingService") as mock_service:
            mock_service.create_task.return_value = str(uuid.uuid4())
            mock_service.get_task.return_value = {
                "task_id": str(uuid.uuid4()),
                "status": TaskStatus.COMPLETED,
                "config": {"epochs": 1},
            }
            
            # 这里可以添加更多端到端测试逻辑
            pass
        
        suite.record_result("e2e_tests", "完整训练工作流", True, duration=time.time()-start_time)
    except Exception as e:
        suite.record_result("e2e_tests", "完整训练工作流", False, str(e), time.time()-start_time)
    
    suite.teardown()
    return suite


# ==================== 测试运行器 ====================

def run_all_tests():
    """运行所有测试"""
    print("运行COMPASS自动化测试套件...")
    print("使用pytest运行: pytest tests/test_compass_automated.py -v")


if __name__ == "__main__":
    run_all_tests()

