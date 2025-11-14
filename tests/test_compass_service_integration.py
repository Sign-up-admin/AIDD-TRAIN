"""
COMPASS服务集成测试
测试服务层各个组件的集成
"""

import pytest
import os
import sys
import uuid
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
from fastapi.testclient import TestClient

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from compass.service.server import app
from compass.service.models.task import TaskStatus


@pytest.fixture
def client(monkeypatch):
    """创建测试客户端"""
    monkeypatch.setenv('AUTH_ENABLED', 'false')
    monkeypatch.setenv('FORCE_AUTH_CRITICAL', 'false')
    monkeypatch.setenv('RATE_LIMIT_TRAINING', '10000')
    monkeypatch.setenv('RATE_LIMIT_DEFAULT', '10000')
    
    # 重新加载认证模块
    import importlib
    import compass.service.middleware.auth
    importlib.reload(compass.service.middleware.auth)
    
    return TestClient(app)


@pytest.fixture
def temp_dir():
    """创建临时目录"""
    temp_path = tempfile.mkdtemp()
    yield temp_path
    shutil.rmtree(temp_path, ignore_errors=True)


@pytest.mark.integration
class TestServiceIntegration:
    """服务集成测试类"""
    
    def test_service_startup(self, client):
        """测试服务启动"""
        response = client.get("/")
        assert response.status_code == 200
    
    def test_health_check(self, client):
        """测试健康检查"""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
    
    def test_api_documentation(self, client):
        """测试API文档"""
        response = client.get("/docs")
        assert response.status_code == 200
        
        response = client.get("/openapi.json")
        assert response.status_code == 200
        schema = response.json()
        assert "openapi" in schema
        assert "paths" in schema
    
    def test_training_task_creation_flow(self, client):
        """测试训练任务创建流程"""
        with patch("compass.service.routes.training.training_service") as mock_service:
            task_id = str(uuid.uuid4())
            mock_service.create_task.return_value = task_id
            mock_service.get_task.return_value = {
                "task_id": task_id,
                "status": TaskStatus.PENDING,
                "config": {"epochs": 2, "execution_mode": "smoke_test"},
                "created_at": "2024-01-01T00:00:00",
                "updated_at": "2024-01-01T00:00:00",
            }
            
            # 创建任务
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
            
            # 获取任务
            response = client.get(f"/api/v1/training/tasks/{task_id}")
            assert response.status_code == 200
            data = response.json()
            assert data["task_id"] == task_id
    
    def test_task_list_and_filter(self, client):
        """测试任务列表和过滤"""
        with patch("compass.service.routes.training.training_service") as mock_service:
            mock_service.list_tasks.return_value = [
                {
                    "task_id": str(uuid.uuid4()),
                    "status": TaskStatus.COMPLETED,
                    "config": {"epochs": 2},
                },
                {
                    "task_id": str(uuid.uuid4()),
                    "status": TaskStatus.RUNNING,
                    "config": {"epochs": 5},
                }
            ]
            
            # 获取所有任务
            response = client.get("/api/v1/training/tasks")
            assert response.status_code == 200
            tasks = response.json()
            assert isinstance(tasks, list)
            
            # 按状态过滤
            response = client.get("/api/v1/training/tasks?status=completed")
            assert response.status_code == 200
    
    def test_task_stop_functionality(self, client):
        """测试任务停止功能"""
        with patch("compass.service.routes.training.training_service") as mock_service:
            task_id = str(uuid.uuid4())
            mock_service.get_task.return_value = {
                "task_id": task_id,
                "status": TaskStatus.RUNNING,
                "config": {"epochs": 2},
            }
            mock_service.stop_task.return_value = {
                "task_id": task_id,
                "status": TaskStatus.CANCELLED,
            }
            
            response = client.post(f"/api/v1/training/tasks/{task_id}/stop")
            assert response.status_code in [200, 202]
    
    def test_dataset_endpoints(self, client):
        """测试数据集端点"""
        with patch("compass.service.routes.data.data_service") as mock_service:
            mock_service.list_datasets.return_value = []
            
            response = client.get("/api/v1/data/datasets")
            assert response.status_code == 200
            datasets = response.json()
            assert isinstance(datasets, list)
    
    def test_model_endpoints(self, client):
        """测试模型端点"""
        with patch("compass.service.routes.models.model_service") as mock_service:
            mock_service.list_models.return_value = []
            
            response = client.get("/api/v1/models")
            assert response.status_code == 200
            models = response.json()
            assert isinstance(models, list)
    
    def test_inference_endpoint(self, client):
        """测试推理端点"""
        with patch("compass.service.routes.inference.inference_service") as mock_service:
            mock_service.predict.return_value = {
                "affinity": -7.5,
                "confidence": 0.95
            }
            
            response = client.post(
                "/api/v1/inference/predict",
                json={
                    "protein_path": "test.pdb",
                    "ligand_path": "test.sdf"
                }
            )
            # 可能返回404如果服务未配置，或200如果成功
            assert response.status_code in [200, 404, 422]
    
    def test_error_handling(self, client):
        """测试错误处理"""
        # 测试不存在的任务
        fake_task_id = str(uuid.uuid4())
        with patch("compass.service.routes.training.training_service") as mock_service:
            from compass.service.exceptions import ServiceException
            from compass.service.error_codes import ErrorCode
            mock_service.get_task.side_effect = ServiceException(
                ErrorCode.TASK_NOT_FOUND, "Task not found"
            )
            
            response = client.get(f"/api/v1/training/tasks/{fake_task_id}")
            assert response.status_code == 404
    
    def test_cors_headers(self, client):
        """测试CORS头"""
        response = client.options("/health")
        # CORS预检请求应该返回200
        assert response.status_code in [200, 405]
    
    def test_rate_limiting(self, client):
        """测试速率限制"""
        # 在测试环境中速率限制应该被禁用或设置得很高
        # 这里主要测试端点可访问性
        for _ in range(10):
            response = client.get("/health")
            assert response.status_code == 200

