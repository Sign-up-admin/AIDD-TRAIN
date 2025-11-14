"""
集成测试：COMPASS服务集成
"""
import pytest
import sys
import requests_mock
from pathlib import Path
from unittest.mock import Mock, patch

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from compass_client import CompassClient, CompassError
from services.common.service_protocol import ServiceInfo


@pytest.fixture
def mock_service_info():
    """创建模拟服务信息"""
    service = Mock(spec=ServiceInfo)
    service.service_id = "compass-1"
    service.host = "localhost"
    service.port = 8080
    service.status = Mock(value="healthy")
    return service


@pytest.fixture
def compass_client(mock_service_info):
    """创建CompassClient实例"""
    with patch('compass_client.FlashDockRegistryClient') as mock_registry_class:
        mock_registry = Mock()
        mock_registry.discover_compass_services.return_value = [mock_service_info]
        mock_registry_class.return_value = mock_registry
        
        with patch('compass_client.LoadBalancer') as mock_lb_class:
            mock_lb = Mock()
            mock_lb.select_service.return_value = mock_service_info
            mock_lb_class.return_value = mock_lb
            
            client = CompassClient()
            yield client


@pytest.mark.integration
class TestCompassIntegration:
    """COMPASS服务集成测试类"""
    
    def test_create_and_list_tasks(self, compass_client):
        """测试创建任务并列出任务"""
        with requests_mock.Mocker() as m:
            # Mock创建任务
            m.post("http://localhost:8080/api/v1/training/tasks", json={
                "task_id": "task-123",
                "status": "pending",
                "config": {"epochs": 10}
            }, status_code=201)
            
            # Mock列出任务
            m.get("http://localhost:8080/api/v1/training/tasks", json={
                "tasks": [
                    {
                        "task_id": "task-123",
                        "status": "pending",
                        "config": {"epochs": 10}
                    }
                ]
            }, status_code=200)
            
            # 创建任务
            task = compass_client.create_training_task(
                config={"epochs": 10, "batch_size": 2}
            )
            assert task["task_id"] == "task-123"
            
            # 列出任务
            tasks = compass_client.list_training_tasks()
            assert len(tasks) == 1
            assert tasks[0]["task_id"] == "task-123"
    
    def test_task_lifecycle(self, compass_client):
        """测试任务完整生命周期"""
        with requests_mock.Mocker() as m:
            task_id = "task-123"
            
            # 创建任务
            m.post("http://localhost:8080/api/v1/training/tasks", json={
                "task_id": task_id,
                "status": "pending",
                "config": {"epochs": 10}
            }, status_code=201)
            
            # 获取任务
            m.get(f"http://localhost:8080/api/v1/training/tasks/{task_id}", json={
                "task_id": task_id,
                "status": "pending",
                "config": {"epochs": 10}
            }, status_code=200)
            
            # 启动任务
            m.post(f"http://localhost:8080/api/v1/training/tasks/{task_id}/start", json={
                "message": "Task started",
                "task_id": task_id
            }, status_code=200)
            
            # 停止任务
            m.post(f"http://localhost:8080/api/v1/training/tasks/{task_id}/stop", json={
                "message": "Task stopped",
                "task_id": task_id
            }, status_code=200)
            
            # 执行生命周期
            task = compass_client.create_training_task(
                config={"epochs": 10, "batch_size": 2}
            )
            assert task["task_id"] == task_id
            
            task = compass_client.get_training_task(task_id)
            assert task["status"] == "pending"
            
            result = compass_client.start_training_task(task_id)
            assert "started" in result["message"].lower()
            
            result = compass_client.stop_training_task(task_id)
            assert "stop" in result["message"].lower()
    
    def test_dataset_upload_and_list(self, compass_client, tmp_path):
        """测试数据集上传和列表"""
        with requests_mock.Mocker() as m:
            # 创建测试文件
            test_file = tmp_path / "test_dataset.zip"
            test_file.write_bytes(b"test data")
            
            # Mock上传
            m.post("http://localhost:8080/api/v1/data/upload", json={
                "dataset_id": "dataset-123",
                "name": "test_dataset.zip"
            }, status_code=200)
            
            # Mock列表
            m.get("http://localhost:8080/api/v1/data/datasets", json={
                "datasets": [
                    {
                        "dataset_id": "dataset-123",
                        "name": "test_dataset.zip",
                        "size": 1024,
                        "file_count": 10,
                        "status": "ready"
                    }
                ]
            }, status_code=200)
            
            # 上传数据集
            dataset_id = compass_client.upload_dataset(
                str(test_file),
                name="Test Dataset"
            )
            assert dataset_id == "dataset-123"
            
            # 列出数据集
            datasets = compass_client.list_datasets()
            assert len(datasets) == 1
            assert datasets[0]["dataset_id"] == "dataset-123"
    
    def test_error_handling_and_retry(self, compass_client):
        """测试错误处理和重试机制"""
        with requests_mock.Mocker() as m:
            # 第一次请求失败
            m.get("http://localhost:8080/api/v1/training/tasks/task-123", [
                {"status_code": 500, "json": {"error": "Internal server error"}},
                {"status_code": 200, "json": {
                    "task_id": "task-123",
                    "status": "running"
                }}
            ])
            
            # 由于重试机制，应该最终成功
            task = compass_client.get_training_task("task-123")
            assert task["task_id"] == "task-123"
            assert task["status"] == "running"
    
    def test_service_discovery_integration(self):
        """测试服务发现集成"""
        with patch('compass_client.FlashDockRegistryClient') as mock_registry_class:
            mock_service1 = Mock(spec=ServiceInfo)
            mock_service1.service_id = "compass-1"
            mock_service1.host = "localhost"
            mock_service1.port = 8080
            mock_service1.status = Mock(value="healthy")
            
            mock_service2 = Mock(spec=ServiceInfo)
            mock_service2.service_id = "compass-2"
            mock_service2.host = "localhost"
            mock_service2.port = 8081
            mock_service2.status = Mock(value="healthy")
            
            mock_registry = Mock()
            mock_registry.discover_compass_services.return_value = [mock_service1, mock_service2]
            mock_registry_class.return_value = mock_registry
            
            with patch('compass_client.LoadBalancer') as mock_lb_class:
                mock_lb = Mock()
                # 第一次返回service1，第二次返回service2
                mock_lb.select_service.side_effect = [mock_service1, mock_service2]
                mock_lb_class.return_value = mock_lb
                
                client = CompassClient()
                
                # 验证服务发现被调用
                assert len(client._cached_services) == 2

