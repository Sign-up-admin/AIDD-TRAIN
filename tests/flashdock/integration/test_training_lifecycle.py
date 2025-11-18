"""
集成测试：训练任务生命周期
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

from compass_client import CompassClient
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
class TestTrainingLifecycle:
    """训练任务生命周期集成测试类"""
    
    def test_complete_training_lifecycle(self, compass_client):
        """测试完整的训练任务生命周期"""
        with requests_mock.Mocker() as m:
            task_id = "lifecycle-task-123"
            base_url = "http://localhost:8080/api/v1/training/tasks"
            
            # 1. 创建任务
            m.post(f"{base_url}", json={
                "task_id": task_id,
                "status": "pending",
                "config": {
                    "execution_mode": "smoke_test",
                    "epochs": 10,
                    "batch_size": 2,
                    "learning_rate": 0.001,
                    "optimizer": "adam"
                },
                "description": "Lifecycle test task",
                "created_at": "2024-01-01T00:00:00Z"
            }, status_code=201)
            
            # 设置所有需要的mock响应（按调用顺序）
            # requests_mock会按顺序匹配请求
            
            # 2. 第一次获取任务（pending状态）- 在创建后调用
            m.get(f"{base_url}/{task_id}", [
                {"json": {
                    "task_id": task_id,
                    "status": "pending",
                    "config": {"epochs": 10}
                }, "status_code": 200},
                # 4. 第二次获取任务（running状态）- 在启动后调用
                {"json": {
                    "task_id": task_id,
                    "status": "running",
                    "config": {"epochs": 10},
                    "progress": {"current_epoch": 5, "total_epochs": 10}
                }, "status_code": 200},
                # 7. 第三次获取任务（cancelled状态）- 在停止后调用
                {"json": {
                    "task_id": task_id,
                    "status": "cancelled",
                    "config": {"epochs": 10}
                }, "status_code": 200}
            ])
            
            # 3. 启动任务
            m.post(f"{base_url}/{task_id}/start", json={
                "message": "Task started successfully",
                "task_id": task_id
            }, status_code=200)
            
            # 5. 获取任务进度 - get_task_progress实际调用的是/metrics端点
            # 实际API返回的格式：{"metrics": {"training": {"current_epoch": ..., "total_epochs": ...}, ...}}
            m.get(f"{base_url}/{task_id}/metrics", json={
                "task_id": task_id,
                "metrics": {
                    "stage": "training",
                    "progress": 0.5,
                    "message": "Training in progress",
                    "cancelled": False,
                    "training": {
                        "current_epoch": 5,
                        "total_epochs": 10,
                        "current_batch": 0,
                        "total_batches": 0,
                        "train_loss": 0.0,
                        "val_loss": 0.0,
                        "epoch_progress": 50.0
                    },
                    "elapsed_time": 100.0
                }
            }, status_code=200)
            
            # 6. 停止任务
            m.post(f"{base_url}/{task_id}/stop", json={
                "message": "Task stopped successfully",
                "task_id": task_id
            }, status_code=200)
            
            # 执行生命周期
            # 1. 创建
            task = compass_client.create_training_task(
                config={
                    "execution_mode": "smoke_test",
                    "epochs": 10,
                    "batch_size": 2,
                    "learning_rate": 0.001,
                    "optimizer": "adam"
                },
                description="Lifecycle test task"
            )
            assert task["task_id"] == task_id
            assert task["status"] == "pending"
            
            # 2. 获取（pending状态）
            task = compass_client.get_training_task(task_id)
            assert task["status"] == "pending"
            
            # 3. 启动
            result = compass_client.start_training_task(task_id)
            assert "started" in result["message"].lower()
            
            # 4. 检查运行状态
            task = compass_client.get_training_task(task_id)
            assert task["status"] == "running"
            
            # 5. 获取进度 - get_task_progress返回metrics字典（包含training子字典）
            progress = compass_client.get_task_progress(task_id)
            training_info = progress.get("training", {})
            assert training_info.get("current_epoch") == 5
            assert training_info.get("total_epochs") == 10
            
            # 6. 停止
            result = compass_client.stop_training_task(task_id)
            assert "stop" in result["message"].lower()
            
            # 7. 验证停止状态
            task = compass_client.get_training_task(task_id)
            # 停止后状态应该是cancelled
            assert task["status"] == "cancelled"
    
    def test_task_with_dataset(self, compass_client, tmp_path):
        """测试使用数据集创建任务"""
        with requests_mock.Mocker() as m:
            # 创建测试数据集文件
            dataset_file = tmp_path / "dataset.zip"
            dataset_file.write_bytes(b"dataset data")
            
            dataset_id = "dataset-123"
            task_id = "task-with-dataset-123"
            
            # Mock上传数据集
            m.post("http://localhost:8080/api/v1/data/upload", json={
                "dataset_id": dataset_id,
                "name": "dataset.zip"
            }, status_code=200)
            
            # Mock创建任务（带数据集）
            m.post("http://localhost:8080/api/v1/training/tasks", json={
                "task_id": task_id,
                "status": "pending",
                "config": {"epochs": 10},
                "dataset_id": dataset_id
            }, status_code=201)
            
            # 上传数据集
            uploaded_dataset_id = compass_client.upload_dataset(str(dataset_file))
            assert uploaded_dataset_id == dataset_id
            
            # 使用数据集创建任务
            task = compass_client.create_training_task(
                config={"epochs": 10, "batch_size": 2},
                dataset_id=dataset_id
            )
            assert task["task_id"] == task_id
    
    def test_task_error_handling(self, compass_client):
        """测试任务错误处理"""
        with requests_mock.Mocker() as m:
            task_id = "error-task-123"
            
            # Mock创建任务失败
            m.post("http://localhost:8080/api/v1/training/tasks", json={
                "error": "Invalid configuration",
                "error_code": "ERR_1000",
                "detail": {"field": "epochs", "message": "Epochs must be between 1 and 10000"}
            }, status_code=400)
            
            # 应该抛出CompassError
            with pytest.raises(Exception):  # CompassError或HTTPError
                compass_client.create_training_task(
                    config={"epochs": 100000, "batch_size": 2}
                )

