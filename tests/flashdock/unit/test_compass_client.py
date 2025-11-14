"""
单元测试：CompassClient COMPASS服务客户端
"""
import pytest
import sys
import json
import requests
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add paths
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "FLASH_DOCK-main" / "services"))

from compass_client import CompassClient, CompassError
from load_balancer import LoadBalanceStrategy
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
def mock_registry_client(mock_service_info):
    """创建模拟注册中心客户端"""
    registry = Mock()
    registry.discover_compass_services.return_value = [mock_service_info]
    registry.registry_url = "http://localhost:8500"
    return registry


class TestCompassClient:
    """CompassClient测试类"""
    
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_init_default(self, mock_load_balancer_class, mock_registry_class):
        """测试默认初始化"""
        mock_registry = Mock()
        mock_registry.discover_compass_services.return_value = []
        mock_registry_class.return_value = mock_registry
        
        client = CompassClient()
        
        assert client.registry_client is not None
        assert client.timeout == 30.0
        mock_registry_class.assert_called_once_with(
            registry_url="http://localhost:8500",
            timeout=5.0
        )
    
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_init_custom(self, mock_load_balancer_class, mock_registry_class):
        """测试自定义初始化"""
        mock_registry = Mock()
        mock_registry.discover_compass_services.return_value = []
        mock_registry_class.return_value = mock_registry
        
        client = CompassClient(
            registry_url="http://test:9000",
            timeout=60.0,
            load_balance_strategy=LoadBalanceStrategy.RANDOM
        )
        
        assert client.timeout == 60.0
        mock_registry_class.assert_called_once_with(
            registry_url="http://test:9000",
            timeout=5.0
        )
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_create_training_task_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功创建训练任务"""
        # 设置mock
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        # Mock响应
        mock_response = Mock()
        mock_response.status_code = 201
        mock_response.json.return_value = {
            "task_id": "test-task-123",
            "status": "pending",
            "config": {"epochs": 10}
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        result = client.create_training_task(
            config={"epochs": 10, "batch_size": 2},
            description="Test task"
        )
        
        assert result["task_id"] == "test-task-123"
        assert result["status"] == "pending"
        mock_request.assert_called_once()
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_create_training_task_missing_task_id(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试创建任务响应缺少task_id"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 201
        mock_response.json.return_value = {"status": "pending"}  # 缺少task_id
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        with pytest.raises(CompassError) as exc_info:
            client.create_training_task(config={"epochs": 10})
        
        assert "missing 'task_id'" in str(exc_info.value)
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_list_training_tasks_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功列出训练任务"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "tasks": [
                {"task_id": "task-1", "status": "running"},
                {"task_id": "task-2", "status": "completed"}
            ]
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        tasks = client.list_training_tasks()
        
        assert len(tasks) == 2
        assert tasks[0]["task_id"] == "task-1"
        assert tasks[1]["task_id"] == "task-2"
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_get_training_task_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功获取训练任务"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "task_id": "task-123",
            "status": "running",
            "config": {"epochs": 10}
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        task = client.get_training_task("task-123")
        
        assert task["task_id"] == "task-123"
        assert task["status"] == "running"
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_start_training_task_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功启动训练任务"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "message": "Task started",
            "task_id": "task-123"
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        result = client.start_training_task("task-123")
        
        assert result["task_id"] == "task-123"
        assert "started" in result["message"].lower()
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_stop_training_task_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功停止训练任务"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "message": "Task stopped",
            "task_id": "task-123"
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        result = client.stop_training_task("task-123")
        
        assert result["task_id"] == "task-123"
        assert "stop" in result["message"].lower()
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_upload_dataset_success(self, mock_load_balancer_class, mock_registry_class, mock_request, tmp_path):
        """测试成功上传数据集"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        # 创建测试文件
        test_file = tmp_path / "test_dataset.zip"
        test_file.write_bytes(b"test data")
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "dataset_id": "dataset-123",
            "name": "test_dataset.zip"
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        dataset_id = client.upload_dataset(str(test_file), name="Test Dataset")
        
        assert dataset_id == "dataset-123"
        mock_request.assert_called_once()
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_upload_dataset_file_not_found(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试上传不存在的文件"""
        mock_registry = Mock()
        mock_registry.discover_compass_services.return_value = []
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer_class.return_value = mock_load_balancer
        
        client = CompassClient()
        with pytest.raises(FileNotFoundError):
            client.upload_dataset("nonexistent_file.zip")
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_list_datasets_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功列出数据集"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "datasets": [
                {
                    "dataset_id": "dataset-1",
                    "name": "Dataset 1",
                    "size": 1024000
                }
            ]
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        datasets = client.list_datasets()
        
        assert len(datasets) == 1
        assert datasets[0]["dataset_id"] == "dataset-1"
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_list_models_success(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试成功列出模型"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "models": [
                {
                    "model_id": "model-1",
                    "name": "Model 1",
                    "version": "1.0.0",
                    "file_size": 50000000
                }
            ]
        }
        mock_response.raise_for_status = Mock()
        mock_request.return_value = mock_response
        
        client = CompassClient()
        models = client.list_models()
        
        assert len(models) == 1
        assert models[0]["model_id"] == "model-1"
    
    @patch('compass_client.requests.request')
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_compass_error_handling(self, mock_load_balancer_class, mock_registry_class, mock_request):
        """测试CompassError错误处理"""
        mock_registry = Mock()
        mock_service = Mock(spec=ServiceInfo)
        mock_service.host = "localhost"
        mock_service.port = 8080
        mock_service.status = Mock(value="healthy")
        mock_registry.discover_compass_services.return_value = [mock_service]
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = mock_service
        mock_load_balancer_class.return_value = mock_load_balancer
        
        # Mock HTTP错误响应
        mock_response = Mock()
        mock_response.status_code = 400
        mock_response.headers = {"content-type": "application/json"}
        mock_response.json.return_value = {
            "error": "Bad request",
            "error_code": "ERR_1000",
            "detail": {"field": "config"}
        }
        mock_response.text = json.dumps(mock_response.json.return_value)
        mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(response=mock_response)
        mock_request.return_value = mock_response
        
        client = CompassClient()
        with pytest.raises(CompassError) as exc_info:
            client.get_training_task("task-123")
        
        assert exc_info.value.status_code == 400
        assert exc_info.value.error_code == "ERR_1000"
        assert "Bad request" in str(exc_info.value)
    
    @patch('compass_client.FlashDockRegistryClient')
    @patch('compass_client.LoadBalancer')
    def test_no_service_available(self, mock_load_balancer_class, mock_registry_class):
        """测试没有可用服务时的错误"""
        mock_registry = Mock()
        mock_registry.discover_compass_services.return_value = []
        mock_registry_class.return_value = mock_registry
        
        mock_load_balancer = Mock()
        mock_load_balancer.select_service.return_value = None
        mock_load_balancer_class.return_value = mock_load_balancer
        
        client = CompassClient()
        with pytest.raises(ConnectionError) as exc_info:
            client.get_training_task("task-123")
        
        assert "没有可用的COMPASS服务" in str(exc_info.value)

