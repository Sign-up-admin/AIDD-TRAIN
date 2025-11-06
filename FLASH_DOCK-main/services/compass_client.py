"""
COMPASS service client for FLASH-DOCK.
"""
import sys
import logging
import requests
from pathlib import Path
from typing import List, Dict, Optional

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from registry_client import FlashDockRegistryClient
from load_balancer import LoadBalancer, LoadBalanceStrategy

logger = logging.getLogger(__name__)


class CompassClient:
    """Client for interacting with COMPASS service."""
    
    def __init__(
        self, 
        registry_url: str = "http://localhost:8500",
        timeout: float = 30.0,
        load_balance_strategy: LoadBalanceStrategy = LoadBalanceStrategy.ROUND_ROBIN
    ):
        """
        Initialize COMPASS client.
        
        Args:
            registry_url: Registry URL for service discovery
            timeout: Request timeout in seconds (default: 30.0)
            load_balance_strategy: Load balancing strategy
        """
        self.registry_client = FlashDockRegistryClient(registry_url=registry_url, timeout=5.0)
        self.load_balancer = LoadBalancer(strategy=load_balance_strategy)
        self.timeout = timeout
        self._cached_services = []
        self._refresh_services()
    
    def _refresh_services(self):
        """Refresh the list of available COMPASS services."""
        try:
            self._cached_services = self.registry_client.discover_compass_services(healthy_only=True)
            logger.debug(f"Refreshed services cache: {len(self._cached_services)} services available")
        except Exception as e:
            logger.warning(f"Failed to refresh services: {e}")
            # Keep existing cache if refresh fails
            if not self._cached_services:
                self._cached_services = []
    
    def _get_service_url(self) -> Optional[str]:
        """
        Get URL of a COMPASS service using load balancer.
        
        Returns:
            Optional[str]: Service URL or None if no service available
        """
        if not self._cached_services:
            self._refresh_services()
        
        if not self._cached_services:
            raise ConnectionError("No COMPASS services available")
        
        service = self.load_balancer.select_service(self._cached_services)
        if not service:
            raise ConnectionError("Failed to select COMPASS service")
        
        return f"http://{service.host}:{service.port}"
    
    def _make_request(
        self, 
        method: str, 
        endpoint: str, 
        **kwargs
    ) -> requests.Response:
        """
        Make HTTP request to COMPASS service with automatic retry and service selection.
        
        Args:
            method: HTTP method (GET, POST, DELETE, etc.)
            endpoint: API endpoint (e.g., "/api/v1/data/datasets")
            **kwargs: Additional arguments to pass to requests
            
        Returns:
            requests.Response: Response object
            
        Raises:
            ConnectionError: If no service is available
            requests.RequestException: If request fails
        """
        service_url = self._get_service_url()
        url = f"{service_url}{endpoint}"
        
        # Set default timeout if not provided
        if 'timeout' not in kwargs:
            kwargs['timeout'] = self.timeout
        
        try:
            response = requests.request(method, url, **kwargs)
            response.raise_for_status()
            return response
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed to {url}: {e}")
            # Try refreshing services and retry once
            self._refresh_services()
            service_url = self._get_service_url()
            url = f"{service_url}{endpoint}"
            response = requests.request(method, url, **kwargs)
            response.raise_for_status()
            return response
    
    def upload_dataset(
        self, 
        file_path: str, 
        name: Optional[str] = None, 
        description: Optional[str] = None
    ) -> str:
        """
        Upload a dataset file.
        
        Args:
            file_path: Path to the dataset file (zip, tar, tar.gz)
            name: Optional dataset name
            description: Optional dataset description
            
        Returns:
            str: Dataset ID
            
        Raises:
            FileNotFoundError: If file doesn't exist
            requests.RequestException: If upload fails
        """
        file_path_obj = Path(file_path)
        if not file_path_obj.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        with open(file_path, 'rb') as f:
            files = {'file': (file_path_obj.name, f, 'application/octet-stream')}
            data = {}
            if name:
                data['name'] = name
            if description:
                data['description'] = description
            
            response = self._make_request('POST', '/api/v1/data/upload', files=files, data=data)
            result = response.json()
            return result['dataset_id']
    
    def list_datasets(self) -> List[Dict]:
        """
        List all datasets.
        
        Returns:
            List[Dict]: List of dataset dictionaries with keys:
                dataset_id, name, description, size, file_count, status, 
                created_at, updated_at, metadata, error
        """
        response = self._make_request('GET', '/api/v1/data/datasets')
        result = response.json()
        
        # Convert datasets to list of dicts with datetime serialization
        datasets = []
        for ds in result.get('datasets', []):
            dataset_dict = {
                'dataset_id': ds['dataset_id'],
                'name': ds['name'],
                'description': ds.get('description'),
                'size': ds.get('size', 0),
                'file_count': ds.get('file_count', 0),
                'status': ds.get('status', 'pending'),
                'created_at': ds.get('created_at'),
                'updated_at': ds.get('updated_at'),
                'metadata': ds.get('metadata', {}),
                'error': ds.get('error')
            }
            datasets.append(dataset_dict)
        
        return datasets
    
    def delete_dataset(self, dataset_id: str) -> None:
        """
        Delete a dataset.
        
        Args:
            dataset_id: ID of the dataset to delete
            
        Raises:
            requests.RequestException: If deletion fails
        """
        self._make_request('DELETE', f'/api/v1/data/datasets/{dataset_id}')
    
    def get_inference_status(self) -> Dict:
        """
        Get inference service status.
        
        Returns:
            Dict: Status dictionary with keys:
                status, loaded_models, max_cache_size, cached_model_ids,
                latest_model, device
        """
        response = self._make_request('GET', '/api/v1/inference/status')
        return response.json()
    
    def list_models(self) -> List[Dict]:
        """
        List all available models.
        
        Returns:
            List[Dict]: List of model dictionaries with keys:
                model_id, name, version, task_id, checkpoint_path,
                file_size, created_at, metadata, metrics
        """
        response = self._make_request('GET', '/api/v1/models')
        result = response.json()
        
        # Convert models to list of dicts
        models = []
        for model in result.get('models', []):
            model_dict = {
                'model_id': model['model_id'],
                'name': model['name'],
                'version': model['version'],
                'task_id': model.get('task_id'),
                'checkpoint_path': model.get('checkpoint_path'),
                'file_size': model.get('file_size', 0),
                'created_at': model.get('created_at'),
                'metadata': model.get('metadata', {}),
                'metrics': model.get('metrics', {})
            }
            models.append(model_dict)
        
        return models
    
    def create_training_task(
        self,
        config: Dict,
        dataset_id: Optional[str] = None,
        description: Optional[str] = None
    ) -> Dict:
        """
        Create a new training task.
        
        Args:
            config: Training configuration dictionary with keys:
                execution_mode: One of ['validation_tuned', 'validation', 'prototyping', 'smoke_test', 'production']
                epochs: Number of training epochs (1-10000)
                batch_size: Batch size (1-128)
                learning_rate: Learning rate (0-1.0)
                optimizer: Optimizer name ('adam', 'adamw', 'sgd', 'rmsprop')
            dataset_id: Optional dataset ID to use
            description: Optional task description
            
        Returns:
            Dict: Task dictionary with keys:
                task_id, status, config, created_at, updated_at, 
                started_at, completed_at, progress, error, description
                
        Raises:
            requests.RequestException: If creation fails
        """
        data = {'config': config}
        if dataset_id:
            data['dataset_id'] = dataset_id
        if description:
            data['description'] = description
        
        response = self._make_request('POST', '/api/v1/training/tasks', json=data)
        return response.json()
    
    def list_training_tasks(self) -> List[Dict]:
        """
        List all training tasks.
        
        Returns:
            List[Dict]: List of task dictionaries with keys:
                task_id, status, config, created_at, updated_at,
                started_at, completed_at, progress, error, description
        """
        response = self._make_request('GET', '/api/v1/training/tasks')
        result = response.json()
        
        # Convert tasks to list of dicts
        tasks = []
        for task in result.get('tasks', []):
            task_dict = {
                'task_id': task['task_id'],
                'status': task['status'],
                'config': task.get('config', {}),
                'created_at': task.get('created_at'),
                'updated_at': task.get('updated_at'),
                'started_at': task.get('started_at'),
                'completed_at': task.get('completed_at'),
                'progress': task.get('progress', {}),
                'error': task.get('error'),
                'description': task.get('description')
            }
            tasks.append(task_dict)
        
        return tasks
    
    def get_training_task(self, task_id: str) -> Dict:
        """
        Get detailed information about a specific training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Task dictionary with all details
            
        Raises:
            requests.RequestException: If task not found or request fails
        """
        response = self._make_request('GET', f'/api/v1/training/tasks/{task_id}')
        task = response.json()
        
        return {
            'task_id': task['task_id'],
            'status': task['status'],
            'config': task.get('config', {}),
            'created_at': task.get('created_at'),
            'updated_at': task.get('updated_at'),
            'started_at': task.get('started_at'),
            'completed_at': task.get('completed_at'),
            'progress': task.get('progress', {}),
            'error': task.get('error'),
            'description': task.get('description')
        }
    
    def start_training_task(self, task_id: str) -> Dict:
        """
        Start a training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Confirmation message with task_id
            
        Raises:
            requests.RequestException: If task cannot be started
        """
        response = self._make_request('POST', f'/api/v1/training/tasks/{task_id}/start')
        return response.json()
    
    def stop_training_task(self, task_id: str) -> Dict:
        """
        Stop a running training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Confirmation message with task_id
            
        Raises:
            requests.RequestException: If task cannot be stopped
        """
        response = self._make_request('POST', f'/api/v1/training/tasks/{task_id}/stop')
        return response.json()
    
    def pause_training_task(self, task_id: str) -> Dict:
        """
        Pause a training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Confirmation message with task_id
            
        Raises:
            requests.RequestException: If task cannot be paused
        """
        response = self._make_request('POST', f'/api/v1/training/tasks/{task_id}/pause')
        return response.json()
    
    def get_task_logs(self, task_id: str, limit: int = 100) -> List[str]:
        """
        Get log messages for a training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            limit: Maximum number of log lines to return (default: 100, max: 1000)
            
        Returns:
            List[str]: List of log messages
            
        Raises:
            requests.RequestException: If task not found or request fails
        """
        params = {'limit': min(limit, 1000)}
        response = self._make_request('GET', f'/api/v1/training/tasks/{task_id}/logs', params=params)
        result = response.json()
        return result.get('logs', [])
    
    def get_task_progress(self, task_id: str) -> Dict:
        """
        Get task progress and metrics.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Progress dictionary with metrics
            
        Raises:
            requests.RequestException: If task not found or request fails
        """
        response = self._make_request('GET', f'/api/v1/training/tasks/{task_id}/metrics')
        result = response.json()
        return result.get('metrics', {})
    
    def delete_training_task(self, task_id: str) -> Dict:
        """
        Delete a training task.
        
        Args:
            task_id: Unique task identifier (UUID)
            
        Returns:
            Dict: Confirmation message with task_id
            
        Raises:
            requests.RequestException: If task not found or request fails
        """
        response = self._make_request('DELETE', f'/api/v1/training/tasks/{task_id}')
        return response.json()