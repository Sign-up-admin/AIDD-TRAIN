"""
COMPASS service client for FLASH-DOCK.
"""
import sys
import logging
import requests
from pathlib import Path
from typing import Dict, Optional, List

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Import from local FLASH_DOCK services directory
sys.path.insert(0, str(Path(__file__).parent))
from service_manager import ServiceManager
from services.common.service_protocol import ServiceInfo

logger = logging.getLogger(__name__)


class CompassClient:
    """Client for interacting with COMPASS service."""
    
    def __init__(self, registry_url: str = "http://localhost:8500"):
        """
        Initialize COMPASS client.
        
        Args:
            registry_url: Service registry URL
        """
        self.service_manager = ServiceManager(registry_url=registry_url)
        self.session = requests.Session()
        self.session.headers.update({'Content-Type': 'application/json'})
    
    def _get_service(self) -> Optional[ServiceInfo]:
        """Get a COMPASS service instance."""
        return self.service_manager.get_compass_service()
    
    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """
        Make HTTP request to COMPASS service with detailed error handling.
        
        Args:
            method: HTTP method
            endpoint: API endpoint
            **kwargs: Additional arguments for requests
            
        Returns:
            requests.Response: Response object
            
        Raises:
            ConnectionError: If no service is available
            requests.exceptions.RequestException: If request fails
        """
        service = None
        service_id = None
        
        try:
        # First attempt
        service = self._get_service()
        if not service:
            error_msg = "No COMPASS service available from registry"
            logger.error(error_msg)
            # Try to refresh and retry once
            logger.info("Attempting to refresh service list...")
            self.service_manager.refresh_services()
            service = self._get_service()
            if not service:
                raise ConnectionError(error_msg)
        
        service_id = service.service_id
        # Increment connection count
        self.service_manager.load_balancer.increment_connections(service_id)
        
        url = f"{service.base_url}{endpoint}"
        logger.debug(f"Making {method} request to {url}")
        
        try:
            response = self.session.request(method, url, timeout=30, **kwargs)
            response.raise_for_status()
            logger.debug(f"Request successful: {response.status_code}")
            return response
            except requests.exceptions.ConnectionError as e:
                logger.error(f"Connection error to {url}: {e}")
                # Try to refresh service list and retry once
                logger.info("Connection failed, refreshing service list...")
                self.service_manager.refresh_services()
                service = self._get_service()
                if service:
                    # Decrement old service, increment new one
                    if service_id:
                        self.service_manager.load_balancer.decrement_connections(service_id)
                    service_id = service.service_id
                    self.service_manager.load_balancer.increment_connections(service_id)
                    
                    url = f"{service.base_url}{endpoint}"
                    logger.info(f"Retrying request to {url}")
                    try:
                        response = self.session.request(method, url, timeout=30, **kwargs)
                        response.raise_for_status()
                        return response
                    except Exception as retry_e:
                        logger.error(f"Retry also failed: {retry_e}")
                        raise
                raise ConnectionError(f"Failed to connect to COMPASS service: {e}")
            except requests.exceptions.Timeout as e:
                logger.error(f"Request timeout to {url}: {e}")
                raise
            except requests.exceptions.HTTPError as e:
                logger.error(f"HTTP error {e.response.status_code} from {url}: {e}")
                raise
            except requests.exceptions.RequestException as e:
                logger.error(f"Request failed to {url}: {e}")
                raise
        finally:
            # Decrement connection count when request completes
            if service_id:
                self.service_manager.load_balancer.decrement_connections(service_id)
    
    # Training API
    
    def create_training_task(self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None) -> str:
        """
        Create a training task.
        
        Args:
            config: Training configuration
            dataset_id: Optional dataset ID
            description: Optional description
            
        Returns:
            str: Task ID
        """
        data = {
            "config": config,
            "dataset_id": dataset_id,
            "description": description
        }
        response = self._make_request("POST", "/api/v1/training/tasks", json=data)
        result = response.json()
        return result['task_id']
    
    def start_training(self, task_id: str):
        """Start a training task."""
        self._make_request("POST", f"/api/v1/training/tasks/{task_id}/start")
    
    def stop_training(self, task_id: str):
        """Stop a training task."""
        self._make_request("POST", f"/api/v1/training/tasks/{task_id}/stop")
    
    def pause_training(self, task_id: str):
        """Pause a training task."""
        self._make_request("POST", f"/api/v1/training/tasks/{task_id}/pause")
    
    def get_task_status(self, task_id: str) -> Dict:
        """Get task status."""
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}")
        return response.json()
    
    def list_tasks(self) -> List[Dict]:
        """List all tasks."""
        response = self._make_request("GET", "/api/v1/training/tasks")
        result = response.json()
        return result['tasks']
    
    def get_task_logs(self, task_id: str, limit: int = 100) -> List[str]:
        """Get task logs."""
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}/logs", params={"limit": limit})
        result = response.json()
        return result['logs']
    
    def delete_task(self, task_id: str):
        """Delete a task."""
        self._make_request("DELETE", f"/api/v1/training/tasks/{task_id}")
    
    # Data API
    
    def upload_dataset(self, file_path: str, name: Optional[str] = None, description: Optional[str] = None) -> str:
        """
        Upload a dataset.
        
        Args:
            file_path: Path to dataset file
            name: Optional dataset name
            description: Optional description
            
        Returns:
            str: Dataset ID
        """
        with open(file_path, 'rb') as f:
            files = {'file': (Path(file_path).name, f)}
            data = {}
            if name:
                data['name'] = name
            if description:
                data['description'] = description
            
            response = self._make_request("POST", "/api/v1/data/upload", files=files, data=data)
            result = response.json()
            return result['dataset_id']
    
    def list_datasets(self) -> List[Dict]:
        """List all datasets."""
        response = self._make_request("GET", "/api/v1/data/datasets")
        result = response.json()
        return result['datasets']
    
    def delete_dataset(self, dataset_id: str):
        """Delete a dataset."""
        self._make_request("DELETE", f"/api/v1/data/datasets/{dataset_id}")
    
    # Model API
    
    def list_models(self) -> List[Dict]:
        """List all models."""
        response = self._make_request("GET", "/api/v1/models")
        result = response.json()
        return result['models']
    
    def download_model(self, model_id: str, save_path: str):
        """
        Download a model checkpoint.
        
        Args:
            model_id: Model ID
            save_path: Path to save the checkpoint
        """
        response = self._make_request("GET", f"/api/v1/models/{model_id}/download", stream=True)
        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
    
    # Inference API
    
    def predict(self, protein_path: str, ligand_path: str, model_id: Optional[str] = None) -> float:
        """
        Predict binding affinity.
        
        Args:
            protein_path: Path to protein PDB file
            ligand_path: Path to ligand SDF file
            model_id: Optional model ID (uses latest if None)
            
        Returns:
            float: Binding affinity
        """
        data = {
            "protein_path": protein_path,
            "ligand_path": ligand_path,
            "model_id": model_id
        }
        response = self._make_request("POST", "/api/v1/inference/predict", json=data)
        result = response.json()
        return result['binding_affinity']
    
    def batch_predict(self, protein_ligand_pairs: List[Dict[str, str]], model_id: Optional[str] = None) -> Dict:
        """
        Perform batch predictions.
        
        Args:
            protein_ligand_pairs: List of {"protein": path, "ligand": path}
            model_id: Optional model ID
            
        Returns:
            Dict: Batch prediction results
        """
        data = {
            "protein_ligand_pairs": protein_ligand_pairs,
            "model_id": model_id
        }
        response = self._make_request("POST", "/api/v1/inference/batch", json=data)
        return response.json()
    
    def get_inference_status(self) -> Dict:
        """Get inference service status."""
        response = self._make_request("GET", "/api/v1/inference/status")
        return response.json()

