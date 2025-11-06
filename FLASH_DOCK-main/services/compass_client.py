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
    
    def __init__(self, registry_url: str = "http://localhost:8500", timeout: float = 5.0):
        """
        Initialize COMPASS client.
        
        Args:
            registry_url: Service registry URL
            timeout: Request timeout in seconds for registry operations (default: 5.0)
        """
        self.registry_url = registry_url
        self.service_manager = ServiceManager(registry_url=registry_url, timeout=timeout)
        self.session = requests.Session()
        self.session.headers.update({'Content-Type': 'application/json'})
        self.registry_timeout = timeout
    
    def _get_service(self) -> Optional[ServiceInfo]:
        """
        Get a COMPASS service instance.
        
        Returns:
            Optional[ServiceInfo]: Service instance or None if unavailable
        """
        return self.service_manager.get_compass_service()
    
    def is_service_available(self) -> bool:
        """
        Check if COMPASS service is available.
        
        Returns:
            bool: True if service is available, False otherwise
        """
        return self.service_manager.is_registry_available() and self._get_service() is not None
    
    def check_connection(self) -> Dict:
        """
        Check connection status to COMPASS service with detailed diagnostics.
        
        Returns:
            Dict: Connection status information including:
                - registry_available: bool
                - compass_service_available: bool
                - compass_service_info: Optional[ServiceInfo]
                - connection_test: bool (whether test request succeeded)
                - error_message: Optional[str]
        """
        result = {
            'registry_available': False,
            'compass_service_available': False,
            'compass_service_info': None,
            'connection_test': False,
            'error_message': None
        }
        
        try:
            # Check registry
            result['registry_available'] = self.service_manager.is_registry_available()
            
            if not result['registry_available']:
                result['error_message'] = f"服务注册中心不可用 ({self.registry_url})"
                return result
            
            # Check compass service
            service = self._get_service()
            if service:
                result['compass_service_available'] = True
                result['compass_service_info'] = {
                    'service_id': service.service_id,
                    'base_url': service.base_url,
                    'status': service.status.value if hasattr(service.status, 'value') else str(service.status)
                }
                
                # Test connection by making a simple request
                try:
                    test_response = self._make_request("GET", "/health", timeout=5)
                    if test_response.status_code == 200:
                        result['connection_test'] = True
                except Exception as e:
                    result['error_message'] = f"连接测试失败: {str(e)}"
            else:
                result['error_message'] = "未找到可用的COMPASS服务"
                
        except Exception as e:
            result['error_message'] = f"检查连接时发生错误: {str(e)}"
            logger.error(f"Error checking connection: {e}", exc_info=True)
        
        return result
    
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
            # Check registry availability first
            if not self.service_manager.is_registry_available():
                error_msg = (
                    f"服务注册中心不可用 (Registry at {self.registry_url} is not available). "
                    "请确保服务注册中心正在运行。"
                )
                logger.error(error_msg)
                raise ConnectionError(error_msg)
            
            # First attempt
            service = self._get_service()
            if not service:
                error_msg = "没有可用的COMPASS服务。请确保COMPASS服务已注册并正在运行。"
                logger.error(error_msg)
                # Try to refresh and retry once
                logger.info("尝试刷新服务列表...")
                refresh_success = self.service_manager.refresh_services()
                if not refresh_success:
                    raise ConnectionError(
                        f"无法刷新服务列表。注册中心可能不可用 ({self.registry_url})"
                    )
                service = self._get_service()
                if not service:
                    raise ConnectionError(
                        "没有可用的COMPASS服务。请检查：\n"
                        "  1. 服务注册中心是否正在运行\n"
                        "  2. COMPASS服务是否已注册并处于健康状态"
                    )
            
            service_id = service.service_id
            # Increment connection count
            self.service_manager.load_balancer.increment_connections(service_id)
            
            url = f"{service.base_url}{endpoint}"
            logger.debug(f"Making {method} request to {url}")
            
            # Extract timeout from kwargs or use default
            timeout = kwargs.pop('timeout', 30)
            try:
                response = self.session.request(method, url, timeout=timeout, **kwargs)
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
                        timeout = kwargs.pop('timeout', 30) if 'timeout' in kwargs else 30
                        response = self.session.request(method, url, timeout=timeout, **kwargs)
                        response.raise_for_status()
                        return response
                    except Exception as retry_e:
                        logger.error(f"Retry also failed: {retry_e}")
                        raise
                raise ConnectionError(
                    f"无法连接到COMPASS服务 ({url})。请检查服务是否正在运行。"
                )
            except requests.exceptions.Timeout as e:
                timeout_val = kwargs.get('timeout', 30)
                logger.error(f"请求超时 ({timeout_val}s) 到 {url}: {e}")
                raise TimeoutError(
                    f"请求超时（超过{timeout_val}秒）。服务可能响应缓慢或不可用。"
                )
            except ConnectionError:
                # Re-raise ConnectionError as-is (we already have good error messages)
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
    
    def create_training_task(self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None, timeout: int = 60) -> str:
        """
        Create a training task with timeout and detailed error handling.
        
        Args:
            config: Training configuration
            dataset_id: Optional dataset ID
            description: Optional description
            timeout: Request timeout in seconds (default: 60)
            
        Returns:
            str: Task ID
            
        Raises:
            ConnectionError: If service is unavailable
            TimeoutError: If request times out
            ValueError: If resource check fails or invalid config
        """
        data = {
            "config": config,
            "dataset_id": dataset_id,
            "description": description
        }
        try:
            # Use longer timeout for task creation (includes resource checks)
            original_timeout = self.session.timeout if hasattr(self.session, 'timeout') else 30
            response = self._make_request("POST", "/api/v1/training/tasks", json=data, timeout=timeout)
            result = response.json()
            
            # Check if task was created successfully
            if 'task_id' not in result:
                raise ValueError("Server response missing task_id")
            
            return result['task_id']
        except requests.exceptions.Timeout as e:
            logger.error(f"Task creation timeout after {timeout}s: {e}")
            raise TimeoutError(f"任务创建超时（超过{timeout}秒）。可能是系统资源不足或服务繁忙，请稍后重试。")
        except requests.exceptions.HTTPError as e:
            # Parse error response for user-friendly messages
            if e.response is not None:
                try:
                    error_data = e.response.json()
                    error_code = error_data.get('error_code', '')
                    error_detail = error_data.get('detail', {})
                    
                    # Map error codes to user-friendly messages
                    if error_code == 'ERR_1003':  # TIMEOUT
                        raise TimeoutError("任务创建超时。请检查系统资源或稍后重试。")
                    elif error_code == 'ERR_1002':  # SERVICE_UNAVAILABLE
                        resource_check = error_detail.get('resource_check', {})
                        issues = resource_check.get('issues', [])
                        if issues:
                            raise ValueError(f"系统资源不足：\n" + "\n".join(f"  • {issue}" for issue in issues))
                        raise ConnectionError("服务暂时不可用，请稍后重试。")
                    elif error_code == 'ERR_2000':  # VALIDATION_ERROR
                        raise ValueError(f"配置验证失败：{error_data.get('detail', {}).get('error_message', str(e))}")
                    else:
                        error_msg = error_data.get('error', error_data.get('detail', {}).get('error_message', str(e)))
                        raise ValueError(f"创建任务失败：{error_msg}")
                except ValueError:
                    raise  # Re-raise our custom ValueError
                except Exception:
                    # If JSON parsing fails, use HTTP error message
                    if e.response.status_code == 503:
                        raise ConnectionError(f"服务暂时不可用：{str(e)}")
                    elif e.response.status_code == 400:
                        raise ValueError(f"请求参数错误：{str(e)}")
                    else:
                        raise ValueError(f"创建任务失败（HTTP {e.response.status_code}）：{str(e)}")
            raise ValueError(f"创建任务失败：{str(e)}")
        except ConnectionError as e:
            # Re-raise with original message
            logger.error(f"Connection error during task creation: {e}")
            raise
        except requests.exceptions.ConnectionError as e:
            logger.error(f"Connection error during task creation: {e}")
            # Check if registry is available
            if not self.service_manager.is_registry_available():
                raise ConnectionError(
                    f"无法连接到服务注册中心 ({self.registry_url})。请确保服务注册中心正在运行。"
                )
            raise ConnectionError(
                "无法连接到COMPASS服务。请检查：\n"
                "  1. COMPASS服务是否正在运行\n"
                "  2. 服务是否已在注册中心注册"
            )
        except Exception as e:
            logger.error(f"Unexpected error during task creation: {e}", exc_info=True)
            raise ValueError(f"创建任务时发生未知错误：{str(e)}")
    
    def start_training(self, task_id: str, timeout: int = 30, verify_start: bool = True, max_retries: int = 3):
        """
        Start a training task with verification and retry mechanism.
        
        Args:
            task_id: Task ID
            timeout: Request timeout in seconds (default: 30)
            verify_start: Whether to verify task actually started by checking status (default: True)
            max_retries: Maximum number of retry attempts if start fails (default: 3)
            
        Returns:
            Dict: Response from server containing message and task_id
            
        Raises:
            ConnectionError: If service is unavailable
            TimeoutError: If request times out
            ValueError: If task cannot be started or verification fails
        """
        import time as time_module
        
        # Check service availability before starting
        if not self.is_service_available():
            raise ConnectionError(
                "COMPASS服务不可用，无法启动任务。请检查：\n"
                "  1. COMPASS服务是否正在运行\n"
                "  2. 服务是否已在注册中心注册\n"
                "  3. 服务注册中心是否可用"
            )
        
        last_exception = None
        
        # Retry logic
        for attempt in range(max_retries):
            try:
                # Get current task status before starting
                initial_status = None
                if verify_start:
                    try:
                        task_status = self.get_task_status(task_id)
                        initial_status = task_status.get('status')
                        logger.debug(f"Task {task_id} initial status: {initial_status}")
                    except Exception as e:
                        logger.warning(f"Could not get initial status for task {task_id}: {e}")
                
                # Send start request
                response = self._make_request("POST", f"/api/v1/training/tasks/{task_id}/start", timeout=timeout)
                
                # Parse and verify response
                try:
                    result = response.json()
                    if not isinstance(result, dict):
                        raise ValueError(f"Invalid response format: expected dict, got {type(result)}")
                    
                    # Verify response contains expected fields
                    if 'task_id' not in result or result['task_id'] != task_id:
                        logger.warning(f"Response task_id mismatch: expected {task_id}, got {result.get('task_id')}")
                    
                    response_message = result.get('message', '')
                    logger.info(f"Start command sent successfully: {response_message}")
                    
                    # Verify task actually started by checking status change (快速验证模式)
                    if verify_start:
                        # 快速验证：减少等待时间，最多检查2次
                        verification_attempts = 2
                        initial_wait = 0.3  # 从1秒减少到0.3秒
                        check_interval = 0.5  # 检查间隔0.5秒
                        max_verification_time = 2.0  # 最大验证时间2秒
                        
                        status_verified = False
                        verification_start_time = time_module.time()
                        
                        # 初始等待
                        time_module.sleep(initial_wait)
                        
                        for v_attempt in range(verification_attempts):
                            # 检查是否超过最大验证时间
                            elapsed = time_module.time() - verification_start_time
                            if elapsed > max_verification_time:
                                logger.info(f"达到最大验证时间 ({max_verification_time}s)，停止验证")
                                break
                            
                            try:
                                task_status = self.get_task_status(task_id)
                                current_status = task_status.get('status')
                                
                                # Status should change from pending/paused to initializing/running
                                if initial_status:
                                    # 如果状态已经是 INITIALIZING 或 RUNNING，立即返回成功
                                    if current_status in ['initializing', 'running']:
                                        logger.info(f"✓ Task {task_id} status verified: {initial_status} -> {current_status} (快速验证成功)")
                                        status_verified = True
                                        break
                                    elif current_status == 'failed':
                                        error_msg = task_status.get('error', 'Unknown error')
                                        raise ValueError(f"任务启动失败，状态变为failed: {error_msg}")
                                    elif current_status != initial_status:
                                        # 状态变化了，但不一定是预期的状态，记录日志
                                        logger.info(f"Task {task_id} status changed: {initial_status} -> {current_status}")
                                else:
                                    # Can't verify without initial status, but check if it's running
                                    if current_status in ['initializing', 'running']:
                                        logger.info(f"✓ Task {task_id} is in {current_status} status (快速验证成功)")
                                        status_verified = True
                                        break
                                
                                # 如果不是最后一次尝试，等待后继续
                                if v_attempt < verification_attempts - 1:
                                    time_module.sleep(check_interval)
                                    
                            except requests.exceptions.HTTPError as e:
                                # 处理速率限制 (HTTP 429) - 不阻塞，记录警告
                                if e.response and e.response.status_code == 429:
                                    logger.warning(f"速率限制 (HTTP 429)，跳过状态验证。启动命令已发送，任务应已开始。")
                                    # 不标记为失败，因为启动命令已经成功发送
                                    status_verified = None  # None表示不确定，但不阻塞
                                    break
                                else:
                                    logger.warning(f"Status verification attempt {v_attempt + 1} failed: {e}")
                                    if v_attempt == verification_attempts - 1:
                                        # 最后一次尝试失败，但不抛出异常，只记录警告
                                        logger.warning(
                                            f"无法验证任务状态变化。启动命令已发送，但状态检查失败: {e}\n"
                                            f"请手动检查任务状态。"
                                        )
                            except Exception as e:
                                logger.warning(f"Status verification attempt {v_attempt + 1} failed: {e}")
                                if v_attempt == verification_attempts - 1:
                                    # 最后一次尝试失败，但不抛出异常，只记录警告
                                    logger.warning(
                                        f"无法验证任务状态变化。启动命令已发送，但状态检查失败: {e}\n"
                                        f"请手动检查任务状态。"
                                    )
                        
                        # 记录验证结果
                        if status_verified is True:
                            logger.info(f"✓ Task {task_id} 启动验证成功")
                        elif status_verified is None:
                            logger.warning(
                                f"⚠ Task {task_id} 启动命令已发送，但无法验证状态（可能因速率限制）。"
                                f"原始状态: {initial_status}，请检查任务状态。"
                            )
                        elif initial_status:
                            logger.warning(
                                f"⚠ Task {task_id} 启动命令已发送，但状态未在 {max_verification_time}s 内变化。"
                                f"原始状态: {initial_status}，请检查任务状态。"
                            )
                    
                    return result
                    
                except ValueError:
                    raise  # Re-raise ValueError as-is
                except Exception as e:
                    logger.error(f"Unexpected error parsing response: {e}")
                    raise ValueError(f"无法解析服务器响应: {e}")
                    
            except requests.exceptions.Timeout as e:
                last_exception = e
                logger.error(f"Start training timeout for task {task_id} (attempt {attempt + 1}/{max_retries}): {e}")
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 2  # Exponential backoff
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time_module.sleep(wait_time)
                    continue
                else:
                    raise TimeoutError(
                        f"启动训练超时（超过{timeout}秒，已重试{max_retries}次）。\n"
                        f"可能原因：\n"
                        f"  1. COMPASS服务响应缓慢\n"
                        f"  2. 任务初始化时间过长\n"
                        f"  3. 网络连接问题\n"
                        f"请稍后检查任务状态。"
                    )
                    
            except requests.exceptions.HTTPError as e:
                last_exception = e
                if e.response is not None:
                    status_code = e.response.status_code
                    
                    # For 400 errors, don't retry (invalid request)
                    if status_code == 400:
                        try:
                            error_data = e.response.json()
                            error_msg = error_data.get('error', error_data.get('detail', str(e)))
                            error_code = error_data.get('error_code', '')
                            
                            # Check if task exists and get current status
                            diagnostic_info = ""
                            try:
                                task_status = self.get_task_status(task_id)
                                current_status = task_status.get('status')
                                diagnostic_info = f"\n当前任务状态: {current_status}"
                                if current_status == 'failed':
                                    diagnostic_info += f"\n错误信息: {task_status.get('error', 'N/A')}"
                            except Exception:
                                diagnostic_info = "\n无法获取任务状态信息"
                            
                            raise ValueError(
                                f"无法启动任务：{error_msg}\n"
                                f"错误代码: {error_code}{diagnostic_info}"
                            )
                        except ValueError:
                            raise  # Re-raise ValueError as-is
                        except Exception:
                            raise ValueError(f"无法启动任务（HTTP 400）：{str(e)}")
                    
                    # For 404, don't retry (task not found)
                    elif status_code == 404:
                        raise ValueError(f"任务不存在（task_id: {task_id}）")
                    
                    # For 503 (service unavailable), retry
                    elif status_code == 503:
                        if attempt < max_retries - 1:
                            wait_time = (attempt + 1) * 2
                            logger.warning(f"Service unavailable (503), retrying in {wait_time} seconds...")
                            time_module.sleep(wait_time)
                            continue
                        else:
                            raise ConnectionError(
                                f"COMPASS服务暂时不可用（已重试{max_retries}次）。"
                                f"请检查服务状态并稍后重试。"
                            )
                    
                    # For other HTTP errors, retry once
                    else:
                        if attempt < max_retries - 1:
                            wait_time = (attempt + 1) * 2
                            logger.warning(f"HTTP {status_code} error, retrying in {wait_time} seconds...")
                            time_module.sleep(wait_time)
                            continue
                        else:
                            raise ValueError(f"启动任务失败（HTTP {status_code}）：{str(e)}")
                else:
                    if attempt < max_retries - 1:
                        wait_time = (attempt + 1) * 2
                        logger.warning(f"HTTP error without response, retrying in {wait_time} seconds...")
                        time_module.sleep(wait_time)
                        continue
                    else:
                        raise ValueError(f"启动任务失败：{str(e)}")
                        
            except ConnectionError as e:
                # Connection errors should not be retried (service discovery issue)
                logger.error(f"Connection error during start_training: {e}")
                raise
                
            except Exception as e:
                last_exception = e
                logger.error(f"Unexpected error during start_training (attempt {attempt + 1}/{max_retries}): {e}", exc_info=True)
                if attempt < max_retries - 1:
                    wait_time = (attempt + 1) * 2
                    logger.info(f"Retrying in {wait_time} seconds...")
                    time_module.sleep(wait_time)
                    continue
                else:
                    raise ValueError(f"启动任务时发生未知错误（已重试{max_retries}次）：{str(e)}")
        
        # If we get here, all retries failed
        if last_exception:
            raise ValueError(f"启动任务失败（已重试{max_retries}次）：{str(last_exception)}")
        else:
            raise ValueError(f"启动任务失败：未知错误")
    
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
    
    def list_tasks(self, timeout: int = 10) -> List[Dict]:
        """
        List all tasks.
        
        Args:
            timeout: Request timeout in seconds (default: 10)
            
        Returns:
            List of task dictionaries
        """
        try:
            response = self._make_request("GET", "/api/v1/training/tasks", timeout=timeout)
            result = response.json()
            return result.get('tasks', [])
        except requests.exceptions.Timeout as e:
            logger.error(f"List tasks timeout: {e}")
            raise TimeoutError(f"获取任务列表超时（超过{timeout}秒）。请稍后重试。")
        except Exception as e:
            logger.error(f"Failed to list tasks: {e}")
            raise
    
    def get_task_logs(self, task_id: str, limit: int = 100) -> List[str]:
        """Get task logs."""
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}/logs", params={"limit": limit})
        result = response.json()
        return result['logs']
    
    def get_task_resources(self, task_id: str) -> Dict:
        """Get task resource usage (GPU, CPU, Memory)."""
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}/resources")
        return response.json()
    
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

