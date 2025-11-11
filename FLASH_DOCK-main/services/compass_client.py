"""
COMPASS service client for FLASH-DOCK.
"""

import sys
import logging
import requests
import json
import threading
import queue
from pathlib import Path
from typing import List, Dict, Optional, Any, Callable
from urllib.parse import urlparse

# Add parent directory to path to import services
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from registry_client import FlashDockRegistryClient
from load_balancer import LoadBalancer, LoadBalanceStrategy

logger = logging.getLogger(__name__)

# Try to import websockets, but make it optional
# Note: This is only for Python WebSocket client (TaskStreamClient)
# Frontend uses browser native WebSocket API, which doesn't require this library
try:
    import websockets

    WEBSOCKETS_AVAILABLE = True
    logger.debug(
        f"websockets library available (version: {getattr(websockets, '__version__', 'unknown')})"
    )
except ImportError:
    WEBSOCKETS_AVAILABLE = False
    logger.warning(
        "websockets library not available. Python WebSocket streaming will be disabled. "
        "Note: Frontend browser WebSocket connections are not affected."
    )


class CompassError(Exception):
    """Custom exception for COMPASS service errors with detailed error information."""

    def __init__(
        self,
        message: str,
        status_code: Optional[int] = None,
        error_code: Optional[str] = None,
        detail: Optional[Dict[str, Any]] = None,
        original_exception: Optional[Exception] = None,
    ):
        """
        Initialize CompassError.

        Args:
            message: Error message
            status_code: HTTP status code
            error_code: COMPASS error code (e.g., "ERR_1000")
            detail: Additional error details
            original_exception: Original exception that caused this error
        """
        super().__init__(message)
        self.message = message
        self.status_code = status_code
        self.error_code = error_code
        self.detail = detail or {}
        self.original_exception = original_exception

    def __str__(self) -> str:
        """Return formatted error message."""
        parts = [self.message]
        if self.error_code:
            parts.append(f"Error Code: {self.error_code}")
        if self.status_code:
            parts.append(f"HTTP Status: {self.status_code}")
        if self.detail:
            detail_str = json.dumps(self.detail, indent=2)
            parts.append(f"Details: {detail_str}")
        return "\n".join(parts)


class CompassClient:
    """Client for interacting with COMPASS service."""

    def __init__(
        self,
        registry_url: str = "http://localhost:8500",
        timeout: float = 30.0,
        load_balance_strategy: LoadBalanceStrategy = LoadBalanceStrategy.ROUND_ROBIN,
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
        self._cached_services: List[Dict[str, Any]] = []
        self._refresh_services()

    def _refresh_services(self):
        """Refresh the list of available COMPASS services."""
        try:
            self._cached_services = self.registry_client.discover_compass_services(
                healthy_only=True
            )
            logger.debug(
                f"Refreshed services cache: {len(self._cached_services)} services available"
            )
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
            raise ConnectionError(
                "没有可用的COMPASS服务。"
                f"请确保COMPASS服务已启动并注册到服务注册中心 ({self.registry_client.registry_url})。"
            )

        service = self.load_balancer.select_service(self._cached_services)
        if not service:
            raise ConnectionError(
                "无法选择COMPASS服务。"
                f"服务列表中有 {len(self._cached_services)} 个服务，但负载均衡器无法选择。"
            )

        return f"http://{service.host}:{service.port}"

    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
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
            CompassError: If request fails with COMPASS error details
            requests.RequestException: If request fails with non-COMPASS error
        """
        service_url = self._get_service_url()
        url = f"{service_url}{endpoint}"

        # Set default timeout if not provided
        if "timeout" not in kwargs:
            kwargs["timeout"] = self.timeout

        def _extract_error_details(response: requests.Response) -> Dict[str, Any]:
            """Extract error details from COMPASS error response."""
            error_details = {
                "status_code": response.status_code,
                "error_code": None,
                "error": None,
                "detail": None,
                "url": url,
                "method": method,
            }
            try:
                # Try to parse JSON error response
                if response.headers.get("content-type", "").startswith("application/json"):
                    error_json = response.json()
                    error_details["error"] = error_json.get("error") or error_json.get("message")
                    error_details["error_code"] = error_json.get("error_code")
                    error_details["detail"] = error_json.get("detail")
                    # Include additional fields if available
                    if "message" in error_json and not error_details["error"]:
                        error_details["error"] = error_json.get("message")
            except (json.JSONDecodeError, ValueError):
                # If JSON parsing fails, use response text
                response_text = response.text[:500] if response.text else "Unknown error"
                error_details["error"] = response_text
                # Try to extract useful information from non-JSON response
                if response_text and len(response_text) > 0:
                    error_details["detail"] = {"raw_response": response_text}
            return error_details

        def _make_single_request() -> requests.Response:
            """Make a single request and handle errors."""
            try:
                response = requests.request(method, url, **kwargs)
                response.raise_for_status()
                return response
            except requests.exceptions.HTTPError as e:
                # Extract error details before raising
                error_details = _extract_error_details(e.response)
                logger.error(
                    f"HTTP error {error_details['status_code']} from {url}: "
                    f"{error_details.get('error', str(e))}"
                )
                raise CompassError(
                    message=error_details.get("error")
                    or f"HTTP {error_details['status_code']} error",
                    status_code=error_details["status_code"],
                    error_code=error_details.get("error_code"),
                    detail=error_details.get("detail"),
                    original_exception=e,
                )
            except requests.exceptions.ConnectionError as e:
                # Connection errors - provide more context
                logger.error(f"Connection failed to {url}: {e}")
                raise ConnectionError(
                    f"无法连接到COMPASS服务: {e}. " f"请确保服务正在运行并可通过 {url} 访问。"
                ) from e
            except requests.exceptions.Timeout as e:
                # Timeout errors - provide more context
                timeout_value = kwargs.get("timeout", self.timeout)
                logger.error(f"Request timeout to {url} (timeout={timeout_value}s): {e}")
                raise requests.exceptions.Timeout(
                    f"请求超时 (超时时间: {timeout_value}秒). " f"服务可能响应缓慢或不可用: {url}"
                ) from e
            except requests.exceptions.RequestException as e:
                logger.error(f"Request failed to {url}: {e}")
                raise

        try:
            return _make_single_request()
        except (CompassError, ConnectionError):
            # Don't retry CompassError or ConnectionError
            raise
        except requests.exceptions.RequestException as e:
            # Try refreshing services and retry once for other request errors
            logger.warning(f"Retrying request after refreshing services: {url}")
            self._refresh_services()
            try:
                service_url = self._get_service_url()
                url = f"{service_url}{endpoint}"
                return _make_single_request()
            except requests.exceptions.RequestException as retry_e:
                # If retry also fails, check if it's an HTTP error with details
                if isinstance(retry_e, requests.exceptions.HTTPError) and retry_e.response:
                    error_details = _extract_error_details(retry_e.response)
                    raise CompassError(
                        message=error_details.get("error")
                        or f"HTTP {error_details['status_code']} error",
                        status_code=error_details["status_code"],
                        error_code=error_details.get("error_code"),
                        detail=error_details.get("detail"),
                        original_exception=retry_e,
                    )
                raise

    def upload_dataset(
        self, file_path: str, name: Optional[str] = None, description: Optional[str] = None
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

        with open(file_path, "rb") as f:
            files = {"file": (file_path_obj.name, f, "application/octet-stream")}
            data = {}
            if name:
                data["name"] = name
            if description:
                data["description"] = description

            response = self._make_request("POST", "/api/v1/data/upload", files=files, data=data)
            result = response.json()
            if not isinstance(result, dict) or "dataset_id" not in result:
                raise CompassError(
                    "Invalid response format: missing 'dataset_id'",
                    detail={"response": result},
                )
            return result["dataset_id"]

    def list_datasets(self) -> List[Dict]:
        """
        List all datasets.

        Returns:
            List[Dict]: List of dataset dictionaries with keys:
                dataset_id, name, description, size, file_count, status,
                created_at, updated_at, metadata, error
        """
        response = self._make_request("GET", "/api/v1/data/datasets")
        result = response.json()

        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )

        # Convert datasets to list of dicts with datetime serialization
        datasets = []
        datasets_list = result.get("datasets", [])
        if not isinstance(datasets_list, list):
            logger.warning("Expected 'datasets' to be a list, got %s", type(datasets_list))
            datasets_list = []

        for ds in datasets_list:
            if not isinstance(ds, dict):
                logger.warning("Skipping invalid dataset entry: %s", ds)
                continue
            # Validate required fields
            if "dataset_id" not in ds:
                logger.warning("Skipping dataset entry missing 'dataset_id': %s", ds)
                continue
            if "name" not in ds:
                logger.warning("Skipping dataset entry missing 'name': %s", ds)
                continue
            dataset_dict = {
                "dataset_id": ds["dataset_id"],
                "name": ds["name"],
                "description": ds.get("description"),
                "size": ds.get("size", 0),
                "file_count": ds.get("file_count", 0),
                "status": ds.get("status", "pending"),
                "created_at": ds.get("created_at"),
                "updated_at": ds.get("updated_at"),
                "metadata": ds.get("metadata", {}),
                "error": ds.get("error"),
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
        self._make_request("DELETE", f"/api/v1/data/datasets/{dataset_id}")

    def get_inference_status(self) -> Dict:
        """
        Get inference service status.

        Returns:
            Dict: Status dictionary with keys:
                status, loaded_models, max_cache_size, cached_model_ids,
                latest_model, device
        """
        response = self._make_request("GET", "/api/v1/inference/status")
        result: Dict[str, Any] = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        return result

    def list_models(self) -> List[Dict]:
        """
        List all available models.

        Returns:
            List[Dict]: List of model dictionaries with keys:
                model_id, name, version, task_id, checkpoint_path,
                file_size, created_at, metadata, metrics
        """
        response = self._make_request("GET", "/api/v1/models")
        result = response.json()

        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )

        # Convert models to list of dicts
        models = []
        models_list = result.get("models", [])
        if not isinstance(models_list, list):
            logger.warning("Expected 'models' to be a list, got %s", type(models_list))
            models_list = []

        for model in models_list:
            if not isinstance(model, dict):
                logger.warning("Skipping invalid model entry: %s", model)
                continue
            # Validate required fields
            if "model_id" not in model:
                logger.warning("Skipping model entry missing 'model_id': %s", model)
                continue
            if "name" not in model:
                logger.warning("Skipping model entry missing 'name': %s", model)
                continue
            if "version" not in model:
                logger.warning("Skipping model entry missing 'version': %s", model)
                continue
            model_dict = {
                "model_id": model["model_id"],
                "name": model["name"],
                "version": model["version"],
                "task_id": model.get("task_id"),
                "checkpoint_path": model.get("checkpoint_path"),
                "file_size": model.get("file_size", 0),
                "created_at": model.get("created_at"),
                "metadata": model.get("metadata", {}),
                "metrics": model.get("metrics", {}),
            }
            models.append(model_dict)

        return models

    def create_training_task(
        self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None
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
        data = {"config": config}
        if dataset_id:
            data["dataset_id"] = dataset_id
        if description:
            data["description"] = description

        response = self._make_request("POST", "/api/v1/training/tasks", json=data)
        result: Dict[str, Any] = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        # Validate required fields
        if "task_id" not in result:
            raise CompassError(
                "Invalid response format: missing 'task_id'",
                detail={"response": result},
            )
        return result

    def list_training_tasks(self) -> List[Dict]:
        """
        List all training tasks.

        Returns:
            List[Dict]: List of task dictionaries with keys:
                task_id, status, config, created_at, updated_at,
                started_at, completed_at, progress, error, description
        """
        response = self._make_request("GET", "/api/v1/training/tasks")
        result = response.json()

        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )

        # Convert tasks to list of dicts
        tasks = []
        tasks_list = result.get("tasks", [])
        if not isinstance(tasks_list, list):
            logger.warning("Expected 'tasks' to be a list, got %s", type(tasks_list))
            tasks_list = []

        for task in tasks_list:
            if not isinstance(task, dict):
                logger.warning("Skipping invalid task entry: %s", task)
                continue
            # Validate required fields
            if "task_id" not in task:
                logger.warning("Skipping task entry missing 'task_id': %s", task)
                continue
            if "status" not in task:
                logger.warning("Skipping task entry missing 'status': %s", task)
                continue
            task_dict = {
                "task_id": task["task_id"],
                "status": task["status"],
                "config": task.get("config", {}),
                "created_at": task.get("created_at"),
                "updated_at": task.get("updated_at"),
                "started_at": task.get("started_at"),
                "completed_at": task.get("completed_at"),
                "progress": task.get("progress", {}),
                "error": task.get("error"),
                "description": task.get("description"),
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
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}")
        task = response.json()

        if not isinstance(task, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": task},
            )

        # Validate required fields
        if "task_id" not in task:
            raise CompassError(
                "Invalid response format: missing 'task_id'",
                detail={"response": task},
            )
        if "status" not in task:
            raise CompassError(
                "Invalid response format: missing 'status'",
                detail={"response": task},
            )

        return {
            "task_id": task["task_id"],
            "status": task["status"],
            "config": task.get("config", {}),
            "created_at": task.get("created_at"),
            "updated_at": task.get("updated_at"),
            "started_at": task.get("started_at"),
            "completed_at": task.get("completed_at"),
            "progress": task.get("progress", {}),
            "error": task.get("error"),
            "description": task.get("description"),
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
        response = self._make_request("POST", f"/api/v1/training/tasks/{task_id}/start")
        result: Dict[str, Any] = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        return result

    def stop_training_task(self, task_id: str, timeout: Optional[float] = None) -> Dict:
        """
        Stop a running training task.

        Args:
            task_id: Unique task identifier (UUID)
            timeout: Request timeout in seconds (default: 10.0, which is sufficient since
                     the backend now returns immediately and monitors cancellation in background)

        Returns:
            Dict: Confirmation message with task_id

        Raises:
            requests.RequestException: If task cannot be stopped
        """
        # Use shorter timeout since backend now returns immediately (within 2 seconds for quick response)
        # The backend monitors cancellation in background, so we don't need to wait 30+ seconds
        request_timeout = timeout if timeout is not None else 10.0
        response = self._make_request("POST", f"/api/v1/training/tasks/{task_id}/stop", timeout=request_timeout)
        result: Dict[str, Any] = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        return result

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
        response = self._make_request("POST", f"/api/v1/training/tasks/{task_id}/pause")
        result = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        return result

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
        params = {"limit": min(limit, 1000)}
        response = self._make_request(
            "GET", f"/api/v1/training/tasks/{task_id}/logs", params=params
        )
        result = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        logs = result.get("logs", [])
        if not isinstance(logs, list):
            logger.warning("Expected 'logs' to be a list, got %s", type(logs))
            logs = []
        return logs

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
        response = self._make_request("GET", f"/api/v1/training/tasks/{task_id}/metrics")
        result = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        metrics = result.get("metrics", {})
        if not isinstance(metrics, dict):
            logger.warning("Expected 'metrics' to be a dictionary, got %s", type(metrics))
            metrics = {}
        return metrics

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
        response = self._make_request("DELETE", f"/api/v1/training/tasks/{task_id}")
        result: Dict[str, Any] = response.json()
        if not isinstance(result, dict):
            raise CompassError(
                "Invalid response format: expected dictionary",
                detail={"response": result},
            )
        return result

    def stream_task_logs(
        self,
        task_id: str,
        on_log: Optional[Callable[[str], None]] = None,
        on_resources: Optional[Callable[[Dict], None]] = None,
        on_error: Optional[Callable[[Exception], None]] = None,
        on_connect: Optional[Callable[[], None]] = None,
    ) -> "TaskStreamClient":
        """
        Create a WebSocket client for streaming task logs and resources.

        Args:
            task_id: Task ID to stream
            on_log: Callback function for log messages (receives log data string)
            on_resources: Callback function for resource updates (receives resource dict)
            on_error: Callback function for errors (receives Exception)
            on_connect: Callback function when connected

        Returns:
            TaskStreamClient: Stream client instance
        """
        if not WEBSOCKETS_AVAILABLE:
            raise ImportError(
                "websockets library is required for streaming. Install it with: pip install websockets"
            )

        service_url = self._get_service_url()
        # Convert HTTP URL to WebSocket URL
        parsed = urlparse(service_url)
        # Use wss:// for HTTPS, ws:// for HTTP
        ws_scheme = "wss" if parsed.scheme == "https" else "ws"
        ws_url = f"{ws_scheme}://{parsed.netloc}/api/v1/training/tasks/{task_id}/stream"

        return TaskStreamClient(
            ws_url=ws_url,
            on_log=on_log,
            on_resources=on_resources,
            on_error=on_error,
            on_connect=on_connect,
        )


class TaskStreamClient:
    """WebSocket client for streaming task logs and resources."""

    def __init__(
        self,
        ws_url: str,
        on_log: Optional[Callable[[str], None]] = None,
        on_resources: Optional[Callable[[Dict], None]] = None,
        on_error: Optional[Callable[[Exception], None]] = None,
        on_connect: Optional[Callable[[], None]] = None,
    ):
        """
        Initialize stream client.

        Args:
            ws_url: WebSocket URL
            on_log: Callback for log messages
            on_resources: Callback for resource updates
            on_error: Callback for errors
            on_connect: Callback when connected
        """
        self.ws_url = ws_url
        self.on_log = on_log
        self.on_resources = on_resources
        self.on_error = on_error
        self.on_connect = on_connect

        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._ws = None
        self._message_queue: queue.Queue = queue.Queue()
        self._reconnect_delay = 1.0
        self._max_reconnect_delay = 60.0

    def start(self):
        """Start the WebSocket client in a background thread."""
        if self._running:
            return

        self._running = True
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        """Stop the WebSocket client."""
        self._running = False
        if self._ws:
            try:
                # Close WebSocket connection
                import asyncio

                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                loop.run_until_complete(self._ws.close())
            except Exception:
                pass

        if self._thread:
            self._thread.join(timeout=5)

    def _run(self):
        """Run WebSocket client in background thread."""
        import asyncio

        # Create new event loop for this thread
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

        while self._running:
            try:
                loop.run_until_complete(self._connect_and_listen())
            except Exception as e:
                if self.on_error:
                    try:
                        self.on_error(e)
                    except Exception:
                        pass

                if self._running:
                    # Exponential backoff for reconnection
                    logger.warning(
                        f"WebSocket connection lost, reconnecting in {self._reconnect_delay}s..."
                    )
                    import time

                    time.sleep(self._reconnect_delay)
                    self._reconnect_delay = min(
                        self._reconnect_delay * 2, self._max_reconnect_delay
                    )
                else:
                    break

        loop.close()

    def _handle_message(self, message: str) -> None:
        """Handle a single WebSocket message."""
        try:
            data = json.loads(message)
            msg_type = data.get("type")
            msg_data = data.get("data")

            if msg_type == "log" and self.on_log:
                self._call_callback(self.on_log, msg_data, "on_log")
            elif msg_type == "resources" and self.on_resources:
                self._call_callback(self.on_resources, msg_data, "on_resources")
            elif msg_type == "connected":
                logger.info("Connected to stream: %s", msg_data)
            elif msg_type == "error":
                logger.error("Stream error: %s", msg_data)
                if self.on_error:
                    self._call_callback(self.on_error, Exception(msg_data), "on_error")

        except json.JSONDecodeError:
            logger.warning("Invalid JSON message: %s", message)
        except Exception as e:
            logger.error("Error processing message: %s", e)
            if self.on_error:
                self._call_callback(self.on_error, e, "on_error")

    def _call_callback(
        self, callback: Optional[Callable[[Any], None]], data: Any, callback_name: str
    ) -> None:
        """Safely call a callback function."""
        if callback is None:
            return
        try:
            if data is None:
                # For callbacks that don't take arguments (like on_connect)
                callback()  # type: ignore
            else:
                callback(data)
        except Exception as e:
            logger.error("Error in %s callback: %s", callback_name, e)

    async def _connect_and_listen(self):
        """Connect to WebSocket and listen for messages."""
        try:
            async with websockets.connect(self.ws_url) as ws:
                self._ws = ws
                self._reconnect_delay = 1.0  # Reset reconnect delay on successful connection

                if self.on_connect:
                    try:
                        self.on_connect()
                    except Exception as e:
                        logger.error("Error in on_connect callback: %s", e)

                # Listen for messages
                async for message in ws:
                    if not self._running:
                        break
                    self._handle_message(message)

        except websockets.exceptions.ConnectionClosed:
            if self._running:
                raise  # Re-raise to trigger reconnection
        except Exception as e:
            if self._running:
                logger.error("WebSocket error: %s", e)
                raise  # Re-raise to trigger reconnection

    def send_command(self, command: str):
        """
        Send a command to the server (for future use).

        Args:
            command: Command string to send
        """
        if not self._ws:
            logger.warning("WebSocket not connected, cannot send command")
            return

        try:
            import asyncio

            # This would need to be called from the async context
            # For now, we'll just log it
            logger.info("Command to send: %s", command)
        except Exception as e:
            logger.error(f"Error sending command: {e}")

    def is_connected(self) -> bool:
        """Check if WebSocket is connected."""
        return self._ws is not None and self._running
