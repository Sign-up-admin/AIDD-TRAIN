"""
Training service implementation.
"""
import os
import sys
import io
import uuid
import threading
import logging
import time
import traceback
from typing import Dict, Optional, List
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

try:
    import psutil
except ImportError:
    psutil = None

try:
    import torch
except ImportError:
    torch = None

from compass.service.models.task import TaskStatus, TrainingTaskResponse
from compass.service.config import SERVICE_CONFIG
from compass.logger import TrainingLogger
from compass.main import main
from compass.config import get_config
from compass.service.services.progress_tracker import ProgressTracker
from compass.service.services.progress_logger import ProgressAwareLogger
from compass.service.services.task_lifecycle_logger import TaskLifecycleLogger
from compass.service.exceptions import ServiceException, ValidationError
from compass.service.error_codes import ErrorCode

logger = logging.getLogger(__name__)

# Log warning about psutil after logger is initialized
if psutil is None:
    logger.warning("psutil not available, resource monitoring will be limited")


def get_resource_usage() -> Dict:
    """
    Get current resource usage (GPU, CPU, Memory).
    
    Returns:
        Dict: Resource usage information
    """
    resources = {
        "cpu_percent": 0.0,
        "memory": {
            "total_gb": 0.0,
            "used_gb": 0.0,
            "percent": 0.0
        },
        "gpu": None
    }
    
    # Get CPU and memory info if psutil is available
    if psutil:
        try:
            resources["cpu_percent"] = psutil.cpu_percent(interval=0.1)
            mem = psutil.virtual_memory()
            resources["memory"] = {
                "total_gb": mem.total / (1024**3),
                "used_gb": mem.used / (1024**3),
                "percent": mem.percent
            }
        except Exception as e:
            logger.warning(f"Failed to get CPU/memory info: {e}")
    
    # Get GPU information if available
    if torch and torch.cuda.is_available():
        try:
            allocated = torch.cuda.memory_allocated() / (1024**2)  # MB
            reserved = torch.cuda.memory_reserved() / (1024**2)  # MB
            total = torch.cuda.get_device_properties(0).total_memory / (1024**2)  # MB
            
            resources["gpu"] = {
                "available": True,
                "device_name": torch.cuda.get_device_name(0),
                "memory": {
                    "allocated_mb": allocated,
                    "reserved_mb": reserved,
                    "total_mb": total,
                    "allocated_percent": (allocated / total * 100) if total > 0 else 0,
                    "reserved_percent": (reserved / total * 100) if total > 0 else 0
                }
            }
        except Exception as e:
            logger.warning(f"Failed to get GPU info: {e}")
            resources["gpu"] = {"available": False, "error": str(e)}
    else:
        resources["gpu"] = {"available": False}
    
    return resources


class TrainingService:
    """Service for managing training tasks."""
    
    def __init__(self):
        """Initialize training service."""
        self.tasks: Dict[str, TrainingTaskResponse] = {}
        self.task_logs: Dict[str, List[str]] = {}
        self.executor = ThreadPoolExecutor(max_workers=SERVICE_CONFIG['max_workers'])
        self.task_threads: Dict[str, threading.Thread] = {}
        self.progress_trackers: Dict[str, ProgressTracker] = {}
        self.task_resources: Dict[str, Dict] = {}  # Store resource usage for each task
        self.resource_monitor_threads: Dict[str, threading.Thread] = {}  # Track resource monitor threads
        self.task_lifecycle_loggers: Dict[str, TaskLifecycleLogger] = {}  # Track lifecycle loggers
        self.lock = threading.Lock()
        
        # Initialize shutdown flag BEFORE starting health monitor thread
        self._shutdown_flag = threading.Event()
        
        # Start health monitoring thread
        self.health_monitor_thread = threading.Thread(target=self._health_monitor_loop, daemon=True)
        self.health_monitor_thread.start()
    
    def _check_resources(self) -> Dict:
        """
        Check if system has enough resources to create a new task.
        
        Returns:
            Dict with check results and any issues
        """
        result = {
            "available": True,
            "issues": [],
            "resources": {}
        }
        
        enable_check = SERVICE_CONFIG.get('enable_resource_check', True)
        if not enable_check:
            return result
        
        # Get current resource usage
        resources = get_resource_usage()
        result["resources"] = resources
        
        # Check concurrent tasks limit
        max_concurrent = SERVICE_CONFIG.get('max_concurrent_tasks', 4)
        running_count = sum(1 for t in self.tasks.values() if t.status == TaskStatus.RUNNING)
        if running_count >= max_concurrent:
            result["available"] = False
            result["issues"].append(f"Maximum concurrent tasks ({max_concurrent}) reached")
        
        # Check CPU usage
        if psutil:
            max_cpu = SERVICE_CONFIG.get('max_cpu_percent', 90.0)
            cpu_percent = resources.get('cpu_percent', 0)
            if cpu_percent > max_cpu:
                result["available"] = False
                result["issues"].append(f"CPU usage ({cpu_percent:.1f}%) exceeds limit ({max_cpu}%)")
        
        # Check memory availability
        if psutil:
            min_memory_gb = SERVICE_CONFIG.get('min_available_memory_gb', 2.0)
            mem_info = resources.get('memory', {})
            available_gb = mem_info.get('total_gb', 0) - mem_info.get('used_gb', 0)
            if available_gb < min_memory_gb:
                result["available"] = False
                result["issues"].append(f"Available memory ({available_gb:.2f} GB) below minimum ({min_memory_gb} GB)")
        
        # Check GPU availability if needed (could be enhanced based on config)
        if resources.get('gpu') and resources['gpu'].get('available'):
            gpu_mem = resources['gpu'].get('memory', {})
            allocated_percent = gpu_mem.get('allocated_percent', 0)
            if allocated_percent > 95:
                result["available"] = False
                result["issues"].append(f"GPU memory usage ({allocated_percent:.1f}%) too high")
        
        return result
    
    def create_task(self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None) -> str:
        """
        Create a new training task with resource validation and lifecycle tracking.
        
        Args:
            config: Training configuration
            dataset_id: Optional dataset ID
            description: Optional task description
            
        Returns:
            str: Task ID
            
        Raises:
            ServiceException: If resource check fails or creation times out
        """
        task_id = str(uuid.uuid4())
        creation_start = datetime.now()
        creation_timeout = SERVICE_CONFIG.get('task_creation_timeout', 30)
        
        # Initialize lifecycle logger
        lifecycle_logger = TaskLifecycleLogger(task_id, SERVICE_CONFIG.get('log_dir', 'logs'))
        self.task_lifecycle_loggers[task_id] = lifecycle_logger
        
        # Log creation start
        resource_check = self._check_resources()
        lifecycle_logger.log_create_start(config, resource_check)
        
        # Check if creation is taking too long
        if (datetime.now() - creation_start).total_seconds() > creation_timeout:
            lifecycle_logger.log_create_failure(
                TimeoutError(f"Task creation timeout after {creation_timeout}s"),
                "Creation timeout"
            )
            raise ServiceException(
                f"Task creation timeout after {creation_timeout} seconds",
                error_code=ErrorCode.TIMEOUT
            )
        
        # Perform resource check
        if not resource_check.get("available", True):
            issues = resource_check.get("issues", [])
            error_msg = "Resource check failed: " + "; ".join(issues)
            lifecycle_logger.log_create_failure(
                ValueError(error_msg),
                "Resource check failed"
            )
            raise ServiceException(
                error_msg,
                error_code=ErrorCode.SERVICE_UNAVAILABLE,
                detail={"resource_check": resource_check}
            )
        
        # Create task with CREATING status
        task = TrainingTaskResponse(
            task_id=task_id,
            status=TaskStatus.CREATING,
            config=config,
            created_at=datetime.now(),
            updated_at=datetime.now(),
            description=description
        )
        
        try:
            with self.lock:
                self.tasks[task_id] = task
                self.task_logs[task_id] = []
                self.progress_trackers[task_id] = ProgressTracker(task_id)
                self.task_resources[task_id] = {}
            
            # Update status to PENDING after successful creation
            with self.lock:
                task.status = TaskStatus.PENDING
                task.updated_at = datetime.now()
            
            lifecycle_logger.log_state_transition(TaskStatus.CREATING.value, TaskStatus.PENDING.value, "Task created")
            lifecycle_logger.log_create_success()
            
            logger.info(f"Created training task: {task_id}")
            return task_id
            
        except Exception as e:
            # Clean up on failure
            with self.lock:
                if task_id in self.tasks:
                    del self.tasks[task_id]
                if task_id in self.task_logs:
                    del self.task_logs[task_id]
                if task_id in self.progress_trackers:
                    del self.progress_trackers[task_id]
                if task_id in self.task_resources:
                    del self.task_resources[task_id]
                if task_id in self.task_lifecycle_loggers:
                    del self.task_lifecycle_loggers[task_id]
            
            lifecycle_logger.log_create_failure(e, "Exception during task creation")
            raise ServiceException(
                f"Failed to create task: {str(e)}",
                error_code=ErrorCode.INTERNAL_ERROR
            ) from e
    
    def get_task(self, task_id: str) -> Optional[TrainingTaskResponse]:
        """
        Get task by ID.
        
        Args:
            task_id: Task ID
            
        Returns:
            Optional[TrainingTaskResponse]: Task or None if not found
        """
        task = self.tasks.get(task_id)
        if task and task_id in self.progress_trackers:
            # Update progress in task
            progress = self.progress_trackers[task_id].get_progress()
            task.progress = progress
        return task
    
    def list_tasks(self) -> List[TrainingTaskResponse]:
        """
        List all tasks.
        
        Returns:
            List[TrainingTaskResponse]: List of tasks with updated progress
        """
        tasks = list(self.tasks.values())
        # Update progress for all tasks
        for task in tasks:
            if task.task_id in self.progress_trackers:
                task.progress = self.progress_trackers[task.task_id].get_progress()
        return tasks
    
    def start_task(self, task_id: str) -> bool:
        """
        Start a training task with initialization tracking and thread management.
        
        Args:
            task_id: Task ID
            
        Returns:
            bool: True if started successfully
        """
        lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
        init_timeout = SERVICE_CONFIG.get('task_initialization_timeout', 60)
        
        with self.lock:
            if task_id not in self.tasks:
                if lifecycle_logger:
                    lifecycle_logger.log_initialize_failure(
                        ValueError("Task not found"),
                        "Task not found"
                    )
                return False
            
            task = self.tasks[task_id]
            if task.status != TaskStatus.PENDING and task.status != TaskStatus.PAUSED:
                if lifecycle_logger:
                    lifecycle_logger.log_initialize_failure(
                        ValueError(f"Task in invalid state: {task.status}"),
                        f"Invalid state: {task.status}"
                    )
                return False
            
            # Update to INITIALIZING status
            task.status = TaskStatus.INITIALIZING
            task.updated_at = datetime.now()
            if lifecycle_logger:
                lifecycle_logger.log_state_transition(TaskStatus.PENDING.value, TaskStatus.INITIALIZING.value, "Starting task initialization")
                lifecycle_logger.log_initialize_start()
        
        # Start training in background thread
        init_start_time = datetime.now()
        
        def run_training():
            try:
                # Check initialization timeout
                if (datetime.now() - init_start_time).total_seconds() > init_timeout:
                    if lifecycle_logger:
                        lifecycle_logger.log_timeout("initializing", init_timeout)
                    raise TimeoutError(f"Task initialization timeout after {init_timeout}s")
                
                # Update to RUNNING status before starting
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.RUNNING
                        self.tasks[task_id].started_at = datetime.now()
                        self.tasks[task_id].updated_at = datetime.now()
                        if lifecycle_logger:
                            lifecycle_logger.log_state_transition(TaskStatus.INITIALIZING.value, TaskStatus.RUNNING.value, "Initialization complete")
                            lifecycle_logger.log_initialize_success()
                            lifecycle_logger.log_start()
                
                self._run_training(task_id)
            except Exception as e:
                logger.error(f"Training task {task_id} failed: {e}", exc_info=True)
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.FAILED
                        self.tasks[task_id].error = str(e)
                        self.tasks[task_id].updated_at = datetime.now()
                        self._log(task_id, f"ERROR: {e}")
                        if lifecycle_logger:
                            lifecycle_logger.log_failure(e, "running", {"error": str(e)})
            finally:
                # Clean up thread references when task completes
                with self.lock:
                    if task_id in self.task_threads:
                        thread = self.task_threads.pop(task_id)
                        logger.debug(f"Cleaned up thread for task {task_id}")
                    # Clean up resource monitor thread
                    if task_id in self.resource_monitor_threads:
                        monitor_thread = self.resource_monitor_threads.pop(task_id)
                        logger.debug(f"Cleaned up resource monitor thread for task {task_id}")
                    # Note: Keep progress tracker for querying completed tasks
                    # Only remove if task is deleted
        
        try:
            thread = threading.Thread(target=run_training, daemon=True, name=f"TrainingThread-{task_id}")
            thread.start()
            self.task_threads[task_id] = thread
            
            logger.info(f"Started training task: {task_id}")
            return True
        except Exception as e:
            # Rollback on thread creation failure
            with self.lock:
                if task_id in self.tasks:
                    self.tasks[task_id].status = TaskStatus.FAILED
                    self.tasks[task_id].error = f"Failed to create thread: {str(e)}"
                    self.tasks[task_id].updated_at = datetime.now()
                if lifecycle_logger:
                    lifecycle_logger.log_initialize_failure(e, "Thread creation failed")
            logger.error(f"Failed to create thread for task {task_id}: {e}", exc_info=True)
            return False
    
    def _run_training(self, task_id: str):
        """Run training task."""
        task = self.tasks[task_id]
        config = task.config
        
        # Get progress tracker
        progress_tracker = self.progress_trackers.get(task_id)
        if not progress_tracker:
            progress_tracker = ProgressTracker(task_id)
            with self.lock:
                self.progress_trackers[task_id] = progress_tracker
        
        progress_tracker.set_stage("initializing", "Initializing training environment")
        
        # Create log directory
        log_dir = os.path.join(SERVICE_CONFIG['log_dir'], f"task_{task_id}")
        os.makedirs(log_dir, exist_ok=True)
        
        # Create progress-aware logger
        logger_instance = ProgressAwareLogger(log_dir=log_dir, progress_tracker=progress_tracker)
        
        # Capture stdout/stderr for console output
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        
        # Create a custom stdout that writes to both console and capture
        class TeeOutput:
            def __init__(self, stream, task_id, original_stream, log_func):
                self.stream = stream
                self.task_id = task_id
                self.original_stream = original_stream
                self.log_func = log_func
                self.buffer = ""
                
            def write(self, text):
                # Write to original stream immediately
                self.original_stream.write(text)
                self.original_stream.flush()
                
                # Buffer the text
                self.buffer += text
                
                # Process complete lines
                while '\n' in self.buffer:
                    line, self.buffer = self.buffer.split('\n', 1)
                    if line.strip():  # Only log non-empty lines
                        self.log_func(self.task_id, line.strip())
                    
            def flush(self):
                self.original_stream.flush()
                # Flush any remaining buffer as a line
                if self.buffer.strip():
                    self.log_func(self.task_id, self.buffer.strip())
                    self.buffer = ""
                
            def __getattr__(self, name):
                return getattr(self.original_stream, name)
        
        # Replace stdout and stderr
        sys.stdout = TeeOutput(stdout_capture, task_id, old_stdout, self._log)
        sys.stderr = TeeOutput(stderr_capture, task_id, old_stderr, self._log)
        
        # Start resource monitoring thread
        lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
        
        def monitor_resources():
            """Monitor resource usage periodically."""
            while task_id in self.tasks and self.tasks[task_id].status == TaskStatus.RUNNING:
                try:
                    resources = get_resource_usage()
                    with self.lock:
                        if task_id in self.task_resources:
                            self.task_resources[task_id] = resources
                    # Log resource usage periodically (every 30 seconds)
                    if lifecycle_logger and int(time.time()) % 30 == 0:
                        lifecycle_logger.log_resource_usage(resources)
                except Exception as e:
                    logger.warning(f"Failed to get resources for task {task_id}: {e}")
                time.sleep(5)  # Update every 5 seconds
        
        resource_monitor_thread = threading.Thread(
            target=monitor_resources, 
            daemon=True, 
            name=f"ResourceMonitor-{task_id}"
        )
        resource_monitor_thread.start()
        with self.lock:
            self.resource_monitor_threads[task_id] = resource_monitor_thread
        
        # Update config with task-specific settings
        training_config = get_config(config.get('execution_mode', 'validation_tuned'))
        training_config.update(config)
        training_config['log_dir'] = log_dir
        training_config['checkpoint_dir'] = os.path.join(
            SERVICE_CONFIG['checkpoint_dir'],
            f"task_{task_id}"
        )
        
        # Use default dataset if no dataset_id provided
        # The default paths are already set in get_config(), but we can override if needed
        if not config.get('dataset_id'):
            # Check if default paths exist (relative to project root)
            project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
            default_index = os.path.join(project_root, 'index', 'INDEX_general_PL.2020R1.lst')
            default_dataset = os.path.join(project_root, 'PDBbind-2025.8.4', 'P-L')
            
            # Use relative paths (as in config.py) if they exist, otherwise try absolute paths
            if not os.path.exists(default_index):
                default_index = training_config.get('index_file', 'index/INDEX_general_PL.2020R1.lst')
            else:
                training_config['index_file'] = default_index
                
            if not os.path.exists(default_dataset):
                default_dataset = training_config.get('dataset_path', 'PDBbind-2025.8.4/P-L/')
            else:
                training_config['dataset_path'] = default_dataset
            
            logger_instance.log(f"Using default PDBbind dataset: {default_dataset}")
            logger_instance.log(f"Using index file: {default_index}")
        
        self._log(task_id, f"Starting training with config: {config}")
        
        # Run training
        try:
            # Update progress before starting
            progress_tracker.set_stage("data_processing", "Starting data processing")
            task.progress = progress_tracker.get_progress()
            
            # Run main training
            lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
            if lifecycle_logger:
                lifecycle_logger.log_running("training", progress_tracker.get_progress())
            
            main(training_config, logger_instance)
            
            with self.lock:
                if task_id in self.tasks:
                    lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
                    # Check if cancelled before marking as completed
                    if progress_tracker.is_cancelled():
                        self.tasks[task_id].status = TaskStatus.CANCELLED
                        self.tasks[task_id].updated_at = datetime.now()
                        self._log(task_id, "Training cancelled")
                        if lifecycle_logger:
                            lifecycle_logger.log_cancel("Training cancelled by user")
                    else:
                        self.tasks[task_id].status = TaskStatus.COMPLETED
                        self.tasks[task_id].completed_at = datetime.now()
                        self.tasks[task_id].updated_at = datetime.now()
                        progress_tracker.set_completed("Training completed successfully")
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        self._log(task_id, "Training completed successfully")
                        if lifecycle_logger:
                            lifecycle_logger.log_complete("Training completed successfully")
        except Exception as e:
            # Check if it's a cancellation exception
            from compass.training.exceptions import TrainingCancelled
            if isinstance(e, TrainingCancelled):
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.CANCELLED
                        self.tasks[task_id].updated_at = datetime.now()
                        progress_tracker.set_stage("cancelled", "Training cancelled by user")
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        self._log(task_id, "Training cancelled")
                logger.info(f"Training task {task_id} was cancelled")
            else:
                with self.lock:
                    if task_id in self.tasks:
                        progress_tracker.set_stage("failed", f"Training failed: {str(e)}")
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
                        if lifecycle_logger:
                            lifecycle_logger.log_failure(e, "training", {"error": str(e), "traceback": traceback.format_exc()})
                raise
        finally:
            # Restore stdout/stderr
            sys.stdout = old_stdout
            sys.stderr = old_stderr
    
    def stop_task(self, task_id: str) -> bool:
        """
        Stop a training task.
        
        Args:
            task_id: Task ID
            
        Returns:
            bool: True if stopped successfully
        """
        with self.lock:
            if task_id not in self.tasks:
                return False
            
            task = self.tasks[task_id]
            if task.status != TaskStatus.RUNNING:
                return False
            
            # Set cancellation flag in progress tracker
            if task_id in self.progress_trackers:
                self.progress_trackers[task_id].cancel()
                logger.info(f"Set cancellation flag for task {task_id}")
            
            task.status = TaskStatus.CANCELLED
            task.updated_at = datetime.now()
        
        self._log(task_id, "Training task cancellation requested")
        logger.info(f"Stopped training task: {task_id}")
        return True
    
    def pause_task(self, task_id: str) -> bool:
        """
        Pause a training task.
        
        Args:
            task_id: Task ID
            
        Returns:
            bool: True if paused successfully
        """
        with self.lock:
            if task_id not in self.tasks:
                return False
            
            task = self.tasks[task_id]
            if task.status != TaskStatus.RUNNING:
                return False
            
            task.status = TaskStatus.PAUSED
            task.updated_at = datetime.now()
        
        self._log(task_id, "Training task paused")
        logger.info(f"Paused training task: {task_id}")
        return True
    
    def delete_task(self, task_id: str) -> bool:
        """
        Delete a training task and clean up all resources.
        
        Args:
            task_id: Task ID
            
        Returns:
            bool: True if deleted successfully
        """
        with self.lock:
            if task_id in self.tasks:
                del self.tasks[task_id]
            if task_id in self.task_logs:
                del self.task_logs[task_id]
            if task_id in self.task_threads:
                del self.task_threads[task_id]
            if task_id in self.resource_monitor_threads:
                del self.resource_monitor_threads[task_id]
            if task_id in self.progress_trackers:
                del self.progress_trackers[task_id]
            if task_id in self.task_resources:
                del self.task_resources[task_id]
            if task_id in self.task_lifecycle_loggers:
                del self.task_lifecycle_loggers[task_id]
        
        logger.info(f"Deleted training task: {task_id}")
        return True
    
    def get_logs(self, task_id: str, limit: int = 100) -> List[str]:
        """
        Get task logs.
        
        Args:
            task_id: Task ID
            limit: Maximum number of log lines
            
        Returns:
            List[str]: Log lines
        """
        logs = self.task_logs.get(task_id, [])
        return logs[-limit:] if limit > 0 else logs
    
    def _log(self, task_id: str, message: str):
        """Add log message to task."""
        if task_id not in self.task_logs:
            self.task_logs[task_id] = []
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Avoid duplicate timestamps if message already has one
        if message.startswith('[') and ']' in message:
            self.task_logs[task_id].append(message)
        else:
            self.task_logs[task_id].append(f"[{timestamp}] {message}")
        
        # Limit log size to prevent memory issues
        if len(self.task_logs[task_id]) > 10000:
            self.task_logs[task_id] = self.task_logs[task_id][-5000:]
    
    def get_resource_usage(self, task_id: str) -> Optional[Dict]:
        """
        Get resource usage for a task.
        
        Args:
            task_id: Task ID
            
        Returns:
            Optional[Dict]: Resource usage information or None if not available
        """
        return self.task_resources.get(task_id)
    
    def _health_monitor_loop(self):
        """
        Background thread to monitor task health and clean up zombie tasks.
        Runs every 30 seconds.
        """
        while not self._shutdown_flag.is_set():
            try:
                time.sleep(30)  # Check every 30 seconds
                
                with self.lock:
                    tasks_to_cleanup = []
                    
                    # Check for zombie tasks (status is RUNNING but thread is dead)
                    for task_id, task in list(self.tasks.items()):
                        thread = self.task_threads.get(task_id)
                        
                        # Check if thread is dead but task is still marked as running
                        if (task.status == TaskStatus.RUNNING and 
                            thread is not None and not thread.is_alive()):
                            logger.warning(f"Detected zombie task {task_id}: thread is dead but status is RUNNING")
                            tasks_to_cleanup.append((task_id, "zombie_thread"))
                        
                        # Check for tasks stuck in INITIALIZING for too long
                        if task.status == TaskStatus.INITIALIZING:
                            init_timeout = SERVICE_CONFIG.get('task_initialization_timeout', 60)
                            if task.started_at:
                                elapsed = (datetime.now() - task.started_at).total_seconds()
                                if elapsed > init_timeout:
                                    logger.warning(f"Task {task_id} stuck in INITIALIZING for {elapsed}s")
                                    tasks_to_cleanup.append((task_id, "init_timeout"))
                            elif task.created_at:
                                elapsed = (datetime.now() - task.created_at).total_seconds()
                                if elapsed > init_timeout:
                                    logger.warning(f"Task {task_id} stuck in INITIALIZING for {elapsed}s")
                                    tasks_to_cleanup.append((task_id, "init_timeout"))
                    
                    # Clean up zombie tasks
                    for task_id, reason in tasks_to_cleanup:
                        try:
                            lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
                            if lifecycle_logger:
                                lifecycle_logger.log_failure(
                                    Exception(f"Zombie task detected: {reason}"),
                                    "health_check",
                                    {"reason": reason}
                                )
                            
                            if task_id in self.tasks:
                                self.tasks[task_id].status = TaskStatus.FAILED
                                self.tasks[task_id].error = f"Task health check failed: {reason}"
                                self.tasks[task_id].updated_at = datetime.now()
                            
                            # Clean up thread references
                            if task_id in self.task_threads:
                                self.task_threads.pop(task_id)
                            if task_id in self.resource_monitor_threads:
                                self.resource_monitor_threads.pop(task_id)
                            
                            logger.info(f"Cleaned up zombie task {task_id}: {reason}")
                        except Exception as e:
                            logger.error(f"Error cleaning up zombie task {task_id}: {e}", exc_info=True)
                    
                    # Clean up dead resource monitor threads
                    for task_id, monitor_thread in list(self.resource_monitor_threads.items()):
                        if not monitor_thread.is_alive():
                            logger.debug(f"Removing dead resource monitor thread for task {task_id}")
                            self.resource_monitor_threads.pop(task_id, None)
                    
            except Exception as e:
                logger.error(f"Error in health monitor loop: {e}", exc_info=True)
                time.sleep(30)  # Wait before retrying
    
    def shutdown(self):
        """Shutdown the service and clean up resources."""
        logger.info("Shutting down training service...")
        self._shutdown_flag.set()
        
        # Wait for health monitor thread
        if self.health_monitor_thread.is_alive():
            self.health_monitor_thread.join(timeout=5)
        
        # Clean up executor
        self.executor.shutdown(wait=True)
        
        logger.info("Training service shutdown complete")
    
    def get_health_status(self) -> Dict:
        """
        Get health status of the training service.
        
        Returns:
            Dict with health information
        """
        with self.lock:
            running_tasks = sum(1 for t in self.tasks.values() if t.status == TaskStatus.RUNNING)
            pending_tasks = sum(1 for t in self.tasks.values() if t.status == TaskStatus.PENDING)
            active_threads = len([t for t in self.task_threads.values() if t.is_alive()])
            active_monitors = len([t for t in self.resource_monitor_threads.values() if t.is_alive()])
            
            return {
                "status": "healthy",
                "total_tasks": len(self.tasks),
                "running_tasks": running_tasks,
                "pending_tasks": pending_tasks,
                "active_threads": active_threads,
                "active_monitor_threads": active_monitors,
                "max_concurrent_tasks": SERVICE_CONFIG.get('max_concurrent_tasks', 4),
                "resource_check_enabled": SERVICE_CONFIG.get('enable_resource_check', True)
            }

