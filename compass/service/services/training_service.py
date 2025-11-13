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
import asyncio
import json
from typing import Dict, Optional, List, Set
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from collections import deque

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
from compass.service.utils.training_helpers import (
    OutputRedirector,
    ResourceMonitor,
    prepare_training_config,
)
from compass.service.utils.task_stop_helpers import (
    validate_task_can_be_stopped,
    set_cancellation_flag,
    wait_for_task_completion,
)

logger = logging.getLogger(__name__)

# Log warning about psutil after logger is initialized
if psutil is None:
    logger.warning("psutil not available, resource monitoring will be limited")


def get_resource_usage() -> Dict:
    """
    Get current resource usage (GPU, CPU, Memory, Storage).

    Returns:
        Dict: Resource usage information
    """
    resources = {
        "cpu_percent": 0.0,
        "memory": {"total_gb": 0.0, "used_gb": 0.0, "percent": 0.0},
        "gpu": None,
        "storage": None,
    }

    # Get CPU and memory info if psutil is available
    if psutil:
        try:
            resources["cpu_percent"] = psutil.cpu_percent(interval=0.1)
            mem = psutil.virtual_memory()
            resources["memory"] = {
                "total_gb": mem.total / (1024**3),
                "used_gb": mem.used / (1024**3),
                "percent": mem.percent,
            }

            # Get disk/storage usage
            try:
                disk = psutil.disk_usage("/")
                resources["storage"] = {
                    "total_gb": disk.total / (1024**3),
                    "used_gb": disk.used / (1024**3),
                    "free_gb": disk.free / (1024**3),
                    "percent": (disk.used / disk.total * 100) if disk.total > 0 else 0,
                }
            except Exception as e:
                logger.warning(f"Failed to get disk info: {e}")
                resources["storage"] = {"available": False, "error": str(e)}
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
                    "reserved_percent": (reserved / total * 100) if total > 0 else 0,
                },
            }
        except Exception as e:
            logger.warning(f"Failed to get GPU info: {e}")
            resources["gpu"] = {"available": False, "error": str(e)}
    else:
        resources["gpu"] = {"available": False}

    if resources["storage"] is None:
        resources["storage"] = {"available": False}

    return resources


class StreamManager:
    """Manages real-time log streams for training tasks."""

    def __init__(self, max_queue_size: int = 10000):
        """
        Initialize stream manager.

        Args:
            max_queue_size: Maximum size of log queue per task
        """
        self.max_queue_size = max_queue_size
        self.log_queues: Dict[str, asyncio.Queue] = {}
        self.resource_queues: Dict[str, asyncio.Queue] = {}
        self.lock = threading.Lock()
        # Event loop for thread-safe queue operations
        self._loop = None
        self._loop_thread = None
        self._start_loop()

    def _start_loop(self):
        """Start event loop in a separate thread for thread-safe queue operations."""

        def run_loop():
            self._loop = asyncio.new_event_loop()
            asyncio.set_event_loop(self._loop)
            self._loop.run_forever()

        self._loop_thread = threading.Thread(target=run_loop, daemon=True)
        self._loop_thread.start()
        # Wait for loop to be ready
        while self._loop is None:
            time.sleep(0.01)

    def _run_in_loop(self, coro):
        """Run coroutine in the event loop thread."""
        if self._loop is None:
            return
        future = asyncio.run_coroutine_threadsafe(coro, self._loop)
        return future

    def create_stream(self, task_id: str):
        """Create log and resource streams for a task."""
        import logging

        stream_logger = logging.getLogger(__name__)

        # 先检查是否已存在（在锁保护下）
        with self.lock:
            if task_id in self.log_queues and task_id in self.resource_queues:
                stream_logger.debug(f"Streams already exist for task {task_id}")
                return

        def _create():
            """在事件循环线程中创建队列"""
            created_log = False
            created_resource = False
            # 在锁保护下检查并创建队列
            with self.lock:
                # 再次检查（可能在等待期间已被其他线程创建）
                if task_id not in self.log_queues:
                    self.log_queues[task_id] = asyncio.Queue(maxsize=self.max_queue_size)
                    created_log = True
                if task_id not in self.resource_queues:
                    self.resource_queues[task_id] = asyncio.Queue(maxsize=100)
                    created_resource = True

            if created_log or created_resource:
                stream_logger.info(
                    f"Created streams for task {task_id}: "
                    f"log_queue={created_log}, resource_queue={created_resource}"
                )
            else:
                stream_logger.debug(f"Streams already exist for task {task_id}")

        # 如果事件循环可用，在事件循环线程中创建队列
        if self._loop:
            # 使用call_soon_threadsafe确保在事件循环线程中创建asyncio.Queue
            self._loop.call_soon_threadsafe(_create)
        else:
            # 如果事件循环不可用，直接执行（这种情况应该很少发生）
            _create()

    def push_log(self, task_id: str, data: str):
        """
        Push log data to stream queue.

        Args:
            task_id: Task ID
            data: Log data (may contain ANSI escape codes)
        """

        def _push():
            if task_id in self.log_queues:
                queue = self.log_queues[task_id]
                try:
                    # Use put_nowait to avoid blocking
                    queue.put_nowait(
                        {"type": "log", "data": data, "timestamp": datetime.now().isoformat()}
                    )
                except asyncio.QueueFull:
                    # Remove oldest item if queue is full
                    try:
                        queue.get_nowait()
                        queue.put_nowait(
                            {"type": "log", "data": data, "timestamp": datetime.now().isoformat()}
                        )
                    except asyncio.QueueEmpty:
                        pass

        with self.lock:
            if task_id in self.log_queues and self._loop:
                self._loop.call_soon_threadsafe(_push)

    def push_resources(self, task_id: str, resources: Dict):
        """
        Push resource data to stream queue.

        Args:
            task_id: Task ID
            resources: Resource usage data
        """

        def _push():
            if task_id in self.resource_queues:
                queue = self.resource_queues[task_id]
                try:
                    queue.put_nowait(
                        {
                            "type": "resources",
                            "data": resources,
                            "timestamp": datetime.now().isoformat(),
                        }
                    )
                except asyncio.QueueFull:
                    # Remove oldest item if queue is full
                    try:
                        queue.get_nowait()
                        queue.put_nowait(
                            {
                                "type": "resources",
                                "data": resources,
                                "timestamp": datetime.now().isoformat(),
                            }
                        )
                    except asyncio.QueueEmpty:
                        pass

        with self.lock:
            if task_id in self.resource_queues and self._loop:
                self._loop.call_soon_threadsafe(_push)

    def get_log_queue(self, task_id: str) -> Optional[asyncio.Queue]:
        """Get log queue for a task."""
        with self.lock:
            return self.log_queues.get(task_id)

    def get_resource_queue(self, task_id: str) -> Optional[asyncio.Queue]:
        """Get resource queue for a task."""
        with self.lock:
            return self.resource_queues.get(task_id)

    def remove_stream(self, task_id: str):
        """Remove streams for a task."""
        with self.lock:
            self.log_queues.pop(task_id, None)
            self.resource_queues.pop(task_id, None)

    def shutdown(self):
        """Shutdown stream manager and clean up resources."""
        with self.lock:
            # Clear all queues first
            self.log_queues.clear()
            self.resource_queues.clear()

            # Stop event loop gracefully
            if self._loop and not self._loop.is_closed():
                try:
                    # Schedule loop stop
                    self._loop.call_soon_threadsafe(self._loop.stop)
                    # Wait a bit for loop to stop
                    if self._loop_thread and self._loop_thread.is_alive():
                        self._loop_thread.join(timeout=2.0)
                        if self._loop_thread.is_alive():
                            logger.warning(
                                "StreamManager event loop thread did not stop within timeout"
                            )
                            # Force close the loop if thread is still alive
                            try:
                                if self._loop and not self._loop.is_closed():
                                    self._loop.close()
                            except Exception as close_error:
                                logger.warning(f"Error force-closing event loop: {close_error}")
                except (RuntimeError, AttributeError) as e:
                    # Loop may already be closed or in invalid state
                    logger.debug(f"Event loop already closed or invalid: {e}")
                except Exception as e:
                    logger.warning(f"Error stopping StreamManager event loop: {e}", exc_info=True)
                finally:
                    # Always clear references
                    self._loop = None
                    self._loop_thread = None

    def __getstate__(self):
        """
        Custom pickle state getter.
        Excludes threading.Lock and asyncio objects which cannot be pickled.
        """
        state = self.__dict__.copy()
        # Remove objects that cannot be pickled
        state.pop("lock", None)
        state.pop("_loop", None)
        state.pop("_loop_thread", None)
        # Note: log_queues and resource_queues contain asyncio.Queue objects
        # which cannot be pickled, so we clear them
        state["log_queues"] = {}
        state["resource_queues"] = {}
        return state

    def __setstate__(self, state):
        """
        Custom pickle state setter.
        Recreates threading.Lock and asyncio objects after unpickling.
        """
        self.__dict__.update(state)
        # Recreate lock object after unpickling
        self.lock = threading.Lock()
        # Recreate event loop and thread
        self._loop = None
        self._loop_thread = None
        self._start_loop()

    def cleanup(self):
        """Cleanup resources: stop event loop and close all queues."""
        if self._loop and not self._loop.is_closed():
            try:
                # Stop the event loop
                self._loop.call_soon_threadsafe(self._loop.stop)
                # Wait for loop thread to finish (with timeout)
                if self._loop_thread and self._loop_thread.is_alive():
                    self._loop_thread.join(timeout=2.0)
                # Close the loop
                if not self._loop.is_closed():
                    self._loop.close()
            except Exception as e:
                logger.warning(f"Error cleaning up StreamManager event loop: {e}")

        # Clear queues
        with self.lock:
            self.log_queues.clear()
            self.resource_queues.clear()


class TrainingService:
    """Service for managing training tasks."""

    def __init__(self):
        """Initialize training service."""
        self.tasks: Dict[str, TrainingTaskResponse] = {}
        self.task_logs: Dict[str, List[str]] = {}
        self.executor = ThreadPoolExecutor(max_workers=SERVICE_CONFIG["max_workers"])
        self.task_threads: Dict[str, threading.Thread] = {}
        self.progress_trackers: Dict[str, ProgressTracker] = {}
        self.task_resources: Dict[str, Dict] = {}  # Store resource usage for each task
        self.resource_monitor_threads: Dict[str, threading.Thread] = (
            {}
        )  # Track resource monitor threads
        self.task_lifecycle_loggers: Dict[str, TaskLifecycleLogger] = {}  # Track lifecycle loggers
        self.lock = threading.Lock()

        # Stream manager for real-time log streaming
        self.stream_manager = StreamManager()

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
        result = {"available": True, "issues": [], "resources": {}}

        enable_check = SERVICE_CONFIG.get("enable_resource_check", True)
        test_mode = SERVICE_CONFIG.get("test_mode", False)
        resource_check_relaxed = SERVICE_CONFIG.get("resource_check_relaxed", False)

        # In test mode, disable resource check
        if test_mode:
            logger.info("Test mode enabled - resource check disabled")
            return result

        if not enable_check:
            return result

        # Get current resource usage with improved CPU sampling
        # Use multiple samples for more accurate CPU usage
        resources = get_resource_usage()
        if psutil and resource_check_relaxed:
            # In relaxed mode, use average of multiple samples for CPU
            cpu_samples = []
            for _ in range(3):
                sample = psutil.cpu_percent(interval=0.1)
                cpu_samples.append(sample)
            if cpu_samples:
                resources["cpu_percent"] = sum(cpu_samples) / len(cpu_samples)
                logger.debug(f"CPU usage (averaged): {resources['cpu_percent']:.1f}%")

        result["resources"] = resources

        # Check concurrent tasks limit
        max_concurrent = SERVICE_CONFIG.get("max_concurrent_tasks", 4)
        running_count = sum(1 for t in self.tasks.values() if t.status == TaskStatus.RUNNING)
        if running_count >= max_concurrent:
            result["available"] = False
            result["issues"].append(f"Maximum concurrent tasks ({max_concurrent}) reached")

        # Check CPU usage
        if psutil:
            max_cpu = SERVICE_CONFIG.get("max_cpu_percent", 90.0)
            # In relaxed mode, allow higher CPU usage (up to 98%)
            if resource_check_relaxed:
                max_cpu = min(max_cpu * 1.1, 98.0)  # Allow 10% more or up to 98%
                logger.debug(f"Relaxed mode: CPU limit adjusted to {max_cpu:.1f}%")

            cpu_percent = resources.get("cpu_percent", 0)
            if cpu_percent > max_cpu:
                result["available"] = False
                result["issues"].append(
                    f"CPU usage ({cpu_percent:.1f}%) exceeds limit ({max_cpu:.1f}%)"
                )

        # Check memory availability
        if psutil:
            min_memory_gb = SERVICE_CONFIG.get("min_available_memory_gb", 2.0)
            # In relaxed mode, reduce memory requirement by 50%
            if resource_check_relaxed:
                min_memory_gb = min_memory_gb * 0.5
                logger.debug(f"Relaxed mode: Memory requirement reduced to {min_memory_gb:.2f} GB")

            mem_info = resources.get("memory", {})
            available_gb = mem_info.get("total_gb", 0) - mem_info.get("used_gb", 0)
            if available_gb < min_memory_gb:
                result["available"] = False
                result["issues"].append(
                    f"Available memory ({available_gb:.2f} GB) below minimum ({min_memory_gb:.2f} GB)"
                )

        # Check GPU availability if needed (could be enhanced based on config)
        if resources.get("gpu") and resources["gpu"].get("available"):
            gpu_mem = resources["gpu"].get("memory", {})
            allocated_percent = gpu_mem.get("allocated_percent", 0)
            # In relaxed mode, allow higher GPU usage (up to 99%)
            gpu_threshold = 99.0 if resource_check_relaxed else 95.0
            if allocated_percent > gpu_threshold:
                result["available"] = False
                result["issues"].append(
                    f"GPU memory usage ({allocated_percent:.1f}%) too high (limit: {gpu_threshold}%)"
                )

        return result

    def create_task(
        self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None
    ) -> str:
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
        creation_timeout = SERVICE_CONFIG.get("task_creation_timeout", 30)

        # Initialize lifecycle logger
        lifecycle_logger = TaskLifecycleLogger(task_id, SERVICE_CONFIG.get("log_dir", "logs"))
        self.task_lifecycle_loggers[task_id] = lifecycle_logger

        # Log creation start
        resource_check = self._check_resources()
        lifecycle_logger.log_create_start(config, resource_check)

        # Check if creation is taking too long
        if (datetime.now() - creation_start).total_seconds() > creation_timeout:
            lifecycle_logger.log_create_failure(
                TimeoutError(f"Task creation timeout after {creation_timeout}s"), "Creation timeout"
            )
            raise ServiceException(
                f"Task creation timeout after {creation_timeout} seconds",
                error_code=ErrorCode.TIMEOUT,
            )

        # Perform resource check
        if not resource_check.get("available", True):
            issues = resource_check.get("issues", [])
            error_msg = "Resource check failed: " + "; ".join(issues)
            lifecycle_logger.log_create_failure(ValueError(error_msg), "Resource check failed")
            raise ServiceException(
                error_msg,
                error_code=ErrorCode.SERVICE_UNAVAILABLE,
                detail={"resource_check": resource_check},
            )

        # Create task with CREATING status
        task = TrainingTaskResponse(
            task_id=task_id,
            status=TaskStatus.CREATING,
            config=config,
            created_at=datetime.now(),
            updated_at=datetime.now(),
            description=description,
        )

        try:
            with self.lock:
                self.tasks[task_id] = task
                self.task_logs[task_id] = []
                self.progress_trackers[task_id] = ProgressTracker(task_id)
                self.task_resources[task_id] = {}
                # Create stream for this task
                logger.info(f"Creating stream for task {task_id}")
                self.stream_manager.create_stream(task_id)
                logger.debug(f"Stream created for task {task_id}")

            # Update status to PENDING after successful creation
            with self.lock:
                task.status = TaskStatus.PENDING
                task.updated_at = datetime.now()

            lifecycle_logger.log_state_transition(
                TaskStatus.CREATING.value, TaskStatus.PENDING.value, "Task created"
            )
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
                # 清理stream
                self.stream_manager.remove_stream(task_id)

            lifecycle_logger.log_create_failure(e, "Exception during task creation")
            raise ServiceException(
                f"Failed to create task: {str(e)}", error_code=ErrorCode.INTERNAL_ERROR
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
        init_timeout = SERVICE_CONFIG.get("task_initialization_timeout", 60)

        with self.lock:
            if task_id not in self.tasks:
                if lifecycle_logger:
                    lifecycle_logger.log_initialize_failure(
                        ValueError("Task not found"), "Task not found"
                    )
                return False

            task = self.tasks[task_id]
            if task.status != TaskStatus.PENDING and task.status != TaskStatus.PAUSED:
                if lifecycle_logger:
                    lifecycle_logger.log_initialize_failure(
                        ValueError(f"Task in invalid state: {task.status}"),
                        f"Invalid state: {task.status}",
                    )
                return False

            # Update to INITIALIZING status
            task.status = TaskStatus.INITIALIZING
            task.updated_at = datetime.now()
            if lifecycle_logger:
                lifecycle_logger.log_state_transition(
                    TaskStatus.PENDING.value,
                    TaskStatus.INITIALIZING.value,
                    "Starting task initialization",
                )
                lifecycle_logger.log_initialize_start()

            # Ensure WebSocket stream is created and send initialization message
            logger.info(f"Ensuring stream exists for task {task_id} before starting")
            self.stream_manager.create_stream(task_id)
            logger.debug(f"Stream verified/created for task {task_id}")

            # Update progress tracker with initializing status
            if task_id in self.progress_trackers:
                self.progress_trackers[task_id].set_stage(
                    "initializing", "Task initialization started"
                )

            # Send initialization message to WebSocket stream
            init_message = f"\n{'='*70}\n"
            init_message += "Training Task Initialization Started\n"
            init_message += f"Task ID: {task_id}\n"
            init_message += f"Status: {task.status.value}\n"
            init_message += f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            init_message += f"{'='*70}\n"
            self.stream_manager.push_log(task_id, init_message)

            # Send task configuration
            task_config = task.config
            config_message = "\n[Configuration]\n"
            config_message += f"Execution Mode: {task_config.get('execution_mode', 'N/A')}\n"
            config_message += f"Epochs: {task_config.get('epochs', 'N/A')}\n"
            config_message += f"Batch Size: {task_config.get('batch_size', 'N/A')}\n"
            config_message += f"Learning Rate: {task_config.get('learning_rate', 'N/A')}\n"
            config_message += f"Optimizer: {task_config.get('optimizer', 'N/A')}\n"
            if task.description:
                config_message += f"Description: {task.description}\n"
            config_message += "\n"
            self.stream_manager.push_log(task_id, config_message)

        # Start training in background thread
        init_start_time = datetime.now()

        def run_training():
            # Set up signal handlers in this thread context
            import signal

            def signal_handler(sig, frame):
                """Handle stop signals (SIGINT, SIGTERM) for graceful shutdown."""
                logger.warning(f"Received signal {sig} in training thread for task {task_id}")
                # Set cancellation flag
                if task_id in self.progress_trackers:
                    self.progress_trackers[task_id].cancel()
                    logger.info(f"Set cancellation flag for task {task_id} via signal")
                # Send message to stream
                signal_message = f"\n[Signal] Received stop signal {sig}\n"
                signal_message += "Training will stop gracefully...\n\n"
                self.stream_manager.push_log(task_id, signal_message)

            # Register signal handlers for this thread (Python signal handlers work in main thread,
            # but we can use them if called from the main process)
            # Note: On Windows, signal handling is limited
            if sys.platform != "win32":
                signal.signal(signal.SIGTERM, signal_handler)
                signal.signal(signal.SIGINT, signal_handler)

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
                            lifecycle_logger.log_state_transition(
                                TaskStatus.INITIALIZING.value,
                                TaskStatus.RUNNING.value,
                                "Initialization complete",
                            )
                            lifecycle_logger.log_initialize_success()
                            lifecycle_logger.log_start()

                        # Send status update to WebSocket stream
                        status_message = "\n[Status Update] Task is now RUNNING\n"
                        status_message += f"Started at: {self.tasks[task_id].started_at.strftime('%Y-%m-%d %H:%M:%S')}\n"
                        status_message += "Preparing training environment...\n\n"
                        self.stream_manager.push_log(task_id, status_message)

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
                        self.task_threads.pop(task_id)
                        logger.debug(f"Cleaned up thread for task {task_id}")
                    # Clean up resource monitor thread
                    if task_id in self.resource_monitor_threads:
                        self.resource_monitor_threads.pop(task_id)
                        logger.debug(f"Cleaned up resource monitor thread for task {task_id}")
                    # Note: Keep progress tracker for querying completed tasks
                    # Only remove if task is deleted

        try:
            thread = threading.Thread(
                target=run_training, daemon=True, name=f"TrainingThread-{task_id}"
            )
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

        # Ensure WebSocket stream is created for this task
        self.stream_manager.create_stream(task_id)

        # Get progress tracker
        progress_tracker = self.progress_trackers.get(task_id)
        if not progress_tracker:
            progress_tracker = ProgressTracker(task_id)
            with self.lock:
                self.progress_trackers[task_id] = progress_tracker

        progress_tracker.set_stage("initializing", "Initializing training environment")

        # Send welcome message to WebSocket stream BEFORE setting up logger
        welcome_message = f"\n{'='*70}\n"
        welcome_message += "Training Environment Setup Started\n"
        welcome_message += "Preparing log directories and logger...\n"
        welcome_message += f"{'='*70}\n\n"
        self.stream_manager.push_log(task_id, welcome_message)

        # Create log directory
        log_dir = os.path.join(SERVICE_CONFIG["log_dir"], f"task_{task_id}")
        os.makedirs(log_dir, exist_ok=True)

        # Send log directory info
        self.stream_manager.push_log(task_id, f"[Setup] Log directory: {log_dir}\n")

        # Create progress-aware logger
        logger_instance = ProgressAwareLogger(log_dir=log_dir, progress_tracker=progress_tracker)

        # Send logger setup confirmation
        self.stream_manager.push_log(task_id, "[Setup] Logger initialized\n")
        self.stream_manager.push_log(task_id, "[Setup] Starting training main function...\n\n")

        # Set up output redirection
        output_redirector = OutputRedirector(task_id, self._log, self.stream_manager)
        self.stream_manager.push_log(task_id, "[Setup] Setting up stdout/stderr redirection...\n")
        output_redirector.setup()

        # Send confirmation that redirection is set up
        self.stream_manager.push_log(
            task_id, "[Setup] Output redirection active. All training logs will appear below.\n"
        )
        self.stream_manager.push_log(task_id, f"{'='*70}\n")
        self.stream_manager.push_log(task_id, "COMPASS Training Main Function Starting\n")
        self.stream_manager.push_log(task_id, f"{'='*70}\n\n")

        # Start resource monitoring
        lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
        resource_monitor = ResourceMonitor(
            task_id,
            self.tasks,
            self.task_resources,
            self.stream_manager,
            lifecycle_logger,
            self.lock,
            get_resource_usage,
        )
        resource_monitor.start()
        with self.lock:
            self.resource_monitor_threads[task_id] = resource_monitor.thread

        # Prepare training configuration
        checkpoint_dir = os.path.join(SERVICE_CONFIG["checkpoint_dir"], f"task_{task_id}")
        training_config = prepare_training_config(config, log_dir, checkpoint_dir)

        if not config.get("dataset_id"):
            logger_instance.log(
                f"Using default PDBbind dataset: {training_config.get('dataset_path')}"
            )
            logger_instance.log(f"Using index file: {training_config.get('index_file')}")

        self._log(task_id, f"Starting training with config: {config}")

        # Send pre-training message
        self.stream_manager.push_log(
            task_id, "\n[Training] About to start main training function...\n"
        )
        self.stream_manager.push_log(task_id, "[Training] Training configuration:\n")
        self.stream_manager.push_log(
            task_id, f"  - Execution Mode: {training_config.get('execution_mode', 'N/A')}\n"
        )
        self.stream_manager.push_log(
            task_id, f"  - Epochs: {training_config.get('epochs', 'N/A')}\n"
        )
        self.stream_manager.push_log(
            task_id, f"  - Batch Size: {training_config.get('batch_size', 'N/A')}\n"
        )
        self.stream_manager.push_log(
            task_id, f"  - Learning Rate: {training_config.get('learning_rate', 'N/A')}\n"
        )
        self.stream_manager.push_log(
            task_id, f"  - Optimizer: {training_config.get('optimizer', 'N/A')}\n"
        )
        self.stream_manager.push_log(task_id, f"\n{'='*70}\n\n")

        # Run training
        try:
            # Update progress before starting
            progress_tracker.set_stage("data_processing", "Starting data processing")
            task.progress = progress_tracker.get_progress()

            # Run main training
            lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
            if lifecycle_logger:
                lifecycle_logger.log_running("training", progress_tracker.get_progress())

            # Send message right before calling main
            self.stream_manager.push_log(task_id, "[Training] Calling main() function now...\n\n")

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
            # Check if it's a cancellation exception or if cancellation was requested
            from compass.training.exceptions import TrainingCancelled

            # Check if cancellation was requested (even if exception is not TrainingCancelled)
            is_cancelled = isinstance(e, TrainingCancelled) or progress_tracker.is_cancelled()

            if is_cancelled:
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.CANCELLED
                        self.tasks[task_id].updated_at = datetime.now()
                        if isinstance(e, TrainingCancelled):
                            cancel_reason = f"Training cancelled: {str(e)}"
                        else:
                            # Cancellation was requested but got a different exception (e.g., file cleanup error)
                            cancel_reason = (
                                f"Training cancelled (encountered error during cleanup: {str(e)})"
                            )
                        progress_tracker.set_stage("cancelled", cancel_reason)
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        self._log(task_id, cancel_reason)
                        # Send cancellation message to WebSocket stream
                        cancel_message = "\n[CANCELLATION] Training cancelled successfully\n"
                        cancel_message += f"Reason: {cancel_reason}\n"
                        cancel_message += (
                            f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
                        )
                        self.stream_manager.push_log(task_id, cancel_message)
                logger.info(f"[STOP] Training task {task_id} was cancelled: {cancel_reason}")
            else:
                with self.lock:
                    if task_id in self.tasks:
                        progress_tracker.set_stage("failed", f"Training failed: {str(e)}")
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
                        if lifecycle_logger:
                            lifecycle_logger.log_failure(
                                e,
                                "training",
                                {"error": str(e), "traceback": traceback.format_exc()},
                            )
                raise
        finally:
            # Restore stdout/stderr
            output_redirector.restore()

    def stop_task(self, task_id: str) -> tuple[bool, str | None]:
        """
        Stop a training task.

        Args:
            task_id: Task ID

        Returns:
            tuple[bool, str | None]: (success, error_message)
                - success: True if stopped successfully
                - error_message: None if successful, error message if failed
        """
        stop_request_time = datetime.now()
        logger.info(
            f"[STOP] Stop request received for task {task_id} at {stop_request_time.isoformat()}"
        )

        with self.lock:
            if task_id not in self.tasks:
                logger.warning(f"[STOP] Task {task_id} not found in tasks dictionary")
                return False, "Task not found"

            task = self.tasks[task_id]
            logger.info(f"[STOP] Task {task_id} current status: {task.status.value}")

            # Validate task can be stopped
            can_stop, error_msg = validate_task_can_be_stopped(task.status)
            if not can_stop:
                logger.warning(f"[STOP] Cannot stop task {task_id}: {error_msg}")
                return False, error_msg

            # Immediately update status to CANCELLING for better user feedback
            task.status = TaskStatus.CANCELLING
            task.updated_at = datetime.now()
            logger.info(
                f"[STOP] Task {task_id} status updated to CANCELLING at {task.updated_at.isoformat()}"
            )

            # Set cancellation flag in progress tracker
            progress_tracker = self.progress_trackers.get(task_id)
            if not progress_tracker:
                logger.warning(f"[STOP] No progress tracker found for task {task_id}")
                # Try to create one if task exists
                if task_id in self.tasks:
                    logger.info(f"[STOP] Creating progress tracker for task {task_id}")
                    progress_tracker = ProgressTracker(task_id)
                    self.progress_trackers[task_id] = progress_tracker

            cancellation_flag_set_time = set_cancellation_flag(
                progress_tracker, task_id, stop_request_time
            )

            # Log thread status
            thread = self.task_threads.get(task_id)
            if thread:
                logger.info(
                    f"[STOP] Training thread for task {task_id} is {'alive' if thread.is_alive() else 'not alive'}"
                )
            else:
                logger.warning(f"[STOP] No training thread found for task {task_id}")

            # Send stop message to WebSocket stream
            stop_message = f"\n{'='*70}\n"
            stop_message += "Training Task Stop Requested\n"
            stop_message += f"Task ID: {task_id}\n"
            stop_message += f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"
            stop_message += f"Current Status: {task.status.value}\n"
            if task_id in self.progress_trackers:
                progress = self.progress_trackers[task_id].get_progress()
                stop_message += f"Current Stage: {progress.get('stage', 'unknown')}\n"
                stop_message += "Waiting for training to respond to cancellation...\n"
            stop_message += f"{'='*70}\n\n"
            self.stream_manager.push_log(task_id, stop_message)

        self._log(task_id, "Training task cancellation requested")
        logger.info(
            f"[STOP] Stop requested for training task: {task_id} (waiting for graceful shutdown)"
        )
        logger.info(f"[STOP] Current task status: {task.status.value}")

        # Quick check: if task is already in final state, return immediately
        with self.lock:
            if task_id in self.tasks:
                task = self.tasks[task_id]
                if task.status in [TaskStatus.CANCELLED, TaskStatus.FAILED, TaskStatus.COMPLETED]:
                    logger.info(
                        f"[STOP] Task {task_id} already in final state: {task.status.value}"
                    )
                    return True, None

        # Do a quick wait (2 seconds) to see if training responds quickly
        # This improves responsiveness for tasks that stop quickly
        quick_wait_time = 2.0  # seconds
        quick_wait_interval = 0.2  # seconds
        quick_waited = 0.0

        while quick_waited < quick_wait_time:
            time.sleep(quick_wait_interval)
            quick_waited += quick_wait_interval

            with self.lock:
                if task_id not in self.tasks:
                    logger.warning(f"[STOP] Task {task_id} no longer exists")
                    return True, None
                task = self.tasks[task_id]
                # Check if task has responded to cancellation
                if task.status in [TaskStatus.CANCELLED, TaskStatus.FAILED, TaskStatus.COMPLETED]:
                    logger.info(
                        f"[STOP] Task {task_id} responded quickly to stop request (status: {task.status.value}, waited: {quick_waited:.1f}s)"
                    )
                    return True, None

                # Check if progress tracker shows cancelled and thread is done
                if task_id in self.progress_trackers:
                    progress = self.progress_trackers[task_id].get_progress()
                    if progress.get("cancelled", False):
                        thread = self.task_threads.get(task_id)
                        if thread is None or not thread.is_alive():
                            if task.status == TaskStatus.CANCELLING:
                                task.status = TaskStatus.CANCELLED
                                task.updated_at = datetime.now()
                                self._log(task_id, "Training task cancelled successfully")
                                response_time = (datetime.now() - stop_request_time).total_seconds()
                                logger.info(
                                    f"[STOP] Task {task_id} cancelled quickly (thread finished, waited: {quick_waited:.1f}s, "
                                    f"total response time: {response_time:.1f}s)"
                                )
                                return True, None
                        else:
                            # Thread still alive but cancellation flag is set
                            # Status should already be CANCELLING, continue to background monitoring
                            logger.info(
                                f"[STOP] Task {task_id} cancellation acknowledged quickly (waited: {quick_waited:.1f}s, "
                                f"stage: {progress.get('stage')}, thread still processing, status: {task.status.value})"
                            )
                            # Continue to background monitoring for thread cleanup

        # If task hasn't responded quickly, start background thread to monitor
        # and return immediately so API doesn't block
        def _wait_for_cancellation_background():
            """Background thread to wait for training to respond to cancellation"""
            max_wait_time = 30  # seconds - timeout for graceful cancellation
            wait_interval = 0.2  # seconds - reduced interval for faster response
            waited = quick_waited  # Start from where quick wait left off
            last_log_time = 0

            logger.info(
                f"[STOP] Background monitoring started for task {task_id}, max wait: {max_wait_time}s"
            )

            while waited < max_wait_time:
                time.sleep(wait_interval)
                waited += wait_interval

                with self.lock:
                    if task_id not in self.tasks:
                        logger.warning(f"[STOP] Task {task_id} no longer exists")
                        return
                    task = self.tasks[task_id]
                    # Check if task has responded to cancellation
                    if task.status in [
                        TaskStatus.CANCELLED,
                        TaskStatus.FAILED,
                        TaskStatus.COMPLETED,
                    ]:
                        response_time = (datetime.now() - stop_request_time).total_seconds()
                        logger.info(
                            f"[STOP] Task {task_id} responded to stop request (status: {task.status.value}, "
                            f"waited: {waited:.1f}s, total response time: {response_time:.1f}s)"
                        )
                        if cancellation_flag_set_time:
                            flag_to_response_time = (
                                datetime.now() - cancellation_flag_set_time
                            ).total_seconds()
                            logger.info(
                                f"[STOP] Time from flag set to response: {flag_to_response_time:.3f}s"
                            )
                        return

                    # Check if progress tracker shows cancelled
                    if task_id in self.progress_trackers:
                        progress = self.progress_trackers[task_id].get_progress()
                        if progress.get("cancelled", False):
                            # Check if training thread is still alive
                            thread = self.task_threads.get(task_id)
                            if thread is None or not thread.is_alive():
                                # Thread has finished, update status
                                if task.status == TaskStatus.CANCELLING:
                                    task.status = TaskStatus.CANCELLED
                                    task.updated_at = datetime.now()
                                    self._log(task_id, "Training task cancelled successfully")
                                    response_time = (
                                        datetime.now() - stop_request_time
                                    ).total_seconds()
                                    logger.info(
                                        f"[STOP] Task {task_id} cancelled (thread finished, waited: {waited:.1f}s, "
                                        f"total response time: {response_time:.1f}s)"
                                    )
                                    return
                            else:
                                # Thread still alive but cancellation flag is set
                                # Status should already be CANCELLING, just log progress
                                if waited - last_log_time >= 2.0:  # Log every 2 seconds
                                    logger.info(
                                        f"[STOP] Task {task_id} cancellation in progress (waited: {waited:.1f}s, "
                                        f"stage: {progress.get('stage')}, thread still processing)"
                                    )
                                    last_log_time = waited

            # Timeout reached - task didn't respond gracefully
            timeout_time = datetime.now()
            total_wait_time = (timeout_time - stop_request_time).total_seconds()
            logger.warning(
                f"[STOP] Timeout reached for task {task_id} after {max_wait_time}s "
                f"(total wait: {total_wait_time:.1f}s)"
            )

            with self.lock:
                if task_id in self.tasks:
                    task = self.tasks[task_id]
                    # Force status to CANCELLED if still cancelling
                    if task.status == TaskStatus.CANCELLING:
                        task.status = TaskStatus.CANCELLED
                        task.updated_at = datetime.now()
                        self._log(
                            task_id,
                            f"Training task force-cancelled after {max_wait_time}s timeout",
                        )
                        logger.warning(
                            f"[STOP] Task {task_id} force-cancelled after {max_wait_time}s timeout"
                        )
                        logger.warning(f"[STOP] Final task status: {task.status.value}")
                        if task_id in self.progress_trackers:
                            progress = self.progress_trackers[task_id].get_progress()
                            logger.warning(
                                f"[STOP] Final progress: stage={progress.get('stage')}, "
                                f"cancelled={progress.get('cancelled')}"
                            )

                        # Check thread status
                        thread = self.task_threads.get(task_id)
                        if thread and thread.is_alive():
                            logger.warning(
                                f"[STOP] Task {task_id} thread still alive after timeout - "
                                f"may need manual cleanup. Thread name: {thread.name}"
                            )
                            # Note: Python doesn't support forcefully killing threads
                            # The thread should eventually detect cancellation and exit
                            logger.info(
                                f"Task {task_id}: Cancellation flag set, waiting for thread to respond"
                            )

                    # Send timeout message to stream
                    timeout_message = f"\n[Warning] Stop request timeout after {max_wait_time}s\n"
                    timeout_message += "Task status forced to CANCELLED\n\n"
                    self.stream_manager.push_log(task_id, timeout_message)

        # Start background thread to monitor cancellation
        cancellation_thread = threading.Thread(
            target=_wait_for_cancellation_background, daemon=True, name=f"stop-task-{task_id}"
        )
        cancellation_thread.start()
        logger.info(f"[STOP] Started background thread to monitor cancellation for task {task_id}")
        logger.info(
            f"[STOP] Cancellation request accepted for task {task_id}, monitoring in background"
        )

        # Return immediately - cancellation is in progress
        # Frontend should poll task status to check if it's cancelled
        return True, None

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
            # Remove stream
            self.stream_manager.remove_stream(task_id)

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
        if message.startswith("[") and "]" in message:
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
                        if (
                            task.status == TaskStatus.RUNNING
                            and thread is not None
                            and not thread.is_alive()
                        ):
                            logger.warning(
                                f"Detected zombie task {task_id}: thread is dead but status is RUNNING"
                            )
                            tasks_to_cleanup.append((task_id, "zombie_thread"))

                        # Check for tasks stuck in INITIALIZING for too long
                        if task.status == TaskStatus.INITIALIZING:
                            init_timeout = SERVICE_CONFIG.get("task_initialization_timeout", 60)
                            if task.started_at:
                                elapsed = (datetime.now() - task.started_at).total_seconds()
                                if elapsed > init_timeout:
                                    logger.warning(
                                        f"Task {task_id} stuck in INITIALIZING for {elapsed}s"
                                    )
                                    tasks_to_cleanup.append((task_id, "init_timeout"))
                            elif task.created_at:
                                elapsed = (datetime.now() - task.created_at).total_seconds()
                                if elapsed > init_timeout:
                                    logger.warning(
                                        f"Task {task_id} stuck in INITIALIZING for {elapsed}s"
                                    )
                                    tasks_to_cleanup.append((task_id, "init_timeout"))

                    # Clean up zombie tasks
                    for task_id, reason in tasks_to_cleanup:
                        try:
                            lifecycle_logger = self.task_lifecycle_loggers.get(task_id)
                            if lifecycle_logger:
                                lifecycle_logger.log_failure(
                                    Exception(f"Zombie task detected: {reason}"),
                                    "health_check",
                                    {"reason": reason},
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
                            logger.error(
                                f"Error cleaning up zombie task {task_id}: {e}", exc_info=True
                            )

                    # Clean up dead resource monitor threads
                    for task_id, monitor_thread in list(self.resource_monitor_threads.items()):
                        if not monitor_thread.is_alive():
                            logger.debug(
                                f"Removing dead resource monitor thread for task {task_id}"
                            )
                            self.resource_monitor_threads.pop(task_id, None)

            except Exception as e:
                logger.error(f"Error in health monitor loop: {e}", exc_info=True)
                time.sleep(30)  # Wait before retrying

    def shutdown(self):
        """Shutdown the service and clean up resources."""
        logger.info("Shutting down training service...")
        self._shutdown_flag.set()

        # Stop all running tasks gracefully
        with self.lock:
            running_task_ids = [
                task_id for task_id, task in self.tasks.items() if task.status == TaskStatus.RUNNING
            ]

        for task_id in running_task_ids:
            try:
                logger.info(f"Stopping task {task_id} during shutdown...")
                self.stop_task(task_id)
            except Exception as e:
                logger.warning(f"Error stopping task {task_id} during shutdown: {e}", exc_info=True)

        # Wait for all task threads to complete (with timeout)
        with self.lock:
            task_threads = list(self.task_threads.values())
            resource_threads = list(self.resource_monitor_threads.values())

        for thread in task_threads:
            if thread.is_alive():
                thread.join(timeout=3.0)
                if thread.is_alive():
                    logger.warning(f"Task thread {thread.name} did not stop within timeout")

        for thread in resource_threads:
            if thread.is_alive():
                thread.join(timeout=2.0)
                if thread.is_alive():
                    logger.warning(
                        f"Resource monitor thread {thread.name} did not stop within timeout"
                    )

        # Shutdown stream manager
        if hasattr(self, "stream_manager") and self.stream_manager:
            try:
                self.stream_manager.shutdown()
            except (RuntimeError, AttributeError) as e:
                logger.debug(f"Stream manager already closed or invalid: {e}")
            except Exception as e:
                logger.warning(f"Error shutting down stream manager: {e}", exc_info=True)

        # Wait for health monitor thread
        if hasattr(self, "health_monitor_thread") and self.health_monitor_thread.is_alive():
            self.health_monitor_thread.join(timeout=5.0)
            if self.health_monitor_thread.is_alive():
                logger.warning("Health monitor thread did not stop within timeout")

        # Clean up executor
        try:
            # ThreadPoolExecutor.shutdown(wait=True) blocks until all tasks complete
            # For timeout control, we use wait=False and manually wait with timeout
            self.executor.shutdown(wait=False)
            # Give threads time to finish (max 10 seconds)
            import time

            start_time = time.time()
            while time.time() - start_time < 10.0:
                if all(
                    not thread.is_alive()
                    for thread in self.executor._threads
                    if hasattr(self.executor, "_threads")
                ):
                    break
                time.sleep(0.1)
        except (RuntimeError, AttributeError) as e:
            logger.debug(f"Executor already shut down or invalid: {e}")
        except Exception as e:
            logger.warning(f"Unexpected error shutting down executor: {e}", exc_info=True)

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
            active_monitors = len(
                [t for t in self.resource_monitor_threads.values() if t.is_alive()]
            )

            return {
                "status": "healthy",
                "total_tasks": len(self.tasks),
                "running_tasks": running_tasks,
                "pending_tasks": pending_tasks,
                "active_threads": active_threads,
                "active_monitor_threads": active_monitors,
                "max_concurrent_tasks": SERVICE_CONFIG.get("max_concurrent_tasks", 4),
                "resource_check_enabled": SERVICE_CONFIG.get("enable_resource_check", True),
            }
