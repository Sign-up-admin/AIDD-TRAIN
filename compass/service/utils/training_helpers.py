"""
Helper classes and functions for training execution.

This module provides utilities for managing training output redirection,
resource monitoring, and configuration setup.
"""

import os
import sys
import io
import time
import threading
import logging
from typing import Optional, Dict, Callable, Any

from compass.service.models.task import TaskStatus

logger = logging.getLogger(__name__)


class TeeOutput:
    """Custom output stream that writes to both console and log stream."""

    def __init__(
        self,
        stream: io.StringIO,
        task_id: str,
        original_stream: Any,
        log_func: Callable[[str, str], None],
        stream_manager: Optional[Any],
    ):
        """
        Initialize TeeOutput.

        Args:
            stream: StringIO capture stream
            task_id: Task ID
            original_stream: Original stdout/stderr stream
            log_func: Function to log lines
            stream_manager: Stream manager for WebSocket streaming
        """
        self.stream = stream
        self.task_id = task_id
        self.original_stream = original_stream
        self.log_func = log_func
        self.stream_manager = stream_manager
        self.buffer = ""
        self._isatty = getattr(original_stream, "isatty", lambda: False)()

    def write(self, text: str) -> None:
        """Write text to both original stream and log stream."""
        if not text:  # Skip empty writes
            return

        # Write to original stream immediately (for console output)
        try:
            self.original_stream.write(text)
            self.original_stream.flush()
        except (OSError, ValueError):
            # Ignore errors writing to original stream (e.g., if it's closed)
            pass

        # Push raw text to WebSocket stream immediately (preserves ANSI escape codes)
        if self.stream_manager:
            try:
                self.stream_manager.push_log(self.task_id, text)
            except Exception as e:
                logger.warning(f"Failed to push log to stream for task {self.task_id}: {e}")

        # Buffer the text for line-based logging
        self.buffer += text

        # Process complete lines for structured logging
        while "\n" in self.buffer:
            line, self.buffer = self.buffer.split("\n", 1)
            if line.strip():  # Only log non-empty lines
                try:
                    self.log_func(self.task_id, line.strip())
                except Exception as e:
                    logger.warning(f"Failed to log line for task {self.task_id}: {e}")

    def flush(self) -> None:
        """Flush the stream."""
        try:
            self.original_stream.flush()
        except (OSError, ValueError):
            pass

        # Flush any remaining buffer as a line
        if self.buffer.strip():
            try:
                self.log_func(self.task_id, self.buffer.strip())
                self.buffer = ""
            except Exception as e:
                logger.warning(f"Failed to flush log for task {self.task_id}: {e}")

    def isatty(self) -> bool:
        """Return True if this is a TTY, needed for tqdm and other libraries."""
        return self._isatty

    def fileno(self) -> int:
        """Return file descriptor if available."""
        if hasattr(self.original_stream, "fileno"):
            try:
                return self.original_stream.fileno()
            except (OSError, ValueError):
                pass
        return -1

    def __getattr__(self, name: str):
        """Delegate other attributes to original stream."""
        return getattr(self.original_stream, name)


class OutputRedirector:
    """Manages stdout/stderr redirection for training tasks."""

    def __init__(
        self,
        task_id: str,
        log_func: Callable[[str, str], None],
        stream_manager: Optional[Any],
    ):
        """
        Initialize output redirector.

        Args:
            task_id: Task ID
            log_func: Function to log lines
            stream_manager: Stream manager for WebSocket streaming
        """
        self.task_id = task_id
        self.log_func = log_func
        self.stream_manager = stream_manager
        self.old_stdout = sys.stdout
        self.old_stderr = sys.stderr
        self.stdout_capture = io.StringIO()
        self.stderr_capture = io.StringIO()
        self.tee_stdout: Optional[TeeOutput] = None
        self.tee_stderr: Optional[TeeOutput] = None

    def setup(self) -> None:
        """Set up output redirection."""
        self.tee_stdout = TeeOutput(
            self.stdout_capture, self.task_id, self.old_stdout, self.log_func, self.stream_manager
        )
        self.tee_stderr = TeeOutput(
            self.stderr_capture, self.task_id, self.old_stderr, self.log_func, self.stream_manager
        )
        sys.stdout = self.tee_stdout
        sys.stderr = self.tee_stderr

    def restore(self) -> None:
        """Restore original stdout/stderr."""
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr


class ResourceMonitor:
    """Monitors resource usage for training tasks."""

    def __init__(
        self,
        task_id: str,
        tasks: Dict[str, Any],
        task_resources: Dict[str, Dict[str, Any]],
        stream_manager: Optional[Any],
        lifecycle_logger: Optional[Any],
        lock: threading.Lock,
        get_resource_usage_func: Callable[[], Dict[str, Any]],
    ):
        """
        Initialize resource monitor.

        Args:
            task_id: Task ID
            tasks: Dictionary of tasks
            task_resources: Dictionary to store resource usage
            stream_manager: Stream manager for WebSocket streaming
            lifecycle_logger: Lifecycle logger instance
            lock: Thread lock
            get_resource_usage_func: Function to get resource usage
        """
        self.task_id = task_id
        self.tasks = tasks
        self.task_resources = task_resources
        self.stream_manager = stream_manager
        self.lifecycle_logger = lifecycle_logger
        self.lock = lock
        self.get_resource_usage = get_resource_usage_func
        self.thread: Optional[threading.Thread] = None

    def start(self) -> None:
        """Start resource monitoring thread."""
        self.thread = threading.Thread(
            target=self._monitor_loop, daemon=True, name=f"ResourceMonitor-{self.task_id}"
        )
        self.thread.start()

    def _monitor_loop(self) -> None:
        """Monitor resource usage periodically."""
        last_resource_push = 0
        resource_push_interval = 2  # Push resources every 2 seconds via WebSocket

        while self.task_id in self.tasks and self.tasks[self.task_id].status == TaskStatus.RUNNING:
            try:
                resources = self.get_resource_usage()
                with self.lock:
                    if self.task_id in self.task_resources:
                        self.task_resources[self.task_id] = resources

                # Push resources to stream queue for WebSocket clients
                current_time = time.time()
                if current_time - last_resource_push >= resource_push_interval:
                    if self.stream_manager:
                        self.stream_manager.push_resources(self.task_id, resources)
                    last_resource_push = current_time

                # Log resource usage periodically (every 30 seconds)
                if self.lifecycle_logger and int(time.time()) % 30 == 0:
                    self.lifecycle_logger.log_resource_usage(resources)
            except Exception as e:
                logger.warning(f"Failed to get resources for task {self.task_id}: {e}")
            time.sleep(1)  # Update every 1 second for more responsive monitoring


def prepare_training_config(
    config: Dict[str, Any], log_dir: str, checkpoint_dir: str
) -> Dict[str, Any]:
    """
    Prepare training configuration with task-specific settings.

    Args:
        config: Original task configuration
        log_dir: Log directory path
        checkpoint_dir: Checkpoint directory path

    Returns:
        Prepared training configuration
    """
    from compass.config import get_config

    training_config = get_config(config.get("execution_mode", "validation_tuned"))
    training_config.update(config)
    training_config["log_dir"] = log_dir
    training_config["checkpoint_dir"] = checkpoint_dir

    # Use default dataset if no dataset_id provided
    if not config.get("dataset_id"):
        project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
        default_index = os.path.join(project_root, "index", "INDEX_general_PL.2020R1.lst")
        default_dataset = os.path.join(project_root, "PDBbind-2025.8.4", "P-L")

        if os.path.exists(default_index):
            training_config["index_file"] = default_index
        else:
            training_config["index_file"] = training_config.get(
                "index_file", "index/INDEX_general_PL.2020R1.lst"
            )

        if os.path.exists(default_dataset):
            training_config["dataset_path"] = default_dataset
        else:
            training_config["dataset_path"] = training_config.get(
                "dataset_path", "PDBbind-2025.8.4/P-L/"
            )

    return training_config
