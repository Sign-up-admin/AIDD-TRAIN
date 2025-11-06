"""
Upload queue manager for concurrent upload control.
"""

import asyncio
import logging
from typing import Optional, Callable, Any
from datetime import datetime
from enum import Enum
import threading

logger = logging.getLogger(__name__)


class UploadStatus(Enum):
    """Upload status."""

    PENDING = "pending"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"


class UploadTask:
    """Upload task information."""

    def __init__(self, task_id: str, dataset_id: str):
        self.task_id = task_id
        self.dataset_id = dataset_id
        self.status = UploadStatus.PENDING
        self.created_at = datetime.now()
        self.started_at: Optional[datetime] = None
        self.completed_at: Optional[datetime] = None
        self.error: Optional[str] = None
        self.lock = threading.Lock()


class UploadQueueManager:
    """Manager for upload queue with concurrency control."""

    def __init__(self, max_concurrent: int = 2):
        """
        Initialize upload queue manager.

        Args:
            max_concurrent: Maximum number of concurrent uploads
        """
        self.max_concurrent = max_concurrent
        self.semaphore = threading.Semaphore(max_concurrent)
        self.active_uploads: dict[str, UploadTask] = {}
        self.lock = threading.Lock()
        self.logger = logging.getLogger(__name__)

    def submit_upload(
        self, task_id: str, dataset_id: str, upload_func: Callable[[str], Any]
    ) -> UploadTask:
        """
        Submit an upload task to the queue.

        Args:
            task_id: Unique task ID
            dataset_id: Dataset ID
            upload_func: Function to execute for upload (synchronous)

        Returns:
            UploadTask: Upload task object
        """
        task = UploadTask(task_id, dataset_id)
        with self.lock:
            self.active_uploads[task_id] = task

        # Start upload in background thread
        thread = threading.Thread(
            target=self._process_upload, args=(task, upload_func), daemon=True
        )
        thread.start()

        return task

    def _process_upload(self, task: UploadTask, upload_func: Callable[[str], Any]):
        """Process upload with concurrency control."""
        with self.semaphore:
            with task.lock:
                task.status = UploadStatus.PROCESSING
                task.started_at = datetime.now()

            try:
                self.logger.info(
                    f"Starting upload task {task.task_id} for dataset {task.dataset_id}"
                )
                # Execute synchronous upload function
                upload_func(task.dataset_id)

                with task.lock:
                    task.status = UploadStatus.COMPLETED
                    task.completed_at = datetime.now()
                self.logger.info(f"Completed upload task {task.task_id}")
            except Exception as e:
                with task.lock:
                    task.status = UploadStatus.FAILED
                    task.error = str(e)
                    task.completed_at = datetime.now()
                self.logger.error(f"Failed upload task {task.task_id}: {e}", exc_info=True)

    def get_task_status(self, task_id: str) -> Optional[UploadTask]:
        """Get upload task status."""
        with self.lock:
            return self.active_uploads.get(task_id)

    def get_queue_size(self) -> int:
        """Get current queue size (pending tasks)."""
        with self.lock:
            return len(
                [t for t in self.active_uploads.values() if t.status == UploadStatus.PENDING]
            )

    def get_active_count(self) -> int:
        """Get number of active uploads."""
        with self.lock:
            return len(
                [t for t in self.active_uploads.values() if t.status == UploadStatus.PROCESSING]
            )
