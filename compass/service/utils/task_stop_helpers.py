"""
Helper functions for task stopping and resource cleanup.

This module provides utilities for managing task cancellation and resource cleanup.
"""

import logging
import time
from typing import Tuple, Optional, Dict, Any
from datetime import datetime

from compass.service.models.task import TaskStatus

logger = logging.getLogger(__name__)


def validate_task_can_be_stopped(task_status: TaskStatus) -> Tuple[bool, Optional[str]]:
    """
    Validate if a task can be stopped.

    Args:
        task_status: Current task status

    Returns:
        Tuple of (can_stop, error_message)
    """
    if task_status not in [TaskStatus.RUNNING, TaskStatus.INITIALIZING]:
        error_msg = (
            f"Task is in '{task_status.value}' status. "
            "Only running or initializing tasks can be stopped."
        )
        return False, error_msg
    return True, None


def set_cancellation_flag(
    progress_tracker, task_id: str, stop_request_time: datetime
) -> Optional[datetime]:
    """
    Set cancellation flag in progress tracker.

    Args:
        progress_tracker: Progress tracker instance
        task_id: Task ID
        stop_request_time: Time when stop was requested

    Returns:
        Cancellation flag set time or None
    """
    if progress_tracker:
        cancellation_flag_set_time = datetime.now()
        progress_tracker.cancel()

        logger.info(
            f"[STOP] Cancellation flag set for task {task_id} at {cancellation_flag_set_time.isoformat()}"
        )
        logger.info(
            f"[STOP] Time from request to flag set: "
            f"{(cancellation_flag_set_time - stop_request_time).total_seconds():.3f}s"
        )

        # Verify cancellation flag was set
        progress_after = progress_tracker.get_progress()
        if not progress_after.get("cancelled", False):
            logger.error(f"[STOP] ERROR: Cancellation flag not set correctly for task {task_id}!")
        else:
            logger.info(
                f"[STOP] Verified cancellation flag set - Task {task_id} cancelled: "
                f"{progress_after.get('cancelled')}"
            )

        return cancellation_flag_set_time
    return None


def wait_for_task_completion(
    task_id: str,
    tasks: Dict[str, Any],
    progress_trackers: Dict[str, Any],
    quick_wait_time: float = 2.0,
    quick_wait_interval: float = 0.2,
) -> bool:
    """
    Wait for task to complete cancellation quickly.

    Args:
        task_id: Task ID
        tasks: Dictionary of tasks
        progress_trackers: Dictionary of progress trackers
        quick_wait_time: Maximum time to wait in seconds
        quick_wait_interval: Check interval in seconds

    Returns:
        True if task completed quickly, False otherwise
    """
    quick_waited = 0.0

    while quick_waited < quick_wait_time:
        time.sleep(quick_wait_interval)
        quick_waited += quick_wait_interval

        if task_id not in tasks:
            logger.warning(f"[STOP] Task {task_id} no longer exists")
            return True

        task = tasks[task_id]
        # Check if task has responded to cancellation
        if task.status in [TaskStatus.CANCELLED, TaskStatus.FAILED, TaskStatus.COMPLETED]:
            logger.info(
                f"[STOP] Task {task_id} responded quickly to stop request "
                f"(status: {task.status.value}, waited: {quick_waited:.1f}s)"
            )
            return True

        # Check if progress tracker shows cancelled and thread is done
        if task_id in progress_trackers:
            progress = progress_trackers[task_id].get_progress()
            if progress.get("cancelled", False):
                # Check if thread is still alive
                # If cancelled flag is set and we've waited a bit, assume it's processing
                if quick_waited >= quick_wait_time * 0.5:  # Wait at least half the quick wait time
                    return False  # Continue to longer wait

    return False  # Didn't complete quickly, need longer wait
