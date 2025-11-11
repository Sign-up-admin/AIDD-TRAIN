"""
Task lifecycle logging system for tracking training task creation, execution, and completion.
"""

import json
import logging
import traceback
from typing import Dict, Optional, Any
from datetime import datetime
from pathlib import Path


class TaskLifecycleLogger:
    """Logger for tracking task lifecycle events with structured logging."""

    def __init__(self, task_id: str, log_dir: str = "logs"):
        """
        Initialize task lifecycle logger.

        Args:
            task_id: Task ID
            log_dir: Directory for log files
        """
        self.task_id = task_id
        self.log_dir = Path(log_dir)
        self.log_dir.mkdir(parents=True, exist_ok=True)

        # Create task-specific log file
        self.lifecycle_log_file = self.log_dir / f"task_{task_id}_lifecycle.jsonl"
        self.error_log_file = self.log_dir / f"task_{task_id}_errors.log"

        # Standard logger for compatibility
        self.logger = logging.getLogger(f"task_lifecycle.{task_id}")

        # Event timestamps tracking
        self.event_timestamps: Dict[str, datetime] = {}
        self.start_time: Optional[datetime] = None

    def _write_lifecycle_event(
        self,
        event_type: str,
        stage: str,
        status: str,
        message: str,
        data: Optional[Dict] = None,
        error: Optional[Exception] = None,
    ):
        """
        Write a structured lifecycle event to the log file.

        Args:
            event_type: Type of event (create, start, run, complete, fail, cancel)
            stage: Current stage (creating, initializing, running, etc.)
            status: Event status (success, failed, timeout)
            message: Event message
            data: Additional data to log
            error: Exception if any
        """
        event = {
            "timestamp": datetime.now().isoformat(),
            "task_id": self.task_id,
            "event_type": event_type,
            "stage": stage,
            "status": status,
            "message": message,
        }

        # Calculate elapsed time if start time is set
        if self.start_time:
            elapsed = (datetime.now() - self.start_time).total_seconds()
            event["elapsed_seconds"] = elapsed

        # Add timing for this event
        if event_type in self.event_timestamps:
            prev_time = self.event_timestamps[event_type]
            event["duration_seconds"] = (datetime.now() - prev_time).total_seconds()

        # Store timestamp for this event type
        self.event_timestamps[event_type] = datetime.now()

        # Add resource usage if available
        if data:
            event["data"] = data

        # Add error information if present
        if error:
            event["error"] = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            }
            # Also write to error log file
            self._write_error_log(event_type, stage, message, error)

        # Write to JSONL file
        try:
            with open(self.lifecycle_log_file, "a", encoding="utf-8") as f:
                f.write(json.dumps(event, ensure_ascii=False) + "\n")
        except Exception as e:
            self.logger.error(f"Failed to write lifecycle event: {e}", exc_info=True)

    def _write_error_log(self, event_type: str, stage: str, message: str, error: Exception):
        """Write error details to error log file."""
        try:
            with open(self.error_log_file, "a", encoding="utf-8") as f:
                f.write(f"[{datetime.now().isoformat()}] {event_type} - {stage}\n")
                f.write(f"Message: {message}\n")
                f.write(f"Error: {type(error).__name__}: {str(error)}\n")
                f.write(f"Traceback:\n{traceback.format_exc()}\n")
                f.write("-" * 80 + "\n")
        except Exception as e:
            self.logger.error(f"Failed to write error log: {e}", exc_info=True)

    def log_create_start(self, config: Dict, resource_check: Optional[Dict] = None):
        """Log task creation start."""
        self.start_time = datetime.now()
        data = {"config": config}
        if resource_check:
            data["resource_check"] = resource_check
        self._write_lifecycle_event(
            event_type="create",
            stage="creating",
            status="start",
            message="Task creation started",
            data=data,
        )
        self.logger.info(f"Task {self.task_id} creation started")

    def log_create_success(self):
        """Log successful task creation."""
        self._write_lifecycle_event(
            event_type="create",
            stage="creating",
            status="success",
            message="Task created successfully",
        )
        self.logger.info(f"Task {self.task_id} created successfully")

    def log_create_failure(self, error: Exception, reason: str):
        """Log task creation failure."""
        self._write_lifecycle_event(
            event_type="create",
            stage="creating",
            status="failed",
            message=f"Task creation failed: {reason}",
            error=error,
        )
        self.logger.error(f"Task {self.task_id} creation failed: {reason}", exc_info=True)

    def log_initialize_start(self):
        """Log task initialization start."""
        self._write_lifecycle_event(
            event_type="initialize",
            stage="initializing",
            status="start",
            message="Task initialization started",
        )
        self.logger.info(f"Task {self.task_id} initialization started")

    def log_initialize_success(self):
        """Log successful task initialization."""
        self._write_lifecycle_event(
            event_type="initialize",
            stage="initializing",
            status="success",
            message="Task initialized successfully",
        )
        self.logger.info(f"Task {self.task_id} initialized successfully")

    def log_initialize_failure(self, error: Exception, reason: str):
        """Log task initialization failure."""
        self._write_lifecycle_event(
            event_type="initialize",
            stage="initializing",
            status="failed",
            message=f"Task initialization failed: {reason}",
            error=error,
        )
        self.logger.error(f"Task {self.task_id} initialization failed: {reason}", exc_info=True)

    def log_start(self):
        """Log task start."""
        self._write_lifecycle_event(
            event_type="start", stage="pending", status="success", message="Task started"
        )
        self.logger.info(f"Task {self.task_id} started")

    def log_running(self, stage: str, progress: Optional[Dict] = None):
        """Log task running status."""
        data = {"progress": progress} if progress else None
        self._write_lifecycle_event(
            event_type="running",
            stage=stage,
            status="success",
            message=f"Task running in stage: {stage}",
            data=data,
        )

    def log_resource_usage(self, resources: Dict):
        """Log resource usage."""
        self._write_lifecycle_event(
            event_type="resource_check",
            stage="running",
            status="success",
            message="Resource usage check",
            data={"resources": resources},
        )

    def log_complete(self, message: str = "Task completed successfully"):
        """Log task completion."""
        self._write_lifecycle_event(
            event_type="complete", stage="completed", status="success", message=message
        )
        self.logger.info(f"Task {self.task_id} completed: {message}")

    def log_failure(self, error: Exception, stage: str, context: Optional[Dict] = None):
        """Log task failure."""
        data = {"context": context} if context else None
        self._write_lifecycle_event(
            event_type="fail",
            stage=stage,
            status="failed",
            message=f"Task failed in stage: {stage}",
            error=error,
            data=data,
        )
        self.logger.error(f"Task {self.task_id} failed in stage {stage}", exc_info=True)

    def log_cancel(self, reason: str = "Task cancelled by user"):
        """Log task cancellation."""
        self._write_lifecycle_event(
            event_type="cancel", stage="cancelled", status="success", message=reason
        )
        self.logger.info(f"Task {self.task_id} cancelled: {reason}")

    def log_timeout(self, stage: str, timeout_seconds: float):
        """Log task timeout."""
        self._write_lifecycle_event(
            event_type="timeout",
            stage=stage,
            status="failed",
            message=f"Task timed out in stage: {stage} after {timeout_seconds}s",
            data={"timeout_seconds": timeout_seconds},
        )
        self.logger.warning(
            f"Task {self.task_id} timed out in stage {stage} after {timeout_seconds}s"
        )

    def log_state_transition(self, from_status: str, to_status: str, reason: str = ""):
        """Log state transition."""
        self._write_lifecycle_event(
            event_type="state_transition",
            stage=to_status,
            status="success",
            message=f"State transition: {from_status} -> {to_status}"
            + (f" ({reason})" if reason else ""),
            data={"from_status": from_status, "to_status": to_status, "reason": reason},
        )
        self.logger.debug(f"Task {self.task_id} state transition: {from_status} -> {to_status}")

    def get_lifecycle_events(self, limit: int = 100) -> list:
        """
        Read lifecycle events from log file.

        Args:
            limit: Maximum number of events to return

        Returns:
            List of event dictionaries
        """
        events = []
        if not self.lifecycle_log_file.exists():
            return events

        try:
            with open(self.lifecycle_log_file, "r", encoding="utf-8") as f:
                lines = f.readlines()
                for line in lines[-limit:]:
                    try:
                        event = json.loads(line.strip())
                        events.append(event)
                    except json.JSONDecodeError:
                        continue
        except Exception as e:
            self.logger.error(f"Failed to read lifecycle events: {e}", exc_info=True)

        return events










