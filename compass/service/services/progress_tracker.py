"""
Progress tracker for training tasks.
"""

import threading
import time
import multiprocessing
from typing import Dict, Optional
from datetime import datetime


class ProgressTracker:
    """Tracks training progress in real-time."""

    def __init__(self, task_id: str):
        """
        Initialize progress tracker.

        Args:
            task_id: Task ID
        """
        self.task_id = task_id
        self.lock = threading.Lock()
        # Use Event for cancellation flag to avoid lock contention
        # Event is thread-safe and allows lock-free reads
        self._cancellation_event = threading.Event()
        # Use multiprocessing.Value for cross-process cancellation flag
        # This allows worker processes in multiprocessing.Pool to check cancellation
        self._mp_cancellation_flag = multiprocessing.Value(
            "i", 0
        )  # 0 = not cancelled, 1 = cancelled
        self.current_stage = (
            "initializing"  # initializing, data_processing, training, validation, completed
        )
        self.stage_progress = 0.0  # 0.0 to 1.0
        self.stage_message = ""

        # Cancellation flag (kept for backward compatibility, but use _cancellation_event for checks)
        self.cancelled = False

        # Data processing progress
        self.data_processing_total = 0
        self.data_processing_completed = 0

        # Training progress
        self.current_epoch = 0
        self.total_epochs = 0
        self.current_batch = 0
        self.total_batches = 0
        self.train_loss = 0.0
        self.val_loss = 0.0

        # Metrics
        self.start_time = None
        self.last_update_time = None

    def set_stage(self, stage: str, message: str = ""):
        """
        Set current stage.

        Args:
            stage: Stage name (initializing, data_processing, training, validation, completed)
            message: Stage message
        """
        with self.lock:
            self.current_stage = stage
            self.stage_message = message
            self.last_update_time = datetime.now()
            if self.start_time is None:
                self.start_time = datetime.now()

    def update_data_processing(self, completed: int, total: int, message: str = ""):
        """
        Update data processing progress.

        Args:
            completed: Number of completed items
            total: Total number of items
            message: Progress message
        """
        with self.lock:
            self.current_stage = "data_processing"
            self.data_processing_total = total
            self.data_processing_completed = completed
            if total > 0:
                self.stage_progress = completed / total
            self.stage_message = message or f"Processing {completed}/{total} data points"
            self.last_update_time = datetime.now()

    def update_training(
        self,
        epoch: int,
        total_epochs: int,
        batch: int,
        total_batches: int,
        train_loss: float = 0.0,
        val_loss: float = 0.0,
        message: str = "",
    ):
        """
        Update training progress.

        Args:
            epoch: Current epoch
            total_epochs: Total epochs
            batch: Current batch
            total_batches: Total batches in current epoch
            train_loss: Current training loss
            val_loss: Current validation loss
            message: Progress message
        """
        with self.lock:
            self.current_stage = "training"
            self.current_epoch = epoch
            self.total_epochs = total_epochs
            self.current_batch = batch
            self.total_batches = total_batches
            self.train_loss = train_loss
            self.val_loss = val_loss

            # Calculate overall progress
            if total_epochs > 0 and total_batches > 0:
                epoch_progress = (epoch - 1) / total_epochs
                batch_progress = batch / total_batches / total_epochs
                self.stage_progress = epoch_progress + batch_progress

            if not message:
                message = f"Epoch {epoch}/{total_epochs}, Batch {batch}/{total_batches}"
                if train_loss > 0:
                    message += f", Train Loss: {train_loss:.4f}"
                if val_loss > 0:
                    message += f", Val Loss: {val_loss:.4f}"
            self.stage_message = message
            self.last_update_time = datetime.now()

    def set_completed(self, message: str = "Training completed"):
        """Mark training as completed."""
        with self.lock:
            self.current_stage = "completed"
            self.stage_progress = 1.0
            self.stage_message = message
            self.last_update_time = datetime.now()

    def cancel(self):
        """Mark task as cancelled."""
        with self.lock:
            self.cancelled = True
            self._cancellation_event.set()  # Set event (thread-safe, no lock needed for readers)
            # Set multiprocessing flag for cross-process cancellation
            with self._mp_cancellation_flag.get_lock():
                self._mp_cancellation_flag.value = 1
            self.stage_message = "Training cancelled by user"
            self.last_update_time = datetime.now()

    def is_cancelled(self) -> bool:
        """
        Check if task is cancelled.

        Uses threading.Event for lock-free reads, improving performance
        when called frequently from training loops.
        """
        # Lock-free read using Event.is_set() - much faster than acquiring lock
        return self._cancellation_event.is_set()

    def is_cancelled_mp(self) -> bool:
        """
        Check if task is cancelled (multiprocessing-safe version).

        This method can be called from worker processes in multiprocessing.Pool.
        Returns True if cancellation has been requested.
        """
        with self._mp_cancellation_flag.get_lock():
            return self._mp_cancellation_flag.value == 1

    def get_mp_cancellation_flag(self):
        """
        Get the multiprocessing cancellation flag for passing to worker processes.

        Returns:
            multiprocessing.Value: The shared cancellation flag
        """
        return self._mp_cancellation_flag

    def get_progress(self) -> Dict:
        """
        Get current progress information.

        Returns:
            Dict: Progress information
        """
        with self.lock:
            elapsed_time = None
            if self.start_time:
                elapsed = datetime.now() - self.start_time
                elapsed_time = elapsed.total_seconds()

            # Use event state for cancelled flag (more reliable)
            cancelled = self._cancellation_event.is_set()
            return {
                "stage": self.current_stage,
                "progress": self.stage_progress,
                "message": self.stage_message,
                "cancelled": cancelled,
                "data_processing": {
                    "completed": self.data_processing_completed,
                    "total": self.data_processing_total,
                    "percentage": (
                        (self.data_processing_completed / self.data_processing_total * 100)
                        if self.data_processing_total > 0
                        else 0
                    ),
                },
                "training": {
                    "current_epoch": self.current_epoch,
                    "total_epochs": self.total_epochs,
                    "current_batch": self.current_batch,
                    "total_batches": self.total_batches,
                    "train_loss": self.train_loss,
                    "val_loss": self.val_loss,
                    "epoch_progress": (
                        (self.current_epoch / self.total_epochs * 100)
                        if self.total_epochs > 0
                        else 0
                    ),
                },
                "elapsed_time": elapsed_time,
                "last_update": self.last_update_time.isoformat() if self.last_update_time else None,
            }

    def __getstate__(self):
        """
        Custom pickle state getter.
        Excludes threading.Lock, threading.Event, and multiprocessing.Value objects which cannot be pickled.
        """
        state = self.__dict__.copy()
        # Remove lock, event, and multiprocessing objects as they cannot be pickled
        state.pop("lock", None)
        state.pop("_cancellation_event", None)
        state.pop("_mp_cancellation_flag", None)
        return state

    def __setstate__(self, state):
        """
        Custom pickle state setter.
        Recreates threading.Lock, threading.Event, and multiprocessing.Value objects after unpickling.
        """
        self.__dict__.update(state)
        # Recreate lock, event, and multiprocessing objects after unpickling
        self.lock = threading.Lock()
        self._cancellation_event = threading.Event()
        self._mp_cancellation_flag = multiprocessing.Value("i", 0)
        # Restore cancellation state if it was set
        if state.get("cancelled", False):
            self._cancellation_event.set()
            with self._mp_cancellation_flag.get_lock():
                self._mp_cancellation_flag.value = 1
