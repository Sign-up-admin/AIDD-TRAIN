"""
Progress-aware logger wrapper for training tasks.
"""

import re
from compass.logger import TrainingLogger
from compass.service.services.progress_tracker import ProgressTracker


class ProgressAwareLogger(TrainingLogger):
    """Logger that updates progress tracker based on log messages."""

    def __init__(self, log_dir, progress_tracker: ProgressTracker):
        """
        Initialize progress-aware logger.

        Args:
            log_dir: Log directory
            progress_tracker: Progress tracker instance
        """
        super().__init__(log_dir)
        self.progress_tracker = progress_tracker

    def log(self, message):
        """Log message and update progress if applicable."""
        super().log(message)
        self._update_progress_from_message(message)

    def log_warning(self, message):
        """Log warning and update progress if applicable."""
        super().log_warning(message)
        self._update_progress_from_message(message)

    def log_error(self, message):
        """Log error and update progress if applicable."""
        super().log_error(message)
        self._update_progress_from_message(message)

    def _update_progress_from_message(self, message: str):
        """Parse message and update progress tracker."""
        # Check if progress_tracker exists (may be None after unpickling)
        if not hasattr(self, "progress_tracker") or self.progress_tracker is None:
            return

        # Data processing progress
        if "Processing" in message and "items" in message:
            match = re.search(r"Processing (\d+)/(\d+)", message)
            if match:
                completed = int(match.group(1))
                total = int(match.group(2))
                self.progress_tracker.update_data_processing(completed, total, message)

        # Training epoch progress - match format: "Epoch 01/200 | Train Loss: 1.2345 | Val Loss: 1.2345 | LR: 0.000100"
        epoch_match = re.search(r"Epoch (\d+)/(\d+)", message)
        train_loss_match = re.search(r"Train Loss[:\s]+([\d.]+)", message, re.IGNORECASE)
        val_loss_match = re.search(r"Val Loss[:\s]+([\d.]+)", message, re.IGNORECASE)

        if epoch_match:
            epoch = int(epoch_match.group(1))
            total_epochs = int(epoch_match.group(2))

            # Get batch info from training state if available
            batch = 0
            total_batches = 0
            # Try to parse from message or use default
            batch_match = re.search(r"batch (\d+)/(\d+)", message, re.IGNORECASE)
            if batch_match:
                batch = int(batch_match.group(1))
                total_batches = int(batch_match.group(2))

            train_loss = 0.0
            if train_loss_match:
                train_loss = float(train_loss_match.group(1))

            val_loss = 0.0
            if val_loss_match:
                val_loss = float(val_loss_match.group(1))

            self.progress_tracker.update_training(
                epoch,
                total_epochs,
                batch,
                total_batches,
                train_loss=train_loss,
                val_loss=val_loss,
                message=message,
            )

        # Stage detection (only if progress_tracker exists)
        if hasattr(self, "progress_tracker") and self.progress_tracker is not None:
            if "Step 1:" in message or "Parsing PDBbind" in message:
                self.progress_tracker.set_stage("data_processing", "Parsing PDBbind index files")
            elif "Step 2:" in message or "Verifying data" in message:
                self.progress_tracker.set_stage("data_processing", "Processing dataset")
            elif "Step 3:" in message or "Splitting data" in message:
                self.progress_tracker.set_stage("data_processing", "Preparing data loaders")
            elif "Step 4:" in message or "Setting up model" in message:
                self.progress_tracker.set_stage("initializing", "Initializing model")
            elif "Step 5:" in message or "Selecting training recipe" in message:
                self.progress_tracker.set_stage("initializing", "Setting up training recipe")
            elif "Step 6:" in message or "Starting training" in message:
                self.progress_tracker.set_stage("training", "Starting training")
            elif "Training Finished" in message or "completed" in message.lower():
                self.progress_tracker.set_completed("Training completed successfully")

    def __getstate__(self):
        """
        Custom pickle state getter.
        Excludes progress_tracker which contains threading.Lock that cannot be pickled.
        """
        state = self.__dict__.copy()
        # Remove progress_tracker as it contains threading.Lock
        # The progress_tracker will be recreated if needed after unpickling
        state.pop("progress_tracker", None)
        return state

    def __setstate__(self, state):
        """
        Custom pickle state setter.
        Progress tracker will be None after unpickling - this is acceptable
        as the logger can function without it for basic logging.
        """
        self.__dict__.update(state)
        # Progress tracker is not restored - this is intentional
        # as it cannot be pickled and is not needed in subprocesses
        # Methods that use progress_tracker will check if it exists before using it
