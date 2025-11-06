"""
Training service implementation.
"""
import os
import uuid
import threading
import logging
from typing import Dict, Optional, List
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor

from compass.service.models.task import TaskStatus, TrainingTaskResponse
from compass.service.config import SERVICE_CONFIG
from compass.logger import TrainingLogger
from compass.main import main
from compass.config import get_config
from compass.service.services.progress_tracker import ProgressTracker
from compass.service.services.progress_logger import ProgressAwareLogger

logger = logging.getLogger(__name__)


class TrainingService:
    """Service for managing training tasks."""
    
    def __init__(self):
        """Initialize training service."""
        self.tasks: Dict[str, TrainingTaskResponse] = {}
        self.task_logs: Dict[str, List[str]] = {}
        self.executor = ThreadPoolExecutor(max_workers=SERVICE_CONFIG['max_workers'])
        self.task_threads: Dict[str, threading.Thread] = {}
        self.progress_trackers: Dict[str, ProgressTracker] = {}
        self.lock = threading.Lock()
    
    def create_task(self, config: Dict, dataset_id: Optional[str] = None, description: Optional[str] = None) -> str:
        """
        Create a new training task.
        
        Args:
            config: Training configuration
            dataset_id: Optional dataset ID
            description: Optional task description
            
        Returns:
            str: Task ID
        """
        task_id = str(uuid.uuid4())
        
        task = TrainingTaskResponse(
            task_id=task_id,
            status=TaskStatus.PENDING,
            config=config,
            created_at=datetime.now(),
            updated_at=datetime.now(),
            description=description
        )
        
        with self.lock:
            self.tasks[task_id] = task
            self.task_logs[task_id] = []
            self.progress_trackers[task_id] = ProgressTracker(task_id)
        
        logger.info(f"Created training task: {task_id}")
        return task_id
    
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
        Start a training task.
        
        Args:
            task_id: Task ID
            
        Returns:
            bool: True if started successfully
        """
        with self.lock:
            if task_id not in self.tasks:
                return False
            
            task = self.tasks[task_id]
            if task.status != TaskStatus.PENDING and task.status != TaskStatus.PAUSED:
                return False
            
            task.status = TaskStatus.RUNNING
            task.started_at = datetime.now()
            task.updated_at = datetime.now()
            # Initialize progress
            if task_id in self.progress_trackers:
                task.progress = self.progress_trackers[task_id].get_progress()
        
        # Start training in background thread
        def run_training():
            try:
                self._run_training(task_id)
            except Exception as e:
                logger.error(f"Training task {task_id} failed: {e}", exc_info=True)
                with self.lock:
                    if task_id in self.tasks:
                        self.tasks[task_id].status = TaskStatus.FAILED
                        self.tasks[task_id].error = str(e)
                        self.tasks[task_id].updated_at = datetime.now()
                        self._log(task_id, f"ERROR: {e}")
            finally:
                # Clean up thread reference and progress tracker when task completes
                with self.lock:
                    if task_id in self.task_threads:
                        thread = self.task_threads.pop(task_id)
                        logger.debug(f"Cleaned up thread for task {task_id}")
                    # Note: Keep progress tracker for querying completed tasks
                    # Only remove if task is deleted
        
        thread = threading.Thread(target=run_training, daemon=True)
        thread.start()
        self.task_threads[task_id] = thread
        
        logger.info(f"Started training task: {task_id}")
        return True
    
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
            
            main(training_config, logger_instance)
            
            with self.lock:
                if task_id in self.tasks:
                    # Check if cancelled before marking as completed
                    if progress_tracker.is_cancelled():
                        self.tasks[task_id].status = TaskStatus.CANCELLED
                        self.tasks[task_id].updated_at = datetime.now()
                        self._log(task_id, "Training cancelled")
                    else:
                        self.tasks[task_id].status = TaskStatus.COMPLETED
                        self.tasks[task_id].completed_at = datetime.now()
                        self.tasks[task_id].updated_at = datetime.now()
                        progress_tracker.set_completed("Training completed successfully")
                        self.tasks[task_id].progress = progress_tracker.get_progress()
                        self._log(task_id, "Training completed successfully")
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
                raise
    
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
        Delete a training task.
        
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
            if task_id in self.progress_trackers:
                del self.progress_trackers[task_id]
        
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
        self.task_logs[task_id].append(f"[{timestamp}] {message}")

