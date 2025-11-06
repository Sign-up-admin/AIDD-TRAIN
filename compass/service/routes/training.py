"""
Training task routes.
"""

from fastapi import APIRouter, HTTPException, status
from typing import Optional
import logging

from compass.service.models.task import (
    TrainingTaskCreate,
    TrainingTaskResponse,
    TaskListResponse,
    TaskLogResponse,
    TaskMetricsResponse,
)
from compass.service.services.training_service import TrainingService
from compass.service.exceptions import ServiceException, NotFoundError, ValidationError
from compass.service.error_codes import ErrorCode

router = APIRouter(tags=["training"])
training_service = TrainingService()
logger = logging.getLogger(__name__)


@router.post(
    "/api/v1/training/tasks",
    response_model=TrainingTaskResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Create Training Task",
    description="Create a new training task with specified configuration",
    response_description="Created training task details",
)
async def create_task(request: TrainingTaskCreate):
    """
    Create a new training task.

    Args:
        request: Training task creation request containing:
            - config: Training configuration dictionary
                - execution_mode: One of ['validation_tuned', 'validation', 'prototyping', 'smoke_test', 'production']
                - epochs: Number of training epochs (1-10000)
                - batch_size: Batch size (1-128)
                - learning_rate: Learning rate (0-1.0)
                - optimizer: Optimizer name ('adam', 'adamw', 'sgd', 'rmsprop')
            - dataset_id: Optional dataset ID to use
            - description: Optional task description

    Returns:
        TrainingTaskResponse: Created task with task_id and initial status

    Raises:
        400: Invalid configuration parameters
        500: Internal server error

    Example request:
        {
            "config": {
                "execution_mode": "validation_tuned",
                "epochs": 100,
                "batch_size": 32,
                "learning_rate": 0.001,
                "optimizer": "adam"
            },
            "description": "Training run for protein-ligand binding"
        }
    """
    try:
        task_id = training_service.create_task(
            config=request.config, dataset_id=request.dataset_id, description=request.description
        )
        task = training_service.get_task(task_id)
        if not task:
            raise NotFoundError("Task", task_id)
        return task
    except ValueError as e:
        raise ValidationError(str(e))
    except ServiceException:
        # Re-raise service exceptions as-is (they already have proper error codes)
        raise
    except Exception as e:
        logger.error(f"Failed to create task: {e}", exc_info=True)
        # Provide more detailed error information
        error_detail = {
            "error_type": type(e).__name__,
            "error_message": str(e),
            "config": request.config,
        }
        raise ServiceException(
            "Failed to create task",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_code=ErrorCode.INTERNAL_ERROR,
            detail=error_detail,
        )


@router.get(
    "/api/v1/training/tasks",
    response_model=TaskListResponse,
    summary="List Training Tasks",
    description="Get a list of all training tasks with their current status",
)
async def list_tasks():
    """
    List all training tasks.

    Returns:
        TaskListResponse: List of all training tasks with their status and progress

    Example response:
        {
            "tasks": [
                {
                    "task_id": "uuid",
                    "status": "running",
                    "progress": {...},
                    ...
                }
            ],
            "count": 1
        }
    """
    tasks = training_service.list_tasks()
    return TaskListResponse(tasks=tasks, count=len(tasks))


@router.get(
    "/api/v1/training/tasks/{task_id}",
    response_model=TrainingTaskResponse,
    summary="Get Training Task",
    description="Get detailed information about a specific training task",
)
async def get_task(task_id: str):
    """
    Get task by ID.

    Args:
        task_id: Unique task identifier (UUID)

    Returns:
        TrainingTaskResponse: Task details including status, progress, and configuration

    Raises:
        404: Task not found

    Example response:
        {
            "task_id": "uuid",
            "status": "running",
            "config": {...},
            "progress": {
                "stage": "training",
                "current_epoch": 5,
                "total_epochs": 100,
                ...
            },
            ...
        }
    """
    task = training_service.get_task(task_id)
    if not task:
        raise NotFoundError("Task", task_id)
    return task


@router.post(
    "/api/v1/training/tasks/{task_id}/start",
    summary="Start Training Task",
    description="Start a pending or paused training task",
    response_description="Task start confirmation",
)
async def start_task(task_id: str):
    """
    Start a training task.

    Args:
        task_id: Unique task identifier (UUID)

    Returns:
        Dict with confirmation message and task_id

    Raises:
        400: Task cannot be started (wrong status or not found)
        500: Internal server error

    Note:
        - Task must be in 'pending' or 'paused' status
        - Starting a task will begin training in a background thread
    """
    success = training_service.start_task(task_id)
    if not success:
        raise ServiceException(
            f"Failed to start task {task_id}", status_code=status.HTTP_400_BAD_REQUEST
        )
    return {"message": f"Task {task_id} started", "task_id": task_id}


@router.post(
    "/api/v1/training/tasks/{task_id}/stop",
    summary="Stop Training Task",
    description="Stop a running training task. The training loop will check for cancellation and stop gracefully.",
    response_description="Task stop confirmation",
)
async def stop_task(task_id: str):
    """
    Stop a training task.

    Args:
        task_id: Unique task identifier (UUID)

    Returns:
        Dict with confirmation message and task_id

    Raises:
        400: Task cannot be stopped (wrong status or not found)

    Note:
        - Task must be in 'running' status
        - Training will stop at the next batch/epoch checkpoint
        - Task status will be set to 'cancelled'
    """
    success = training_service.stop_task(task_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"Failed to stop task {task_id}"
        )
    return {"message": f"Task {task_id} stopped", "task_id": task_id}


@router.post("/api/v1/training/tasks/{task_id}/pause")
async def pause_task(task_id: str):
    """Pause a training task."""
    success = training_service.pause_task(task_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"Failed to pause task {task_id}"
        )
    return {"message": f"Task {task_id} paused", "task_id": task_id}


@router.delete("/api/v1/training/tasks/{task_id}")
async def delete_task(task_id: str):
    """Delete a training task."""
    success = training_service.delete_task(task_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Task {task_id} not found"
        )
    return {"message": f"Task {task_id} deleted", "task_id": task_id}


@router.get(
    "/api/v1/training/tasks/{task_id}/logs",
    response_model=TaskLogResponse,
    summary="Get Task Logs",
    description="Get log messages for a training task",
)
async def get_task_logs(task_id: str, limit: int = 100):
    """
    Get task logs.

    Args:
        task_id: Unique task identifier (UUID)
        limit: Maximum number of log lines to return (default: 100, max: 1000)

    Returns:
        TaskLogResponse: List of log messages with timestamps

    Raises:
        404: Task not found

    Example response:
        {
            "logs": [
                "[2025-01-15 10:30:00] Starting training...",
                "[2025-01-15 10:31:00] Epoch 1/100 completed"
            ],
            "total_lines": 2
        }
    """
    task = training_service.get_task(task_id)
    if not task:
        raise NotFoundError("Task", task_id)

    logs = training_service.get_logs(task_id, limit=limit)
    return TaskLogResponse(logs=logs, total_lines=len(logs))


@router.get("/api/v1/training/tasks/{task_id}/metrics", response_model=TaskMetricsResponse)
async def get_task_metrics(task_id: str):
    """Get task metrics."""
    task = training_service.get_task(task_id)
    if not task:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND, detail=f"Task {task_id} not found"
        )

    return TaskMetricsResponse(metrics=task.progress)


@router.get(
    "/api/v1/training/tasks/{task_id}/resources",
    summary="Get Task Resource Usage",
    description="Get resource usage information (GPU, CPU, Memory) for a training task",
)
async def get_task_resources(task_id: str):
    """
    Get resource usage for a training task.

    Args:
        task_id: Unique task identifier (UUID)

    Returns:
        Dict: Resource usage information including GPU, CPU, and memory usage

    Raises:
        404: Task not found
    """
    task = training_service.get_task(task_id)
    if not task:
        raise NotFoundError("Task", task_id)

    resources = training_service.get_resource_usage(task_id)
    if resources is None:
        # Return empty resources if not available yet
        return {
            "cpu_percent": 0.0,
            "memory": {"total_gb": 0.0, "used_gb": 0.0, "percent": 0.0},
            "gpu": {"available": False},
        }

    return resources
