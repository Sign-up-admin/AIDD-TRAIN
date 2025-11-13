"""
Training task routes.
"""

from fastapi import APIRouter, HTTPException, status, WebSocket, WebSocketDisconnect
from typing import Optional
import logging
import json
import asyncio

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
from compass.service.utils.input_sanitizer import sanitize_task_id
from compass.service.utils.websocket_manager import (
    WebSocketConnectionState,
    wait_for_stream_queues,
    send_log_messages,
    send_resource_updates,
    receive_client_messages,
    send_heartbeat,
)

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
    except (KeyError, TypeError, AttributeError) as e:
        # Configuration-related errors
        logger.error(f"Invalid configuration when creating task: {e}", exc_info=True)
        raise ValidationError(f"Invalid configuration: {str(e)}")
    except (OSError, IOError, PermissionError) as e:
        # File system or permission errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(
            f"File system error when creating task: {e}",
            exc_info=True,
            extra={"error_type": type(e).__name__},
        )
        raise ServiceException(
            f"Failed to create task: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_code=ErrorCode.INTERNAL_ERROR,
        )
    except (ConnectionError, TimeoutError) as e:
        # Network or timeout errors
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(
            f"Network error when creating task: {e}",
            exc_info=True,
            extra={"error_type": type(e).__name__},
        )
        raise ServiceException(
            f"Failed to create task: {error_message}",
            status_code=status.HTTP_503_SERVICE_UNAVAILABLE,
            error_code=ErrorCode.SERVICE_UNAVAILABLE,
        )
    except Exception as e:
        # Unexpected errors - sanitize error message before exposing
        from compass.service.exceptions import sanitize_error_message

        error_message = sanitize_error_message(e, include_details=False)
        logger.error(
            f"Unexpected error when creating task: {e}",
            exc_info=True,
            extra={
                "error_type": type(e).__name__,
                "config_keys": list(request.config.keys()) if request.config else None,
            },
        )
        raise ServiceException(
            f"Failed to create task: {error_message}",
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            error_code=ErrorCode.INTERNAL_ERROR,
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
        400: Invalid task ID format

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
    # Validate and sanitize task_id to prevent injection attacks
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

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
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

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
        400: Task cannot be stopped (wrong status or not found) or invalid task ID format
        404: Task not found

    Note:
        - Task must be in 'running' status
        - Training will stop at the next batch/epoch checkpoint
        - Task status will be set to 'cancelled'
    """
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

    success, error_message = training_service.stop_task(task_id)
    if not success:
        # Check if task exists to determine appropriate status code
        task = training_service.get_task(task_id)
        if not task:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND, detail=f"Task {task_id} not found"
            )
        else:
            # Task exists but cannot be stopped (wrong status)
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=error_message or f"Failed to stop task {task_id}",
            )
    return {"message": f"Task {task_id} stopped", "task_id": task_id}


@router.post("/api/v1/training/tasks/{task_id}/pause")
async def pause_task(task_id: str):
    """Pause a training task."""
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

    success = training_service.pause_task(task_id)
    if not success:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST, detail=f"Failed to pause task {task_id}"
        )
    return {"message": f"Task {task_id} paused", "task_id": task_id}


@router.delete("/api/v1/training/tasks/{task_id}")
async def delete_task(task_id: str):
    """Delete a training task."""
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

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
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

    # Validate limit
    if limit < 1 or limit > 1000:
        raise ValidationError("Limit must be between 1 and 1000")

    task = training_service.get_task(task_id)
    if not task:
        raise NotFoundError("Task", task_id)

    logs = training_service.get_logs(task_id, limit=limit)
    return TaskLogResponse(logs=logs, total_lines=len(logs))


@router.get("/api/v1/training/tasks/{task_id}/metrics", response_model=TaskMetricsResponse)
async def get_task_metrics(task_id: str):
    """Get task metrics."""
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

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
        400: Invalid task ID format
    """
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        raise ValidationError(f"Invalid task ID format: {str(e)}")

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


@router.websocket("/api/v1/training/tasks/{task_id}/stream")
async def stream_task_logs(websocket: WebSocket, task_id: str):
    """
    WebSocket endpoint for real-time log streaming and resource monitoring.

    Args:
        websocket: WebSocket connection
        task_id: Task ID to stream logs for

    Message format:
        - Outgoing: {"type": "log"|"resources", "data": ..., "timestamp": ...}
        - Incoming: {"type": "command"|"ping", "data": ...} (optional, for future use)
    """
    # Validate and sanitize task_id
    try:
        task_id = sanitize_task_id(task_id)
    except ValueError as e:
        await websocket.close(code=1008, reason=f"Invalid task ID format: {str(e)}")
        return

    # Accept WebSocket connection first (required by FastAPI)
    await websocket.accept()
    logger.info(f"WebSocket connection accepted for task {task_id}")

    # Verify task exists
    logger.info(f"WebSocket connection attempt for task {task_id}")
    task = training_service.get_task(task_id)
    if not task:
        logger.warning(f"WebSocket connection rejected: Task {task_id} not found")
        await websocket.send_json(
            {
                "type": "error",
                "data": f"Task {task_id} not found",
                "timestamp": None,
            }
        )
        await websocket.close(code=status.WS_1008_POLICY_VIOLATION, reason="Task not found")
        return

    logger.info(f"Task {task_id} found with status: {task.status.value}")
    logger.info(f"WebSocket connection established for task {task_id}")

    # Wait for stream queues to be available
    log_queue, resource_queue = await wait_for_stream_queues(
        training_service.stream_manager, task_id, websocket
    )

    # Final check - if still not available, send error and close
    if not log_queue or not resource_queue:
        error_msg = f"Stream not available for task {task_id}. Task may not be started yet."
        logger.error(error_msg)
        await websocket.send_json({"type": "error", "data": error_msg, "timestamp": None})
        await websocket.send_json(
            {
                "type": "error",
                "data": "Please ensure the task is in 'running' or 'initializing' status before connecting.",
                "timestamp": None,
            }
        )
        await websocket.close(code=status.WS_1008_POLICY_VIOLATION, reason="Stream not available")
        return

    # Send initial connection success message
    await websocket.send_json(
        {
            "type": "connected",
            "data": f"Connected to task {task_id} stream. Stream is ready.",
            "timestamp": None,
        }
    )
    logger.info(f"Stream ready for task {task_id}, starting log and resource streaming")

    # Create connection state manager
    connection_state = WebSocketConnectionState(task_id, ping_interval=30.0, ping_timeout=60.0)

    # Start background tasks
    log_task = asyncio.create_task(send_log_messages(websocket, log_queue, connection_state))
    resource_task = asyncio.create_task(
        send_resource_updates(websocket, resource_queue, connection_state)
    )
    receive_task = asyncio.create_task(receive_client_messages(websocket, connection_state))
    heartbeat_task = asyncio.create_task(send_heartbeat(websocket, connection_state))

    try:
        # Wait for any task to complete (usually receive_messages when client disconnects)
        done, pending = await asyncio.wait(
            [log_task, resource_task, receive_task, heartbeat_task],
            return_when=asyncio.FIRST_COMPLETED,
        )

        # Cancel remaining tasks
        for task in pending:
            task.cancel()
            try:
                await task
            except asyncio.CancelledError:
                pass
    except (asyncio.CancelledError, KeyboardInterrupt):
        logger.debug(f"WebSocket stream cancelled for task {task_id}")
    except (ConnectionError, OSError) as e:
        logger.warning(f"Connection error in WebSocket stream for task {task_id}: {e}")
    except (WebSocketDisconnect, RuntimeError) as e:
        logger.info(f"WebSocket disconnected for task {task_id}: {e}")
    except Exception as e:
        logger.error(
            f"Unexpected error in WebSocket stream for task {task_id}: {e}",
            exc_info=True,
            extra={"task_id": task_id, "error_type": type(e).__name__},
        )
    finally:
        connection_state.mark_dead()
        logger.info(f"WebSocket connection closed for task {task_id}")
        try:
            await websocket.close()
        except (ConnectionError, OSError, RuntimeError) as e:
            logger.debug(f"Error closing WebSocket for task {task_id}: {e}")
        except Exception as e:
            logger.warning(f"Unexpected error closing WebSocket for task {task_id}: {e}")
