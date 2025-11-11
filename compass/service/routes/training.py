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
        404: Task not found

    Note:
        - Task must be in 'running' status
        - Training will stop at the next batch/epoch checkpoint
        - Task status will be set to 'cancelled'
    """
    success, error_message = training_service.stop_task(task_id)
    if not success:
        # Check if task exists to determine appropriate status code
        task = training_service.get_task(task_id)
        if not task:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Task {task_id} not found"
            )
        else:
            # Task exists but cannot be stopped (wrong status)
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=error_message or f"Failed to stop task {task_id}"
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
    # Accept WebSocket connection first (required by FastAPI)
    await websocket.accept()
    logger.info(f"WebSocket connection accepted for task {task_id}")
    
    # Verify task exists
    logger.info(f"WebSocket connection attempt for task {task_id}")
    task = training_service.get_task(task_id)
    if not task:
        logger.warning(f"WebSocket connection rejected: Task {task_id} not found")
        await websocket.send_json({
            "type": "error",
            "data": f"Task {task_id} not found",
            "timestamp": None,
        })
        await websocket.close(code=status.WS_1008_POLICY_VIOLATION, reason="Task not found")
        return

    logger.info(f"Task {task_id} found with status: {task.status.value}")
    logger.info(f"WebSocket connection established for task {task_id}")

    # Wait for stream to be created (with timeout)
    # This handles the race condition where WebSocket connects before stream is created
    max_wait_time = 10.0  # Maximum wait time in seconds
    wait_interval = 0.1  # Check interval in seconds
    waited_time = 0.0

    log_queue = None
    resource_queue = None

    while waited_time < max_wait_time:
        log_queue = training_service.stream_manager.get_log_queue(task_id)
        resource_queue = training_service.stream_manager.get_resource_queue(task_id)

        if log_queue and resource_queue:
            logger.info(f"Stream queues found for task {task_id} after {waited_time:.2f}s")
            break

        # Send waiting message to client
        if waited_time == 0:
            await websocket.send_json(
                {
                    "type": "connected",
                    "data": f"Waiting for stream to be created for task {task_id}...",
                    "timestamp": None,
                }
            )

        await asyncio.sleep(wait_interval)
        waited_time += wait_interval

    # If streams still not available, try creating them
    if not log_queue or not resource_queue:
        logger.warning(
            f"Stream queues not found for task {task_id} after {max_wait_time}s, attempting to create"
        )
        training_service.stream_manager.create_stream(task_id)

        # Wait for creation to complete (with retry)
        max_retries = 10
        retry_count = 0
        while retry_count < max_retries and (not log_queue or not resource_queue):
            await asyncio.sleep(0.2)  # Wait 200ms between checks
            log_queue = training_service.stream_manager.get_log_queue(task_id)
            resource_queue = training_service.stream_manager.get_resource_queue(task_id)
            retry_count += 1
        
        if log_queue and resource_queue:
            logger.info(f"Stream queues created successfully for task {task_id} after {retry_count * 0.2:.1f}s")

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

    # Connection state
    connection_alive = True
    last_ping_time = asyncio.get_event_loop().time()
    ping_interval = 30.0  # Send ping every 30 seconds
    ping_timeout = 60.0  # Close connection if no pong received in 60 seconds

    def check_connection_alive():
        """Check if connection is still alive."""
        return connection_alive

    # Start background tasks
    async def send_logs():
        """Send log messages from queue."""
        nonlocal connection_alive
        try:
            while connection_alive:
                try:
                    # Wait for log message with timeout
                    message = await asyncio.wait_for(log_queue.get(), timeout=0.1)
                    if not connection_alive:
                        break

                    # Send log data directly as text to preserve ANSI escape codes
                    # Use JSON only for metadata, send data as raw text
                    log_data = message.get("data", "")
                    # Send as JSON but ensure ANSI codes are preserved (JSON will escape them properly)
                    await websocket.send_json(message)
                except asyncio.TimeoutError:
                    # Check if connection is still alive
                    if not connection_alive:
                        break
                    continue
                except Exception as e:
                    logger.error(
                        f"Error sending log message for task {task_id}: {e}", exc_info=True
                    )
                    connection_alive = False
                    break
        except asyncio.CancelledError:
            connection_alive = False
        except Exception as e:
            logger.error(f"Unexpected error in send_logs for task {task_id}: {e}", exc_info=True)
            connection_alive = False

    async def send_resources():
        """Send resource updates from queue."""
        nonlocal connection_alive
        try:
            while connection_alive:
                try:
                    # Wait for resource message with timeout
                    message = await asyncio.wait_for(resource_queue.get(), timeout=0.1)
                    if not connection_alive:
                        break
                    await websocket.send_json(message)
                except asyncio.TimeoutError:
                    if not connection_alive:
                        break
                    continue
                except Exception as e:
                    logger.error(
                        f"Error sending resource message for task {task_id}: {e}", exc_info=True
                    )
                    connection_alive = False
                    break
        except asyncio.CancelledError:
            connection_alive = False
        except Exception as e:
            logger.error(
                f"Unexpected error in send_resources for task {task_id}: {e}", exc_info=True
            )
            connection_alive = False

    async def receive_messages():
        """Receive messages from client (for ping/pong and future command support)."""
        nonlocal connection_alive, last_ping_time
        try:
            while connection_alive:
                try:
                    data = await websocket.receive_text()
                    try:
                        message = json.loads(data)
                        msg_type = message.get("type")

                        # Handle ping/pong for connection keepalive
                        if msg_type == "ping":
                            last_ping_time = asyncio.get_event_loop().time()
                            await websocket.send_json(
                                {"type": "pong", "data": "pong", "timestamp": None}
                            )
                        # Handle command messages (for future use)
                        elif msg_type == "command":
                            logger.info(
                                f"Received command from client for task {task_id}: {message.get('data')}"
                            )
                            # Future: Execute command and send output
                        else:
                            logger.debug(f"Received unknown message type: {msg_type}")
                    except json.JSONDecodeError:
                        # Handle plain text pings
                        if data.strip() == "ping":
                            last_ping_time = asyncio.get_event_loop().time()
                            await websocket.send_text("pong")
                        else:
                            logger.warning(
                                f"Invalid JSON received from client for task {task_id}: {data}"
                            )
                except WebSocketDisconnect:
                    logger.info(f"WebSocket disconnected by client for task {task_id}")
                    connection_alive = False
                    break
                except Exception as e:
                    logger.error(f"Error receiving message for task {task_id}: {e}", exc_info=True)
                    connection_alive = False
                    break
        except asyncio.CancelledError:
            connection_alive = False
        except Exception as e:
            logger.error(
                f"Unexpected error in receive_messages for task {task_id}: {e}", exc_info=True
            )
            connection_alive = False

    async def heartbeat():
        """Send periodic ping to check connection health."""
        nonlocal connection_alive, last_ping_time
        try:
            while connection_alive:
                await asyncio.sleep(ping_interval)
                if not connection_alive:
                    break

                current_time = asyncio.get_event_loop().time()
                # Check if we haven't received a pong in too long
                if current_time - last_ping_time > ping_timeout:
                    logger.warning(f"WebSocket connection timeout for task {task_id}, closing")
                    connection_alive = False
                    break

                # Send ping
                try:
                    await websocket.send_json({"type": "ping", "data": "ping", "timestamp": None})
                except Exception as e:
                    logger.debug(f"Error sending ping for task {task_id}: {e}")
                    connection_alive = False
                    break
        except asyncio.CancelledError:
            connection_alive = False
        except Exception as e:
            logger.error(f"Unexpected error in heartbeat for task {task_id}: {e}", exc_info=True)
            connection_alive = False

    # Run all tasks concurrently
    log_task = asyncio.create_task(send_logs())
    resource_task = asyncio.create_task(send_resources())
    receive_task = asyncio.create_task(receive_messages())
    heartbeat_task = asyncio.create_task(heartbeat())

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
    except Exception as e:
        logger.error(f"Error in WebSocket stream for task {task_id}: {e}", exc_info=True)
    finally:
        connection_alive = False
        logger.info(f"WebSocket connection closed for task {task_id}")
        try:
            await websocket.close()
        except Exception:
            pass
