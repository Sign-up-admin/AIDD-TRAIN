"""
WebSocket connection manager for task log streaming.

This module provides classes and functions to manage WebSocket connections
for real-time log streaming and resource monitoring.
"""

import asyncio
import json
import logging
from typing import Optional, Dict, Any, Tuple
from fastapi import WebSocket, WebSocketDisconnect, status

logger = logging.getLogger(__name__)


class WebSocketConnectionState:
    """Manages WebSocket connection state."""

    def __init__(self, task_id: str, ping_interval: float = 30.0, ping_timeout: float = 60.0):
        """
        Initialize connection state.

        Args:
            task_id: Task ID for this connection
            ping_interval: Interval between pings in seconds
            ping_timeout: Timeout for pong response in seconds
        """
        self.task_id = task_id
        self.connection_alive = True
        self.last_ping_time = asyncio.get_event_loop().time()
        self.ping_interval = ping_interval
        self.ping_timeout = ping_timeout

    def is_alive(self) -> bool:
        """Check if connection is alive."""
        return self.connection_alive

    def mark_dead(self) -> None:
        """Mark connection as dead."""
        self.connection_alive = False

    def update_ping_time(self) -> None:
        """Update last ping time."""
        self.last_ping_time = asyncio.get_event_loop().time()

    def check_timeout(self) -> bool:
        """Check if connection has timed out."""
        current_time = asyncio.get_event_loop().time()
        return (current_time - self.last_ping_time) > self.ping_timeout


async def wait_for_stream_queues(
    stream_manager: Any, task_id: str, websocket: WebSocket, max_wait_time: float = 10.0
) -> Tuple[Optional[Any], Optional[Any]]:
    """
    Wait for stream queues to be available.

    Args:
        stream_manager: Stream manager instance
        task_id: Task ID
        websocket: WebSocket connection
        max_wait_time: Maximum wait time in seconds

    Returns:
        Tuple of (log_queue, resource_queue) or (None, None) if not available
    """
    wait_interval = 0.1
    waited_time = 0.0

    while waited_time < max_wait_time:
        log_queue = stream_manager.get_log_queue(task_id)
        resource_queue = stream_manager.get_resource_queue(task_id)

        if log_queue and resource_queue:
            logger.info(f"Stream queues found for task {task_id} after {waited_time:.2f}s")
            return log_queue, resource_queue

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

    # Try creating streams if not available
    if not log_queue or not resource_queue:
        logger.warning(
            f"Stream queues not found for task {task_id} after {max_wait_time}s, attempting to create"
        )
        stream_manager.create_stream(task_id)

        # Wait for creation to complete (with retry)
        max_retries = 10
        retry_count = 0
        while retry_count < max_retries:
            await asyncio.sleep(0.2)
            log_queue = stream_manager.get_log_queue(task_id)
            resource_queue = stream_manager.get_resource_queue(task_id)
            if log_queue and resource_queue:
                logger.info(
                    f"Stream queues created successfully for task {task_id} after {retry_count * 0.2:.1f}s"
                )
                return log_queue, resource_queue
            retry_count += 1

    return None, None


async def send_log_messages(
    websocket: WebSocket, log_queue: Any, connection_state: WebSocketConnectionState
) -> None:
    """
    Send log messages from queue to WebSocket.

    Args:
        websocket: WebSocket connection
        log_queue: Queue containing log messages
        connection_state: Connection state manager
    """
    task_id = connection_state.task_id
    try:
        while connection_state.is_alive():
            try:
                message = await asyncio.wait_for(log_queue.get(), timeout=0.1)
                if not connection_state.is_alive():
                    break

                await websocket.send_json(message)
            except asyncio.TimeoutError:
                if not connection_state.is_alive():
                    break
                continue
            except Exception as e:
                logger.error(f"Error sending log message for task {task_id}: {e}", exc_info=True)
                connection_state.mark_dead()
                break
    except asyncio.CancelledError:
        connection_state.mark_dead()
    except Exception as e:
        logger.error(
            f"Unexpected error in send_log_messages for task {task_id}: {e}", exc_info=True
        )
        connection_state.mark_dead()


async def send_resource_updates(
    websocket: WebSocket, resource_queue: Any, connection_state: WebSocketConnectionState
) -> None:
    """
    Send resource updates from queue to WebSocket.

    Args:
        websocket: WebSocket connection
        resource_queue: Queue containing resource updates
        connection_state: Connection state manager
    """
    task_id = connection_state.task_id
    try:
        while connection_state.is_alive():
            try:
                message = await asyncio.wait_for(resource_queue.get(), timeout=0.1)
                if not connection_state.is_alive():
                    break
                await websocket.send_json(message)
            except asyncio.TimeoutError:
                if not connection_state.is_alive():
                    break
                continue
            except Exception as e:
                logger.error(
                    f"Error sending resource message for task {task_id}: {e}", exc_info=True
                )
                connection_state.mark_dead()
                break
    except asyncio.CancelledError:
        connection_state.mark_dead()
    except Exception as e:
        logger.error(
            f"Unexpected error in send_resource_updates for task {task_id}: {e}", exc_info=True
        )
        connection_state.mark_dead()


async def receive_client_messages(
    websocket: WebSocket, connection_state: WebSocketConnectionState
) -> None:
    """
    Receive messages from client (ping/pong and commands).

    Args:
        websocket: WebSocket connection
        connection_state: Connection state manager
    """
    task_id = connection_state.task_id
    try:
        while connection_state.is_alive():
            try:
                data = await websocket.receive_text()
                try:
                    message = json.loads(data)
                    msg_type = message.get("type")

                    # Handle ping/pong for connection keepalive
                    if msg_type == "ping":
                        connection_state.update_ping_time()
                        await websocket.send_json(
                            {"type": "pong", "data": "pong", "timestamp": None}
                        )
                    # Handle command messages (for future use)
                    elif msg_type == "command":
                        logger.info(
                            f"Received command from client for task {task_id}: {message.get('data')}"
                        )
                    else:
                        logger.debug(f"Received unknown message type: {msg_type}")
                except json.JSONDecodeError:
                    # Handle plain text pings
                    if data.strip() == "ping":
                        connection_state.update_ping_time()
                        await websocket.send_text("pong")
                    else:
                        logger.warning(
                            f"Invalid JSON received from client for task {task_id}: {data}"
                        )
            except WebSocketDisconnect:
                logger.info(f"WebSocket disconnected by client for task {task_id}")
                connection_state.mark_dead()
                break
            except (json.JSONDecodeError, UnicodeDecodeError) as e:
                logger.warning(f"Invalid message format for task {task_id}: {e}")
                try:
                    await websocket.send_text(json.dumps({"error": "Invalid message format"}))
                except (ConnectionError, OSError):
                    connection_state.mark_dead()
                    break
            except (ConnectionError, OSError, RuntimeError) as e:
                logger.warning(f"Connection error for task {task_id}: {e}")
                connection_state.mark_dead()
                break
            except Exception as e:
                logger.error(
                    f"Unexpected error receiving message for task {task_id}: {e}",
                    exc_info=True,
                    extra={"error_type": type(e).__name__, "task_id": task_id},
                )
                connection_state.mark_dead()
                break
    except asyncio.CancelledError:
        connection_state.mark_dead()
    except (WebSocketDisconnect, ConnectionError, OSError) as e:
        logger.info(f"Connection closed in receive_client_messages for task {task_id}: {e}")
        connection_state.mark_dead()
    except (RuntimeError, ValueError, TypeError) as e:
        logger.error(
            f"Runtime error in WebSocket handler for task {task_id}: {e}",
            exc_info=True,
            extra={"error_type": type(e).__name__, "task_id": task_id},
        )
        connection_state.mark_dead()
    except Exception as e:
        logger.error(
            f"Unexpected error in receive_client_messages for task {task_id}: {e}",
            exc_info=True,
            extra={"task_id": task_id, "error_type": type(e).__name__},
        )
        connection_state.mark_dead()


async def send_heartbeat(websocket: WebSocket, connection_state: WebSocketConnectionState) -> None:
    """
    Send periodic ping to check connection health.

    Args:
        websocket: WebSocket connection
        connection_state: Connection state manager
    """
    task_id = connection_state.task_id
    try:
        while connection_state.is_alive():
            await asyncio.sleep(connection_state.ping_interval)
            if not connection_state.is_alive():
                break

            # Check if we haven't received a pong in too long
            if connection_state.check_timeout():
                logger.warning(f"WebSocket connection timeout for task {task_id}, closing")
                connection_state.mark_dead()
                break

            # Send ping
            try:
                await websocket.send_json({"type": "ping", "data": "ping", "timestamp": None})
            except Exception as e:
                logger.debug(f"Error sending ping for task {task_id}: {e}")
                connection_state.mark_dead()
                break
    except asyncio.CancelledError:
        connection_state.mark_dead()
    except Exception as e:
        logger.error(f"Unexpected error in send_heartbeat for task {task_id}: {e}", exc_info=True)
        connection_state.mark_dead()
