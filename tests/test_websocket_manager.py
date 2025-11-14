"""
Tests for WebSocket connection manager.
"""

import pytest
import asyncio
from unittest.mock import Mock, AsyncMock, patch
from fastapi import WebSocket

from compass.service.utils.websocket_manager import (
    WebSocketConnectionState,
    wait_for_stream_queues,
    send_log_messages,
    send_resource_updates,
    receive_client_messages,
    send_heartbeat,
)


@pytest.fixture
def mock_websocket():
    """Create a mock WebSocket."""
    websocket = AsyncMock(spec=WebSocket)
    websocket.accept = AsyncMock()
    websocket.send_json = AsyncMock()
    websocket.send_text = AsyncMock()
    websocket.receive_text = AsyncMock()
    websocket.close = AsyncMock()
    return websocket


@pytest.fixture
def mock_stream_manager():
    """Create a mock stream manager."""
    manager = Mock()
    log_queue = asyncio.Queue()
    resource_queue = asyncio.Queue()
    manager.get_log_queue = Mock(return_value=log_queue)
    manager.get_resource_queue = Mock(return_value=resource_queue)
    manager.create_stream = Mock()
    return manager


class TestWebSocketConnectionState:
    """Test WebSocketConnectionState class."""

    def test_initialization(self):
        """Test connection state initialization."""
        state = WebSocketConnectionState("test_task", ping_interval=30.0, ping_timeout=60.0)
        assert state.task_id == "test_task"
        assert state.connection_alive is True
        assert state.ping_interval == 30.0
        assert state.ping_timeout == 60.0

    def test_is_alive(self):
        """Test is_alive method."""
        state = WebSocketConnectionState("test_task")
        assert state.is_alive() is True
        state.mark_dead()
        assert state.is_alive() is False

    def test_mark_dead(self):
        """Test mark_dead method."""
        state = WebSocketConnectionState("test_task")
        state.mark_dead()
        assert state.connection_alive is False

    def test_update_ping_time(self):
        """Test update_ping_time method."""
        state = WebSocketConnectionState("test_task")
        initial_time = state.last_ping_time
        state.update_ping_time()
        assert state.last_ping_time >= initial_time

    def test_check_timeout(self):
        """Test check_timeout method."""
        state = WebSocketConnectionState("test_task", ping_timeout=0.1)
        assert state.check_timeout() is False
        # Simulate timeout by setting last_ping_time far in the past
        import time
        state.last_ping_time = time.time() - 1.0
        assert state.check_timeout() is True


@pytest.mark.asyncio
class TestWaitForStreamQueues:
    """Test wait_for_stream_queues function."""

    async def test_queues_available_immediately(self, mock_stream_manager, mock_websocket):
        """Test when queues are available immediately."""
        log_queue, resource_queue = await wait_for_stream_queues(
            mock_stream_manager, "test_task", mock_websocket, max_wait_time=1.0
        )
        assert log_queue is not None
        assert resource_queue is not None

    async def test_queues_created_after_wait(self, mock_stream_manager, mock_websocket):
        """Test when queues are created after waiting."""
        # Initially return None
        mock_stream_manager.get_log_queue = Mock(side_effect=[None, None, asyncio.Queue()])
        mock_stream_manager.get_resource_queue = Mock(side_effect=[None, None, asyncio.Queue()])

        log_queue, resource_queue = await wait_for_stream_queues(
            mock_stream_manager, "test_task", mock_websocket, max_wait_time=1.0
        )
        # Should eventually get queues
        assert log_queue is not None or resource_queue is not None


@pytest.mark.asyncio
class TestSendLogMessages:
    """Test send_log_messages function."""

    async def test_send_logs(self, mock_websocket):
        """Test sending log messages."""
        log_queue = asyncio.Queue()
        connection_state = WebSocketConnectionState("test_task")

        # Add a message to queue
        await log_queue.put({"type": "log", "data": "test message"})

        # Start sending (will run until connection dies)
        task = asyncio.create_task(
            send_log_messages(mock_websocket, log_queue, connection_state)
        )

        # Wait a bit for message to be sent
        await asyncio.sleep(0.1)

        # Mark connection as dead to stop the loop
        connection_state.mark_dead()

        # Wait for task to complete
        await asyncio.sleep(0.1)
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass

        # Verify message was sent
        mock_websocket.send_json.assert_called()


@pytest.mark.asyncio
class TestSendResourceUpdates:
    """Test send_resource_updates function."""

    async def test_send_resources(self, mock_websocket):
        """Test sending resource updates."""
        resource_queue = asyncio.Queue()
        connection_state = WebSocketConnectionState("test_task")

        # Add a resource update to queue
        await resource_queue.put({"type": "resources", "data": {"cpu": 50.0}})

        # Start sending
        task = asyncio.create_task(
            send_resource_updates(mock_websocket, resource_queue, connection_state)
        )

        # Wait a bit
        await asyncio.sleep(0.1)

        # Mark connection as dead
        connection_state.mark_dead()

        # Wait for task to complete
        await asyncio.sleep(0.1)
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass

        # Verify resource update was sent
        mock_websocket.send_json.assert_called()


@pytest.mark.asyncio
class TestReceiveClientMessages:
    """Test receive_client_messages function."""

    async def test_receive_ping(self, mock_websocket):
        """Test receiving ping message."""
        connection_state = WebSocketConnectionState("test_task")
        mock_websocket.receive_text = AsyncMock(return_value='{"type": "ping", "data": "ping"}')

        task = asyncio.create_task(receive_client_messages(mock_websocket, connection_state))

        # Wait a bit
        await asyncio.sleep(0.1)

        # Mark connection as dead
        connection_state.mark_dead()

        # Wait for task to complete
        await asyncio.sleep(0.1)
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass

        # Verify pong was sent
        mock_websocket.send_json.assert_called()


@pytest.mark.asyncio
class TestSendHeartbeat:
    """Test send_heartbeat function."""

    async def test_send_heartbeat(self, mock_websocket):
        """Test sending heartbeat."""
        connection_state = WebSocketConnectionState("test_task", ping_interval=0.1)

        task = asyncio.create_task(send_heartbeat(mock_websocket, connection_state))

        # Wait for at least one ping
        await asyncio.sleep(0.2)

        # Mark connection as dead
        connection_state.mark_dead()

        # Wait for task to complete
        await asyncio.sleep(0.1)
        task.cancel()
        try:
            await task
        except asyncio.CancelledError:
            pass

        # Verify ping was sent
        assert mock_websocket.send_json.called




