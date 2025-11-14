"""
Tests for training helper classes and functions.
"""

import pytest
import sys
import io
from unittest.mock import Mock

from compass.service.utils.training_helpers import (
    TeeOutput,
    OutputRedirector,
    ResourceMonitor,
    prepare_training_config,
)


class TestTeeOutput:
    """Test TeeOutput class."""

    def test_write(self):
        """Test writing to TeeOutput."""
        stream = io.StringIO()
        original_stream = sys.stdout
        log_func = Mock()
        stream_manager = Mock()

        tee = TeeOutput(stream, "test_task", original_stream, log_func, stream_manager)
        tee.write("test message\n")

        # Check that log_func was called
        log_func.assert_called()

    def test_flush(self):
        """Test flushing TeeOutput."""
        stream = io.StringIO()
        original_stream = sys.stdout
        log_func = Mock()
        stream_manager = Mock()

        tee = TeeOutput(stream, "test_task", original_stream, log_func, stream_manager)
        tee.buffer = "test"
        tee.flush()

        # Check that log_func was called with buffered content
        log_func.assert_called()

    def test_isatty(self):
        """Test isatty method."""
        stream = io.StringIO()
        original_stream = sys.stdout
        log_func = Mock()
        stream_manager = Mock()

        tee = TeeOutput(stream, "test_task", original_stream, log_func, stream_manager)
        # Should return the original stream's isatty value
        assert isinstance(tee.isatty(), bool)


class TestOutputRedirector:
    """Test OutputRedirector class."""

    def test_setup_and_restore(self):
        """Test setting up and restoring output redirection."""
        original_stdout = sys.stdout
        original_stderr = sys.stderr

        redirector = OutputRedirector("test_task", Mock(), Mock())
        redirector.setup()

        # Check that stdout/stderr were replaced
        assert sys.stdout != original_stdout
        assert sys.stderr != original_stderr

        redirector.restore()

        # Check that stdout/stderr were restored
        assert sys.stdout == original_stdout
        assert sys.stderr == original_stderr


class TestResourceMonitor:
    """Test ResourceMonitor class."""

    def test_initialization(self):
        """Test ResourceMonitor initialization."""
        monitor = ResourceMonitor(
            "test_task",
            {},
            {},
            Mock(),
            None,
            Mock(),
            Mock(return_value={"cpu": 50.0}),
        )
        assert monitor.task_id == "test_task"
        assert monitor.thread is None

    def test_start(self):
        """Test starting resource monitor."""
        monitor = ResourceMonitor(
            "test_task",
            {},
            {},
            Mock(),
            None,
            Mock(),
            Mock(return_value={"cpu": 50.0}),
        )
        monitor.start()
        assert monitor.thread is not None
        assert monitor.thread.is_alive() is True

        # Cleanup
        monitor.thread.join(timeout=1.0)


class TestPrepareTrainingConfig:
    """Test prepare_training_config function."""

    def test_prepare_config(self):
        """Test preparing training configuration."""
        config = {"execution_mode": "validation_tuned", "epochs": 100}
        log_dir = "/test/logs"
        checkpoint_dir = "/test/checkpoints"

        training_config = prepare_training_config(config, log_dir, checkpoint_dir)

        assert training_config["log_dir"] == log_dir
        assert training_config["checkpoint_dir"] == checkpoint_dir
        assert training_config["epochs"] == 100




