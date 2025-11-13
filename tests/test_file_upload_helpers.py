"""
Tests for file upload helper functions.
"""

import pytest
import tempfile
import os
from pathlib import Path
from unittest.mock import Mock, AsyncMock
from fastapi import UploadFile

from compass.service.utils.file_upload_helpers import (
    validate_file_extension,
    read_and_validate_file_size,
    sanitize_dataset_metadata,
    create_temp_file,
    cleanup_temp_file,
    check_zip_bomb,
)


class TestValidateFileExtension:
    """Test validate_file_extension function."""

    def test_valid_zip_extension(self):
        """Test valid .zip extension."""
        ext = validate_file_extension("test.zip")
        assert ext == ".zip"

    def test_valid_tar_extension(self):
        """Test valid .tar extension."""
        ext = validate_file_extension("test.tar")
        assert ext == ".tar"

    def test_valid_tar_gz_extension(self):
        """Test valid .tar.gz extension."""
        ext = validate_file_extension("test.tar.gz")
        assert ext == ".gz"

    def test_missing_filename(self):
        """Test missing filename raises exception."""
        with pytest.raises(Exception):
            validate_file_extension(None)

    def test_invalid_extension(self):
        """Test invalid extension raises exception."""
        with pytest.raises(Exception):
            validate_file_extension("test.txt")


@pytest.mark.asyncio
class TestReadAndValidateFileSize:
    """Test read_and_validate_file_size function."""

    async def test_read_small_file(self):
        """Test reading a small file."""
        file = AsyncMock(spec=UploadFile)
        file.read = AsyncMock(side_effect=[b"test content", b""])
        file.seek = AsyncMock()

        chunks, size = await read_and_validate_file_size(file, max_size=1024)
        assert size == len(b"test content")
        assert len(chunks) == 2
        file.seek.assert_called_once()

    async def test_file_too_large(self):
        """Test file exceeding max size raises exception."""
        file = AsyncMock(spec=UploadFile)
        large_chunk = b"x" * 1024
        file.read = AsyncMock(side_effect=[large_chunk, large_chunk, b""])
        file.seek = AsyncMock()

        with pytest.raises(Exception):
            await read_and_validate_file_size(file, max_size=1500)


class TestSanitizeDatasetMetadata:
    """Test sanitize_dataset_metadata function."""

    def test_with_name(self):
        """Test with provided name."""
        name, desc = sanitize_dataset_metadata("test_dataset", None, None)
        assert name == "test_dataset"
        assert desc is None

    def test_with_filename(self):
        """Test with filename."""
        name, desc = sanitize_dataset_metadata(None, None, "test_file.zip")
        assert name is not None
        assert desc is None

    def test_with_description(self):
        """Test with description."""
        name, desc = sanitize_dataset_metadata("test", "test description", None)
        assert name == "test"
        assert desc == "test description"


class TestCreateTempFile:
    """Test create_temp_file function."""

    def test_create_temp_file(self):
        """Test creating temporary file."""
        content = [b"test content"]
        file_path = create_temp_file(content, ".txt")
        assert os.path.exists(file_path)
        assert Path(file_path).suffix == ".txt"

        # Cleanup
        cleanup_temp_file(file_path)


class TestCleanupTempFile:
    """Test cleanup_temp_file function."""

    def test_cleanup_existing_file(self):
        """Test cleaning up existing file."""
        # Create a temp file
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
            f.write(b"test")

        # Cleanup
        result = cleanup_temp_file(temp_path)
        assert result is True
        assert not os.path.exists(temp_path)

    def test_cleanup_nonexistent_file(self):
        """Test cleaning up non-existent file."""
        result = cleanup_temp_file("/nonexistent/file.txt")
        assert result is True  # Should return True even if file doesn't exist

    def test_cleanup_none_path(self):
        """Test cleaning up None path."""
        result = cleanup_temp_file(None)
        assert result is True


class TestCheckZipBomb:
    """Test check_zip_bomb function."""

    def test_normal_zip_not_bomb(self):
        """Test normal zip file is not detected as bomb."""
        # Create a simple zip file
        import zipfile
        with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as f:
            zip_path = f.name

        try:
            with zipfile.ZipFile(zip_path, "w") as zf:
                zf.writestr("test.txt", "test content")

            result = check_zip_bomb(zip_path)
            assert result is False
        finally:
            cleanup_temp_file(zip_path)


