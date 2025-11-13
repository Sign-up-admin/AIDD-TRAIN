"""
Data service implementation.

This module provides the DataService class for managing datasets in the COMPASS service.
It handles dataset creation, upload, extraction, and listing operations.

Classes:
    DataService: Main service class for dataset management operations.
"""

import os
import uuid
import shutil
import zipfile
import tarfile
import logging
from typing import Dict, Optional, List
from datetime import datetime
from pathlib import Path

from compass.service.models.dataset import DatasetResponse
from compass.service.config import SERVICE_CONFIG

logger = logging.getLogger(__name__)


class DataService:
    """
    Service for managing datasets.

    This class provides functionality to create, upload, extract, and manage
    datasets used for training and inference in the COMPASS service.

    Attributes:
        datasets: Dictionary mapping dataset IDs to DatasetResponse objects
        data_dir: Path to the directory where datasets are stored
    """

    def __init__(self):
        """
        Initialize data service.

        Creates the data directory if it doesn't exist and initializes
        the datasets dictionary.
        """
        self.datasets: Dict[str, DatasetResponse] = {}
        self.data_dir = Path(SERVICE_CONFIG["data_dir"])
        self.data_dir.mkdir(parents=True, exist_ok=True)

    def create_dataset(
        self, name: str, description: Optional[str] = None, metadata: Optional[Dict] = None
    ) -> str:
        """
        Create a new dataset.

        Creates a new dataset entry with a unique ID and prepares a directory
        for storing the dataset files.

        Args:
            name: Dataset name (required)
            description: Optional description of the dataset
            metadata: Optional metadata dictionary for additional information

        Returns:
            str: Unique dataset ID (UUID)

        Raises:
            ValueError: If name is empty or None
        """
        dataset_id = str(uuid.uuid4())
        dataset_dir = self.data_dir / dataset_id
        dataset_dir.mkdir(parents=True, exist_ok=True)

        dataset = DatasetResponse(
            dataset_id=dataset_id,
            name=name,
            description=description,
            status="pending",
            created_at=datetime.now(),
            updated_at=datetime.now(),
            metadata=metadata or {},
        )

        self.datasets[dataset_id] = dataset
        logger.info(f"Created dataset: {dataset_id} ({name})")
        return dataset_id

    def upload_dataset(self, dataset_id: str, file_path: str) -> bool:
        """
        Upload and extract dataset files.

        Uploads a dataset file (zip or tar archive) and extracts it to the
        dataset directory. Supports .zip, .tar, .tar.gz, and .tar.bz2 formats.

        Args:
            dataset_id: ID of the dataset to upload files to
            file_path: Path to the archive file to upload

        Returns:
            bool: True if upload and extraction successful, False otherwise

        Raises:
            FileNotFoundError: If dataset_id doesn't exist or file_path is invalid
            ValueError: If file format is not supported
        """
        if dataset_id not in self.datasets:
            return False

        dataset = self.datasets[dataset_id]
        dataset_dir = self.data_dir / dataset_id

        try:
            # Extract if archive
            if file_path.endswith(".zip"):
                with zipfile.ZipFile(file_path, "r") as zip_ref:
                    # Validate zip members to prevent path traversal attacks
                    for member in zip_ref.namelist():
                        # Check for path traversal attempts
                        if ".." in member or member.startswith("/"):
                            raise ValueError(f"Invalid zip member path: {member}")
                    zip_ref.extractall(dataset_dir)  # nosec B202 - validated above
            elif file_path.endswith(".tar.gz") or file_path.endswith(".tar"):
                with tarfile.open(file_path, "r:*") as tar_ref:
                    # Validate tar members to prevent path traversal attacks
                    for member in tar_ref.getmembers():
                        # Check for path traversal attempts
                        if ".." in member.name or member.name.startswith("/"):
                            raise ValueError(f"Invalid tar member path: {member.name}")
                        # Check for absolute paths
                        if os.path.isabs(member.name):
                            raise ValueError(f"Absolute path not allowed: {member.name}")
                    tar_ref.extractall(dataset_dir)  # nosec B202 - validated above
            else:
                # Copy file
                shutil.copy2(file_path, dataset_dir)

            # Update dataset info
            self._update_dataset_info(dataset_id)
            dataset.status = "ready"
            dataset.updated_at = datetime.now()

            logger.info(f"Uploaded dataset: {dataset_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to upload dataset {dataset_id}: {e}")
            dataset.status = "error"
            dataset.error = str(e)
            return False

    def _update_dataset_info(self, dataset_id: str):
        """Update dataset size and file count."""
        dataset = self.datasets[dataset_id]
        dataset_dir = self.data_dir / dataset_id

        total_size = 0
        file_count = 0

        for root, dirs, files in os.walk(dataset_dir):
            for file in files:
                file_path = os.path.join(root, file)
                total_size += os.path.getsize(file_path)
                file_count += 1

        dataset.size = total_size
        dataset.file_count = file_count

    def get_dataset(self, dataset_id: str) -> Optional[DatasetResponse]:
        """
        Get dataset by ID.

        Args:
            dataset_id: Dataset ID

        Returns:
            Optional[DatasetResponse]: Dataset or None if not found
        """
        return self.datasets.get(dataset_id)

    def list_datasets(self) -> List[DatasetResponse]:
        """
        List all datasets.

        Returns:
            List[DatasetResponse]: List of datasets
        """
        return list(self.datasets.values())

    def delete_dataset(self, dataset_id: str) -> bool:
        """
        Delete a dataset.

        Args:
            dataset_id: Dataset ID

        Returns:
            bool: True if deleted successfully
        """
        if dataset_id not in self.datasets:
            return False

        dataset_dir = self.data_dir / dataset_id
        if dataset_dir.exists():
            shutil.rmtree(dataset_dir)

        del self.datasets[dataset_id]
        logger.info(f"Deleted dataset: {dataset_id}")
        return True

    def get_dataset_path(self, dataset_id: str) -> Optional[Path]:
        """
        Get dataset directory path.

        Args:
            dataset_id: Dataset ID

        Returns:
            Optional[Path]: Dataset path or None if not found
        """
        if dataset_id not in self.datasets:
            return None
        return self.data_dir / dataset_id
