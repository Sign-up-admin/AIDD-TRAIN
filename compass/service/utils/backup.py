"""
Data backup utilities for COMPASS service.
"""

import os
import shutil
import logging
import tarfile
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Dict
from dataclasses import dataclass
import json

logger = logging.getLogger(__name__)


@dataclass
class BackupConfig:
    """Backup configuration."""

    enabled: bool = True
    backup_dir: str = "./backups"
    max_backups: int = 10
    backup_interval_hours: int = 24
    backup_on_startup: bool = False
    backup_datasets: bool = True
    backup_models: bool = True
    backup_config: bool = True
    backup_database: bool = True
    compress: bool = True


class BackupManager:
    """Manages data backups for COMPASS service."""

    def __init__(self, config: BackupConfig):
        """
        Initialize backup manager.

        Args:
            config: Backup configuration
        """
        self.config = config
        self.backup_dir = Path(config.backup_dir)
        self.backup_dir.mkdir(parents=True, exist_ok=True)
        self.last_backup_time: Optional[datetime] = None

    def create_backup(
        self,
        datasets_dir: Optional[str] = None,
        models_dir: Optional[str] = None,
        config_file: Optional[str] = None,
        database_file: Optional[str] = None,
        backup_name: Optional[str] = None,
    ) -> Optional[str]:
        """
        Create a backup of service data.

        Args:
            datasets_dir: Path to datasets directory
            models_dir: Path to models/checkpoints directory
            config_file: Path to configuration file
            database_file: Path to database file
            backup_name: Optional backup name (defaults to timestamp)

        Returns:
            Path to created backup file or None if failed
        """
        if not self.config.enabled:
            logger.info("Backup disabled, skipping")
            return None

        try:
            # Generate backup name
            if backup_name is None:
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                backup_name = f"compass_backup_{timestamp}"

            backup_path = self.backup_dir / backup_name
            backup_path.mkdir(parents=True, exist_ok=True)

            # Backup datasets
            if self.config.backup_datasets and datasets_dir:
                self._backup_directory(datasets_dir, backup_path / "datasets")

            # Backup models/checkpoints
            if self.config.backup_models and models_dir:
                self._backup_directory(models_dir, backup_path / "checkpoints")

            # Backup configuration
            if self.config.backup_config and config_file:
                self._backup_file(config_file, backup_path / "config.json")

            # Backup database
            if self.config.backup_database and database_file:
                self._backup_file(database_file, backup_path / Path(database_file).name)

            # Create backup metadata
            metadata = {
                "backup_name": backup_name,
                "created_at": datetime.now().isoformat(),
                "backup_config": {
                    "datasets": self.config.backup_datasets,
                    "models": self.config.backup_models,
                    "config": self.config.backup_config,
                    "database": self.config.backup_database,
                },
            }
            with open(backup_path / "backup_metadata.json", "w") as f:
                json.dump(metadata, f, indent=2)

            # Compress if enabled
            if self.config.compress:
                archive_path = self._compress_backup(backup_path)
                # Remove uncompressed directory
                shutil.rmtree(backup_path)
                backup_path = archive_path
                logger.info(f"Created compressed backup: {backup_path}")
            else:
                logger.info(f"Created backup: {backup_path}")

            # Update last backup time
            self.last_backup_time = datetime.now()

            # Cleanup old backups
            self._cleanup_old_backups()

            return str(backup_path)

        except Exception as e:
            logger.error(f"Failed to create backup: {e}", exc_info=True)
            return None

    def _backup_directory(self, source_dir: str, target_dir: Path):
        """Backup a directory."""
        source = Path(source_dir)
        if not source.exists():
            logger.warning(f"Source directory does not exist: {source_dir}")
            return

        target_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, target_dir, dirs_exist_ok=True)
        logger.debug(f"Backed up directory: {source_dir} -> {target_dir}")

    def _backup_file(self, source_file: str, target_file: Path):
        """Backup a file."""
        source = Path(source_file)
        if not source.exists():
            logger.warning(f"Source file does not exist: {source_file}")
            return

        target_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target_file)
        logger.debug(f"Backed up file: {source_file} -> {target_file}")

    def _compress_backup(self, backup_dir: Path) -> Path:
        """Compress backup directory to tar.gz."""
        archive_name = f"{backup_dir.name}.tar.gz"
        archive_path = backup_dir.parent / archive_name

        with tarfile.open(archive_path, "w:gz") as tar:
            tar.add(backup_dir, arcname=backup_dir.name)

        return archive_path

    def _cleanup_old_backups(self):
        """Remove old backups beyond max_backups limit."""
        try:
            # Get all backup files/directories
            backups = []
            for item in self.backup_dir.iterdir():
                if item.is_file() and (item.suffix == ".gz" or item.suffix == ".tar"):
                    backups.append((item.stat().st_mtime, item))
                elif item.is_dir() and item.name.startswith("compass_backup_"):
                    backups.append((item.stat().st_mtime, item))

            # Sort by modification time (newest first)
            backups.sort(key=lambda x: x[0], reverse=True)

            # Remove old backups
            if len(backups) > self.config.max_backups:
                for _, backup_path in backups[self.config.max_backups :]:
                    try:
                        if backup_path.is_file():
                            backup_path.unlink()
                        else:
                            shutil.rmtree(backup_path)
                        logger.info(f"Removed old backup: {backup_path}")
                    except Exception as e:
                        logger.warning(f"Failed to remove old backup {backup_path}: {e}")

        except Exception as e:
            logger.error(f"Failed to cleanup old backups: {e}", exc_info=True)

    def list_backups(self) -> List[Dict]:
        """
        List all available backups.

        Returns:
            List of backup information dictionaries
        """
        backups = []

        try:
            for item in self.backup_dir.iterdir():
                backup_info = {
                    "name": item.name,
                    "path": str(item),
                    "size": (
                        item.stat().st_size
                        if item.is_file()
                        else sum(f.stat().st_size for f in item.rglob("*") if f.is_file())
                    ),
                    "created_at": datetime.fromtimestamp(item.stat().st_mtime).isoformat(),
                    "is_compressed": item.suffix == ".gz" or item.suffix == ".tar",
                }

                # Try to load metadata
                if item.is_dir():
                    metadata_file = item / "backup_metadata.json"
                    if metadata_file.exists():
                        try:
                            with open(metadata_file, "r") as f:
                                backup_info["metadata"] = json.load(f)
                        except Exception:
                            pass

                backups.append(backup_info)

            # Sort by creation time (newest first)
            backups.sort(key=lambda x: x["created_at"], reverse=True)

        except Exception as e:
            logger.error(f"Failed to list backups: {e}", exc_info=True)

        return backups

    def restore_backup(self, backup_name: str, target_dir: Optional[str] = None) -> bool:
        """
        Restore a backup.

        Args:
            backup_name: Name of backup to restore
            target_dir: Target directory (defaults to original locations)

        Returns:
            True if restoration successful
        """
        backup_path = self.backup_dir / backup_name

        if not backup_path.exists():
            logger.error(f"Backup not found: {backup_name}")
            return False

        try:
            # Extract if compressed
            if backup_path.suffix == ".gz" or backup_path.suffix == ".tar":
                extract_dir = self.backup_dir / backup_path.stem
                with tarfile.open(backup_path, "r:gz") as tar:
                    tar.extractall(extract_dir)
                backup_path = extract_dir / backup_path.stem.rstrip(".tar.gz")

            # Load metadata
            metadata_file = backup_path / "backup_metadata.json"
            metadata = {}
            if metadata_file.exists():
                with open(metadata_file, "r") as f:
                    metadata = json.load(f)

            # Restore components
            if target_dir:
                target = Path(target_dir)
                target.mkdir(parents=True, exist_ok=True)
                shutil.copytree(backup_path, target, dirs_exist_ok=True)
                logger.info(f"Restored backup to: {target_dir}")
            else:
                # Restore to original locations (if metadata available)
                logger.warning("Restoring to original locations requires metadata")
                # TODO: Implement restoration to original locations

            return True

        except Exception as e:
            logger.error(f"Failed to restore backup: {e}", exc_info=True)
            return False


def get_backup_manager() -> Optional[BackupManager]:
    """Get backup manager instance."""
    import os
    from compass.service.config_manager import get_config_manager

    manager = get_config_manager()
    backup_config = BackupConfig(
        enabled=os.getenv("BACKUP_ENABLED", "true").lower() == "true",
        backup_dir=os.getenv("BACKUP_DIR", "./backups"),
        max_backups=int(os.getenv("BACKUP_MAX_COUNT", "10")),
        backup_interval_hours=int(os.getenv("BACKUP_INTERVAL_HOURS", "24")),
        backup_on_startup=os.getenv("BACKUP_ON_STARTUP", "false").lower() == "true",
        backup_datasets=os.getenv("BACKUP_DATASETS", "true").lower() == "true",
        backup_models=os.getenv("BACKUP_MODELS", "true").lower() == "true",
        backup_config=os.getenv("BACKUP_CONFIG", "true").lower() == "true",
        backup_database=os.getenv("BACKUP_DATABASE", "true").lower() == "true",
        compress=os.getenv("BACKUP_COMPRESS", "true").lower() == "true",
    )

    return BackupManager(backup_config) if backup_config.enabled else None










