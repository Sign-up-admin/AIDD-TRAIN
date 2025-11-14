"""
Persistent storage for service registry using SQLite.
"""

import sqlite3
import json
import logging
from contextlib import contextmanager
from pathlib import Path
from typing import Dict, Optional, List
from datetime import datetime

from services.common.service_protocol import ServiceInfo, ServiceStatus

logger = logging.getLogger(__name__)


class ServiceRegistryStorage:
    """Persistent storage for service registry using SQLite."""

    def __init__(self, db_path: str = "registry.db"):
        """
        Initialize service registry storage.

        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()
        logger.info(f"Initialized service registry storage at {self.db_path}")

    def _init_db(self):
        """Initialize database schema."""
        with self._get_connection() as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS services (
                    service_id TEXT PRIMARY KEY,
                    service_name TEXT NOT NULL,
                    host TEXT NOT NULL,
                    port INTEGER NOT NULL,
                    base_url TEXT NOT NULL,
                    status TEXT NOT NULL,
                    metadata TEXT,
                    version TEXT,
                    registered_at TEXT NOT NULL,
                    last_heartbeat TEXT NOT NULL
                )
            """
            )

            # Create index for faster queries
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_service_name
                ON services(service_name)
            """
            )
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_status
                ON services(status)
            """
            )

    @contextmanager
    def _get_connection(self):
        """
        Get database connection with proper transaction handling.

        Uses connection timeout and ensures proper cleanup.
        Optimized for concurrent access with WAL mode and connection pooling settings.
        """
        conn = None
        try:
            # Connect with timeout (default 10 seconds, configurable via env)
            import os

            timeout = float(os.getenv("DB_CONNECTION_TIMEOUT", "10.0"))
            # Ensure timeout is reasonable (between 1 and 60 seconds)
            timeout = max(1.0, min(60.0, timeout))

            conn = sqlite3.connect(self.db_path, timeout=timeout)
            conn.row_factory = sqlite3.Row

            # Enable WAL mode for better concurrency
            conn.execute("PRAGMA journal_mode=WAL")

            # Optimize for concurrent reads
            # Set busy timeout (milliseconds) - allows retries on locked database
            busy_timeout = int(os.getenv("DB_BUSY_TIMEOUT", "5000"))  # 5 seconds default
            conn.execute(f"PRAGMA busy_timeout = {busy_timeout}")

            # Set synchronous mode for better performance (NORMAL is safe with WAL)
            conn.execute("PRAGMA synchronous = NORMAL")

            # Set cache size (in pages, default 2000 pages = ~8MB)
            cache_size = int(os.getenv("DB_CACHE_SIZE", "-2000"))  # Negative = KB
            conn.execute(f"PRAGMA cache_size = {cache_size}")

            yield conn
            conn.commit()
        except sqlite3.OperationalError as e:
            if conn:
                try:
                    conn.rollback()
                except Exception:
                    pass
            logger.error(f"Database operational error: {e}", exc_info=True)
            raise
        except sqlite3.DatabaseError as e:
            if conn:
                try:
                    conn.rollback()
                except Exception:
                    pass
            logger.error(f"Database error: {e}", exc_info=True)
            raise
        except sqlite3.IntegrityError as e:
            if conn:
                try:
                    conn.rollback()
                except Exception:
                    pass
            logger.error(f"Database integrity error: {e}", exc_info=True)
            raise
        except Exception as e:
            if conn:
                try:
                    conn.rollback()
                except Exception:
                    pass
            logger.error(f"Unexpected database error: {e}", exc_info=True)
            raise
        finally:
            if conn:
                try:
                    conn.close()
                except sqlite3.Error as e:
                    logger.warning(f"Error closing database connection: {e}")
                except Exception as e:
                    logger.warning(f"Unexpected error closing database connection: {e}")

    def load_all_services(self) -> Dict[str, ServiceInfo]:
        """
        Load all services from database.

        Returns:
            Dict[str, ServiceInfo]: Dictionary of service_id -> ServiceInfo
        """
        services = {}
        try:
            with self._get_connection() as conn:
                rows = conn.execute("SELECT * FROM services").fetchall()
                for row in rows:
                    try:
                        # Convert Row to dict
                        service_dict = dict(row)

                        # Parse metadata JSON
                        if service_dict.get("metadata"):
                            try:
                                service_dict["metadata"] = json.loads(service_dict["metadata"])
                            except (json.JSONDecodeError, TypeError):
                                service_dict["metadata"] = {}
                        else:
                            service_dict["metadata"] = {}

                        # Parse datetime strings - ensure they are strings for from_dict
                        # from_dict will handle the conversion
                        if service_dict.get("registered_at"):
                            if not isinstance(service_dict["registered_at"], str):
                                # If it's already a datetime, convert to string
                                if isinstance(service_dict["registered_at"], datetime):
                                    service_dict["registered_at"] = service_dict[
                                        "registered_at"
                                    ].isoformat()
                                else:
                                    # Try to convert to string if it's not already
                                    service_dict["registered_at"] = str(service_dict["registered_at"])
                        if service_dict.get("last_heartbeat"):
                            if not isinstance(service_dict["last_heartbeat"], str):
                                if isinstance(service_dict["last_heartbeat"], datetime):
                                    service_dict["last_heartbeat"] = service_dict[
                                        "last_heartbeat"
                                    ].isoformat()
                                else:
                                    # Try to convert to string if it's not already
                                    service_dict["last_heartbeat"] = str(service_dict["last_heartbeat"])

                        service = ServiceInfo.from_dict(service_dict)
                        services[service.service_id] = service
                    except Exception as e:
                        # Use dict access instead of .get() for Row object
                        service_id = dict(row).get("service_id", "unknown")
                        logger.warning(f"Failed to load service {service_id}: {e}", exc_info=True)
                        continue

            logger.info(f"Loaded {len(services)} services from database")
        except Exception as e:
            logger.error(f"Failed to load services from database: {e}", exc_info=True)

        return services

    def register_service(self, service_info: ServiceInfo) -> None:
        """
        Register and persist a service.

        Args:
            service_info: Service information to register
        """
        try:
            with self._get_connection() as conn:
                conn.execute(
                    """
                    INSERT OR REPLACE INTO services
                    (service_id, service_name, host, port, base_url, status,
                     metadata, version, registered_at, last_heartbeat)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                    (
                        service_info.service_id,
                        service_info.service_name,
                        service_info.host,
                        service_info.port,
                        service_info.base_url,
                        service_info.status.value,
                        json.dumps(service_info.metadata),
                        service_info.version,
                        (
                            service_info.registered_at.isoformat()
                            if service_info.registered_at
                            else datetime.now().isoformat()
                        ),
                        (
                            service_info.last_heartbeat.isoformat()
                            if service_info.last_heartbeat
                            else datetime.now().isoformat()
                        ),
                    ),
                )
            logger.debug(f"Persisted service {service_info.service_id} to database")
        except Exception as e:
            logger.error(f"Failed to persist service {service_info.service_id}: {e}", exc_info=True)
            raise

    def update_service_status(self, service_id: str, status: ServiceStatus) -> bool:
        """
        Update service status in database.

        Args:
            service_id: Service ID
            status: New status

        Returns:
            bool: True if updated successfully
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute(
                    "UPDATE services SET status = ? WHERE service_id = ?",
                    (status.value, service_id),
                )
                return bool(cursor.rowcount > 0)
        except Exception as e:
            logger.error(f"Failed to update service status {service_id}: {e}", exc_info=True)
            return False

    def update_heartbeat(self, service_id: str, heartbeat_time: datetime) -> bool:
        """
        Update service heartbeat timestamp.

        Args:
            service_id: Service ID
            heartbeat_time: Heartbeat timestamp

        Returns:
            bool: True if updated successfully
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute(
                    "UPDATE services SET last_heartbeat = ? WHERE service_id = ?",
                    (heartbeat_time.isoformat(), service_id),
                )
                return bool(cursor.rowcount > 0)
        except Exception as e:
            logger.error(f"Failed to update heartbeat {service_id}: {e}", exc_info=True)
            return False

    def delete_service(self, service_id: str) -> bool:
        """
        Delete a service from database.

        Args:
            service_id: Service ID

        Returns:
            bool: True if deleted successfully
        """
        try:
            with self._get_connection() as conn:
                cursor = conn.execute("DELETE FROM services WHERE service_id = ?", (service_id,))
                return bool(cursor.rowcount > 0)
        except Exception as e:
            logger.error(f"Failed to delete service {service_id}: {e}", exc_info=True)
            return False

    def cleanup_stale_services(self, max_age_seconds: int = 300) -> int:
        """
        Clean up services that haven't sent heartbeat for too long.

        Args:
            max_age_seconds: Maximum age in seconds before considering service stale

        Returns:
            int: Number of services cleaned up
        """
        try:
            cutoff_time = datetime.now().timestamp() - max_age_seconds
            cutoff_iso = datetime.fromtimestamp(cutoff_time).isoformat()

            with self._get_connection() as conn:
                cursor = conn.execute(
                    "DELETE FROM services WHERE last_heartbeat < ?", (cutoff_iso,)
                )
                deleted_count = int(cursor.rowcount)
                if deleted_count > 0:
                    logger.info(f"Cleaned up {deleted_count} stale services")
                return deleted_count
        except Exception as e:
            logger.error(f"Failed to cleanup stale services: {e}", exc_info=True)
            return 0
