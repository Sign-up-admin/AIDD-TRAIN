"""
Monitoring and alerting manager for COMPASS service.
"""

import logging
import time
from typing import Dict, List, Optional, Callable
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from threading import Lock

logger = logging.getLogger(__name__)


class AlertLevel(str, Enum):
    """Alert severity levels."""

    INFO = "info"
    WARNING = "warning"
    ERROR = "error"
    CRITICAL = "critical"


@dataclass
class Alert:
    """Alert information."""

    level: AlertLevel
    message: str
    source: str
    timestamp: datetime = field(default_factory=datetime.now)
    metadata: Dict = field(default_factory=dict)


class AlertRule:
    """Alert rule definition."""

    def __init__(
        self,
        name: str,
        condition: Callable[[Dict], bool],
        level: AlertLevel,
        message_template: str,
        cooldown_seconds: int = 300,
    ):
        """
        Initialize alert rule.

        Args:
            name: Rule name
            condition: Function that returns True if alert should trigger
            level: Alert level
            message_template: Message template (can use {metadata} placeholders)
            cooldown_seconds: Cooldown period between alerts (seconds)
        """
        self.name = name
        self.condition = condition
        self.level = level
        self.message_template = message_template
        self.cooldown_seconds = cooldown_seconds
        self.last_triggered: Optional[float] = None

    def check(self, metrics: Dict) -> Optional[Alert]:
        """
        Check if rule should trigger.

        Args:
            metrics: Current metrics dictionary

        Returns:
            Alert if should trigger, None otherwise
        """
        # Check cooldown
        if self.last_triggered:
            elapsed = time.time() - self.last_triggered
            if elapsed < self.cooldown_seconds:
                return None

        # Check condition
        if self.condition(metrics):
            self.last_triggered = time.time()
            message = self.message_template.format(**metrics)
            return Alert(level=self.level, message=message, source=self.name, metadata=metrics)

        return None


class AlertManager:
    """Manages alerts and monitoring rules."""

    def __init__(self):
        """Initialize alert manager."""
        self.rules: List[AlertRule] = []
        self.alerts: List[Alert] = []
        self.alert_handlers: List[Callable[[Alert], None]] = []
        self.lock = Lock()
        self.max_alerts = 1000  # Keep last 1000 alerts

    def add_rule(self, rule: AlertRule):
        """Add an alert rule."""
        with self.lock:
            self.rules.append(rule)

    def add_handler(self, handler: Callable[[Alert], None]):
        """Add an alert handler (e.g., email, webhook, logging)."""
        with self.lock:
            self.alert_handlers.append(handler)

    def check_metrics(self, metrics: Dict):
        """
        Check metrics against all rules.

        Args:
            metrics: Current metrics dictionary
        """
        with self.lock:
            for rule in self.rules:
                alert = rule.check(metrics)
                if alert:
                    self._trigger_alert(alert)

    def _trigger_alert(self, alert: Alert):
        """Trigger an alert."""
        self.alerts.append(alert)

        # Keep only recent alerts
        if len(self.alerts) > self.max_alerts:
            self.alerts = self.alerts[-self.max_alerts :]

        # Log alert
        log_level = {
            AlertLevel.INFO: logging.INFO,
            AlertLevel.WARNING: logging.WARNING,
            AlertLevel.ERROR: logging.ERROR,
            AlertLevel.CRITICAL: logging.CRITICAL,
        }.get(alert.level, logging.INFO)

        logger.log(log_level, f"ALERT [{alert.level.value.upper()}]: {alert.message}")

        # Call handlers
        for handler in self.alert_handlers:
            try:
                handler(alert)
            except Exception as e:
                logger.error(f"Alert handler failed: {e}", exc_info=True)

    def get_alerts(self, level: Optional[AlertLevel] = None, limit: int = 100) -> List[Dict]:
        """
        Get recent alerts.

        Args:
            level: Filter by alert level
            limit: Maximum number of alerts to return

        Returns:
            List of alert dictionaries
        """
        with self.lock:
            alerts = self.alerts
            if level:
                alerts = [a for a in alerts if a.level == level]

            # Sort by timestamp (newest first)
            alerts.sort(key=lambda x: x.timestamp, reverse=True)
            return [
                {
                    "level": a.level.value,
                    "message": a.message,
                    "source": a.source,
                    "timestamp": a.timestamp.isoformat(),
                    "metadata": a.metadata,
                }
                for a in alerts[:limit]
            ]

    def get_alert_summary(self) -> Dict:
        """Get alert summary statistics."""
        with self.lock:
            summary = {"total_alerts": len(self.alerts), "by_level": {}, "recent_alerts": []}

            # Count by level
            for level in AlertLevel:
                count = len([a for a in self.alerts if a.level == level])
                summary["by_level"][level.value] = count

            # Recent alerts (last 10)
            recent = sorted(self.alerts, key=lambda x: x.timestamp, reverse=True)[:10]
            summary["recent_alerts"] = [
                {"level": a.level.value, "message": a.message, "timestamp": a.timestamp.isoformat()}
                for a in recent
            ]

            return summary


# Global alert manager instance
_alert_manager: Optional[AlertManager] = None


def get_alert_manager() -> AlertManager:
    """Get global alert manager instance."""
    global _alert_manager
    if _alert_manager is None:
        _alert_manager = AlertManager()
        _setup_default_rules(_alert_manager)
    return _alert_manager


def _setup_default_rules(manager: AlertManager):
    """Setup default alert rules."""
    # High error rate
    manager.add_rule(
        AlertRule(
            name="high_error_rate",
            condition=lambda m: m.get("error_rate", 0) > 0.1,  # >10% error rate
            level=AlertLevel.WARNING,
            message_template="High error rate detected: {error_rate:.2%}",
            cooldown_seconds=300,
        )
    )

    # Very high error rate
    manager.add_rule(
        AlertRule(
            name="very_high_error_rate",
            condition=lambda m: m.get("error_rate", 0) > 0.5,  # >50% error rate
            level=AlertLevel.CRITICAL,
            message_template="Critical: Very high error rate: {error_rate:.2%}",
            cooldown_seconds=60,
        )
    )

    # Slow response time
    manager.add_rule(
        AlertRule(
            name="slow_response_time",
            condition=lambda m: m.get("response_time", {}).get("p95", 0) > 5.0,  # P95 > 5s
            level=AlertLevel.WARNING,
            message_template="Slow response time detected: P95 = {response_time[p95]:.2f}s",
            cooldown_seconds=300,
        )
    )

    # Very slow response time
    manager.add_rule(
        AlertRule(
            name="very_slow_response_time",
            condition=lambda m: m.get("response_time", {}).get("p95", 0) > 10.0,  # P95 > 10s
            level=AlertLevel.ERROR,
            message_template="Critical: Very slow response time: P95 = {response_time[p95]:.2f}s",
            cooldown_seconds=60,
        )
    )


def setup_file_alert_handler(log_file: str = "alerts.log"):
    """
    Setup file-based alert handler.

    Args:
        log_file: Path to alert log file
    """
    alert_manager = get_alert_manager()

    def file_handler(alert: Alert):
        from pathlib import Path

        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)

        with open(log_path, "a", encoding="utf-8") as f:
            f.write(
                f"[{alert.timestamp.isoformat()}] "
                f"[{alert.level.value.upper()}] "
                f"[{alert.source}] {alert.message}\n"
            )

    alert_manager.add_handler(file_handler)


def setup_webhook_alert_handler(webhook_url: str):
    """
    Setup webhook-based alert handler.

    Args:
        webhook_url: Webhook URL to send alerts to
    """
    alert_manager = get_alert_manager()

    def webhook_handler(alert: Alert):
        import requests

        try:
            requests.post(
                webhook_url,
                json={
                    "level": alert.level.value,
                    "message": alert.message,
                    "source": alert.source,
                    "timestamp": alert.timestamp.isoformat(),
                    "metadata": alert.metadata,
                },
                timeout=5,
            )
        except Exception as e:
            logger.warning(f"Failed to send alert to webhook: {e}")

    alert_manager.add_handler(webhook_handler)
