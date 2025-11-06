"""
Performance metrics middleware for COMPASS service.
"""

import time
import logging
from typing import Dict, Optional
from collections import defaultdict, deque
from fastapi import Request
from starlette.middleware.base import BaseHTTPMiddleware
from starlette.responses import Response

logger = logging.getLogger(__name__)


class MetricsCollector:
    """Collects performance metrics."""

    def __init__(self, window_size: int = 100):
        """
        Initialize metrics collector.

        Args:
            window_size: Number of recent requests to keep in sliding window
        """
        self.window_size = window_size

        # Request metrics
        self.request_count = 0
        self.request_times: deque = deque(maxlen=window_size)
        self.endpoint_times: Dict[str, deque] = defaultdict(lambda: deque(maxlen=window_size))
        self.endpoint_counts: Dict[str, int] = defaultdict(int)

        # Error metrics
        self.error_count = 0
        self.error_counts_by_code: Dict[int, int] = defaultdict(int)
        self.error_counts_by_endpoint: Dict[str, int] = defaultdict(int)

        # Status code metrics
        self.status_counts: Dict[int, int] = defaultdict(int)

        # Response size metrics (approximate)
        self.response_sizes: deque = deque(maxlen=window_size)

    def record_request(
        self,
        endpoint: str,
        method: str,
        duration: float,
        status_code: int,
        response_size: Optional[int] = None,
    ):
        """
        Record a request metric.

        Args:
            endpoint: Request endpoint
            method: HTTP method
            status_code: HTTP status code
            duration: Request duration in seconds
            response_size: Response size in bytes (optional)
        """
        self.request_count += 1
        self.request_times.append(duration)

        endpoint_key = f"{method} {endpoint}"
        self.endpoint_times[endpoint_key].append(duration)
        self.endpoint_counts[endpoint_key] += 1
        self.status_counts[status_code] += 1

        if response_size:
            self.response_sizes.append(response_size)

        # Track errors
        if status_code >= 400:
            self.error_count += 1
            self.error_counts_by_code[status_code] += 1
            self.error_counts_by_endpoint[endpoint_key] += 1

    def get_metrics(self) -> Dict:
        """
        Get current metrics summary.

        Returns:
            Dictionary of metrics
        """
        # Calculate statistics
        request_times_list = list(self.request_times)
        avg_response_time = (
            sum(request_times_list) / len(request_times_list) if request_times_list else 0
        )
        min_response_time = min(request_times_list) if request_times_list else 0
        max_response_time = max(request_times_list) if request_times_list else 0

        # Calculate percentiles
        sorted_times = sorted(request_times_list)
        p50 = sorted_times[len(sorted_times) // 2] if sorted_times else 0
        p95 = sorted_times[int(len(sorted_times) * 0.95)] if sorted_times else 0
        p99 = sorted_times[int(len(sorted_times) * 0.99)] if sorted_times else 0

        # Per-endpoint statistics
        endpoint_stats = {}
        for endpoint, times in self.endpoint_times.items():
            times_list = list(times)
            if times_list:
                endpoint_stats[endpoint] = {
                    "count": self.endpoint_counts[endpoint],
                    "avg_time": sum(times_list) / len(times_list),
                    "min_time": min(times_list),
                    "max_time": max(times_list),
                }

        # Response size statistics
        response_sizes_list = list(self.response_sizes)
        avg_response_size = (
            sum(response_sizes_list) / len(response_sizes_list) if response_sizes_list else 0
        )

        return {
            "total_requests": self.request_count,
            "total_errors": self.error_count,
            "error_rate": self.error_count / self.request_count if self.request_count > 0 else 0,
            "response_time": {
                "avg": avg_response_time,
                "min": min_response_time,
                "max": max_response_time,
                "p50": p50,
                "p95": p95,
                "p99": p99,
            },
            "response_size": {
                "avg": avg_response_size,
            },
            "status_codes": dict(self.status_counts),
            "error_counts_by_code": dict(self.error_counts_by_code),
            "endpoints": endpoint_stats,
        }

    def reset(self):
        """Reset all metrics."""
        self.request_count = 0
        self.request_times.clear()
        self.endpoint_times.clear()
        self.endpoint_counts.clear()
        self.error_count = 0
        self.error_counts_by_code.clear()
        self.error_counts_by_endpoint.clear()
        self.status_counts.clear()
        self.response_sizes.clear()


# Global metrics collector
_metrics_collector: Optional[MetricsCollector] = None


def get_metrics_collector() -> MetricsCollector:
    """Get global metrics collector instance."""
    global _metrics_collector
    if _metrics_collector is None:
        _metrics_collector = MetricsCollector()
    return _metrics_collector


class MetricsMiddleware(BaseHTTPMiddleware):
    """Middleware to collect performance metrics."""

    def __init__(self, app, metrics_collector: Optional[MetricsCollector] = None):
        """
        Initialize metrics middleware.

        Args:
            app: FastAPI application
            metrics_collector: Optional metrics collector instance
        """
        super().__init__(app)
        self.metrics = metrics_collector or get_metrics_collector()

    async def dispatch(self, request: Request, call_next):
        """Process request and collect metrics."""
        # Skip metrics for health checks and docs
        if request.url.path in [
            "/health",
            "/health/ready",
            "/docs",
            "/openapi.json",
            "/redoc",
            "/metrics",
        ]:
            return await call_next(request)

        # Start timer
        start_time = time.time()

        # Process request
        response = await call_next(request)

        # Calculate duration
        duration = time.time() - start_time

        # Get response size (approximate from headers)
        response_size = None
        if hasattr(response, "headers"):
            content_length = response.headers.get("content-length")
            if content_length:
                try:
                    response_size = int(content_length)
                except ValueError:
                    pass

        # Record metrics
        self.metrics.record_request(
            endpoint=request.url.path,
            method=request.method,
            duration=duration,
            status_code=response.status_code,
            response_size=response_size,
        )

        # Check alerts
        try:
            from compass.service.monitoring.alert_manager import get_alert_manager

            alert_manager = get_alert_manager()
            metrics_dict = self.metrics.get_metrics()
            alert_manager.check_metrics(metrics_dict)
        except Exception as e:
            logger.debug(f"Failed to check alerts: {e}")

        # Add metrics headers
        response.headers["X-Response-Time"] = f"{duration:.4f}s"

        return response
