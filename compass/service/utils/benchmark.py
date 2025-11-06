"""
Performance benchmarking utilities for COMPASS service.
"""

import time
import logging
import statistics
from typing import List, Dict, Optional, Callable
from dataclasses import dataclass, asdict
from datetime import datetime
import json

logger = logging.getLogger(__name__)


@dataclass
class BenchmarkResult:
    """Benchmark result."""

    name: str
    iterations: int
    total_time: float
    avg_time: float
    min_time: float
    max_time: float
    median_time: float
    p95_time: float
    p99_time: float
    std_dev: float
    success_count: int
    error_count: int
    errors: List[str]
    timestamp: str


class BenchmarkRunner:
    """Runs performance benchmarks."""

    def __init__(self):
        """Initialize benchmark runner."""
        self.results: List[BenchmarkResult] = []

    def run_benchmark(
        self, name: str, func: Callable, iterations: int = 10, warmup: int = 2, *args, **kwargs
    ) -> BenchmarkResult:
        """
        Run a benchmark.

        Args:
            name: Benchmark name
            func: Function to benchmark
            iterations: Number of iterations
            warmup: Number of warmup iterations
            *args: Arguments for function
            **kwargs: Keyword arguments for function

        Returns:
            BenchmarkResult
        """
        logger.info(f"Running benchmark: {name} ({iterations} iterations, {warmup} warmup)")

        # Warmup
        for _ in range(warmup):
            try:
                func(*args, **kwargs)
            except Exception:
                pass

        # Run benchmark
        times = []
        errors = []
        success_count = 0
        error_count = 0

        for i in range(iterations):
            start_time = time.time()
            try:
                func(*args, **kwargs)
                elapsed = time.time() - start_time
                times.append(elapsed)
                success_count += 1
            except Exception as e:
                error_count += 1
                errors.append(str(e))
                elapsed = time.time() - start_time
                times.append(elapsed)  # Include failed attempts in timing

        # Calculate statistics
        if not times:
            raise ValueError("No successful iterations")

        sorted_times = sorted(times)
        total_time = sum(times)
        avg_time = total_time / len(times)
        min_time = min(times)
        max_time = max(times)
        median_time = statistics.median(times)
        p95_time = (
            sorted_times[int(len(sorted_times) * 0.95)] if len(sorted_times) > 1 else times[0]
        )
        p99_time = (
            sorted_times[int(len(sorted_times) * 0.99)] if len(sorted_times) > 1 else times[0]
        )
        std_dev = statistics.stdev(times) if len(times) > 1 else 0

        result = BenchmarkResult(
            name=name,
            iterations=iterations,
            total_time=total_time,
            avg_time=avg_time,
            min_time=min_time,
            max_time=max_time,
            median_time=median_time,
            p95_time=p95_time,
            p99_time=p99_time,
            std_dev=std_dev,
            success_count=success_count,
            error_count=error_count,
            errors=errors[:10],  # Limit to 10 errors
            timestamp=datetime.now().isoformat(),
        )

        self.results.append(result)
        logger.info(f"Benchmark completed: {name} - Avg: {avg_time:.4f}s, P95: {p95_time:.4f}s")

        return result

    def get_results(self) -> List[Dict]:
        """Get all benchmark results as dictionaries."""
        return [asdict(result) for result in self.results]

    def save_results(self, file_path: str):
        """Save benchmark results to file."""
        results = self.get_results()
        with open(file_path, "w") as f:
            json.dump(results, f, indent=2)
        logger.info(f"Benchmark results saved to: {file_path}")

    def print_summary(self):
        """Print benchmark summary."""
        print("\n" + "=" * 80)
        print("BENCHMARK SUMMARY")
        print("=" * 80)

        for result in self.results:
            print(f"\n{result.name}")
            print(
                f"  Iterations: {result.iterations} (Success: {result.success_count}, Errors: {result.error_count})"
            )
            print(f"  Total Time: {result.total_time:.4f}s")
            print(f"  Avg Time: {result.avg_time:.4f}s")
            print(f"  Min Time: {result.min_time:.4f}s")
            print(f"  Max Time: {result.max_time:.4f}s")
            print(f"  Median: {result.median_time:.4f}s")
            print(f"  P95: {result.p95_time:.4f}s")
            print(f"  P99: {result.p99_time:.4f}s")
            print(f"  Std Dev: {result.std_dev:.4f}s")

        print("\n" + "=" * 80)


def benchmark_inference_service(
    inference_service,
    num_iterations: int = 10,
    protein_path: Optional[str] = None,
    ligand_path: Optional[str] = None,
) -> BenchmarkResult:
    """
    Benchmark inference service.

    Args:
        inference_service: InferenceService instance
        num_iterations: Number of iterations
        protein_path: Path to protein file
        ligand_path: Path to ligand file

    Returns:
        BenchmarkResult
    """
    from compass.service.models.model import InferenceRequest

    if not protein_path or not ligand_path:
        # Use dummy paths for testing
        protein_path = "test_protein.pdb"
        ligand_path = "test_ligand.sdf"

    runner = BenchmarkRunner()

    def inference_func():
        request = InferenceRequest(protein_path=protein_path, ligand_path=ligand_path)
        try:
            inference_service.predict(request)
        except Exception:
            # Expected to fail if files don't exist, but we're benchmarking the code path
            pass

    return runner.run_benchmark("Inference Service", inference_func, iterations=num_iterations)
