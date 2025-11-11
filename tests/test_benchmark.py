"""
Tests for benchmark runner.
"""
import pytest
import time
from compass.service.utils.benchmark import BenchmarkRunner, BenchmarkResult


def test_benchmark_runner_init():
    """Test BenchmarkRunner initialization."""
    runner = BenchmarkRunner()
    assert len(runner.results) == 0


def test_benchmark_runner_simple():
    """Test running a simple benchmark."""
    runner = BenchmarkRunner()
    
    def test_func():
        time.sleep(0.01)
        return "result"
    
    result = runner.run_benchmark(
        "Test Benchmark",
        test_func,
        iterations=5,
        warmup=1
    )
    
    assert result.name == "Test Benchmark"
    assert result.iterations == 5
    assert result.success_count == 5
    assert result.error_count == 0
    assert result.avg_time > 0


def test_benchmark_runner_with_errors():
    """Test benchmark with errors."""
    runner = BenchmarkRunner()
    
    call_count = 0
    def test_func():
        nonlocal call_count
        call_count += 1
        if call_count == 2:
            raise ValueError("Test error")
        return "result"
    
    result = runner.run_benchmark(
        "Test Benchmark with Errors",
        test_func,
        iterations=5,
        warmup=0
    )
    
    assert result.error_count > 0
    assert len(result.errors) > 0


def test_benchmark_runner_get_results():
    """Test getting benchmark results."""
    runner = BenchmarkRunner()
    
    def test_func():
        time.sleep(0.01)
    
    runner.run_benchmark("Test 1", test_func, iterations=3)
    runner.run_benchmark("Test 2", test_func, iterations=3)
    
    results = runner.get_results()
    assert len(results) == 2
    assert results[0]['name'] == "Test 1"
    assert results[1]['name'] == "Test 2"












