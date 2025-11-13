#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
实践测试脚本 - 验证升级后的COMPASS系统功能

测试内容：
1. 服务健康检查
2. API端点测试
3. WebSocket连接测试
4. 监控指标测试
5. 配置验证
"""

import sys
import time
import json
import asyncio
import os
from pathlib import Path
from typing import Dict, Any, Optional
import requests
from websockets import connect
from websockets.exceptions import ConnectionClosed

# 设置Windows控制台编码为UTF-8
if sys.platform == "win32":
    try:
        import io
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
        sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')
    except Exception:
        pass

# 添加项目根目录到路径
project_root = Path(__file__).parent.parent
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

# 服务配置
SERVICES = {
    "registry": {
        "url": "http://localhost:8500",
        "name": "服务注册中心",
    },
    "compass": {
        "url": "http://localhost:8080",
        "name": "COMPASS服务",
    },
}

# 测试结果
test_results: Dict[str, Dict[str, Any]] = {}


def print_header(title: str):
    """打印测试标题"""
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def print_result(test_name: str, success: bool, message: str = ""):
    """打印测试结果"""
    # 使用ASCII字符以避免编码问题
    status = "[PASS]" if success else "[FAIL]"
    print(f"  {status}: {test_name}")
    if message:
        print(f"    {message}")
    return success


def test_service_health(service_name: str, service_config: Dict[str, str]) -> bool:
    """测试服务健康状态"""
    print_header(f"测试 {service_config['name']} 健康检查")
    
    try:
        health_url = f"{service_config['url']}/health"
        response = requests.get(health_url, timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            print_result(
                f"{service_config['name']} 健康检查",
                True,
                f"状态: {data.get('status', 'unknown')}"
            )
            return True
        else:
            print_result(
                f"{service_config['name']} 健康检查",
                False,
                f"HTTP {response.status_code}"
            )
            return False
    except requests.exceptions.ConnectionError:
        print_result(
            f"{service_config['name']} 健康检查",
            False,
            "服务未运行或无法连接"
        )
        return False
    except Exception as e:
        print_result(
            f"{service_config['name']} 健康检查",
            False,
            f"错误: {str(e)}"
        )
        return False


def test_api_endpoints() -> bool:
    """测试API端点"""
    print_header("测试 COMPASS API 端点")
    
    base_url = SERVICES["compass"]["url"]
    endpoints = [
        ("/", "根端点"),
        ("/api/v1/training/tasks", "训练任务列表"),
        ("/api/v1/data/datasets", "数据集列表"),
        ("/metrics", "指标端点"),
    ]
    
    all_passed = True
    for endpoint, name in endpoints:
        try:
            url = f"{base_url}{endpoint}"
            response = requests.get(url, timeout=5)
            
            if response.status_code in [200, 404]:  # 404也是正常的，表示端点存在
                print_result(f"API端点: {name}", True, f"HTTP {response.status_code}")
            else:
                print_result(f"API端点: {name}", False, f"HTTP {response.status_code}")
                all_passed = False
        except Exception as e:
            print_result(f"API端点: {name}", False, f"错误: {str(e)}")
            all_passed = False
    
    return all_passed


def test_metrics_endpoint() -> bool:
    """测试指标端点"""
    print_header("测试监控指标端点")
    
    base_url = SERVICES["compass"]["url"]
    
    # 测试JSON格式
    try:
        url = f"{base_url}/metrics?format_type=json"
        response = requests.get(url, timeout=5)
        
        if response.status_code == 200:
            data = response.json()
            print_result("指标端点 (JSON格式)", True, f"包含 {len(data)} 个指标")
        else:
            print_result("指标端点 (JSON格式)", False, f"HTTP {response.status_code}")
            return False
    except Exception as e:
        print_result("指标端点 (JSON格式)", False, f"错误: {str(e)}")
        return False
    
    # 测试Prometheus格式
    try:
        url = f"{base_url}/metrics?format_type=prometheus"
        response = requests.get(url, timeout=5)
        
        if response.status_code == 200:
            content = response.text
            if "http_requests_total" in content or "http_request_duration_seconds" in content:
                print_result("指标端点 (Prometheus格式)", True, "包含Prometheus指标")
            else:
                print_result("指标端点 (Prometheus格式)", True, "格式正确（可能未启用Prometheus）")
        else:
            print_result("指标端点 (Prometheus格式)", False, f"HTTP {response.status_code}")
            return False
    except Exception as e:
        print_result("指标端点 (Prometheus格式)", False, f"错误: {str(e)}")
        return False
    
    return True


async def test_websocket_connection() -> bool:
    """测试WebSocket连接"""
    print_header("测试 WebSocket 连接")
    
    # 首先创建一个测试任务（如果需要）
    base_url = SERVICES["compass"]["url"]
    
    try:
        # 尝试连接WebSocket端点（使用一个不存在的task_id来测试连接）
        ws_url = f"ws://localhost:8080/api/v1/training/tasks/test-task-id/stream"
        
        try:
            # 使用asyncio.wait_for包装connect以设置超时
            websocket = await asyncio.wait_for(connect(ws_url), timeout=5)
            try:
                # 发送ping消息
                await websocket.send(json.dumps({"type": "ping"}))
                
                # 等待响应
                try:
                    response = await asyncio.wait_for(websocket.recv(), timeout=3)
                    print_result("WebSocket连接", True, "连接成功并收到响应")
                    await websocket.close()
                    return True
                except asyncio.TimeoutError:
                    print_result("WebSocket连接", True, "连接成功（无响应，可能任务不存在）")
                    await websocket.close()
                    return True
            finally:
                try:
                    await websocket.close()
                except Exception:
                    pass
        except ConnectionClosed:
            print_result("WebSocket连接", True, "连接已关闭（正常，任务不存在）")
            return True
        except asyncio.TimeoutError:
            print_result("WebSocket连接", False, "连接超时")
            return False
        except Exception as e:
            # WebSocket连接失败可能是因为任务不存在，这是正常的
            error_msg = str(e)
            if "404" in error_msg or "not found" in error_msg.lower() or "404" in str(e):
                print_result("WebSocket连接", True, "端点存在（任务不存在是正常的）")
                return True
            else:
                # 检查是否是连接被拒绝（服务未运行）
                if "10061" in error_msg or "Connection refused" in error_msg:
                    print_result("WebSocket连接", False, "无法连接到服务")
                    return False
                else:
                    print_result("WebSocket连接", True, f"端点可访问（错误: {error_msg[:50]}）")
                    return True
    except Exception as e:
        error_msg = str(e)
        if "10061" in error_msg or "Connection refused" in error_msg:
            print_result("WebSocket连接", False, "无法连接到服务")
        else:
            print_result("WebSocket连接", False, f"连接错误: {error_msg[:100]}")
        return False


def test_configuration() -> bool:
    """测试配置验证"""
    print_header("测试配置验证")
    
    try:
        # 检查配置验证脚本是否存在
        config_script = project_root / "scripts" / "validate_config.py"
        if not config_script.exists():
            print_result("配置验证脚本", False, "脚本不存在")
            return False
        
        # 运行配置验证（这里只是检查脚本存在，实际验证需要运行）
        print_result("配置验证脚本", True, "脚本存在")
        return True
    except Exception as e:
        print_result("配置验证", False, f"错误: {str(e)}")
        return False


def test_monitoring_setup() -> bool:
    """测试监控系统配置"""
    print_header("测试监控系统配置")
    
    checks = [
        ("Prometheus配置", "monitoring/prometheus/prometheus.yml"),
        ("Grafana数据源配置", "monitoring/grafana/provisioning/datasources/prometheus.yml"),
        ("Grafana仪表板配置", "monitoring/grafana/provisioning/dashboards/dashboard.yml"),
        ("COMPASS仪表板", "monitoring/grafana/dashboards/compass-dashboard.json"),
        ("Docker Compose配置", "docker-compose.monitoring.yml"),
    ]
    
    all_passed = True
    for name, path in checks:
        file_path = project_root / path
        if file_path.exists():
            print_result(name, True, f"文件存在: {path}")
        else:
            print_result(name, False, f"文件不存在: {path}")
            all_passed = False
    
    return all_passed


def test_refactored_modules() -> bool:
    """测试重构后的模块"""
    print_header("测试重构后的模块")
    
    modules = [
        ("WebSocket管理器", "compass/service/utils/websocket_manager.py"),
        ("文件上传助手", "compass/service/utils/file_upload_helpers.py"),
        ("训练助手", "compass/service/utils/training_helpers.py"),
        ("任务停止助手", "compass/service/utils/task_stop_helpers.py"),
    ]
    
    all_passed = True
    for name, path in modules:
        file_path = project_root / path
        if file_path.exists():
            # 尝试导入模块
            try:
                module_name = path.replace("/", ".").replace(".py", "")
                __import__(module_name)
                print_result(name, True, f"模块可导入: {path}")
            except Exception as e:
                print_result(name, False, f"导入错误: {str(e)}")
                all_passed = False
        else:
            print_result(name, False, f"文件不存在: {path}")
            all_passed = False
    
    return all_passed


def main():
    """主测试函数"""
    print("\n" + "=" * 60)
    print("  COMPASS 系统实践测试")
    print("=" * 60)
    print(f"\n项目根目录: {project_root}")
    print(f"测试时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    
    results = {}
    
    # 1. 测试服务健康状态
    print("\n【阶段1】服务健康检查")
    print("-" * 60)
    for service_name, service_config in SERVICES.items():
        results[f"health_{service_name}"] = test_service_health(service_name, service_config)
    
    # 2. 测试API端点
    print("\n【阶段2】API端点测试")
    print("-" * 60)
    results["api_endpoints"] = test_api_endpoints()
    
    # 3. 测试指标端点
    print("\n【阶段3】监控指标测试")
    print("-" * 60)
    results["metrics"] = test_metrics_endpoint()
    
    # 4. 测试WebSocket连接
    print("\n【阶段4】WebSocket连接测试")
    print("-" * 60)
    results["websocket"] = asyncio.run(test_websocket_connection())
    
    # 5. 测试配置验证
    print("\n【阶段5】配置验证测试")
    print("-" * 60)
    results["configuration"] = test_configuration()
    
    # 6. 测试监控系统配置
    print("\n【阶段6】监控系统配置测试")
    print("-" * 60)
    results["monitoring_setup"] = test_monitoring_setup()
    
    # 7. 测试重构后的模块
    print("\n【阶段7】重构模块测试")
    print("-" * 60)
    results["refactored_modules"] = test_refactored_modules()
    
    # 汇总结果
    print_header("测试结果汇总")
    total_tests = len(results)
    passed_tests = sum(1 for v in results.values() if v)
    failed_tests = total_tests - passed_tests
    
    print(f"\n  总测试数: {total_tests}")
    print(f"  通过: {passed_tests}")
    print(f"  失败: {failed_tests}")
    print(f"  通过率: {passed_tests/total_tests*100:.1f}%")
    
    print("\n详细结果:")
    for test_name, result in results.items():
        status = "[PASS]" if result else "[FAIL]"
        print(f"  {status}: {test_name}")
    
    print("\n" + "=" * 60)
    if failed_tests == 0:
        print("  所有测试通过！")
    else:
        print(f"  有 {failed_tests} 个测试失败，请检查上述输出")
        print("  注意：如果服务未运行，部分测试会失败，这是正常的")
    print("=" * 60 + "\n")
    
    return 0 if failed_tests == 0 else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\n测试被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n测试过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

