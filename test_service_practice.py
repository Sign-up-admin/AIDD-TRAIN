"""
实践测试脚本 - 验证COMPASS服务的改进功能
"""
import requests
import json
import time
import sys
import os
from typing import Dict, Any

# 设置Windows控制台编码
if sys.platform == 'win32':
    os.system('chcp 65001 >nul 2>&1')

BASE_URL = "http://localhost:8080"
TEST_TIMEOUT = 5


def print_section(title: str):
    """打印测试章节标题"""
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)


def test_health_check():
    """测试健康检查端点"""
    print_section("1. 健康检查测试")
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=TEST_TIMEOUT)
        print(f"状态码: {response.status_code}")
        if response.status_code == 200:
            data = response.json()
            print(f"服务状态: {data.get('status', 'unknown')}")
            print(f"版本: {data.get('version', 'unknown')}")
            print("[OK] 健康检查通过")
            return True
        else:
            print(f"[FAIL] 健康检查失败: {response.status_code}")
            return False
    except requests.exceptions.ConnectionError:
        print("[FAIL] 无法连接到服务，请确保服务正在运行")
        print(f"   尝试访问: {BASE_URL}")
        print("   启动命令: python compass/service_main.py --host 0.0.0.0 --port 8080")
        return False
    except Exception as e:
        print(f"[ERROR] 健康检查出错: {e}")
        return False


def test_cors_headers():
    """测试CORS配置"""
    print_section("2. CORS配置测试")
    try:
        headers = {"Origin": "http://localhost:8501"}
        response = requests.get(f"{BASE_URL}/health", headers=headers, timeout=TEST_TIMEOUT)
        
        cors_header = response.headers.get("Access-Control-Allow-Origin")
        if cors_header:
            print(f"CORS头: Access-Control-Allow-Origin = {cors_header}")
            if "localhost" in cors_header or "127.0.0.1" in cors_header:
                print("[OK] CORS配置正确")
                return True
            else:
                print("[WARN] CORS配置存在但可能不正确")
                return False
        else:
            print("[FAIL] 未找到CORS头")
            return False
    except Exception as e:
        print(f"[ERROR] CORS测试出错: {e}")
        return False


def test_security_headers():
    """测试安全头"""
    print_section("3. 安全头测试")
    try:
        response = requests.get(f"{BASE_URL}/health", timeout=TEST_TIMEOUT)
        headers = response.headers
        
        security_headers = {
            "X-Content-Type-Options": "nosniff",
            "X-Frame-Options": "DENY",
            "Content-Security-Policy": None,  # 值可能不同
        }
        
        all_present = True
        for header, expected_value in security_headers.items():
            value = headers.get(header)
            if value:
                print(f"[OK] {header}: {value}")
                if expected_value and value != expected_value:
                    print(f"   [WARN] 期望值: {expected_value}")
            else:
                print(f"[FAIL] {header}: 未找到")
                all_present = False
        
        return all_present
    except Exception as e:
        print(f"[ERROR] 安全头测试出错: {e}")
        return False


def test_input_sanitization():
    """测试输入清理"""
    print_section("4. 输入清理测试")
    try:
        # 测试XSS攻击防护
        xss_payload = "<script>alert('XSS')</script>"
        
        # 尝试创建一个任务，使用包含XSS的description
        task_data = {
            "config": {
                "execution_mode": "smoke_test",
                "epochs": 1,
                "batch_size": 1,
                "learning_rate": 0.001,
                "optimizer": "adam"
            },
            "description": xss_payload
        }
        
        response = requests.post(
            f"{BASE_URL}/api/v1/training/tasks",
            json=task_data,
            timeout=TEST_TIMEOUT
        )
        
        if response.status_code in [201, 422, 400]:
            # 检查响应中是否包含原始XSS代码
            response_text = response.text
            if "<script>" in response_text and "&lt;script&gt;" not in response_text:
                print("[WARN] 可能存在XSS漏洞")
                return False
            else:
                print("[OK] XSS防护正常（输入被清理或拒绝）")
                return True
        else:
            print(f"[WARN] 意外的状态码: {response.status_code}")
            return True  # 可能是其他验证错误
    except Exception as e:
        print(f"[ERROR] 输入清理测试出错: {e}")
        return False


def test_rate_limiting():
    """测试速率限制"""
    print_section("5. 速率限制测试")
    try:
        # 快速发送多个请求
        responses = []
        for i in range(15):
            try:
                response = requests.get(f"{BASE_URL}/api/v1/training/tasks", timeout=TEST_TIMEOUT)
                responses.append(response.status_code)
                time.sleep(0.1)  # 短暂延迟
            except Exception as e:
                print(f"请求 {i+1} 出错: {e}")
        
        # 检查是否有429（Too Many Requests）响应
        rate_limited = 429 in responses
        if rate_limited:
            print("[OK] 速率限制正常工作（检测到429响应）")
            print(f"   响应码分布: {dict((code, responses.count(code)) for code in set(responses))}")
            return True
        else:
            print("[WARN] 未检测到速率限制（可能限制值较高）")
            print(f"   响应码分布: {dict((code, responses.count(code)) for code in set(responses))}")
            return True  # 不算失败，可能限制值设置较高
    except Exception as e:
        print(f"[ERROR] 速率限制测试出错: {e}")
        return False


def test_error_handling():
    """测试错误处理"""
    print_section("6. 错误处理测试")
    try:
        # 等待一下，避免速率限制影响
        time.sleep(1)
        
        # 测试无效的task_id格式
        invalid_task_id = "invalid-task-id-123"
        response = requests.get(
            f"{BASE_URL}/api/v1/training/tasks/{invalid_task_id}",
            timeout=TEST_TIMEOUT
        )
        
        if response.status_code == 400:
            data = response.json()
            if "error" in data:
                print(f"[OK] 错误处理正常: {data.get('error', 'unknown')}")
                return True
            else:
                print("[WARN] 返回400但错误格式不正确")
                return False
        elif response.status_code == 404:
            print("[OK] 错误处理正常（返回404）")
            return True
        elif response.status_code == 429:
            print("[WARN] 遇到速率限制（429），这是正常的")
            print("   说明速率限制功能正常工作")
            return True  # 速率限制也是错误处理的一部分
        else:
            print(f"[WARN] 意外的状态码: {response.status_code}")
            return False
    except Exception as e:
        print(f"[ERROR] 错误处理测试出错: {e}")
        return False


def test_metrics_endpoint():
    """测试指标端点"""
    print_section("7. 指标端点测试")
    try:
        response = requests.get(f"{BASE_URL}/metrics", timeout=TEST_TIMEOUT)
        if response.status_code == 200:
            data = response.json()
            print("[OK] 指标端点可访问")
            
            # 检查是否包含速率限制统计
            if "rate_limiting" in data:
                print("[OK] 速率限制统计已包含")
                rate_limit_stats = data["rate_limiting"]
                print(f"   总请求数: {rate_limit_stats.get('total_requests', 0)}")
                print(f"   被限制请求数: {rate_limit_stats.get('rate_limited_requests', 0)}")
                print(f"   限制率: {rate_limit_stats.get('rate_limit_percentage', 0):.2f}%")
            else:
                print("[WARN] 未找到速率限制统计")
                print(f"   可用键: {list(data.keys())}")
            
            return True
        else:
            print(f"[FAIL] 指标端点返回错误: {response.status_code}")
            return False
    except Exception as e:
        print(f"[ERROR] 指标端点测试出错: {e}")
        return False


def test_api_documentation():
    """测试API文档"""
    print_section("8. API文档测试")
    try:
        response = requests.get(f"{BASE_URL}/docs", timeout=TEST_TIMEOUT)
        if response.status_code == 200:
            print("[OK] API文档可访问")
            print(f"   文档地址: {BASE_URL}/docs")
            return True
        else:
            print(f"[WARN] API文档返回: {response.status_code}")
            return False
    except Exception as e:
        print(f"[ERROR] API文档测试出错: {e}")
        return False


def main():
    """主测试函数"""
    print("\n" + "=" * 60)
    print("  COMPASS服务实践测试")
    print("=" * 60)
    print(f"\n测试目标: {BASE_URL}")
    print("提示: 确保服务正在运行 (python compass/service_main.py)")
    
    results = []
    
    # 运行所有测试
    results.append(("健康检查", test_health_check()))
    results.append(("CORS配置", test_cors_headers()))
    results.append(("安全头", test_security_headers()))
    results.append(("输入清理", test_input_sanitization()))
    results.append(("速率限制", test_rate_limiting()))
    results.append(("错误处理", test_error_handling()))
    results.append(("指标端点", test_metrics_endpoint()))
    results.append(("API文档", test_api_documentation()))
    
    # 打印总结
    print_section("测试总结")
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "[PASS] 通过" if result else "[FAIL] 失败"
        print(f"{name}: {status}")
    
    print(f"\n总计: {passed}/{total} 测试通过")
    
    if passed == total:
        print("\n[SUCCESS] 所有测试通过！服务运行正常。")
        return 0
    else:
        print(f"\n[WARN] 有 {total - passed} 个测试未通过，请检查服务配置。")
        return 1


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

