"""测试FLASH-DOCK访问"""
import requests
import sys

def test_access():
    print("=" * 60)
    print("测试 FLASH-DOCK 访问")
    print("=" * 60)
    print()
    
    url = "http://localhost:8501"
    
    try:
        print(f"1. 测试访问: {url}")
        response = requests.get(url, timeout=5)
        print(f"   状态码: {response.status_code}")
        print(f"   响应长度: {len(response.text)} 字节")
        
        if response.status_code == 200:
            print("   [OK] 服务可以访问")
            
            # 检查响应内容
            if "streamlit" in response.text.lower() or "flashdock" in response.text.lower():
                print("   [OK] 响应内容正常")
            else:
                print("   [WARNING] 响应内容可能异常")
                print(f"   前200字符: {response.text[:200]}")
        else:
            print(f"   [WARNING] 状态码异常: {response.status_code}")
            
    except requests.exceptions.ConnectionError:
        print("   [ERROR] 无法连接到服务")
        print("   可能原因:")
        print("     - 服务未启动")
        print("     - 端口被占用")
        print("     - 防火墙阻止")
        return False
    except requests.exceptions.Timeout:
        print("   [ERROR] 连接超时")
        return False
    except Exception as e:
        print(f"   [ERROR] 访问失败: {e}")
        return False
    
    # 测试健康检查端点
    print()
    print("2. 测试健康检查端点...")
    health_urls = [
        "http://localhost:8501/_stcore/health",
        "http://localhost:8501/healthz",
    ]
    
    for health_url in health_urls:
        try:
            response = requests.get(health_url, timeout=2)
            if response.status_code == 200:
                print(f"   [OK] {health_url} 正常")
            else:
                print(f"   [INFO] {health_url} 状态码: {response.status_code}")
        except:
            print(f"   [INFO] {health_url} 不可用")
    
    return True

if __name__ == "__main__":
    success = test_access()
    print()
    print("=" * 60)
    if success:
        print("测试完成")
    else:
        print("测试失败")
    print("=" * 60)
    sys.exit(0 if success else 1)

