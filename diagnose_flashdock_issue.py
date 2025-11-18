"""诊断FLASH-DOCK问题"""
import requests
import subprocess
import sys

def diagnose():
    print("=" * 60)
    print("FLASH-DOCK 问题诊断")
    print("=" * 60)
    print()
    
    # 1. 检查服务是否运行
    print("1. 检查服务状态...")
    try:
        response = requests.get("http://localhost:8501", timeout=3)
        if response.status_code == 200:
            print("   [OK] 服务正在运行")
            print(f"   响应大小: {len(response.text)} 字节")
        else:
            print(f"   [WARNING] 状态码: {response.status_code}")
    except Exception as e:
        print(f"   [ERROR] 无法访问服务: {e}")
        return
    
    # 2. 检查进程
    print("\n2. 检查进程...")
    try:
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", "ps aux | grep streamlit | grep -v grep"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0 and result.stdout.strip():
            print("   [OK] Streamlit进程正在运行")
            lines = result.stdout.strip().splitlines()
            for line in lines[:2]:
                print(f"     {line[:80]}")
        else:
            print("   [WARNING] 未找到Streamlit进程")
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    # 3. 检查端口
    print("\n3. 检查端口...")
    try:
        result = subprocess.run(
            ["netstat", "-ano"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=5
        )
        if ":8501" in result.stdout and "LISTENING" in result.stdout:
            print("   [OK] 端口8501正在监听")
        else:
            print("   [WARNING] 端口8501未监听")
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    # 4. 测试访问
    print("\n4. 测试访问...")
    test_urls = [
        "http://localhost:8501",
        "http://127.0.0.1:8501",
    ]
    
    for url in test_urls:
        try:
            response = requests.get(url, timeout=3)
            print(f"   [OK] {url} - 状态码: {response.status_code}")
        except Exception as e:
            print(f"   [ERROR] {url} - {e}")
    
    # 5. 检查响应内容
    print("\n5. 检查响应内容...")
    try:
        response = requests.get("http://localhost:8501", timeout=3)
        content = response.text.lower()
        
        checks = {
            "streamlit": "streamlit" in content,
            "flashdock": "flashdock" in content or "flash-dock" in content,
            "html": "<html" in content or "<!doctype" in content,
        }
        
        for check, result in checks.items():
            status = "[OK]" if result else "[WARNING]"
            print(f"   {status} 包含 {check}: {result}")
        
        if not any(checks.values()):
            print("   [WARNING] 响应内容可能异常")
            print(f"   前500字符: {response.text[:500]}")
            
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    print("\n" + "=" * 60)
    print("诊断完成")
    print("=" * 60)
    print("\n如果服务可以访问但功能不正常，可能的原因:")
    print("  1. 页面JavaScript错误")
    print("  2. 后端API连接问题")
    print("  3. 依赖缺失")
    print("  4. 配置文件问题")
    print("\n建议:")
    print("  - 打开浏览器开发者工具查看控制台错误")
    print("  - 检查浏览器控制台的Network标签")
    print("  - 查看WSL窗口中的错误信息")

if __name__ == "__main__":
    diagnose()

