"""
一键实践测试脚本
"""
import subprocess
import sys
import time
import requests
import os

def check_service():
    """检查服务是否运行"""
    try:
        response = requests.get("http://localhost:8080/health", timeout=2)
        return response.status_code == 200
    except:
        return False

def start_service():
    """启动服务"""
    print("=" * 60)
    print("正在启动COMPASS服务...")
    print("=" * 60)
    
    # 设置环境变量
    env = os.environ.copy()
    project_root = os.getcwd()
    env["PYTHONPATH"] = project_root
    
    # 启动服务
    cmd = [
        sys.executable,
        "compass/service_main.py",
        "--host", "127.0.0.1",
        "--port", "8080",
        "--registry-url", "http://localhost:8500"
    ]
    
    print(f"启动命令: {' '.join(cmd)}")
    print("\n提示: 服务将在后台启动")
    print("按 Ctrl+C 可以停止服务")
    print("=" * 60)
    
    # 启动服务进程
    process = subprocess.Popen(
        cmd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd=project_root
    )
    
    # 等待服务启动
    print("\n等待服务启动...")
    for i in range(10):
        time.sleep(1)
        if check_service():
            print(f"\n[OK] 服务已启动！(等待了 {i+1} 秒)")
            return process
        print(f".", end="", flush=True)
    
    print("\n[WARN] 服务启动可能较慢，继续测试...")
    return process

def main():
    """主函数"""
    print("\n" + "=" * 60)
    print("  COMPASS服务实践测试")
    print("=" * 60)
    
    # 检查服务是否已运行
    if check_service():
        print("[INFO] 服务已在运行")
        service_process = None
    else:
        print("[INFO] 服务未运行，尝试启动...")
        try:
            service_process = start_service()
        except Exception as e:
            print(f"[ERROR] 启动服务失败: {e}")
            print("\n请手动启动服务:")
            print("  python compass/service_main.py --host 127.0.0.1 --port 8080")
            return 1
    
    # 等待一下确保服务就绪
    time.sleep(2)
    
    # 运行测试脚本
    print("\n" + "=" * 60)
    print("  运行测试脚本")
    print("=" * 60 + "\n")
    
    try:
        # 导入并运行测试
        import test_service_practice
        result = test_service_practice.main()
        
        # 如果服务是我们启动的，询问是否停止
        if service_process:
            print("\n" + "=" * 60)
            print("提示: 服务仍在后台运行")
            print("要停止服务，请按 Ctrl+C 或关闭终端窗口")
            print("=" * 60)
        
        return result
    except KeyboardInterrupt:
        print("\n\n测试被中断")
        if service_process:
            print("正在停止服务...")
            service_process.terminate()
        return 1
    except Exception as e:
        print(f"\n[ERROR] 运行测试时出错: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n\n操作被用户取消")
        sys.exit(1)


