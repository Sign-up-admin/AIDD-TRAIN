"""测试FLASH-DOCK启动并捕获错误"""
import subprocess
import sys
import time

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def test_flashdock_start():
    print("=" * 60)
    print("测试 FLASH-DOCK 启动")
    print("=" * 60)
    print()
    
    wsl_distro = "Ubuntu-24.04"
    env_name = "flash_dock"
    project_root = "/mnt/e/Qinchaojun/AIDD-TRAIN"
    flashdock_dir = f"{project_root}/FLASH_DOCK-main"
    port = "8501"
    
    # 构建启动命令
    wsl_cmd = (
        f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
        f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
        f"conda activate {env_name} && "
        f"export PYTHONPATH={project_root} && "
        f"cd {flashdock_dir} && "
        f"python -c 'import streamlit; print(\"streamlit OK\")' && "
        f"echo 'Starting streamlit...' && "
        f"streamlit run FlashDock.py --server.port {port} --server.address 0.0.0.0 2>&1 | head -20"
    )
    
    print("执行测试命令...")
    print(f"WSL发行版: {wsl_distro}")
    print(f"环境: {env_name}")
    print(f"目录: {flashdock_dir}")
    print()
    
    try:
        # 先测试环境
        print("1. 测试conda环境...")
        test_env_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {env_name} && "
            f"python -c 'import streamlit; print(\"streamlit OK\")'"
        )
        
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", test_env_cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=15
        )
        
        if result.returncode == 0 and "streamlit OK" in result.stdout:
            print("   [OK] Streamlit可以导入")
        else:
            print("   [FAIL] Streamlit导入失败")
            print(f"   错误: {result.stderr}")
            return False
        
        # 测试文件是否存在
        print("\n2. 测试文件是否存在...")
        test_file_cmd = f"test -f '{flashdock_dir}/FlashDock.py' && echo 'FILE_EXISTS' || echo 'FILE_NOT_FOUND'"
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", test_file_cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        
        if "FILE_EXISTS" in result.stdout:
            print("   [OK] FlashDock.py文件存在")
        else:
            print("   [FAIL] FlashDock.py文件不存在")
            return False
        
        # 尝试启动（短时间测试）
        print("\n3. 尝试启动Streamlit（5秒测试）...")
        start_cmd = (
            f"source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            f"source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {env_name} && "
            f"export PYTHONPATH={project_root} && "
            f"cd {flashdock_dir} && "
            f"timeout 5 streamlit run FlashDock.py --server.port {port} --server.address 0.0.0.0 2>&1 || true"
        )
        
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", start_cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        
        print("   输出:")
        output_lines = result.stdout.splitlines()[:30]  # 只显示前30行
        for line in output_lines:
            print(f"     {line}")
        
        if result.stderr:
            print("\n   错误:")
            error_lines = result.stderr.splitlines()[:30]
            for line in error_lines:
                print(f"     {line}")
        
        # 检查是否有明显错误
        output = result.stdout + result.stderr
        if "error" in output.lower() or "exception" in output.lower() or "traceback" in output.lower():
            print("\n   [WARNING] 检测到错误信息")
        else:
            print("\n   [INFO] 未检测到明显错误")
        
        return True
        
    except subprocess.TimeoutExpired:
        print("[ERROR] 命令执行超时")
        return False
    except Exception as e:
        print(f"[ERROR] 测试失败: {e}")
        return False

if __name__ == "__main__":
    success = test_flashdock_start()
    print("\n" + "=" * 60)
    if success:
        print("测试完成")
    else:
        print("测试失败，请检查上述错误信息")
    print("=" * 60)
    sys.exit(0 if success else 1)

