"""
FLASH-DOCK 启动诊断工具
检查WSL环境、conda环境、依赖和启动问题
"""

import sys
import os
import subprocess
from pathlib import Path

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def run_wsl_command(cmd, wsl_distro="Ubuntu-24.04"):
    """在WSL中运行命令"""
    try:
        result = subprocess.run(
            ["wsl", "-d", wsl_distro, "bash", "-c", cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=30
        )
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def check_wsl_status():
    """检查WSL状态"""
    print("=" * 60)
    print("1. 检查WSL状态")
    print("=" * 60)
    
    try:
        result = subprocess.run(
            ["wsl", "--status"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0:
            print("[OK] WSL已安装并启用")
            print(f"   状态信息: {result.stdout.strip()}")
            return True
        else:
            print("[FAIL] WSL未安装或未启用")
            return False
    except Exception as e:
        print(f"[ERROR] 检查WSL状态失败: {e}")
        return False

def check_wsl_distro(wsl_distro="Ubuntu-24.04"):
    """检查WSL发行版"""
    print("\n" + "=" * 60)
    print(f"2. 检查WSL发行版: {wsl_distro}")
    print("=" * 60)
    
    success, stdout, stderr = run_wsl_command("echo 'WSL发行版可用'", wsl_distro)
    if success:
        print(f"[OK] WSL发行版 '{wsl_distro}' 可用")
        return True
    else:
        print(f"[FAIL] WSL发行版 '{wsl_distro}' 不可用")
        print(f"   错误: {stderr}")
        
        # 列出可用的发行版
        try:
            result = subprocess.run(
                ["wsl", "--list", "--verbose"],
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='ignore',
                timeout=10
            )
            if result.returncode == 0:
                print("\n   可用的WSL发行版:")
                print(result.stdout)
        except:
            pass
        return False

def check_conda_environments(wsl_distro="Ubuntu-24.04"):
    """检查conda环境"""
    print("\n" + "=" * 60)
    print("3. 检查conda环境")
    print("=" * 60)
    
    # 检查conda是否可用
    check_conda_cmd = (
        "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
        "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
        "conda --version 2>/dev/null || echo 'conda not found'"
    )
    
    success, stdout, stderr = run_wsl_command(check_conda_cmd, wsl_distro)
    if "conda not found" in stdout or not success:
        print("[FAIL] Conda未找到或未配置")
        print("   请确保在WSL中安装了miniconda或anaconda")
        return False, []
    
    print(f"[OK] Conda可用: {stdout.strip()}")
    
    # 列出所有conda环境
    list_envs_cmd = (
        "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
        "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
        "conda env list"
    )
    
    success, stdout, stderr = run_wsl_command(list_envs_cmd, wsl_distro)
    if success:
        print("\n   可用的conda环境:")
        print(stdout)
        
        # 检查可能的环境名称
        envs = []
        for line in stdout.splitlines():
            if line.strip() and not line.startswith('#'):
                env_name = line.split()[0]
                envs.append(env_name)
        
        return True, envs
    else:
        print("[WARNING] 无法列出conda环境")
        return True, []

def check_specific_env(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
    """检查特定的conda环境"""
    print("\n" + "=" * 60)
    print(f"4. 检查conda环境: {env_name}")
    print("=" * 60)
    
    check_cmd = (
        "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
        "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
        f"conda env list | grep -q '^{env_name} ' && echo 'EXISTS' || echo 'NOT_FOUND'"
    )
    
    success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
    if "EXISTS" in stdout:
        print(f"[OK] Conda环境 '{env_name}' 存在")
        
        # 检查Python版本
        python_cmd = (
            "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {env_name} && python --version"
        )
        success, stdout, stderr = run_wsl_command(python_cmd, wsl_distro)
        if success:
            print(f"   Python版本: {stdout.strip()}")
        
        return True
    else:
        print(f"[FAIL] Conda环境 '{env_name}' 不存在")
        return False

def check_dependencies(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
    """检查依赖"""
    print("\n" + "=" * 60)
    print(f"5. 检查依赖 (环境: {env_name})")
    print("=" * 60)
    
    deps = {
        "streamlit": "streamlit",
        "rdkit": "rdkit",
        "sh": "sh",
        "streamlit_molstar": "streamlit_molstar",
    }
    
    all_ok = True
    for dep_name, import_name in deps.items():
        check_cmd = (
            "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
            "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
            f"conda activate {env_name} && "
            f"python -c 'import {import_name}' 2>&1"
        )
        
        success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
        if success and not stderr:
            print(f"  [OK] {dep_name}")
        else:
            print(f"  [FAIL] {dep_name}")
            if stderr:
                print(f"    错误: {stderr.strip()[:100]}")
            all_ok = False
    
    return all_ok

def check_project_files(wsl_distro="Ubuntu-24.04"):
    """检查项目文件"""
    print("\n" + "=" * 60)
    print("6. 检查项目文件")
    print("=" * 60)
    
    project_root = "/mnt/e/Qinchaojun/AIDD-TRAIN"
    flashdock_dir = f"{project_root}/FLASH_DOCK-main"
    flashdock_py = f"{flashdock_dir}/FlashDock.py"
    
    check_cmd = f"test -d '{project_root}' && echo 'ROOT_EXISTS' || echo 'ROOT_NOT_FOUND'"
    success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
    if "ROOT_EXISTS" in stdout:
        print(f"[OK] 项目根目录存在: {project_root}")
    else:
        print(f"[FAIL] 项目根目录不存在: {project_root}")
        return False
    
    check_cmd = f"test -d '{flashdock_dir}' && echo 'DIR_EXISTS' || echo 'DIR_NOT_FOUND'"
    success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
    if "DIR_EXISTS" in stdout:
        print(f"[OK] FlashDock目录存在: {flashdock_dir}")
    else:
        print(f"[FAIL] FlashDock目录不存在: {flashdock_dir}")
        return False
    
    check_cmd = f"test -f '{flashdock_py}' && echo 'FILE_EXISTS' || echo 'FILE_NOT_FOUND'"
    success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
    if "FILE_EXISTS" in stdout:
        print(f"[OK] FlashDock.py文件存在")
    else:
        print(f"[FAIL] FlashDock.py文件不存在: {flashdock_py}")
        return False
    
    return True

def check_port(wsl_distro="Ubuntu-24.04", port=8501):
    """检查端口占用"""
    print("\n" + "=" * 60)
    print(f"7. 检查端口占用: {port}")
    print("=" * 60)
    
    check_cmd = f"netstat -tuln 2>/dev/null | grep -q ':{port} ' && echo 'IN_USE' || echo 'AVAILABLE'"
    success, stdout, stderr = run_wsl_command(check_cmd, wsl_distro)
    
    if "IN_USE" in stdout:
        print(f"[WARNING] 端口 {port} 已被占用")
        # 查找占用进程
        find_cmd = f"lsof -i :{port} 2>/dev/null || netstat -tulpn 2>/dev/null | grep ':{port} '"
        success, stdout, stderr = run_wsl_command(find_cmd, wsl_distro)
        if stdout:
            print(f"   占用信息:\n{stdout}")
        return False
    else:
        print(f"[OK] 端口 {port} 可用")
        return True

def test_startup(wsl_distro="Ubuntu-24.04", env_name="flash_dock"):
    """测试启动命令"""
    print("\n" + "=" * 60)
    print(f"8. 测试启动命令 (环境: {env_name})")
    print("=" * 60)
    
    project_root = "/mnt/e/Qinchaojun/AIDD-TRAIN"
    flashdock_dir = f"{project_root}/FLASH_DOCK-main"
    
    test_cmd = (
        "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || "
        "source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; "
        f"conda activate {env_name} && "
        f"export PYTHONPATH={project_root} && "
        f"cd {flashdock_dir} && "
        "python -c 'import sys; print(f\"Python: {sys.executable}\"); import streamlit; print(\"Streamlit OK\")'"
    )
    
    success, stdout, stderr = run_wsl_command(test_cmd, wsl_distro)
    if success:
        print("[OK] 启动命令测试通过")
        print(f"   输出:\n{stdout}")
        return True
    else:
        print("[FAIL] 启动命令测试失败")
        if stderr:
            print(f"   错误:\n{stderr}")
        return False

def main():
    print("=" * 60)
    print("FLASH-DOCK 启动诊断工具")
    print("=" * 60)
    print()
    
    wsl_distro = "Ubuntu-24.04"
    
    # 检查WSL状态
    wsl_ok = check_wsl_status()
    if not wsl_ok:
        print("\n[ERROR] WSL未安装或未启用，无法继续诊断")
        return 1
    
    # 检查WSL发行版
    distro_ok = check_wsl_distro(wsl_distro)
    if not distro_ok:
        print("\n[ERROR] WSL发行版不可用，请检查发行版名称")
        return 1
    
    # 检查conda环境
    conda_ok, envs = check_conda_environments(wsl_distro)
    if not conda_ok:
        print("\n[ERROR] Conda未配置，无法继续")
        return 1
    
    # 检查可能的环境名称
    possible_envs = ["flash_dock", "flash_dock_wsl", "flashdocker"]
    found_env = None
    
    for env_name in possible_envs:
        if check_specific_env(wsl_distro, env_name):
            found_env = env_name
            break
    
    if not found_env:
        print("\n[WARNING] 未找到预期的conda环境")
        print("   请检查环境名称，或创建新的环境")
        if envs:
            print(f"\n   可用的环境: {', '.join(envs)}")
        return 1
    
    # 检查依赖
    deps_ok = check_dependencies(wsl_distro, found_env)
    
    # 检查项目文件
    files_ok = check_project_files(wsl_distro)
    
    # 检查端口
    port_ok = check_port(wsl_distro, 8501)
    
    # 测试启动
    startup_ok = test_startup(wsl_distro, found_env)
    
    # 总结
    print("\n" + "=" * 60)
    print("诊断总结")
    print("=" * 60)
    print(f"WSL状态: {'[OK]' if wsl_ok else '[FAIL]'}")
    print(f"WSL发行版: {'[OK]' if distro_ok else '[FAIL]'}")
    print(f"Conda环境: {'[OK]' if found_env else '[FAIL]'}")
    if found_env:
        print(f"  环境名称: {found_env}")
    print(f"依赖检查: {'[OK]' if deps_ok else '[FAIL]'}")
    print(f"项目文件: {'[OK]' if files_ok else '[FAIL]'}")
    print(f"端口状态: {'[OK]' if port_ok else '[WARNING]'}")
    print(f"启动测试: {'[OK]' if startup_ok else '[FAIL]'}")
    
    if wsl_ok and distro_ok and found_env and deps_ok and files_ok and startup_ok:
        print("\n[SUCCESS] 所有检查通过，FLASH-DOCK应该可以正常启动")
        print(f"\n建议使用以下命令启动:")
        print(f"  start_flashdock_wsl.bat")
        print(f"  或修改脚本中的环境名称为: {found_env}")
        return 0
    else:
        print("\n[WARNING] 部分检查未通过，请根据上述信息修复问题")
        return 1

if __name__ == "__main__":
    sys.exit(main())

