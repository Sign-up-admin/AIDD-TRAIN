#!/usr/bin/env python3
"""
Uni-Core 安装脚本
用于解决 ModuleNotFoundError: No module named 'unicore' 错误
支持 Windows 和 Linux 系统
"""

import subprocess
import sys
import os
import tempfile
import shutil

def run_command(cmd, check=True, shell=False):
    """运行命令并返回结果"""
    print(f"执行: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    try:
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, check=check, 
                                  capture_output=True, text=True)
        else:
            result = subprocess.run(cmd, check=check, 
                                  capture_output=True, text=True, shell=shell)
        if result.stdout:
            print(result.stdout)
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        if e.stdout:
            print(e.stdout)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        return False, e.stdout if hasattr(e, 'stdout') else "", e.stderr if hasattr(e, 'stderr') else ""

def check_environment():
    """检查环境"""
    print("=" * 50)
    print("检查环境...")
    print("=" * 50)
    
    # 检查 Python
    print(f"\nPython 版本: {sys.version}")
    
    # 检查 pip
    success, _, _ = run_command([sys.executable, "-m", "pip", "--version"], check=False)
    if not success:
        print("错误: pip 未安装或不可用")
        return False
    
    # 检查 Git
    success, output, _ = run_command(["git", "--version"], check=False)
    if not success:
        print("错误: Git 未安装，请先安装 Git")
        print("Windows: 从 https://git-scm.com/download/win 下载安装")
        print("Linux: sudo apt-get install git 或 sudo yum install git")
        return False
    print(f"Git: {output.strip()}")
    
    # 检查 conda 环境（可选）
    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    if conda_env:
        print(f"当前 conda 环境: {conda_env}")
    else:
        print("未检测到激活的 conda 环境（可选）")
    
    return True

def method1_install():
    """方法 1: 从 GitHub 直接安装"""
    print("\n" + "=" * 50)
    print("方法 1: 从 GitHub 直接安装 Uni-Core")
    print("=" * 50)
    
    cmd = [sys.executable, "-m", "pip", "install", 
           "git+https://github.com/dptech-corp/Uni-Core.git@stable"]
    success, _, _ = run_command(cmd, check=False)
    return success

def method2_install():
    """方法 2: 使用 git:// 协议"""
    print("\n" + "=" * 50)
    print("方法 2: 使用 git:// 协议安装")
    print("=" * 50)
    
    cmd = [sys.executable, "-m", "pip", "install", 
           "git+git://github.com/dptech-corp/Uni-Core.git@stable"]
    success, _, _ = run_command(cmd, check=False)
    return success

def method3_install():
    """方法 3: 克隆后安装"""
    print("\n" + "=" * 50)
    print("方法 3: 克隆后安装")
    print("=" * 50)
    
    temp_dir = tempfile.mkdtemp(prefix="unicore_install_")
    print(f"临时目录: {temp_dir}")
    
    try:
        # 克隆仓库
        repo_dir = os.path.join(temp_dir, "Uni-Core")
        success, _, _ = run_command(["git", "clone", 
                                    "https://github.com/dptech-corp/Uni-Core.git", 
                                    repo_dir], check=False)
        if not success:
            print("错误: Git clone 失败")
            return False
        
        # 切换到 stable 分支（如果存在）
        os.chdir(repo_dir)
        run_command(["git", "checkout", "stable"], check=False)
        
        # 安装
        success, _, _ = run_command([sys.executable, "-m", "pip", "install", "-e", "."], check=False)
        return success
    finally:
        # 清理临时目录
        try:
            shutil.rmtree(temp_dir)
            print(f"已清理临时目录: {temp_dir}")
        except Exception as e:
            print(f"警告: 清理临时目录失败: {e}")

def verify_installation():
    """验证安装"""
    print("\n" + "=" * 50)
    print("验证安装...")
    print("=" * 50)
    
    # 检查基础模块
    try:
        import unicore
        print("✓ unicore 模块导入成功")
    except ImportError as e:
        print(f"✗ unicore 模块导入失败: {e}")
        return False
    
    # 检查核心模块
    try:
        from unicore import checkpoint_utils, distributed_utils, options, utils
        print("✓ 核心模块导入成功")
        return True
    except ImportError as e:
        print(f"✗ 核心模块导入失败: {e}")
        return False

def main():
    """主函数"""
    print("=" * 50)
    print("Uni-Core 安装脚本")
    print("=" * 50)
    print()
    
    # 检查环境
    if not check_environment():
        print("\n环境检查失败，请解决上述问题后重试")
        sys.exit(1)
    
    # 尝试安装方法
    methods = [
        ("方法 1 (GitHub HTTPS)", method1_install),
        ("方法 2 (Git 协议)", method2_install),
        ("方法 3 (克隆后安装)", method3_install),
    ]
    
    for method_name, method_func in methods:
        if method_func():
            print(f"\n✓ {method_name} 成功！")
            break
        else:
            print(f"\n✗ {method_name} 失败，尝试下一个方法...")
    else:
        print("\n" + "=" * 50)
        print("所有安装方法都失败了")
        print("=" * 50)
        print("\n可能的原因:")
        print("  1. 网络连接问题（无法访问 GitHub）")
        print("  2. Git 未正确安装或配置")
        print("  3. 需要配置代理或使用镜像")
        print("\n建议:")
        print("  1. 检查网络连接")
        print("  2. 尝试手动克隆: git clone https://github.com/dptech-corp/Uni-Core.git")
        print("  3. 然后运行: pip install -e Uni-Core/")
        sys.exit(1)
    
    # 验证安装
    if verify_installation():
        print("\n" + "=" * 50)
        print("✓ Uni-Core 安装并验证成功！")
        print("=" * 50)
        sys.exit(0)
    else:
        print("\n" + "=" * 50)
        print("警告: 安装完成但验证失败")
        print("=" * 50)
        print("请检查错误信息，可能需要:")
        print("  1. 重新启动 Python 环境")
        print("  2. 检查 Python 路径")
        print("  3. 手动验证: python -c 'import unicore'")
        sys.exit(1)

if __name__ == "__main__":
    main()

