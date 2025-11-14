"""
检查WSL状态和配置
"""

import subprocess
import sys
import os

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def run_command(cmd, capture_output=True):
    """执行命令并返回结果"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=capture_output,
            text=True,
            encoding="utf-8",
            errors="replace"
        )
        return result.returncode == 0, result.stdout.strip(), result.stderr.strip()
    except Exception as e:
        return False, "", str(e)


def main():
    print("=" * 60)
    print("WSL 状态检查")
    print("=" * 60)
    print()
    
    # 1. 检查WSL版本
    print("[1] WSL版本信息")
    print("-" * 60)
    success, output, error = run_command("wsl --version")
    if success:
        lines = output.split("\n")
        for line in lines:
            if line.strip():
                print(f"  {line}")
    else:
        print(f"  [ERROR] 无法获取WSL版本: {error}")
    print()
    
    # 2. 检查WSL状态
    print("[2] WSL状态")
    print("-" * 60)
    success, output, error = run_command("wsl --status")
    if success:
        lines = output.split("\n")
        for line in lines:
            if line.strip():
                print(f"  {line}")
    else:
        print(f"  [ERROR] 无法获取WSL状态: {error}")
    print()
    
    # 3. 列出所有WSL发行版
    print("[3] 已安装的WSL发行版")
    print("-" * 60)
    success, output, error = run_command("wsl --list --verbose")
    if success:
        lines = output.split("\n")
        for line in lines:
            if line.strip():
                # 标记默认发行版
                if line.strip().startswith("*"):
                    print(f"  {line} (默认)")
                else:
                    print(f"  {line}")
    else:
        print(f"  [ERROR] 无法列出WSL发行版: {error}")
    print()
    
    # 4. 检查Ubuntu-24.04详细信息
    print("[4] Ubuntu-24.04 详细信息")
    print("-" * 60)
    
    # 检查是否可以访问
    success, output, error = run_command("wsl -d Ubuntu-24.04 -- echo 'OK'")
    if success and "OK" in output:
        print("  [OK] Ubuntu-24.04 可以访问")
        
        # 获取系统信息
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- uname -a")
        if success:
            print(f"  内核信息: {output}")
        
        # 获取发行版信息
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- lsb_release -a")
        if success:
            lines = output.split("\n")
            for line in lines:
                if line.strip() and not line.startswith("wsl:"):
                    print(f"  {line}")
        
        # 获取当前用户
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- whoami")
        if success:
            print(f"  当前用户: {output}")
        
        # 获取当前目录
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- pwd")
        if success:
            print(f"  当前目录: {output}")
        
        # 检查conda是否安装
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- which conda")
        if success and output:
            print(f"  [OK] Conda已安装: {output}")
            # 获取conda版本
            success, version, _ = run_command("wsl -d Ubuntu-24.04 -- conda --version")
            if success:
                print(f"  Conda版本: {version}")
        else:
            print("  [INFO] Conda未安装")
        
        # 检查Python版本
        success, output, error = run_command("wsl -d Ubuntu-24.04 -- python3 --version")
        if success:
            print(f"  Python版本: {output}")
        
        # 检查项目目录是否存在
        project_path = "/mnt/e/Qinchaojun/AIDD-TRAIN"
        success, output, error = run_command(f"wsl -d Ubuntu-24.04 -- test -d {project_path} && echo 'EXISTS'")
        if success and "EXISTS" in output:
            print(f"  [OK] 项目目录存在: {project_path}")
        else:
            print(f"  [WARNING] 项目目录不存在: {project_path}")
    else:
        print("  [WARNING] Ubuntu-24.04 无法访问或未初始化")
        print("  提示: 运行 'wsl -d Ubuntu-24.04' 进行首次设置")
    print()
    
    # 5. 检查默认发行版
    print("[5] 默认WSL发行版")
    print("-" * 60)
    success, output, error = run_command("wsl --list --verbose")
    if success:
        lines = output.split("\n")
        for line in lines:
            if line.strip().startswith("*"):
                print(f"  默认发行版: {line.strip().lstrip('*').strip()}")
                break
    print()
    
    # 6. 总结
    print("=" * 60)
    print("检查完成")
    print("=" * 60)
    print()
    print("提示:")
    print("  - 启动Ubuntu-24.04: wsl -d Ubuntu-24.04")
    print("  - 设置默认发行版: wsl --set-default Ubuntu-24.04")
    print("  - 查看WSL帮助: wsl --help")


if __name__ == "__main__":
    main()

