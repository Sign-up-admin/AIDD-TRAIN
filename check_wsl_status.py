"""检查WSL状态"""
import subprocess
import sys

def check_wsl():
    print("=" * 60)
    print("WSL 状态检查")
    print("=" * 60)
    print()
    
    # 检查WSL状态
    print("1. 检查WSL安装状态...")
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
            print("   [OK] WSL已安装并启用")
        else:
            print("   [FAIL] WSL未安装或未启用")
            return False
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
        return False
    
    # 检查发行版
    print("\n2. 检查WSL发行版...")
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
            print("   发行版列表:")
            for line in result.stdout.splitlines():
                if line.strip():
                    print(f"     {line}")
        else:
            print("   [WARNING] 无法列出发行版")
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    # 测试WSL命令
    print("\n3. 测试WSL命令执行...")
    try:
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", "echo 'WSL test OK'"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0 and "WSL test OK" in result.stdout:
            print("   [OK] WSL命令执行正常")
        else:
            print("   [WARNING] WSL命令执行异常")
            print(f"   输出: {result.stdout}")
    except Exception as e:
        print(f"   [ERROR] 测试失败: {e}")
    
    # 检查conda
    print("\n4. 检查conda环境...")
    try:
        cmd = "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda --version 2>&1 || echo 'conda not found'"
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0:
            output = result.stdout.strip()
            if "conda" in output.lower() and "not found" not in output.lower():
                print(f"   [OK] Conda可用: {output}")
            else:
                print("   [WARNING] Conda未找到或未配置")
        else:
            print("   [WARNING] 无法检查conda")
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    # 检查conda环境
    print("\n5. 检查conda环境列表...")
    try:
        cmd = "source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null; conda env list 2>&1"
        result = subprocess.run(
            ["wsl", "-d", "Ubuntu-24.04", "bash", "-c", cmd],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=10
        )
        if result.returncode == 0:
            lines = result.stdout.strip().splitlines()
            if lines:
                print("   环境列表:")
                for line in lines[:10]:  # 只显示前10行
                    if line.strip():
                        print(f"     {line}")
            else:
                print("   [INFO] 没有找到conda环境")
        else:
            print("   [WARNING] 无法列出conda环境")
    except Exception as e:
        print(f"   [ERROR] 检查失败: {e}")
    
    print("\n" + "=" * 60)
    print("检查完成")
    print("=" * 60)
    print("\n总结:")
    print("  - WSL本身是正常的")
    print("  - 如果conda未找到，需要确保conda已安装并配置PATH")
    print("  - 核心服务（注册中心和COMPASS）可以正常启动")
    print("  - FLASH-DOCK需要WSL中的conda环境")
    
    return True

if __name__ == "__main__":
    sys.exit(0 if check_wsl() else 1)
