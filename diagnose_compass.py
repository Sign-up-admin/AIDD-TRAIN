"""
COMPASS服务诊断脚本
检查环境、依赖、端口等配置问题
"""

import sys
import os
import socket
import subprocess
from pathlib import Path
from typing import List, Tuple, Dict

# 颜色输出（Windows支持）
try:
    import colorama
    colorama.init()
    GREEN = colorama.Fore.GREEN
    RED = colorama.Fore.RED
    YELLOW = colorama.Fore.YELLOW
    RESET = colorama.Fore.RESET
except ImportError:
    GREEN = RED = YELLOW = RESET = ""


def print_header(text: str):
    """打印标题"""
    print(f"\n{'='*60}")
    print(f"  {text}")
    print(f"{'='*60}")


def print_check(name: str, status: bool, message: str = ""):
    """打印检查结果"""
    status_str = f"{GREEN}[OK]{RESET}" if status else f"{RED}[FAIL]{RESET}"
    print(f"{status_str} {name}")
    if message:
        print(f"      {message}")


def check_python_version() -> Tuple[bool, str]:
    """检查Python版本"""
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    is_valid = version.major == 3 and version.minor >= 8
    message = f"Python {version_str} (需要 >= 3.8)"
    return is_valid, message


def check_conda_env() -> Tuple[bool, str]:
    """检查AIDDTRAIN conda环境"""
    env_paths = [
        Path("D:/conda_envs/AIDDTRAIN"),
        Path("C:/ProgramData/Anaconda3/envs/AIDDTRAIN"),
        Path(os.path.expanduser("~/anaconda3/envs/AIDDTRAIN")),
        Path(os.path.expanduser("~/miniconda3/envs/AIDDTRAIN")),
    ]
    
    for env_path in env_paths:
        python_exe = env_path / "python.exe"
        if python_exe.exists():
            return True, f"找到环境: {env_path}"
    
    # 尝试通过conda命令查找
    try:
        result = subprocess.run(
            ["conda", "env", "list"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if "AIDDTRAIN" in result.stdout:
            return True, "环境存在（通过conda命令）"
    except:
        pass
    
    return False, "未找到AIDDTRAIN环境"


def check_dependencies() -> List[Tuple[str, bool, str]]:
    """检查关键依赖包"""
    required_packages = [
        ("fastapi", "FastAPI框架"),
        ("uvicorn", "ASGI服务器"),
        ("pydantic", "数据验证"),
        ("torch", "PyTorch"),
        ("numpy", "NumPy"),
        ("torch_geometric", "PyTorch Geometric"),
        ("rdkit", "RDKit化学信息学库"),
        ("tqdm", "进度条"),
    ]
    
    results = []
    for package, description in required_packages:
        try:
            if package == "rdkit":
                from rdkit import Chem
                version = "已安装"
            elif package == "torch_geometric":
                import torch_geometric
                version = torch_geometric.__version__
            else:
                module = __import__(package)
                version = getattr(module, "__version__", "已安装")
            results.append((package, True, f"{description}: {version}"))
        except ImportError:
            results.append((package, False, f"{description}: 未安装"))
        except Exception as e:
            results.append((package, False, f"{description}: 错误 - {str(e)}"))
    
    return results


def check_pythonpath() -> Tuple[bool, str]:
    """检查PYTHONPATH配置"""
    project_root = Path(__file__).parent.resolve()
    pythonpath = os.environ.get("PYTHONPATH", "")
    
    paths = pythonpath.split(os.pathsep) if pythonpath else []
    paths = [p for p in paths if p]
    
    project_in_path = str(project_root) in paths or any(
        Path(p).resolve() == project_root for p in paths if p
    )
    
    if project_in_path:
        return True, f"PYTHONPATH包含项目根目录: {project_root}"
    else:
        return False, f"PYTHONPATH未包含项目根目录: {project_root}"


def check_port(port: int) -> Tuple[bool, str]:
    """检查端口是否被占用"""
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(1)
        result = sock.connect_ex(('localhost', port))
        sock.close()
        
        if result == 0:
            return False, f"端口 {port} 已被占用"
        else:
            return True, f"端口 {port} 可用"
    except Exception as e:
        return False, f"检查端口 {port} 时出错: {str(e)}"


def check_project_structure() -> List[Tuple[str, bool, str]]:
    """检查项目文件结构"""
    project_root = Path(__file__).parent.resolve()
    
    required_files = [
        ("compass/service_main.py", "COMPASS服务入口"),
        ("compass/service/server.py", "FastAPI服务器"),
        ("compass/service/config.py", "服务配置"),
        ("services/common/utils.py", "通用工具"),
        ("services/registry/server.py", "服务注册中心"),
        ("requirements.txt", "依赖文件"),
    ]
    
    results = []
    for rel_path, description in required_files:
        full_path = project_root / rel_path
        exists = full_path.exists()
        results.append((
            rel_path,
            exists,
            f"{description}: {'存在' if exists else '缺失'}"
        ))
    
    return results


def check_imports() -> List[Tuple[str, bool, str]]:
    """检查关键模块导入"""
    project_root = Path(__file__).parent.resolve()
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    
    test_imports = [
        ("compass.service_main", "COMPASS服务入口"),
        ("compass.service.server", "FastAPI服务器"),
        ("services.common.utils", "通用工具"),
        ("compass.service.config", "服务配置"),
    ]
    
    results = []
    for module_name, description in test_imports:
        try:
            __import__(module_name)
            results.append((module_name, True, f"{description}: 导入成功"))
        except ImportError as e:
            results.append((module_name, False, f"{description}: 导入失败 - {str(e)}"))
        except Exception as e:
            results.append((module_name, False, f"{description}: 错误 - {str(e)}"))
    
    return results


def main():
    """主函数"""
    print_header("COMPASS服务诊断工具")
    print(f"项目根目录: {Path(__file__).parent.resolve()}")
    
    all_ok = True
    
    # 1. Python版本检查
    print_header("1. Python环境检查")
    is_valid, message = check_python_version()
    print_check("Python版本", is_valid, message)
    if not is_valid:
        all_ok = False
    
    # 2. Conda环境检查
    print_header("2. Conda环境检查")
    env_ok, env_msg = check_conda_env()
    print_check("AIDDTRAIN环境", env_ok, env_msg)
    if not env_ok:
        all_ok = False
        print(f"\n{YELLOW}提示: 如果环境存在但未找到，请确保：{RESET}")
        print("  1. 环境名称正确：AIDDTRAIN")
        print("  2. 环境路径在以下位置之一：")
        print("     - D:\\conda_envs\\AIDDTRAIN")
        print("     - C:\\ProgramData\\Anaconda3\\envs\\AIDDTRAIN")
        print("     - %USERPROFILE%\\anaconda3\\envs\\AIDDTRAIN")
    
    # 3. 依赖包检查
    print_header("3. 依赖包检查")
    deps = check_dependencies()
    for package, status, message in deps:
        print_check(package, status, message)
        if not status:
            all_ok = False
    
    # 4. PYTHONPATH检查
    print_header("4. PYTHONPATH检查")
    path_ok, path_msg = check_pythonpath()
    print_check("PYTHONPATH配置", path_ok, path_msg)
    if not path_ok:
        all_ok = False
        project_root = Path(__file__).parent.resolve()
        print(f"\n{YELLOW}修复方法: 设置PYTHONPATH环境变量{RESET}")
        print(f"  set PYTHONPATH={project_root}")
    
    # 5. 端口检查
    print_header("5. 端口检查")
    port8080_ok, port8080_msg = check_port(8080)
    print_check("端口8080 (COMPASS)", port8080_ok, port8080_msg)
    if not port8080_ok:
        all_ok = False
        print(f"\n{YELLOW}提示: 端口8080被占用，可能需要停止现有服务{RESET}")
    
    port8500_ok, port8500_msg = check_port(8500)
    print_check("端口8500 (Registry)", port8500_ok, port8500_msg)
    if not port8500_ok:
        print(f"\n{YELLOW}提示: 端口8500被占用，可能是服务注册中心正在运行{RESET}")
    
    # 6. 项目结构检查
    print_header("6. 项目文件结构检查")
    files = check_project_structure()
    for file_path, status, message in files:
        print_check(file_path, status, message)
        if not status:
            all_ok = False
    
    # 7. 模块导入检查
    print_header("7. 模块导入检查")
    imports = check_imports()
    for module, status, message in imports:
        print_check(module, status, message)
        if not status:
            all_ok = False
    
    # 总结
    print_header("诊断总结")
    if all_ok:
        print(f"{GREEN}[OK] 所有检查通过！COMPASS服务应该可以正常启动。{RESET}")
        return 0
    else:
        print(f"{RED}[FAIL] 发现一些问题，请根据上述提示进行修复。{RESET}")
        print(f"\n{YELLOW}建议的修复步骤：{RESET}")
        print("1. 确保AIDDTRAIN conda环境已创建并激活")
        print("2. 安装缺失的依赖包: pip install -r requirements.txt")
        print("3. 设置PYTHONPATH环境变量指向项目根目录")
        print("4. 检查并释放被占用的端口")
        print("5. 运行此诊断脚本再次检查")
        return 1


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n诊断被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n{RED}诊断过程中发生错误: {str(e)}{RESET}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

