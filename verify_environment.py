"""
环境验证脚本
验证COMPASS服务所需的所有依赖和配置
"""

import sys
import os
from pathlib import Path
from typing import List, Tuple

# 颜色输出
try:
    import colorama
    colorama.init()
    GREEN = colorama.Fore.GREEN
    RED = colorama.Fore.RED
    YELLOW = colorama.Fore.YELLOW
    BLUE = colorama.Fore.BLUE
    RESET = colorama.Fore.RESET
except ImportError:
    GREEN = RED = YELLOW = BLUE = RESET = ""


def print_header(text: str):
    """打印标题"""
    print(f"\n{BLUE}{'='*60}{RESET}")
    print(f"{BLUE}  {text}{RESET}")
    print(f"{BLUE}{'='*60}{RESET}")


def print_check(name: str, status: bool, message: str = ""):
    """打印检查结果"""
    status_str = f"{GREEN}[OK]{RESET}" if status else f"{RED}[FAIL]{RESET}"
    print(f"{status_str} {name}")
    if message:
        indent = "   " if status else "   "
        print(f"{indent}{message}")


def verify_python_version() -> Tuple[bool, str]:
    """验证Python版本"""
    version = sys.version_info
    version_str = f"{version.major}.{version.minor}.{version.micro}"
    is_valid = version.major == 3 and version.minor >= 8
    message = f"Python {version_str} (需要 >= 3.8)"
    return is_valid, message


def verify_dependencies() -> List[Tuple[str, bool, str]]:
    """验证所有必需的依赖包"""
    required_packages = [
        # Web框架
        ("fastapi", "FastAPI框架", ">=0.104.0"),
        ("uvicorn", "ASGI服务器", ">=0.24.0"),
        ("pydantic", "数据验证", ">=2.0.0"),
        ("python-multipart", "文件上传支持", ">=0.0.6"),
        
        # 深度学习
        ("torch", "PyTorch", None),
        ("torch_geometric", "PyTorch Geometric", None),
        ("numpy", "NumPy", None),
        
        # 化学信息学
        ("rdkit", "RDKit化学信息学库", None),
        
        # 工具库
        ("tqdm", "进度条", None),
        
        # 服务依赖（可选但推荐）
        ("websockets", "WebSocket支持", None),
        ("aiofiles", "异步文件操作", None),
        ("aiosqlite", "异步SQLite", None),
        ("python-dotenv", "环境变量管理", None),
        ("pyyaml", "YAML解析", None),
    ]
    
    results = []
    for package, description, version_req in required_packages:
        try:
            if package == "rdkit":
                from rdkit import Chem
                version = "已安装"
            else:
                module = __import__(package)
                version = getattr(module, "__version__", "已安装")
            
            # 检查版本要求（简单检查）
            if version_req and version != "已安装":
                # 这里可以添加更详细的版本检查
                pass
            
            results.append((package, True, f"{description}: {version}"))
        except ImportError:
            is_optional = package in ["websockets", "aiofiles", "aiosqlite", "python-dotenv", "pyyaml"]
            status_msg = "未安装（可选）" if is_optional else "未安装（必需）"
            results.append((package, is_optional, f"{description}: {status_msg}"))
        except Exception as e:
            results.append((package, False, f"{description}: 错误 - {str(e)}"))
    
    return results


def verify_project_structure() -> List[Tuple[str, bool, str]]:
    """验证项目文件结构"""
    project_root = Path(__file__).parent.resolve()
    
    required_files = [
        ("compass/__init__.py", "COMPASS包初始化"),
        ("compass/service_main.py", "服务入口点"),
        ("compass/service/server.py", "FastAPI服务器"),
        ("compass/service/config.py", "服务配置"),
        ("compass/service/routes/__init__.py", "路由模块"),
        ("services/__init__.py", "服务包初始化"),
        ("services/common/utils.py", "通用工具"),
        ("services/registry/server.py", "服务注册中心"),
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


def verify_imports() -> List[Tuple[str, bool, str]]:
    """验证关键模块可以导入"""
    project_root = Path(__file__).parent.resolve()
    project_root_str = str(project_root)
    
    # 确保项目根目录在路径中
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)
    
    # 确保PYTHONPATH包含项目根目录
    pythonpath = os.environ.get("PYTHONPATH", "")
    if project_root_str not in pythonpath:
        if pythonpath:
            os.environ["PYTHONPATH"] = os.pathsep.join([project_root_str, pythonpath])
        else:
            os.environ["PYTHONPATH"] = project_root_str
    
    test_imports = [
        ("compass", "COMPASS主包"),
        ("compass.service_main", "服务入口点"),
        ("compass.service.server", "FastAPI服务器"),
        ("compass.service.config", "服务配置"),
        ("compass.service.routes.health", "健康检查路由"),
        ("compass.service.routes.training", "训练路由"),
        ("services.common.utils", "通用工具"),
        ("services.registry.server", "服务注册中心"),
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


def verify_config_files() -> List[Tuple[str, bool, str]]:
    """验证配置文件"""
    project_root = Path(__file__).parent.resolve()
    
    config_files = [
        ("requirements.txt", "Python依赖列表"),
        ("requirements-service.txt", "服务依赖列表"),
    ]
    
    results = []
    for filename, description in config_files:
        file_path = project_root / filename
        exists = file_path.exists()
        if exists:
            try:
                # 尝试读取文件
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                    line_count = len([l for l in content.split('\n') if l.strip()])
                results.append((filename, True, f"{description}: 存在 ({line_count} 行)"))
            except Exception as e:
                results.append((filename, False, f"{description}: 无法读取 - {str(e)}"))
        else:
            results.append((filename, False, f"{description}: 不存在"))
    
    return results


def verify_pythonpath() -> Tuple[bool, str]:
    """验证PYTHONPATH配置"""
    project_root = Path(__file__).parent.resolve()
    project_root_str = str(project_root)
    
    pythonpath = os.environ.get("PYTHONPATH", "")
    paths = pythonpath.split(os.pathsep) if pythonpath else []
    paths = [p for p in paths if p]
    
    project_in_path = (
        project_root_str in paths or
        any(Path(p).resolve() == project_root for p in paths if p) or
        project_root_str in sys.path
    )
    
    if project_in_path:
        return True, f"PYTHONPATH包含项目根目录: {project_root_str}"
    else:
        return False, f"PYTHONPATH未包含项目根目录: {project_root_str}"


def main():
    """主函数"""
    print_header("COMPASS环境验证工具")
    print(f"项目根目录: {Path(__file__).parent.resolve()}")
    print(f"Python版本: {sys.version}")
    print(f"Python路径: {sys.executable}")
    
    all_critical_ok = True
    warnings = []
    
    # 1. Python版本验证
    print_header("1. Python版本验证")
    is_valid, message = verify_python_version()
    print_check("Python版本", is_valid, message)
    if not is_valid:
        all_critical_ok = False
    
    # 2. 依赖包验证
    print_header("2. 依赖包验证")
    deps = verify_dependencies()
    for package, status, message in deps:
        print_check(package, status, message)
        if not status and package not in ["websockets", "aiofiles", "aiosqlite", "python-dotenv", "pyyaml"]:
            all_critical_ok = False
        elif not status:
            warnings.append(f"可选依赖 {package} 未安装")
    
    # 3. 项目结构验证
    print_header("3. 项目文件结构验证")
    files = verify_project_structure()
    for file_path, status, message in files:
        print_check(file_path, status, message)
        if not status:
            all_critical_ok = False
    
    # 4. 模块导入验证
    print_header("4. 模块导入验证")
    imports = verify_imports()
    for module, status, message in imports:
        print_check(module, status, message)
        if not status:
            all_critical_ok = False
    
    # 5. 配置文件验证
    print_header("5. 配置文件验证")
    configs = verify_config_files()
    for config_file, status, message in configs:
        print_check(config_file, status, message)
        if not status:
            warnings.append(f"配置文件 {config_file} 缺失或无法读取")
    
    # 6. PYTHONPATH验证
    print_header("6. PYTHONPATH验证")
    path_ok, path_msg = verify_pythonpath()
    print_check("PYTHONPATH配置", path_ok, path_msg)
    if not path_ok:
        warnings.append("PYTHONPATH未正确配置")
    
    # 总结
    print_header("验证总结")
    
    if all_critical_ok:
        print(f"{GREEN}[OK] 所有关键检查通过！环境配置正确。{RESET}")
        if warnings:
            print(f"\n{YELLOW}警告：{RESET}")
            for warning in warnings:
                print(f"  - {warning}")
        return 0
    else:
        print(f"{RED}[FAIL] 发现关键问题，请修复后重试。{RESET}")
        if warnings:
            print(f"\n{YELLOW}额外警告：{RESET}")
            for warning in warnings:
                print(f"  - {warning}")
        print(f"\n{YELLOW}修复建议：{RESET}")
        print("1. 安装缺失的依赖: pip install -r requirements.txt")
        print("2. 确保项目文件结构完整")
        print("3. 设置PYTHONPATH环境变量指向项目根目录")
        print("4. 运行 diagnose_compass.py 获取详细诊断信息")
        return 1


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n验证被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n{RED}验证过程中发生错误: {str(e)}{RESET}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

