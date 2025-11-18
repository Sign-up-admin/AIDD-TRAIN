"""
诊断和修复 torch_geometric 依赖库问题
检查 conda 环境中相关库的安装状态
增强版：包含详细的版本兼容性、DLL依赖、文件完整性检查
"""

import sys
import subprocess
import os
import json
import platform
from pathlib import Path
from typing import Dict, List, Tuple, Optional

def run_command(cmd, capture_output=True):
    """运行命令并返回结果"""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=capture_output,
            text=True,
            timeout=60
        )
        return result.returncode == 0, result.stdout, result.stderr
    except Exception as e:
        return False, "", str(e)

def check_conda_env():
    """检查当前conda环境"""
    success, stdout, stderr = run_command("conda info --envs")
    if success:
        print("=" * 60)
        print("当前 Conda 环境列表:")
        print(stdout)
        print("=" * 60)
    
    # 检查当前激活的环境
    env_name = os.environ.get("CONDA_DEFAULT_ENV", "未激活")
    print(f"\n当前激活的环境: {env_name}")
    
    # 检查环境路径
    env_path = os.environ.get("CONDA_PREFIX", "")
    if env_path:
        print(f"环境路径: {env_path}")
    return env_name, env_path

def check_package_installed(package_name):
    """检查包是否已安装"""
    # 使用conda list
    success, stdout, stderr = run_command(f"conda list {package_name}")
    if success and package_name.lower() in stdout.lower():
        lines = stdout.strip().split('\n')
        for line in lines:
            if package_name.lower() in line.lower() and not line.startswith('#'):
                return True, line
        return False, "未找到"
    return False, "未找到"

def check_pip_package(package_name):
    """使用pip检查包"""
    success, stdout, stderr = run_command(f"pip show {package_name}")
    if success:
        return True, stdout
    return False, stderr

def check_torch_info():
    """检查PyTorch信息（增强版）"""
    print("\n" + "=" * 60)
    print("PyTorch 信息（详细）:")
    print("=" * 60)
    
    torch_info = {}
    
    try:
        import torch
        torch_info['version'] = torch.__version__
        torch_info['cuda_available'] = torch.cuda.is_available()
        torch_info['cuda_version'] = torch.version.cuda if hasattr(torch.version, 'cuda') else None
        torch_info['cudnn_version'] = torch.backends.cudnn.version() if torch.cuda.is_available() else None
        
        print(f"PyTorch 版本: {torch_info['version']}")
        print(f"CUDA 可用: {torch_info['cuda_available']}")
        
        if torch_info['cuda_available']:
            print(f"CUDA 版本: {torch_info['cuda_version']}")
            if torch_info['cudnn_version']:
                print(f"cuDNN 版本: {torch_info['cudnn_version']}")
        else:
            print("⚠️  警告: CUDA 不可用，这是 CPU 版本")
            # 检查版本字符串中是否包含 +cpu
            if '+cpu' in torch_info['version']:
                print("   确认: PyTorch 版本字符串包含 '+cpu'，确认为 CPU 版本")
            else:
                print("   注意: PyTorch 版本字符串不包含 '+cpu'，但 CUDA 不可用")
        
        # 检查编译配置
        print(f"\n编译信息:")
        print(f"  Python 版本: {sys.version.split()[0]}")
        print(f"  平台架构: {platform.machine()}")
        print(f"  操作系统: {platform.system()} {platform.release()}")
        
    except ImportError:
        print("❌ 错误: PyTorch 未安装")
        return False, None
    
    return True, torch_info

def check_torch_geometric_deps():
    """检查torch_geometric依赖库（增强版）"""
    print("\n" + "=" * 60)
    print("torch_geometric 依赖库检查（详细）:")
    print("=" * 60)
    
    deps = [
        "pyg-lib",
        "torch-scatter",
        "torch-cluster",
        "torch-spline-conv",
        "torch-sparse",
        "torch-geometric"
    ]
    
    results = {}
    
    for dep in deps:
        print(f"\n检查 {dep}:")
        # 检查conda安装
        conda_installed, conda_info = check_package_installed(dep)
        if conda_installed:
            print(f"  [CONDA] 已安装: {conda_info}")
            results[dep] = {"installer": "conda", "info": conda_info, "version": None}
        else:
            # 检查pip安装
            pip_installed, pip_info = check_pip_package(dep)
            if pip_installed:
                print(f"  [PIP] 已安装")
                version = None
                location = None
                # 提取详细信息
                for line in pip_info.split('\n'):
                    if 'Version:' in line:
                        version = line.split(':', 1)[1].strip()
                        print(f"    版本: {version}")
                        results[dep] = {"installer": "pip", "info": pip_info, "version": version}
                    elif 'Location:' in line:
                        location = line.split(':', 1)[1].strip()
                        print(f"    位置: {location}")
                
                # 检查 wheel 文件名以确定 CUDA/CPU 版本
                if location:
                    try:
                        dist_info = Path(location) / f"{dep.replace('-', '_')}-{version}.dist-info" if version else None
                        if dist_info and dist_info.exists():
                            # 查找 RECORD 文件来获取 wheel 文件名
                            record_file = dist_info / "RECORD"
                            if record_file.exists():
                                with open(record_file, 'r', encoding='utf-8') as f:
                                    for line in f:
                                        if '.whl' in line:
                                            wheel_name = line.split(',')[0]
                                            print(f"    Wheel: {wheel_name}")
                                            # 检查是否包含 CUDA 标识
                                            if '+cu' in wheel_name or '_cu' in wheel_name:
                                                print(f"    ⚠️  检测到 CUDA 版本的 wheel 文件")
                                                results[dep]["wheel_type"] = "cuda"
                                            elif '+cpu' in wheel_name or '_cpu' in wheel_name:
                                                print(f"    ✓ 检测到 CPU 版本的 wheel 文件")
                                                results[dep]["wheel_type"] = "cpu"
                                            break
                    except Exception as e:
                        pass
            else:
                print(f"  ❌ [未安装]")
                results[dep] = {"installer": "none", "info": None, "version": None}
    
    return results

def check_pyd_files(env_path):
    """检查.pyd文件是否存在（增强版）"""
    print("\n" + "=" * 60)
    print("检查 .pyd 文件（详细）:")
    print("=" * 60)
    
    if not env_path:
        print("无法确定环境路径")
        return {}
    
    env_path = Path(env_path)
    site_packages = env_path / "Lib" / "site-packages"
    
    pyd_files = {
        "libpyg.pyd": site_packages / "libpyg.pyd",
        "torch_scatter/_version_cuda.pyd": site_packages / "torch_scatter" / "_version_cuda.pyd",
        "torch_cluster/_version_cuda.pyd": site_packages / "torch_cluster" / "_version_cuda.pyd",
        "torch_spline_conv/_version_cuda.pyd": site_packages / "torch_spline_conv" / "_version_cuda.pyd",
        "torch_sparse/_version_cuda.pyd": site_packages / "torch_sparse" / "_version_cuda.pyd",
    }
    
    # 也检查 CPU 版本的文件
    cpu_pyd_files = {
        "torch_scatter/_version_cpu.pyd": site_packages / "torch_scatter" / "_version_cpu.pyd",
        "torch_cluster/_version_cpu.pyd": site_packages / "torch_cluster" / "_version_cpu.pyd",
        "torch_spline_conv/_version_cpu.pyd": site_packages / "torch_spline_conv" / "_version_cpu.pyd",
        "torch_sparse/_version_cpu.pyd": site_packages / "torch_sparse" / "_version_cpu.pyd",
    }
    
    results = {}
    
    print("\nCUDA 版本的 .pyd 文件:")
    for name, path in pyd_files.items():
        if path.exists():
            stat = path.stat()
            size = stat.st_size
            # 检查文件大小是否合理（通常 .pyd 文件应该 > 1KB）
            size_ok = size > 1024
            readable = os.access(path, os.R_OK)
            print(f"  {'✓' if size_ok and readable else '⚠️ '} {name}:")
            print(f"      存在: 是")
            print(f"      大小: {size} 字节 ({size/1024:.2f} KB)")
            print(f"      可读: {'是' if readable else '否'}")
            print(f"      路径: {path}")
            if not size_ok:
                print(f"      ⚠️  警告: 文件大小异常小，可能已损坏")
            results[name] = {"exists": True, "size": size, "readable": readable, "path": str(path)}
        else:
            print(f"  ✗ {name}: 不存在")
            print(f"      路径: {path}")
            results[name] = {"exists": False, "size": 0, "readable": False, "path": str(path)}
    
    print("\nCPU 版本的 .pyd 文件:")
    for name, path in cpu_pyd_files.items():
        if path.exists():
            stat = path.stat()
            size = stat.st_size
            readable = os.access(path, os.R_OK)
            print(f"  ✓ {name}: 存在 ({size} 字节)")
            results[name] = {"exists": True, "size": size, "readable": readable, "path": str(path)}
        else:
            print(f"  - {name}: 不存在（正常，如果使用 CUDA 版本）")
    
    return results

def check_dll_dependencies(pyd_path: Path) -> Dict:
    """检查 .pyd 文件的 DLL 依赖关系"""
    dependencies = {
        "cuda_dlls": [],
        "other_dlls": [],
        "missing_dlls": [],
        "error": None
    }
    
    if not pyd_path.exists():
        dependencies["error"] = "文件不存在"
        return dependencies
    
    # 在 Windows 上使用 dumpbin 或依赖检查工具
    # 首先尝试使用 dumpbin（Visual Studio 工具）
    dumpbin_paths = [
        r"C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\*\bin\Hostx64\x64\dumpbin.exe",
        r"C:\Program Files (x86)\Microsoft Visual Studio\*\VC\bin\dumpbin.exe",
    ]
    
    dumpbin_found = None
    for pattern in dumpbin_paths:
        import glob
        matches = glob.glob(pattern)
        if matches:
            dumpbin_found = matches[0]
            break
    
    if dumpbin_found:
        try:
            cmd = f'"{dumpbin_found}" /dependents "{pyd_path}"'
            success, stdout, stderr = run_command(cmd)
            if success:
                for line in stdout.split('\n'):
                    line = line.strip()
                    if line.endswith('.dll'):
                        dll_name = line.lower()
                        # 检查是否是 CUDA DLL
                        if any(cuda_dll in dll_name for cuda_dll in ['cudart', 'cublas', 'curand', 'cusparse', 'cusolver', 'cufft', 'cudnn']):
                            dependencies["cuda_dlls"].append(line)
                        else:
                            dependencies["other_dlls"].append(line)
        except Exception as e:
            dependencies["error"] = f"dumpbin 执行失败: {str(e)}"
    else:
        # 如果没有 dumpbin，尝试使用 Python 的 ctypes 来加载并检查
        try:
            import ctypes
            # 尝试加载 .pyd 文件来触发依赖检查
            # 注意：这可能会失败，但可以捕获错误信息
            pass  # 暂时跳过，因为可能会影响诊断
        except Exception as e:
            dependencies["error"] = f"无法检查 DLL 依赖: {str(e)}"
    
    return dependencies

def check_cuda_dlls_in_path() -> Dict:
    """检查 PATH 环境变量中是否包含 CUDA DLL"""
    path_dirs = os.environ.get('PATH', '').split(os.pathsep)
    cuda_paths = []
    cuda_dlls_found = []
    
    # 常见的 CUDA 安装路径
    common_cuda_paths = [
        r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA",
        r"C:\Program Files (x86)\NVIDIA GPU Computing Toolkit\CUDA",
    ]
    
    for path_dir in path_dirs:
        if 'cuda' in path_dir.lower():
            cuda_paths.append(path_dir)
            # 检查该路径下是否有 CUDA DLL
            path_obj = Path(path_dir)
            if path_obj.exists():
                for dll_file in path_obj.glob('*.dll'):
                    if any(cuda_name in dll_file.name.lower() for cuda_name in ['cudart', 'cublas', 'curand']):
                        cuda_dlls_found.append(str(dll_file))
    
    return {
        "cuda_paths": cuda_paths,
        "cuda_dlls_found": cuda_dlls_found,
        "common_cuda_paths_exist": [str(Path(p)) for p in common_cuda_paths if Path(p).exists()]
    }

def analyze_version_compatibility(torch_info: Dict, deps_info: Dict) -> Dict:
    """分析版本兼容性"""
    analysis = {
        "issues": [],
        "warnings": [],
        "recommendations": []
    }
    
    if not torch_info:
        return analysis
    
    torch_version = torch_info.get('version', '')
    torch_is_cpu = not torch_info.get('cuda_available', False) or '+cpu' in torch_version
    
    # 检查每个扩展库的版本兼容性
    for dep_name, dep_info in deps_info.items():
        if dep_info.get('installer') == 'none':
            continue
        
        wheel_type = dep_info.get('wheel_type')
        
        if torch_is_cpu:
            # PyTorch 是 CPU 版本
            if wheel_type == 'cuda':
                analysis["issues"].append(
                    f"{dep_name}: 安装了 CUDA 版本的扩展库，但 PyTorch 是 CPU 版本。"
                    f"这会导致无法加载 .pyd 文件，因为缺少 CUDA runtime DLL。"
                )
            elif wheel_type == 'cpu':
                analysis["warnings"].append(
                    f"{dep_name}: 已正确安装 CPU 版本"
                )
        else:
            # PyTorch 是 CUDA 版本
            if wheel_type == 'cpu':
                analysis["warnings"].append(
                    f"{dep_name}: 安装了 CPU 版本的扩展库，但 PyTorch 是 CUDA 版本。"
                    f"建议使用 CUDA 版本以获得更好的性能。"
                )
            elif wheel_type == 'cuda':
                # 检查 CUDA 版本是否匹配
                cuda_version = torch_info.get('cuda_version')
                if cuda_version:
                    analysis["warnings"].append(
                        f"{dep_name}: 已安装 CUDA 版本，PyTorch CUDA 版本: {cuda_version}"
                    )
    
    # 检查 PyTorch 版本与扩展库的兼容性
    # 提取主版本号
    try:
        torch_major_minor = '.'.join(torch_version.split('+')[0].split('.')[:2])
        analysis["recommendations"].append(
            f"PyTorch 版本: {torch_version}，建议使用兼容的扩展库版本"
        )
    except:
        pass
    
    return analysis

def check_installation_history() -> Dict:
    """检查安装历史和 requirements.txt"""
    history = {
        "requirements_files": [],
        "installation_commands": [],
        "issues": []
    }
    
    # 检查 requirements.txt
    req_files = [
        "requirements.txt",
        "requirements-dev.txt",
        "requirements-service.txt",
    ]
    
    for req_file in req_files:
        req_path = Path(req_file)
        if req_path.exists():
            history["requirements_files"].append(str(req_path))
            try:
                with open(req_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                    # 检查是否包含 torch-geometric 相关依赖
                    if 'torch-geometric' in content.lower() or 'torch_geometric' in content.lower():
                        # 检查是否指定了 wheel 源
                        if 'data.pyg.org' in content:
                            history["installation_commands"].append(
                                f"{req_file} 中使用了 PyG wheel 源"
                            )
                            # 检查是否使用了 CUDA 版本
                            if '+cu' in content or 'cu128' in content or 'cu118' in content:
                                history["issues"].append(
                                    f"{req_file} 中指定了 CUDA 版本的 wheel 源，"
                                    f"但当前 PyTorch 可能是 CPU 版本"
                                )
            except Exception as e:
                history["issues"].append(f"无法读取 {req_file}: {str(e)}")
    
    return history

def check_environment_variables() -> Dict:
    """检查环境变量"""
    env_vars = {
        "pythonpath": os.environ.get('PYTHONPATH', ''),
        "path": os.environ.get('PATH', ''),
        "conda_prefix": os.environ.get('CONDA_PREFIX', ''),
        "conda_default_env": os.environ.get('CONDA_DEFAULT_ENV', ''),
        "cuda_path": os.environ.get('CUDA_PATH', ''),
        "cuda_path_v11_8": os.environ.get('CUDA_PATH_V11_8', ''),
        "cuda_path_v12_1": os.environ.get('CUDA_PATH_V12_1', ''),
    }
    
    # 检查 PATH 中是否包含 CUDA
    path_has_cuda = 'cuda' in env_vars["path"].lower()
    
    return {
        "variables": env_vars,
        "path_has_cuda": path_has_cuda,
        "cuda_env_vars_set": any(env_vars[k] for k in ['cuda_path', 'cuda_path_v11_8', 'cuda_path_v12_1'])
    }

def get_installation_commands(torch_version, cuda_available):
    """生成安装命令"""
    print("\n" + "=" * 60)
    print("推荐的修复命令:")
    print("=" * 60)
    
    # 获取PyTorch版本
    try:
        import torch
        torch_version_str = torch.__version__
        cuda_version = torch.version.cuda if torch.cuda.is_available() else None
    except:
        torch_version_str = "未知"
        cuda_version = None
    
    print(f"\n检测到的 PyTorch 版本: {torch_version_str}")
    if cuda_version:
        print(f"检测到的 CUDA 版本: {cuda_version}")
    else:
        print("CUDA: 不可用（CPU版本）")
    
    print("\n方案1: 使用 pip 重新安装（推荐）")
    print("-" * 60)
    
    if cuda_version:
        # CUDA版本
        print("# 先卸载现有版本")
        print("pip uninstall -y pyg-lib torch-scatter torch-cluster torch-spline-conv torch-sparse torch-geometric")
        print("\n# 重新安装（使用CUDA版本）")
        print("pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-{}.html".format(torch_version_str.replace("+cu", "+cu")))
        print("pip install pyg-lib -f https://data.pyg.org/whl/torch-{}.html".format(torch_version_str.replace("+cu", "+cu")))
        print("pip install torch-geometric")
    else:
        # CPU版本
        print("# 先卸载现有版本")
        print("pip uninstall -y pyg-lib torch-scatter torch-cluster torch-spline-conv torch-sparse torch-geometric")
        print("\n# 重新安装（CPU版本）")
        print("pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-{}+cpu.html".format(torch_version_str))
        print("pip install pyg-lib -f https://data.pyg.org/whl/torch-{}+cpu.html".format(torch_version_str))
        print("pip install torch-geometric")
    
    print("\n方案2: 使用 conda-forge 安装")
    print("-" * 60)
    print("conda install -c conda-forge pyg")
    
    print("\n方案3: 如果上述方法都不行，尝试从源码编译")
    print("-" * 60)
    print("pip install --no-build-isolation torch-scatter torch-sparse torch-cluster torch-spline-conv")
    print("pip install pyg-lib --no-build-isolation")
    print("pip install torch-geometric")

def generate_diagnostic_report(torch_info, deps_info, pyd_info, compatibility_analysis, 
                                installation_history, env_vars, cuda_dll_info) -> str:
    """生成完整的诊断报告"""
    report = []
    report.append("\n" + "=" * 80)
    report.append("完整诊断报告")
    report.append("=" * 80)
    
    # 根本原因分析
    report.append("\n## 根本原因分析")
    report.append("-" * 80)
    
    if torch_info and not torch_info.get('cuda_available', False):
        report.append("\n### 主要问题：版本不匹配")
        report.append("PyTorch 是 CPU 版本，但扩展库可能安装了 CUDA 版本")
        
        # 检查是否有 CUDA 版本的 wheel
        cuda_wheels = [dep for dep, info in deps_info.items() 
                      if info.get('wheel_type') == 'cuda']
        if cuda_wheels:
            report.append(f"\n检测到以下扩展库安装了 CUDA 版本：")
            for dep in cuda_wheels:
                report.append(f"  - {dep}")
            report.append("\n这些 CUDA 版本的扩展库需要 CUDA runtime DLL，")
            report.append("但在 CPU 版本的 PyTorch 环境中这些 DLL 不可用，")
            report.append("导致无法加载 .pyd 文件。")
    
    # 兼容性问题
    if compatibility_analysis.get('issues'):
        report.append("\n### 兼容性问题：")
        for issue in compatibility_analysis['issues']:
            report.append(f"  ❌ {issue}")
    
    # 警告
    if compatibility_analysis.get('warnings'):
        report.append("\n### 警告：")
        for warning in compatibility_analysis['warnings']:
            report.append(f"  ⚠️  {warning}")
    
    # DLL 依赖问题
    if cuda_dll_info:
        if not cuda_dll_info.get('cuda_dlls_found'):
            report.append("\n### DLL 依赖问题：")
            report.append("PATH 环境变量中未找到 CUDA DLL")
            report.append("如果使用 CUDA 版本的扩展库，需要确保 CUDA runtime DLL 在 PATH 中")
    
    # 修复建议
    report.append("\n## 修复建议")
    report.append("-" * 80)
    
    if torch_info and not torch_info.get('cuda_available', False):
        report.append("\n### 方案 1：重新安装 CPU 版本的扩展库（推荐）")
        torch_version = torch_info.get('version', '').split('+')[0]
        report.append(f"\n# 1. 卸载现有版本")
        report.append("pip uninstall -y pyg-lib torch-scatter torch-cluster torch-spline-conv torch-sparse torch-geometric")
        report.append(f"\n# 2. 安装 CPU 版本（匹配 PyTorch {torch_version}）")
        report.append(f"pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-{torch_version}+cpu.html")
        report.append(f"pip install pyg-lib -f https://data.pyg.org/whl/torch-{torch_version}+cpu.html")
        report.append("pip install torch-geometric")
        
        report.append("\n### 方案 2：如果确实需要 CUDA 支持")
        report.append("需要先安装 CUDA 版本的 PyTorch，然后再安装 CUDA 版本的扩展库")
    
    # 安装历史问题
    if installation_history.get('issues'):
        report.append("\n### 安装配置问题：")
        for issue in installation_history['issues']:
            report.append(f"  ⚠️  {issue}")
    
    return "\n".join(report)

def main():
    """主函数（增强版）"""
    print("=" * 80)
    print("torch_geometric 依赖库诊断工具（增强版）")
    print("=" * 80)
    
    # 收集所有诊断信息
    all_results = {}
    
    # 1. 检查conda环境
    env_name, env_path = check_conda_env()
    all_results['env'] = {"name": env_name, "path": env_path}
    
    # 2. 检查PyTorch（增强版）
    torch_ok, torch_info = check_torch_info()
    all_results['torch'] = torch_info
    
    if not torch_ok:
        print("\n❌ PyTorch 未安装，无法继续诊断")
        return
    
    # 3. 检查依赖库（增强版）
    deps_status = check_torch_geometric_deps()
    all_results['dependencies'] = deps_status
    
    # 4. 检查.pyd文件（增强版）
    pyd_info = {}
    if env_path:
        pyd_info = check_pyd_files(env_path)
        all_results['pyd_files'] = pyd_info
        
        # 5. 检查 DLL 依赖（针对存在的 CUDA 版本 .pyd 文件）
        print("\n" + "=" * 60)
        print("检查 DLL 依赖关系:")
        print("=" * 60)
        
        cuda_pyd_files = [name for name, info in pyd_info.items() 
                         if info.get('exists') and '_cuda.pyd' in name]
        
        if cuda_pyd_files:
            print("\n检查 CUDA 版本的 .pyd 文件的 DLL 依赖:")
            for pyd_name in cuda_pyd_files[:2]:  # 只检查前两个作为示例
                pyd_path = Path(pyd_info[pyd_name]['path'])
                dll_deps = check_dll_dependencies(pyd_path)
                print(f"\n  {pyd_name}:")
                if dll_deps.get('error'):
                    print(f"    ⚠️  {dll_deps['error']}")
                else:
                    if dll_deps.get('cuda_dlls'):
                        print(f"    需要 CUDA DLL: {', '.join(dll_deps['cuda_dlls'][:3])}")
                    if dll_deps.get('missing_dlls'):
                        print(f"    缺失 DLL: {', '.join(dll_deps['missing_dlls'])}")
        
        # 检查 PATH 中的 CUDA DLL
        print("\n检查 PATH 环境变量中的 CUDA DLL:")
        cuda_dll_info = check_cuda_dlls_in_path()
        all_results['cuda_dlls'] = cuda_dll_info
        
        if cuda_dll_info.get('cuda_paths'):
            print(f"  找到 CUDA 路径: {len(cuda_dll_info['cuda_paths'])} 个")
        else:
            print("  ⚠️  未在 PATH 中找到 CUDA 路径")
        
        if cuda_dll_info.get('cuda_dlls_found'):
            print(f"  找到 CUDA DLL: {len(cuda_dll_info['cuda_dlls_found'])} 个")
        else:
            print("  ⚠️  未找到 CUDA runtime DLL")
    
    # 6. 版本兼容性分析
    print("\n" + "=" * 60)
    print("版本兼容性分析:")
    print("=" * 60)
    compatibility_analysis = analyze_version_compatibility(torch_info, deps_status)
    all_results['compatibility'] = compatibility_analysis
    
    if compatibility_analysis.get('issues'):
        print("\n发现兼容性问题:")
        for issue in compatibility_analysis['issues']:
            print(f"  ❌ {issue}")
    
    if compatibility_analysis.get('warnings'):
        print("\n警告:")
        for warning in compatibility_analysis['warnings']:
            print(f"  ⚠️  {warning}")
    
    # 7. 检查安装历史
    print("\n" + "=" * 60)
    print("检查安装历史和配置文件:")
    print("=" * 60)
    installation_history = check_installation_history()
    all_results['installation_history'] = installation_history
    
    if installation_history.get('requirements_files'):
        print(f"\n找到 requirements 文件: {len(installation_history['requirements_files'])} 个")
        for req_file in installation_history['requirements_files']:
            print(f"  - {req_file}")
    
    if installation_history.get('issues'):
        print("\n发现配置问题:")
        for issue in installation_history['issues']:
            print(f"  ⚠️  {issue}")
    
    # 8. 检查环境变量
    print("\n" + "=" * 60)
    print("检查环境变量:")
    print("=" * 60)
    env_vars = check_environment_variables()
    all_results['environment_variables'] = env_vars
    
    print(f"\nCONDA_PREFIX: {env_vars['variables'].get('conda_prefix', '未设置')}")
    print(f"CONDA_DEFAULT_ENV: {env_vars['variables'].get('conda_default_env', '未设置')}")
    print(f"CUDA_PATH: {env_vars['variables'].get('cuda_path', '未设置')}")
    print(f"PATH 中包含 CUDA: {'是' if env_vars['path_has_cuda'] else '否'}")
    
    # 9. 生成修复建议
    try:
        import torch
        cuda_available = torch.cuda.is_available()
        torch_version = torch.__version__
    except:
        cuda_available = False
        torch_version = "未知"
    
    get_installation_commands(torch_version, cuda_available)
    
    # 10. 生成完整诊断报告
    diagnostic_report = generate_diagnostic_report(
        torch_info, deps_status, pyd_info, compatibility_analysis,
        installation_history, env_vars, cuda_dll_info
    )
    print(diagnostic_report)
    
    print("\n" + "=" * 80)
    print("诊断完成")
    print("=" * 80)
    
    # 保存诊断结果到文件
    try:
        import json
        report_file = Path("torch_geometric_diagnostic_report.json")
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(all_results, f, indent=2, ensure_ascii=False, default=str)
        print(f"\n诊断结果已保存到: {report_file}")
    except Exception as e:
        print(f"\n无法保存诊断报告: {str(e)}")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n诊断被用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"\n诊断过程中发生错误: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

