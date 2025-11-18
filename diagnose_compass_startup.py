"""
COMPASS服务启动失败诊断脚本
"""
import sys
import os
import subprocess
import socket
from pathlib import Path

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

def print_section(title):
    """打印章节标题"""
    print("\n" + "=" * 60)
    print(f" {title}")
    print("=" * 60)

def check_port(port):
    """检查端口是否被占用"""
    try:
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.settimeout(1)
        result = sock.connect_ex(('localhost', port))
        sock.close()
        return result == 0
    except Exception as e:
        print(f"  检查端口时出错: {e}")
        return False

def check_port_process(port):
    """检查端口占用的进程"""
    try:
        result = subprocess.run(
            ["netstat", "-ano"],
            capture_output=True,
            text=True,
            encoding='utf-8',
            errors='ignore',
            timeout=5
        )
        processes = []
        for line in result.stdout.splitlines():
            if f":{port}" in line and "LISTENING" in line:
                parts = line.split()
                if len(parts) >= 5:
                    pid = parts[-1]
                    processes.append(pid)
        return processes
    except Exception as e:
        print(f"  检查端口进程时出错: {e}")
        return []

def test_import(module_name, description):
    """测试模块导入"""
    try:
        __import__(module_name)
        print(f"  [OK] {description}")
        return True
    except ImportError as e:
        print(f"  [FAIL] {description}: {e}")
        return False
    except Exception as e:
        print(f"  [ERROR] {description}: {e}")
        return False

def check_file_exists(filepath, description):
    """检查文件是否存在"""
    path = Path(filepath)
    if path.exists():
        print(f"  [OK] {description}: {path}")
        return True
    else:
        print(f"  [FAIL] {description}: {path} 不存在")
        return False

def check_directory_exists(dirpath, description):
    """检查目录是否存在"""
    path = Path(dirpath)
    if path.exists() and path.is_dir():
        print(f"  [OK] {description}: {path}")
        # 检查写权限
        try:
            test_file = path / ".write_test"
            test_file.touch()
            test_file.unlink()
            print(f"  [OK] {description} 可写")
        except Exception as e:
            print(f"  [WARNING] {description} 不可写: {e}")
        return True
    else:
        print(f"  [FAIL] {description}: {path} 不存在")
        return False

def main():
    print_section("COMPASS服务启动失败诊断")
    
    project_root = Path(__file__).parent.resolve()
    print(f"项目根目录: {project_root}")
    
    # 1. 检查端口状态
    print_section("1. 检查端口状态")
    
    print("\n检查8080端口（COMPASS服务）:")
    port_8080_in_use = check_port(8080)
    if port_8080_in_use:
        print("  [INFO] 8080端口已被占用")
        pids = check_port_process(8080)
        if pids:
            print(f"  占用进程PID: {', '.join(pids)}")
    else:
        print("  [INFO] 8080端口未被占用")
    
    print("\n检查8500端口（服务注册中心）:")
    port_8500_in_use = check_port(8500)
    if port_8500_in_use:
        print("  [OK] 8500端口正在监听（服务注册中心可能正在运行）")
        pids = check_port_process(8500)
        if pids:
            print(f"  进程PID: {', '.join(pids)}")
    else:
        print("  [WARNING] 8500端口未被占用（服务注册中心未运行）")
    
    # 2. 检查关键文件
    print_section("2. 检查关键文件")
    
    files_to_check = [
        ("compass/service_main.py", "服务入口文件"),
        ("compass/service/server.py", "服务器主文件"),
        ("compass/service/config.py", "配置文件"),
        ("compass/service/config_manager.py", "配置管理器"),
        ("services/common/utils.py", "通用工具模块"),
    ]
    
    all_files_ok = True
    for filepath, description in files_to_check:
        if not check_file_exists(filepath, description):
            all_files_ok = False
    
    # 3. 检查目录
    print_section("3. 检查目录")
    
    dirs_to_check = [
        ("logs", "日志目录"),
        ("data", "数据目录"),
        ("checkpoints", "检查点目录"),
    ]
    
    for dirpath, description in dirs_to_check:
        check_directory_exists(dirpath, description)
    
    # 4. 检查Python路径
    print_section("4. 检查Python路径")
    
    print(f"  Python版本: {sys.version}")
    print(f"  Python可执行文件: {sys.executable}")
    print(f"  当前工作目录: {os.getcwd()}")
    print(f"  PYTHONPATH: {os.environ.get('PYTHONPATH', '未设置')}")
    print(f"  sys.path前3项: {sys.path[:3]}")
    
    # 添加项目根目录到路径
    project_root_str = str(project_root)
    if project_root_str not in sys.path:
        sys.path.insert(0, project_root_str)
        print(f"  [INFO] 已添加项目根目录到sys.path: {project_root_str}")
    
    # 5. 检查依赖导入
    print_section("5. 检查依赖导入")
    
    dependencies = [
        ("fastapi", "FastAPI框架"),
        ("uvicorn", "Uvicorn ASGI服务器"),
        ("pydantic", "Pydantic数据验证"),
        ("requests", "Requests HTTP库"),
        ("torch", "PyTorch"),
        ("torch_geometric", "PyTorch Geometric（可选）"),
    ]
    
    all_deps_ok = True
    for module, description in dependencies:
        if not test_import(module, description):
            all_deps_ok = False
    
    # 6. 检查COMPASS模块导入
    print_section("6. 检查COMPASS模块导入")
    
    compass_modules = [
        ("compass.service.server", "服务器模块"),
        ("compass.service.config", "配置模块"),
        ("compass.service.registry.client", "注册中心客户端"),
        ("services.common.utils", "通用工具模块"),
    ]
    
    all_compass_ok = True
    for module, description in compass_modules:
        if not test_import(module, description):
            all_compass_ok = False
    
    # 7. 检查配置文件加载
    print_section("7. 检查配置文件加载")
    
    try:
        from compass.service.config import SERVICE_CONFIG
        print("  [OK] 配置文件加载成功")
        print(f"  主机: {SERVICE_CONFIG.get('host', 'N/A')}")
        print(f"  端口: {SERVICE_CONFIG.get('port', 'N/A')}")
        print(f"  注册中心URL: {SERVICE_CONFIG.get('registry_url', 'N/A')}")
        print(f"  日志目录: {SERVICE_CONFIG.get('log_dir', 'N/A')}")
    except Exception as e:
        print(f"  [FAIL] 配置文件加载失败: {e}")
        import traceback
        traceback.print_exc()
    
    # 8. 检查日志文件
    print_section("8. 检查日志文件")
    
    log_files = [
        ("logs/compass-service.log", "主日志文件"),
        ("logs/compass-service_errors.log", "错误日志文件"),
    ]
    
    for log_file, description in log_files:
        log_path = Path(log_file)
        if log_path.exists():
            size = log_path.stat().st_size
            print(f"  [OK] {description}: {log_path} ({size} 字节)")
            # 读取最后几行
            try:
                with open(log_path, 'r', encoding='utf-8', errors='ignore') as f:
                    lines = f.readlines()
                    if lines:
                        print(f"    最后3行:")
                        for line in lines[-3:]:
                            print(f"      {line.strip()}")
            except Exception as e:
                print(f"    读取日志文件时出错: {e}")
        else:
            print(f"  [INFO] {description}: {log_path} 不存在（首次运行）")
    
    # 9. 尝试手动启动测试
    print_section("9. 启动测试建议")
    
    print("\n根据诊断结果，建议执行以下步骤：")
    
    issues = []
    
    if not port_8500_in_use:
        issues.append("服务注册中心未运行（8500端口）")
        print("\n  [问题] 服务注册中心未运行")
        print("  [解决] 先启动服务注册中心:")
        print("    python services/registry/server.py --host 0.0.0.0 --port 8500")
    
    if not all_deps_ok:
        issues.append("依赖缺失")
        print("\n  [问题] 部分依赖未安装")
        print("  [解决] 安装缺失的依赖:")
        print("    pip install fastapi uvicorn pydantic requests")
        if not test_import("torch", ""):
            print("    pip install torch")
    
    if not all_compass_ok:
        issues.append("COMPASS模块导入失败")
        print("\n  [问题] COMPASS模块导入失败")
        print("  [解决] 检查PYTHONPATH和模块路径")
    
    if not all_files_ok:
        issues.append("关键文件缺失")
        print("\n  [问题] 关键文件缺失")
        print("  [解决] 检查项目完整性")
    
    if not issues:
        print("\n  [OK] 未发现明显问题，可以尝试启动服务:")
        print("    python compass/service_main.py --host 0.0.0.0 --port 8080")
    else:
        print(f"\n  [总结] 发现 {len(issues)} 个问题需要解决")
        for i, issue in enumerate(issues, 1):
            print(f"    {i}. {issue}")
    
    print("\n" + "=" * 60)
    print("诊断完成")
    print("=" * 60)

if __name__ == "__main__":
    main()

