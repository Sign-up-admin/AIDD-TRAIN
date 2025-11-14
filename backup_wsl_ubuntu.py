"""
备份 WSL Ubuntu 发行版
"""
import subprocess
import sys
import os
from pathlib import Path
import datetime

# 设置UTF-8编码
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")


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


def stop_wsl_distro(distro_name="Ubuntu-24.04"):
    """停止 WSL 发行版"""
    print(f"\n[步骤] 停止 WSL 发行版: {distro_name}")
    
    # 终止指定发行版
    print(f"  终止 {distro_name}...")
    success, output, error = run_command(f"wsl --terminate {distro_name}")
    if success:
        print(f"  [OK] {distro_name} 已终止")
    else:
        print(f"  [INFO] {error}")
    
    # 等待进程完全停止
    import time
    time.sleep(2)
    
    return True


def backup_wsl_distro(distro_name="Ubuntu-24.04", backup_dir=None):
    """备份 WSL 发行版"""
    print(f"\n[步骤] 备份 WSL 发行版: {distro_name}")
    
    # 确定备份目录
    if backup_dir is None:
        # 默认备份到项目父目录的 backups 文件夹
        script_dir = Path(__file__).parent.absolute()
        backup_dir = script_dir.parent / "backups"
    else:
        backup_dir = Path(backup_dir)
    
    backup_dir.mkdir(parents=True, exist_ok=True)
    print(f"  备份目录: {backup_dir}")
    
    # 生成备份文件名
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file = backup_dir / f"{distro_name}_backup_{timestamp}.tar"
    
    print(f"  备份文件: {backup_file}")
    print(f"  开始导出（这可能需要几分钟，请耐心等待）...")
    
    # 执行 WSL 导出
    cmd = f'wsl --export {distro_name} "{backup_file}"'
    success, output, error = run_command(cmd, capture_output=False)
    
    if success:
        # 检查文件是否存在并获取大小
        if backup_file.exists():
            file_size_gb = backup_file.stat().st_size / (1024**3)
            print(f"\n[OK] 备份完成！")
            print(f"  备份文件: {backup_file}")
            print(f"  文件大小: {file_size_gb:.2f} GB")
            return backup_file
        else:
            print(f"\n[ERROR] 备份文件未创建")
            return None
    else:
        print(f"\n[ERROR] 备份失败: {error}")
        return None


def main():
    """主函数"""
    print("=" * 60)
    print("WSL Ubuntu 备份工具")
    print("=" * 60)
    print(f"开始时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    distro_name = "Ubuntu-24.04"
    
    try:
        # 步骤 1: 停止 WSL 发行版
        stop_wsl_distro(distro_name)
        
        # 步骤 2: 备份 WSL 发行版
        backup_file = backup_wsl_distro(distro_name)
        
        if backup_file:
            print("\n" + "=" * 60)
            print("备份流程完成！")
            print("=" * 60)
            print(f"备份文件: {backup_file}")
            print(f"完成时间: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("\n提示: 要恢复备份，使用以下命令:")
            print(f'  wsl --import {distro_name} <安装路径> "{backup_file}"')
            print("=" * 60)
            return 0
        else:
            print("\n[ERROR] 备份失败")
            return 1
        
    except KeyboardInterrupt:
        print("\n\n[ERROR] 用户中断操作")
        return 1
    except Exception as e:
        print(f"\n\n[ERROR] 备份过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())

