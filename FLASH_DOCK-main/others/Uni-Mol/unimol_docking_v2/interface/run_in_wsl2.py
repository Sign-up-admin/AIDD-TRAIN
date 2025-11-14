#!/usr/bin/env python3
"""
WSL2 运行包装脚本
从 Windows 调用 WSL2 运行 Uni-Mol 对接预测程序
"""

import os
import sys
import subprocess
import platform
import argparse
from pathlib import Path


def convert_to_wsl_path(windows_path):
    """将 Windows 路径转换为 WSL2 路径"""
    if not windows_path:
        return windows_path
    
    # 如果已经是 Linux 路径格式，直接返回
    if windows_path.startswith('/'):
        return windows_path
    
    # 转换 Windows 路径到 WSL2 路径
    # E:\path\to\file -> /mnt/e/path/to/file
    path = os.path.abspath(windows_path)
    
    # 提取驱动器字母
    if len(path) >= 2 and path[1] == ':':
        drive_letter = path[0].lower()
        rest_path = path[2:].replace('\\', '/')
        if not rest_path.startswith('/'):
            rest_path = '/' + rest_path
        wsl_path = f'/mnt/{drive_letter}{rest_path}'
        return wsl_path
    
    return path


def main():
    parser = argparse.ArgumentParser(
        description='从 Windows 调用 WSL2 运行 Uni-Mol 对接预测程序',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  python run_in_wsl2.py --mode single \\
    --input-protein E:\\data\\protein.pdb \\
    --input-ligand E:\\data\\ligand.sdf \\
    --input-docking-grid E:\\data\\grid.json \\
    --output-ligand-name result \\
    --output-ligand-dir E:\\output
        """
    )
    
    # 添加所有 demo.py 的参数
    parser.add_argument('--mode', type=str, default='single',
                       choices=['single', 'batch_one2one', 'batch_one2many'],
                       help='运行模式')
    parser.add_argument('--model-dir', type=str,
                       help='模型文件路径')
    parser.add_argument('--input-protein', type=str,
                       help='输入蛋白质文件路径')
    parser.add_argument('--input-ligand', type=str,
                       help='输入配体文件路径')
    parser.add_argument('--input-docking-grid', type=str,
                       help='对接网格文件路径')
    parser.add_argument('--input-batch-file', type=str,
                       help='批量输入文件（CSV）')
    parser.add_argument('--output-ligand-name', type=str, default='ligand_predict',
                       help='输出配体名称')
    parser.add_argument('--output-ligand-dir', type=str, default='./predict_sdf',
                       help='输出目录')
    parser.add_argument('--batch-size', type=int, default=4,
                       help='批处理大小')
    parser.add_argument('--nthreads', type=int, default=8,
                       help='线程数')
    parser.add_argument('--conf-size', type=int, default=10,
                       help='构象数量')
    parser.add_argument('--cluster', action='store_true',
                       help='是否进行构象聚类')
    parser.add_argument('--use-current-ligand-conf', action='store_true',
                       help='使用当前配体构象')
    parser.add_argument('--steric-clash-fix', action='store_true',
                       help='是否进行空间冲突修复')
    parser.add_argument('--wsl-distro', type=str, default='Ubuntu-24.04',
                       help='WSL2 发行版名称')
    
    args = parser.parse_args()
    
    # 检查是否在 Windows 上
    if platform.system() != 'Windows':
        print("错误: 此脚本只能在 Windows 上运行")
        print("如果已在 WSL2 中，请直接运行: python3 demo.py [参数]")
        sys.exit(1)
    
    # 构建 WSL2 命令
    wsl_cmd = ['wsl', '-d', args.wsl_distro, '--']
    
    # 获取脚本所在目录（在 WSL2 中的路径）
    script_dir = Path(__file__).parent
    wsl_script_dir = convert_to_wsl_path(str(script_dir))
    demo_script = f'{wsl_script_dir}/demo.py'
    
    # 构建 Python 命令
    python_cmd = ['python3', demo_script]
    
    # 转换所有路径参数
    path_args = [
        '--model-dir', '--input-protein', '--input-ligand',
        '--input-docking-grid', '--input-batch-file', '--output-ligand-dir'
    ]
    
    for arg_name in dir(args):
        if arg_name.startswith('_'):
            continue
        
        arg_value = getattr(args, arg_name)
        if arg_value is None:
            continue
        
        # 跳过 wsl_distro，这是本脚本的参数
        if arg_name == 'wsl_distro':
            continue
        
        # 转换参数名
        arg_flag = '--' + arg_name.replace('_', '-')
        
        # 如果是路径参数，转换路径
        if arg_flag in path_args and isinstance(arg_value, str):
            wsl_path = convert_to_wsl_path(arg_value)
            python_cmd.extend([arg_flag, wsl_path])
        elif isinstance(arg_value, bool):
            if arg_value:
                python_cmd.append(arg_flag)
        else:
            python_cmd.extend([arg_flag, str(arg_value)])
    
    # 组合完整命令
    full_cmd = wsl_cmd + python_cmd
    
    print("=" * 60)
    print("从 Windows 调用 WSL2 运行 Uni-Mol 对接预测")
    print("=" * 60)
    print(f"WSL2 发行版: {args.wsl_distro}")
    print(f"命令: {' '.join(full_cmd)}")
    print("=" * 60)
    print()
    
    # 执行命令
    try:
        result = subprocess.run(full_cmd, check=False)
        sys.exit(result.returncode)
    except KeyboardInterrupt:
        print("\n用户中断")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()



