# Copyright (c) DP Technology.
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import os
import platform
import shutil
import subprocess
import sys
import logging
import torch
from .processor import Processor

logger = logging.getLogger(__name__)


def is_wsl2():
    """检测是否在 WSL2 环境中运行"""
    try:
        # 检查 /proc/version 文件
        if os.path.exists('/proc/version'):
            with open('/proc/version', 'r') as f:
                version_info = f.read().lower()
                if 'microsoft' in version_info or 'wsl' in version_info:
                    return True
        # 检查环境变量
        if os.environ.get('WSL_DISTRO_NAME') or os.environ.get('WSL_INTEROP'):
            return True
        # 检查 /proc/sys/kernel/osrelease（WSL2 特有）
        if os.path.exists('/proc/sys/kernel/osrelease'):
            with open('/proc/sys/kernel/osrelease', 'r') as f:
                osrelease = f.read().lower()
                if 'microsoft' in osrelease or 'wsl' in osrelease:
                    return True
    except:
        pass
    return False


def convert_to_wsl_path(windows_path):
    """将 Windows 路径转换为 WSL2 路径"""
    if not windows_path:
        return windows_path
    
    # 如果已经是 Linux 路径格式，直接返回
    if windows_path.startswith('/'):
        return windows_path
    
    # 转换 Windows 路径到 WSL2 路径
    # E:\path\to\file -> /mnt/e/path/to/file
    # 处理相对路径和绝对路径
    if os.path.isabs(windows_path):
        path = windows_path
    else:
        path = os.path.abspath(windows_path)
    
    # 提取驱动器字母
    if len(path) >= 2 and path[1] == ':':
        drive_letter = path[0].lower()
        rest_path = path[2:].replace('\\', '/')
        # 确保路径以 / 开头
        if not rest_path.startswith('/'):
            rest_path = '/' + rest_path
        wsl_path = f'/mnt/{drive_letter}{rest_path}'
        return wsl_path
    
    return path


def convert_to_windows_path(wsl_path):
    """将 WSL2 路径转换为 Windows 路径"""
    if not wsl_path:
        return wsl_path
    
    # 如果是 Windows 路径格式，直接返回
    if ':' in wsl_path and (wsl_path[1] == ':' or wsl_path.startswith('\\\\')):
        return wsl_path
    
    # 转换 WSL2 路径到 Windows 路径
    # /mnt/e/path/to/file -> E:\path\to\file
    if wsl_path.startswith('/mnt/'):
        parts = wsl_path[5:].split('/', 1)
        if len(parts) == 2:
            drive_letter = parts[0].upper()
            rest_path = parts[1].replace('/', '\\')
            windows_path = f'{drive_letter}:\\{rest_path}'
            return windows_path
    
    return wsl_path


class UnimolPredictor:
    def __init__(self, model_dir, mode='single', nthreads=8, conf_size=10, cluster=False, use_current_ligand_conf=False, steric_clash_fix=False):
        self.model_dir = model_dir
        self.mode = mode
        self.nthreads = nthreads
        self.use_current_ligand_conf = use_current_ligand_conf
        self.cluster = cluster
        self.steric_clash_fix = steric_clash_fix
        if self.use_current_ligand_conf:
            self.conf_size = 1
        else:
            self.conf_size = conf_size

    def preprocess(self, input_protein, input_ligand, input_docking_grid, output_ligand_name, output_ligand_dir):
        # process the input pocket.pdb and ligand.sdf, store in LMDB.
        preprocessor = Processor.build_processors(
            self.mode, self.nthreads, conf_size=self.conf_size, cluster=self.cluster,
            use_current_ligand_conf=self.use_current_ligand_conf)
        processed_data = preprocessor.preprocess(input_protein, input_ligand, input_docking_grid, output_ligand_name, output_ligand_dir)

        # return lmdb path
        return processed_data

    def predict(self, input_protein:str, 
                input_ligand:str, 
                input_docking_grid:str, 
                output_ligand_name:str, 
                output_ligand_dir:str, 
                batch_size:int):
        # Check CUDA availability and environment
        is_windows = platform.system() == 'Windows'
        is_wsl = is_wsl2()
        is_linux = platform.system() == 'Linux' and not is_wsl
        cuda_available = torch.cuda.is_available()
        
        # Log environment information
        if is_wsl:
            logger.info("[WSL2] Running in WSL2 environment")
            logger.info(f"[WSL2] WSL_DISTRO_NAME: {os.environ.get('WSL_DISTRO_NAME', 'N/A')}")
        elif is_linux:
            logger.info("[Linux] Running in native Linux environment")
        elif is_windows:
            logger.info("[Windows] Running in Windows environment")
        
        # 在 WSL2 中，如果输入路径是 Windows 格式，自动转换
        if is_wsl:
            if input_protein and (':' in input_protein or input_protein.startswith('\\')):
                input_protein = convert_to_wsl_path(input_protein)
                logger.info(f"[WSL2] Converted input_protein path: {input_protein}")
            if input_ligand and (':' in input_ligand or input_ligand.startswith('\\')):
                input_ligand = convert_to_wsl_path(input_ligand)
                logger.info(f"[WSL2] Converted input_ligand path: {input_ligand}")
            if input_docking_grid and (':' in input_docking_grid or input_docking_grid.startswith('\\')):
                input_docking_grid = convert_to_wsl_path(input_docking_grid)
                logger.info(f"[WSL2] Converted input_docking_grid path: {input_docking_grid}")
            if output_ligand_dir and (':' in output_ligand_dir or output_ligand_dir.startswith('\\')):
                output_ligand_dir = convert_to_wsl_path(output_ligand_dir)
                logger.info(f"[WSL2] Converted output_ligand_dir path: {output_ligand_dir}")
            if self.model_dir and (':' in self.model_dir or self.model_dir.startswith('\\')):
                self.model_dir = convert_to_wsl_path(self.model_dir)
                logger.info(f"[WSL2] Converted model_dir path: {self.model_dir}")
        
        if is_windows:
            if cuda_available:
                print(f"CUDA is available. Using GPU: {torch.cuda.get_device_name(0)}")
                print(f"CUDA version: {torch.version.cuda}")
            else:
                print("WARNING: CUDA is not available. The model will run on CPU, which may be very slow.")
                print("If you have a GPU, please ensure:")
                print("  1. CUDA toolkit is installed")
                print("  2. PyTorch with CUDA support is installed (not CPU-only version)")
                print("  3. GPU drivers are up to date")
                print("  4. Verify with: python -c 'import torch; print(torch.cuda.is_available())'")
        
        lmdb_name = self.preprocess(input_protein, input_ligand, input_docking_grid, output_ligand_name, output_ligand_dir)
        
        pkt_data_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "example_data", "dict_pkt.txt")
        mol_data_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "example_data", "dict_mol.txt")
        script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "unimol", "infer.py")
        user_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "unimol")
        
        # Copy required files to output directory (cross-platform compatible)
        output_dir_abs = os.path.abspath(output_ligand_dir)
        os.makedirs(output_dir_abs, exist_ok=True)
        shutil.copy(pkt_data_path, output_dir_abs)
        shutil.copy(mol_data_path, output_dir_abs)
        
        # Convert model path to absolute path
        # If model_dir is relative, resolve it relative to the project root (FLASH_DOCK-main)
        if not os.path.isabs(self.model_dir):
            # Get the unimol_docking_v2 directory
            # __file__ is: FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface/predictor/unimol_predictor.py
            # unimol_docking_v2_dir should be: FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
            unimol_docking_v2_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            # Get the FLASH_DOCK-main directory
            # unimol_docking_v2_dir = FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2
            # Need to go up 3 levels: unimol_docking_v2 -> Uni-Mol -> others -> FLASH_DOCK-main
            flash_dock_main_dir = os.path.dirname(os.path.dirname(os.path.dirname(unimol_docking_v2_dir)))
            
            # Clean the model path (remove ./ prefix if present)
            model_path_clean = self.model_dir.lstrip('./').replace('\\', '/')
            
            # If path starts with 'others/Uni-Mol/unimol_docking_v2/', resolve from FLASH_DOCK-main
            if model_path_clean.startswith('others/Uni-Mol/unimol_docking_v2/'):
                model_path = os.path.join(flash_dock_main_dir, model_path_clean)
            # If path is just a filename, assume it's in unimol_docking_v2 directory
            elif '/' not in model_path_clean and '\\' not in model_path_clean:
                model_path = os.path.join(unimol_docking_v2_dir, model_path_clean)
            # Otherwise, try to resolve from FLASH_DOCK-main
            else:
                model_path = os.path.join(flash_dock_main_dir, model_path_clean)
            
            model_path = os.path.abspath(model_path)
        else:
            model_path = os.path.abspath(self.model_dir)
        
        # Check if model file exists
        if not os.path.exists(model_path):
            # Suggest the expected location
            unimol_docking_v2_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            expected_path = os.path.join(unimol_docking_v2_dir, 'unimol_docking_v2_240517.pt')
            raise FileNotFoundError(
                f"Model file not found: {model_path}\n"
                f"Please ensure the model file exists at the specified path.\n"
                f"Expected location: {expected_path}\n"
                f"Download from: https://www.dropbox.com/scl/fi/sfhrtx1tjprce18wbvmdr/unimol_docking_v2_240517.pt"
            )
        
        # inference - use subprocess for cross-platform compatibility
        # Set CUDA_VISIBLE_DEVICES environment variable
        env = os.environ.copy()
        
        # Only set CUDA_VISIBLE_DEVICES if CUDA is available
        if cuda_available:
            env['CUDA_VISIBLE_DEVICES'] = '0'
        else:
            # Remove CUDA_VISIBLE_DEVICES if set, to allow CPU fallback
            env.pop('CUDA_VISIBLE_DEVICES', None)
            if is_windows:
                print("Running in CPU mode (CUDA not available)")
        
        # Set num_workers based on operating system
        # Windows has issues with multiprocessing DataLoader, use single process
        # Linux/WSL2 can use multiple workers for better performance
        if is_windows:
            num_workers = 0
        else:
            # Linux/WSL2: use configured number of threads
            num_workers = self.nthreads
            logger.info(f"[Linux/WSL2] Using {num_workers} workers for data loading")
        
        # Build command arguments
        cmd_args = [
            sys.executable,  # Use current Python interpreter
            script_path,
            '--user-dir', user_dir,
            output_dir_abs,
            '--valid-subset', lmdb_name,
            '--results-path', output_dir_abs,
            '--num-workers', str(num_workers),
        ]
        
        # Only add ddp-backend for non-Windows systems
        # Windows single GPU doesn't need distributed backend
        if not is_windows:
            cmd_args.extend(['--ddp-backend', 'c10d'])
        
        # Windows上自动降低batch_size以避免崩溃
        # Linux/WSL2 可以使用更大的 batch_size
        if is_windows and batch_size > 1:
            logger.warning(f"[WIN] Reducing batch_size from {batch_size} to 1 for Windows compatibility")
            batch_size = 1
        elif (is_linux or is_wsl) and batch_size < 4:
            # Linux/WSL2 上建议使用更大的 batch_size 以提高性能
            logger.info(f"[Linux/WSL2] Using batch_size={batch_size} (can be increased for better performance)")
        
        cmd_args.extend([
            '--batch-size', str(batch_size),
            '--task', 'docking_pose_v2',
            '--loss', 'docking_pose_v2',
            '--arch', 'docking_pose_v2',
            '--conf-size', str(self.conf_size),
            '--dist-threshold', '8.0',
            '--recycling', '4',
            '--path', model_path,
        ])
        
        # Only use FP16 if CUDA is available (FP16 requires CUDA)
        if cuda_available:
            cmd_args.extend([
                '--fp16',
                '--fp16-init-scale', '4',
                '--fp16-scale-window', '256',
            ])
        else:
            # Use CPU mode if CUDA is not available
            cmd_args.append('--cpu')
        
        cmd_args.extend([
            '--log-interval', '50',
            '--log-format', 'simple',
            '--required-batch-size-multiple', '1',
        ])
        
        # Execute inference command with error capture
        try:
            result = subprocess.run(
                cmd_args, 
                env=env, 
                check=True,
                capture_output=True,
                text=True,
                encoding='utf-8',
                errors='replace'  # Handle encoding errors gracefully
            )
            # Log stdout if available
            if result.stdout:
                print("Inference output:")
                print(result.stdout)
        except subprocess.CalledProcessError as e:
            error_msg = f"Subprocess failed with return code {e.returncode}\n"
            error_msg += f"Command: {' '.join(cmd_args)}\n"
            
            if e.stdout:
                error_msg += f"\nSTDOUT:\n{e.stdout}\n"
            if e.stderr:
                error_msg += f"\nSTDERR:\n{e.stderr}\n"
            
            # Additional Windows-specific error information
            if is_windows:
                error_msg += "\nWindows-specific troubleshooting:\n"
                error_msg += "1. Check if CUDA is properly installed and accessible\n"
                error_msg += "2. Verify PyTorch CUDA support: python -c 'import torch; print(torch.cuda.is_available())'\n"
                error_msg += "3. Check GPU drivers are up to date\n"
                error_msg += "4. Try reducing batch_size if memory issues occur\n"
                error_msg += "5. Check Windows Event Viewer for system-level errors\n"
            
            print(error_msg, file=sys.stderr)
            raise RuntimeError(f"Inference subprocess failed: {error_msg}") from e
        except FileNotFoundError as e:
            error_msg = f"Failed to execute command: {cmd_args[0]}\n"
            error_msg += f"Python executable not found: {sys.executable}\n"
            if is_windows:
                error_msg += "Make sure Python is in your PATH\n"
            raise RuntimeError(error_msg) from e
        except Exception as e:
            error_msg = f"Unexpected error during subprocess execution: {e}\n"
            error_msg += f"Command: {' '.join(cmd_args)}\n"
            if is_windows:
                error_msg += "This may be a Windows-specific issue with subprocess execution.\n"
            raise RuntimeError(error_msg) from e
        
        # return results path and lmdb path
        pkl_file = os.path.join(os.path.abspath(output_ligand_dir), lmdb_name + '.pkl')
        lmdb_file = os.path.join(os.path.abspath(output_ligand_dir), lmdb_name+'.lmdb')
        return pkl_file, lmdb_file

    def postprocess(self, output_pkl, output_lmdb, output_ligand_name, output_ligand_dir, input_ligand, input_protein):
        # output the inference results to SDF file
        postprocessor = Processor.build_processors(self.mode, conf_size=self.conf_size)
        mol_list, smi_list, coords_predict_list, holo_coords_list, holo_center_coords_list, prmsd_score_list = postprocessor.postprocess_data_pre(output_pkl, output_lmdb)
        output_ligand_sdf = postprocessor.get_sdf(mol_list, smi_list, coords_predict_list, holo_center_coords_list, prmsd_score_list, output_ligand_name, output_ligand_dir, tta_times=self.conf_size)
        if self.steric_clash_fix:
            output_ligand_sdf = postprocessor.clash_fix(output_ligand_sdf, input_protein, input_ligand)
        return output_ligand_sdf

    def predict_sdf(self, input_protein:str, 
                    input_ligand:str, input_docking_grid:str, 
                    output_ligand_name:str, output_ligand_dir:str, 
                    batch_size:int = 4):
        output_pkl, output_lmdb = self.predict(input_protein, 
                                               input_ligand, 
                                               input_docking_grid, 
                                               output_ligand_name, 
                                               output_ligand_dir, 
                                               batch_size)
        output_sdf = self.postprocess(output_pkl, 
                                      output_lmdb, 
                                      output_ligand_name, 
                                      output_ligand_dir,
                                      input_ligand,
                                      input_protein)

        # return sdf path
        return input_protein, input_ligand, input_docking_grid, output_sdf

    @classmethod
    def build_predictors(cls, model_dir, mode = 'single', 
                         nthreads = 8, conf_size =10, 
                         cluster=False, use_current_ligand_conf=False, steric_clash_fix=False):
        return cls(model_dir, mode, nthreads, conf_size, 
                   cluster, use_current_ligand_conf=use_current_ligand_conf, steric_clash_fix=steric_clash_fix)    
