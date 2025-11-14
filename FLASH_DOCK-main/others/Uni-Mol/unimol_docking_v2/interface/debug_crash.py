#!/usr/bin/env python3
"""
调试脚本：逐步测试以找出崩溃位置
增强版：添加内存监控、详细错误捕获和逐步测试
"""
import os
import sys
import platform
import logging
import traceback
import torch
import gc
import psutil
import threading
from contextlib import contextmanager
from functools import wraps

# 设置日志
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    stream=sys.stdout,
)
logger = logging.getLogger("debug")

# 内存监控
def get_memory_info():
    """获取当前内存使用情况"""
    try:
        process = psutil.Process(os.getpid())
        mem_info = process.memory_info()
        return {
            'rss_mb': mem_info.rss / 1024 / 1024,
            'vms_mb': mem_info.vms / 1024 / 1024,
        }
    except Exception:
        return {'rss_mb': 0, 'vms_mb': 0}

@contextmanager
def memory_monitor(step_name):
    """内存监控上下文管理器"""
    mem_before = get_memory_info()
    logger.debug(f"[MEM] {step_name} - Before: RSS={mem_before['rss_mb']:.1f}MB, VMS={mem_before['vms_mb']:.1f}MB")
    try:
        yield
    finally:
        gc.collect()
        mem_after = get_memory_info()
        mem_diff = mem_after['rss_mb'] - mem_before['rss_mb']
        logger.debug(f"[MEM] {step_name} - After: RSS={mem_after['rss_mb']:.1f}MB, VMS={mem_after['vms_mb']:.1f}MB, Diff={mem_diff:+.1f}MB")

def test_step(step_name, func, *args, **kwargs):
    """测试单个步骤，带内存监控和详细错误捕获"""
    try:
        logger.info(f"=== Testing: {step_name} ===")
        with memory_monitor(step_name):
            result = func(*args, **kwargs)
        logger.info(f"[OK] {step_name} succeeded")
        return result
    except KeyboardInterrupt:
        logger.warning(f"[INTERRUPTED] {step_name} interrupted by user")
        raise
    except SystemExit as e:
        logger.error(f"[SYSTEM_EXIT] {step_name} caused system exit: {e}")
        raise
    except Exception as e:
        logger.error(f"[FAIL] {step_name} failed: {e}")
        logger.error(f"Error type: {type(e).__name__}")
        logger.error(f"Error args: {e.args}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        
        # 额外的Windows特定错误信息
        if platform.system() == 'Windows':
            if isinstance(e, OSError):
                logger.error(f"[WIN] OSError - errno: {e.errno}, winerror: {getattr(e, 'winerror', 'N/A')}")
            elif isinstance(e, MemoryError):
                logger.error(f"[WIN] MemoryError - 可能是内存不足或内存碎片")
            elif isinstance(e, RuntimeError):
                logger.error(f"[WIN] RuntimeError - 可能是多线程/多进程问题")
        
        raise

def main():
    logger.info("=" * 60)
    logger.info("Debug Script: Investigating Windows Crash (Enhanced)")
    logger.info("=" * 60)
    logger.info(f"Python version: {sys.version}")
    logger.info(f"Platform: {platform.system()} {platform.release()}")
    logger.info(f"Platform details: {platform.platform()}")
    logger.info(f"PyTorch version: {torch.__version__}")
    logger.info(f"CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        logger.info(f"CUDA version: {torch.version.cuda}")
    logger.info(f"Initial memory: {get_memory_info()}")
    
    # Windows特定信息
    if platform.system() == 'Windows':
        logger.info(f"[WIN] Windows version: {platform.win32_ver()}")
        try:
            import ctypes
            kernel32 = ctypes.windll.kernel32
            logger.info(f"[WIN] Maximum path length: {kernel32.GetMaximumPathLength() if hasattr(kernel32, 'GetMaximumPathLength') else 'N/A'}")
        except Exception:
            pass
    
    # 测试路径
    test_output_dir = "./test_output"
    lmdb_path = os.path.join(test_output_dir, "ligand_predict.lmdb")
    
    # 步骤1: 检查LMDB文件
    logger.info("\n" + "=" * 60)
    test_step("Check LMDB file exists", os.path.exists, lmdb_path)
    
    # 步骤2: 测试LMDB读取（增强版）
    logger.info("\n" + "=" * 60)
    try:
        import lmdb
        import pickle
        logger.info(f"LMDB version: {lmdb.__version__ if hasattr(lmdb, '__version__') else 'unknown'}")
        
        def test_lmdb_read():
            """测试LMDB读取，包括环境管理和pickle反序列化"""
            env = None
            try:
                logger.debug("Opening LMDB environment...")
                env = lmdb.open(
                    lmdb_path,
                    subdir=False,
                    readonly=True,
                    lock=False,
                    readahead=False,
                    meminit=False,
                    max_readers=256,
                )
                logger.debug(f"LMDB environment opened: {env}")
                logger.debug(f"LMDB info: map_size={env.info()['map_size']}, max_readers={env.max_readers()}")
                
                with env.begin() as txn:
                    keys = list(txn.cursor().iternext(values=False))
                    logger.info(f"Found {len(keys)} keys in LMDB")
                    
                    if len(keys) > 0:
                        # 尝试读取第一个键
                        first_key = keys[0]
                        logger.debug(f"Reading first key: {first_key}")
                        data = txn.get(first_key)
                        logger.info(f"Successfully read first key: {first_key}, data size: {len(data) if data else 0} bytes")
                        
                        # 测试pickle反序列化
                        if data:
                            try:
                                logger.debug("Testing pickle deserialization...")
                                unpickled = pickle.loads(data)
                                logger.info(f"[OK] Pickle deserialization successful, type: {type(unpickled)}")
                                if isinstance(unpickled, dict):
                                    logger.debug(f"Unpickled dict keys: {list(unpickled.keys())}")
                            except Exception as e:
                                logger.error(f"[FAIL] Pickle deserialization failed: {e}")
                                logger.error(f"Pickle error type: {type(e).__name__}")
                                raise
                
                return len(keys)
            finally:
                if env is not None:
                    try:
                        logger.debug("Closing LMDB environment...")
                        env.close()
                        logger.debug("LMDB environment closed")
                    except Exception as e:
                        logger.warning(f"Error closing LMDB environment: {e}")
        
        num_keys = test_step("Read LMDB file (with pickle test)", test_lmdb_read)
    except Exception as e:
        logger.error(f"LMDB test failed: {e}")
        traceback.print_exc()
        return
    
    # 步骤3: 测试导入unicore模块
    logger.info("\n" + "=" * 60)
    try:
        from unicore import checkpoint_utils, distributed_utils, options, utils
        from unicore.logging import progress_bar
        from unicore import tasks
        logger.info("[OK] Unicore modules imported successfully")
    except Exception as e:
        logger.error(f"✗ Failed to import unicore modules: {e}")
        traceback.print_exc()
        return
    
    # 步骤4: 测试模型加载
    logger.info("\n" + "=" * 60)
    model_path = "../unimol_docking_v2_240517.pt"
    if not os.path.exists(model_path):
        logger.error(f"Model file not found: {model_path}")
        return
    
    try:
        # 创建模拟的args对象
        class MockArgs:
            def __init__(self):
                self.data = test_output_dir
                self.task = 'docking_pose_v2'
                self.arch = 'docking_pose_v2'
                self.loss = 'docking_pose_v2'
                self.user_dir = "../unimol"
                self.distributed_world_size = 1
                self.device_id = 0
                self.cpu = True
                self.fp16 = False
                self.batch_size = 4
                self.num_workers = 0
                self.required_batch_size_multiple = 1
                self.data_buffer_size = 10
                self.seed = 1
                self.valid_subset = "ligand_predict"
                self.results_path = test_output_dir
                self.log_interval = 50
                self.log_format = "simple"
                self.no_progress_bar = False
        
        args = MockArgs()
        
        # 需要先导入任务模块
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'unimol'))
        from unicore import options
        
        # 使用options来解析参数，这样任务会被正确注册
        parser = options.get_validation_parser()
        options.add_model_args(parser)
        
        # 构建命令行参数
        cmd_args = [
            '--user-dir', args.user_dir,
            test_output_dir,
            '--valid-subset', 'ligand_predict',
            '--task', 'docking_pose_v2',
            '--arch', 'docking_pose_v2',
            '--loss', 'docking_pose_v2',
            '--batch-size', str(args.batch_size),
            '--num-workers', str(args.num_workers),
            '--path', model_path,
            '--cpu',
        ]
        
        parsed_args = options.parse_args_and_arch(parser, cmd_args)
        
        # 测试任务设置
        logger.info("Testing task setup...")
        task = tasks.setup_task(parsed_args)
        logger.info("[OK] Task setup successful")
        
        # 测试模型构建
        logger.info("Testing model building...")
        state = checkpoint_utils.load_checkpoint_to_cpu(model_path)
        model = task.build_model(parsed_args)
        model.load_state_dict(state["model"], strict=False)
        logger.info("[OK] Model loaded successfully")
        
        # 测试数据集加载（增强版：逐步测试）
        logger.info("\n" + "=" * 60)
        logger.info("Testing dataset loading (step by step)...")
        
        def test_dataset_loading():
            """逐步测试数据集加载"""
            # 步骤1: 加载数据集
            logger.debug("Step 1: Calling task.load_dataset...")
            task.load_dataset("ligand_predict", combine=False, epoch=1)
            logger.debug("Step 1: task.load_dataset completed")
            
            # 步骤2: 获取数据集对象
            logger.debug("Step 2: Getting dataset object...")
            dataset = task.dataset("ligand_predict")
            logger.debug(f"Step 2: Dataset object obtained: {type(dataset)}")
            
            # 步骤3: 检查数据集长度
            if hasattr(dataset, '__len__'):
                dataset_len = len(dataset)
                logger.info(f"Dataset size: {dataset_len}")
            else:
                logger.warning("Dataset does not have __len__ method")
                dataset_len = None
            
            # 步骤4: 尝试访问第一个数据项（绕过迭代器）
            if dataset_len and dataset_len > 0:
                logger.debug("Step 4: Testing direct dataset access (index 0)...")
                try:
                    first_item = dataset[0]
                    logger.info(f"[OK] Direct dataset access successful, item type: {type(first_item)}")
                    if isinstance(first_item, dict):
                        logger.debug(f"First item keys: {list(first_item.keys())[:10]}...")  # 只显示前10个键
                except Exception as e:
                    logger.error(f"[FAIL] Direct dataset access failed: {e}")
                    logger.error(f"Error type: {type(e).__name__}")
                    logger.error(f"Traceback:\n{traceback.format_exc()}")
                    # 不抛出异常，继续测试
            
            return dataset
        
        dataset = test_step("Dataset loading (step by step)", test_dataset_loading)
        
        # 测试批处理迭代器创建
        logger.info("Testing batch iterator creation...")
        batch_iterator = task.get_batch_iterator(
            dataset=dataset,
            batch_size=parsed_args.batch_size,
            ignore_invalid_inputs=True,
            required_batch_size_multiple=parsed_args.required_batch_size_multiple,
            seed=parsed_args.seed,
            num_shards=1,
            shard_id=0,
            num_workers=parsed_args.num_workers,
            data_buffer_size=parsed_args.data_buffer_size,
        )
        logger.info("[OK] Batch iterator created")
        
        # 测试迭代器初始化
        logger.info("Testing iterator initialization...")
        itr = batch_iterator.next_epoch_itr(shuffle=False)
        logger.info("[OK] Iterator initialized")
        
        # 测试progress_bar创建
        logger.info("Testing progress bar creation...")
        progress = progress_bar.progress_bar(
            itr,
            log_format=parsed_args.log_format,
            log_interval=parsed_args.log_interval,
            prefix="test",
            default_log_format="simple",
        )
        logger.info("[OK] Progress bar created")
        
        # 测试第一次迭代 - 这是崩溃发生的地方（增强版）
        logger.info("\n" + "=" * 60)
        logger.info("Testing first iteration (this is where crash occurs)...")
        logger.info("=" * 60)
        
        try:
            # 步骤1: 尝试创建迭代器
            logger.info("Step 1: Creating iterator from progress...")
            try:
                progress_iter = iter(progress)
                logger.info("[OK] Iterator created successfully")
            except Exception as e:
                logger.error(f"[FAIL] Failed to create iterator: {e}")
                logger.error(f"Error type: {type(e).__name__}")
                logger.error(f"Traceback:\n{traceback.format_exc()}")
                raise
            
            # 步骤2: 尝试获取第一个样本
            logger.info("Step 2: Attempting to get first sample from iterator...")
            logger.debug("Memory before next():")
            mem_before = get_memory_info()
            logger.debug(f"  RSS: {mem_before['rss_mb']:.1f}MB, VMS: {mem_before['vms_mb']:.1f}MB")
            
            try:
                first_sample = next(progress_iter)
                logger.info(f"[OK] First sample retrieved successfully")
                
                mem_after = get_memory_info()
                logger.debug(f"Memory after next(): RSS: {mem_after['rss_mb']:.1f}MB, VMS: {mem_after['vms_mb']:.1f}MB")
                
                if isinstance(first_sample, dict):
                    logger.info(f"Sample keys: {list(first_sample.keys())}")
                    # 检查每个键的数据类型和形状
                    for key, value in list(first_sample.items())[:5]:  # 只检查前5个键
                        if hasattr(value, 'shape'):
                            logger.debug(f"  {key}: shape={value.shape}, dtype={getattr(value, 'dtype', 'N/A')}")
                        elif isinstance(value, (list, tuple)):
                            logger.debug(f"  {key}: type={type(value)}, len={len(value)}")
                        else:
                            logger.debug(f"  {key}: type={type(value)}")
                else:
                    logger.info(f"Sample type: {type(first_sample)}")
            except StopIteration:
                logger.warning("Iterator is empty (StopIteration)")
            except Exception as e:
                logger.error(f"[FAIL] Failed to get first sample: {e}")
                logger.error(f"Error type: {type(e).__name__}")
                logger.error(f"Error args: {e.args}")
                
                # Windows特定错误分析
                if platform.system() == 'Windows':
                    if isinstance(e, RuntimeError):
                        logger.error("[WIN] RuntimeError - 可能是多线程问题或资源访问冲突")
                    elif isinstance(e, OSError):
                        logger.error(f"[WIN] OSError - errno: {e.errno}, winerror: {getattr(e, 'winerror', 'N/A')}")
                    elif isinstance(e, MemoryError):
                        logger.error("[WIN] MemoryError - 内存不足或内存访问错误")
                    elif "pickle" in str(e).lower():
                        logger.error("[WIN] Pickle相关错误 - 可能是序列化/反序列化问题")
                    elif "lmdb" in str(e).lower():
                        logger.error("[WIN] LMDB相关错误 - 可能是数据库访问问题")
                
                logger.error(f"Full traceback:\n{traceback.format_exc()}")
                raise
        except Exception as e:
            logger.error(f"First iteration test failed: {e}")
            raise
        
    except Exception as e:
        logger.error(f"Test failed: {e}")
        traceback.print_exc()
        return
    
    logger.info("\n" + "=" * 60)
    logger.info("All tests passed!")
    logger.info("=" * 60)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        traceback.print_exc()
        sys.exit(1)

