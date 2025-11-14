#!/usr/bin/env python3 -u
# Copyright (c) DP Techonology, Inc. and its affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import logging
import os
import platform
import sys
import pickle
import torch
from unicore import checkpoint_utils, distributed_utils, options, utils
from unicore.logging import progress_bar
from unicore import tasks

logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=os.environ.get("LOGLEVEL", "INFO").upper(),
    stream=sys.stdout,
)
logger = logging.getLogger("unimol.inference")


def main(args):

    assert (
        args.batch_size is not None
    ), "Must specify batch size either with --batch-size"

    # Check CUDA availability
    cuda_available = torch.cuda.is_available()
    is_windows = platform.system() == 'Windows'
    
    if is_windows and not cuda_available:
        logger.warning("CUDA is not available. The model will run on CPU, which may be very slow.")
        logger.warning("If you have a GPU, please ensure CUDA and PyTorch with CUDA support are properly installed.")
    
    use_cuda = cuda_available and not args.cpu
    use_fp16 = args.fp16 and use_cuda

    # Safely set CUDA device with default value
    if use_cuda:
        try:
            # Set default device_id to 0 if not provided
            device_id = getattr(args, 'device_id', 0)
            if device_id is None:
                device_id = 0
            
            # Verify device is available
            if device_id >= torch.cuda.device_count():
                logger.error(f"CUDA device {device_id} is not available. Available devices: {torch.cuda.device_count()}")
                raise RuntimeError(f"CUDA device {device_id} is not available")
            
            torch.cuda.set_device(device_id)
            logger.info(f"Using CUDA device {device_id}: {torch.cuda.get_device_name(device_id)}")
        except Exception as e:
            logger.error(f"Failed to initialize CUDA device: {e}")
            if is_windows:
                logger.error("This may be a Windows-specific CUDA initialization issue.")
                logger.error("Try running with --cpu flag to use CPU instead.")
            raise

    if args.distributed_world_size > 1:
        data_parallel_world_size = distributed_utils.get_data_parallel_world_size()
        data_parallel_rank = distributed_utils.get_data_parallel_rank()
    else:
        data_parallel_world_size = 1
        data_parallel_rank = 0

    # Load model
    try:
        logger.info("loading model(s) from {}".format(args.path))
        if not os.path.exists(args.path):
            raise FileNotFoundError(f"Model file not found: {args.path}")
        
        state = checkpoint_utils.load_checkpoint_to_cpu(args.path)
        logger.info("Model checkpoint loaded successfully")
        
        task = tasks.setup_task(args)
        model = task.build_model(args)
        model.load_state_dict(state["model"], strict=False)
        logger.info("Model state dict loaded successfully")
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        logger.error(f"Model path: {args.path}")
        raise

    # Move models to GPU
    try:
        if use_fp16:
            logger.info("Converting model to FP16")
            model.half()
        if use_cuda:
            logger.info("Moving model to CUDA device")
            model.cuda()
            # Verify model is on GPU
            if next(model.parameters()).is_cuda:
                logger.info("Model successfully moved to CUDA")
            else:
                logger.warning("Model may not be on CUDA despite cuda() call")
    except Exception as e:
        logger.error(f"Failed to move model to GPU: {e}")
        if is_windows:
            logger.error("This may be a Windows-specific CUDA issue. Try running with --cpu flag.")
        raise

    # Print args
    logger.info(args)

    # Build loss
    loss = task.build_loss(args)
    loss.eval()

    for subset in args.valid_subset.split(","):
        try:
            task.load_dataset(subset, combine=False, epoch=1)
            dataset = task.dataset(subset)
        except KeyError:
            raise Exception("Cannot find dataset: " + subset)

        if not os.path.exists(args.results_path):
            os.makedirs(args.results_path)
        save_path = os.path.join(args.results_path, subset + ".pkl")
        
        # Initialize data iterator with error handling
        try:
            logger.info(f"Initializing data iterator for subset: {subset}")
            logger.info(f"Dataset size: {len(dataset) if hasattr(dataset, '__len__') else 'unknown'}")
            logger.info(f"Batch size: {args.batch_size}, num_workers: {args.num_workers}")
            
            batch_iterator = task.get_batch_iterator(
                dataset=dataset,
                batch_size=args.batch_size,
                ignore_invalid_inputs=True,
                required_batch_size_multiple=args.required_batch_size_multiple,
                seed=args.seed,
                num_shards=data_parallel_world_size,
                shard_id=data_parallel_rank,
                num_workers=args.num_workers,
                data_buffer_size=args.data_buffer_size,
            )
            logger.info("Batch iterator created successfully")
            
            itr = batch_iterator.next_epoch_itr(shuffle=False)
            logger.info("Data iterator initialized successfully")
        except Exception as e:
            logger.error(f"Failed to initialize data iterator: {e}")
            logger.error(f"Error type: {type(e).__name__}")
            if is_windows:
                logger.error("This may be a Windows-specific issue with data loading.")
                logger.error("Try reducing batch_size or check LMDB file integrity.")
            raise
        try:
            logger.info("Creating progress bar...")
            progress = progress_bar.progress_bar(
                itr,
                log_format=args.log_format,
                log_interval=args.log_interval,
                prefix=f"valid on '{subset}' subset",
                default_log_format=("tqdm" if not args.no_progress_bar else "simple"),
            )
            logger.info("Progress bar created successfully")
        except Exception as e:
            logger.error(f"Failed to create progress bar: {e}")
            raise
        
        log_outputs = []
        try:
            logger.info("Starting inference loop...")
            
            # Windows特定：跳过预测试，直接进入循环
            # 预测试可能导致崩溃，因为访问违规发生在C扩展层面，Python无法捕获
            # Linux/WSL2 上不需要跳过预测试，可以正常使用
            if is_windows:
                logger.info("[WIN] Skipping pre-test (may cause crashes on Windows)")
                logger.info("[WIN] Proceeding directly to inference loop...")
                import gc
                gc.collect()  # 清理内存
            else:
                logger.info("[Linux/WSL2] Using standard inference loop (no pre-test skip needed)")
            
            for i, sample in enumerate(progress):
                try:
                    # 第一次迭代的特殊处理（这是崩溃常发生的地方）
                    if i == 0:
                        logger.info("=" * 60)
                        logger.info("Processing first sample (this is where crashes often occur)...")
                        logger.info("=" * 60)
                        if is_windows:
                            logger.info("[WIN] Windows-specific checks:")
                            logger.info(f"  - num_workers: {args.num_workers}")
                            logger.info(f"  - batch_size: {args.batch_size}")
                            logger.info(f"  - CUDA: {use_cuda}")
                            import gc
                            gc.collect()  # 强制垃圾回收
                            logger.info("  - Garbage collection performed")
                        else:
                            logger.info("[Linux/WSL2] First sample processing:")
                            logger.info(f"  - num_workers: {args.num_workers}")
                            logger.info(f"  - batch_size: {args.batch_size}")
                            logger.info(f"  - CUDA: {use_cuda}")
                    
                    logger.debug(f"Processing sample {i}...")
                    sample = utils.move_to_cuda(sample) if use_cuda else sample
                    if len(sample) == 0:
                        logger.warning(f"Sample {i} is empty, skipping")
                        continue
                    logger.debug(f"Running valid_step for sample {i}...")
                    _, _, log_output = task.valid_step(sample, model, loss, test=True)
                    progress.log({}, step=i)
                    log_outputs.append(log_output)
                    logger.debug(f"Sample {i} processed successfully")
                    
                    if i == 0:
                        logger.info("[OK] First sample processed successfully!")
                        logger.info("=" * 60)
                except Exception as e:
                    logger.error(f"Error processing sample {i}: {e}")
                    logger.error(f"Error type: {type(e).__name__}")
                    import traceback
                    logger.error(f"Traceback: {traceback.format_exc()}")
                    if is_windows:
                        logger.error("=" * 60)
                        logger.error("[WIN] Windows-specific error analysis:")
                        logger.error("This may be a Windows-specific error.")
                        
                        # 详细的错误分析
                        if isinstance(e, RuntimeError):
                            error_msg = str(e).lower()
                            if "pickle" in error_msg:
                                logger.error("  - Pickle-related error detected")
                                logger.error("  - Possible causes: data serialization issue, encoding problem")
                            elif "lmdb" in error_msg:
                                logger.error("  - LMDB-related error detected")
                                logger.error("  - Possible causes: database access conflict, file locking issue")
                            elif "memory" in error_msg or "out of memory" in error_msg:
                                logger.error("  - Memory-related error detected")
                                logger.error("  - Possible causes: insufficient memory, memory fragmentation")
                            else:
                                logger.error("  - RuntimeError: possible causes include:")
                                logger.error("    * Multi-threading issues")
                                logger.error("    * Resource access conflicts")
                                logger.error("    * Data loading problems")
                        elif isinstance(e, OSError):
                            logger.error(f"  - OSError: errno={e.errno}, winerror={getattr(e, 'winerror', 'N/A')}")
                            logger.error("  - Possible causes: file access issue, permission problem")
                        elif isinstance(e, MemoryError):
                            logger.error("  - MemoryError: insufficient memory or memory access violation")
                        elif isinstance(e, KeyError):
                            logger.error("  - KeyError: missing data key, possible data format issue")
                        
                        if use_cuda:
                            logger.error("  - CUDA is enabled, this may be a CUDA-related error on Windows")
                            logger.error("  - Try running with --cpu flag to use CPU instead")
                        
                        logger.error("Suggested solutions:")
                        logger.error("  1. Reduce batch_size (e.g., --batch-size 1)")
                        logger.error("  2. Ensure num_workers=0 (already set)")
                        logger.error("  3. Check LMDB file integrity")
                        logger.error("  4. Try running with --cpu flag")
                        logger.error("  5. Check Windows Event Viewer for system-level errors")
                        logger.error("=" * 60)
                    raise
            
            # Save results
            try:
                pickle.dump(log_outputs, open(save_path, "wb"))
                logger.info(f"Results saved to {save_path}")
            except Exception as e:
                logger.error(f"Failed to save results to {save_path}: {e}")
                raise
        except Exception as e:
            logger.error(f"Error during inference: {e}")
            raise
        
        logger.info("Done inference! ")
    return None


def cli_main():
    try:
        parser = options.get_validation_parser()
        options.add_model_args(parser)
        args = options.parse_args_and_arch(parser)

        # On Windows, directly call main to avoid distributed_utils issues
        # distributed_utils.call_main may cause crashes on Windows
        is_windows = platform.system() == 'Windows'
        if is_windows:
            # Ensure distributed_world_size is 1 for single GPU on Windows
            args.distributed_world_size = 1
            # Set default device_id if not provided
            if not hasattr(args, 'device_id') or args.device_id is None:
                args.device_id = 0
            main(args)
        else:
            distributed_utils.call_main(args, main)
    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Fatal error: {e}", exc_info=True)
        if platform.system() == 'Windows':
            logger.error("If this is a CUDA-related error on Windows, try:")
            logger.error("1. Ensure CUDA and PyTorch with CUDA support are properly installed")
            logger.error("2. Check that your GPU drivers are up to date")
            logger.error("3. Try running with --cpu flag to use CPU instead")
        sys.exit(1)


if __name__ == "__main__":
    cli_main()
