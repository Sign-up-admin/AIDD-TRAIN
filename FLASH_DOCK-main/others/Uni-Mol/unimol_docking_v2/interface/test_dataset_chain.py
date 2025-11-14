#!/usr/bin/env python3
"""
逐步测试数据集处理链，定位问题所在的包装器
"""
import os
import sys
import platform
import logging
import traceback

# 设置日志
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG,
    stream=sys.stdout,
)
logger = logging.getLogger("test_dataset_chain")

def test_dataset_chain_step_by_step():
    """逐步测试数据集处理链"""
    logger.info("=" * 60)
    logger.info("Step-by-Step Dataset Chain Test")
    logger.info("=" * 60)
    
    # 设置路径
    test_output_dir = "./test_output"
    model_path = "../unimol_docking_v2_240517.pt"
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'unimol'))
    
    try:
        # 导入必要的模块
        from unicore import checkpoint_utils, options
        from unicore import tasks
        
        # 解析参数
        parser = options.get_validation_parser()
        options.add_model_args(parser)
        cmd_args = [
            '--user-dir', "../unimol",
            test_output_dir,
            '--valid-subset', 'ligand_predict',
            '--task', 'docking_pose_v2',
            '--arch', 'docking_pose_v2',
            '--loss', 'docking_pose_v2',
            '--batch-size', '1',  # 使用小batch size
            '--num-workers', '0',
            '--path', model_path,
            '--cpu',
        ]
        parsed_args = options.parse_args_and_arch(parser, cmd_args)
        
        # 设置任务
        logger.info("Step 0: Setting up task...")
        task = tasks.setup_task(parsed_args)
        logger.info("[OK] Task setup successful")
        
        # 开始逐步测试数据集处理链
        logger.info("\n" + "=" * 60)
        logger.info("Step 1: Loading base LMDBDataset...")
        logger.info("=" * 60)
        
        from unimol.data import LMDBDataset
        data_path = os.path.join(parsed_args.data, "ligand_predict.lmdb")
        base_dataset = LMDBDataset(data_path)
        logger.info(f"[OK] Base dataset loaded, size: {len(base_dataset)}")
        
        # 测试基础数据集访问
        try:
            item = base_dataset[0]
            logger.info(f"[OK] Base dataset access successful, type: {type(item)}")
        except Exception as e:
            logger.error(f"[FAIL] Base dataset access failed: {e}")
            raise
        
        # 步骤2: TTADockingPoseDataset
        logger.info("\n" + "=" * 60)
        logger.info("Step 2: Adding TTADockingPoseDataset wrapper...")
        logger.info("=" * 60)
        
        from unimol.data import TTADockingPoseDataset
        try:
            tta_dataset = TTADockingPoseDataset(
                base_dataset,
                'atoms', 'coordinates', 'pocket_atoms', 'pocket_coordinates',
                'holo_coordinates', 'holo_pocket_coordinates',
                True, parsed_args.conf_size
            )
            logger.info(f"[OK] TTADockingPoseDataset created, size: {len(tta_dataset)}")
            
            # 测试访问
            item = tta_dataset[0]
            logger.info(f"[OK] TTADockingPoseDataset access successful")
        except Exception as e:
            logger.error(f"[FAIL] TTADockingPoseDataset failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤3: RemoveHydrogenPocketDataset
        logger.info("\n" + "=" * 60)
        logger.info("Step 3: Adding RemoveHydrogenPocketDataset wrapper...")
        logger.info("=" * 60)
        
        from unimol.data import RemoveHydrogenPocketDataset
        try:
            dataset = RemoveHydrogenPocketDataset(
                tta_dataset, 'pocket_atoms', 'pocket_coordinates',
                'holo_pocket_coordinates', True, True
            )
            logger.info(f"[OK] RemoveHydrogenPocketDataset (1) created")
            
            item = dataset[0]
            logger.info(f"[OK] RemoveHydrogenPocketDataset (1) access successful")
        except Exception as e:
            logger.error(f"[FAIL] RemoveHydrogenPocketDataset (1) failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤4: CroppingPocketDataset
        logger.info("\n" + "=" * 60)
        logger.info("Step 4: Adding CroppingPocketDataset wrapper...")
        logger.info("=" * 60)
        
        from unimol.data import CroppingPocketDataset
        try:
            dataset = CroppingPocketDataset(
                dataset, task.seed, 'pocket_atoms', 'pocket_coordinates',
                'holo_pocket_coordinates', parsed_args.max_pocket_atoms
            )
            logger.info(f"[OK] CroppingPocketDataset created")
            
            item = dataset[0]
            logger.info(f"[OK] CroppingPocketDataset access successful")
        except Exception as e:
            logger.error(f"[FAIL] CroppingPocketDataset failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤5: RemoveHydrogenPocketDataset (2)
        logger.info("\n" + "=" * 60)
        logger.info("Step 5: Adding RemoveHydrogenPocketDataset wrapper (2)...")
        logger.info("=" * 60)
        
        try:
            dataset = RemoveHydrogenPocketDataset(
                dataset, 'atoms', 'coordinates', 'holo_coordinates', True, True
            )
            logger.info(f"[OK] RemoveHydrogenPocketDataset (2) created")
            
            item = dataset[0]
            logger.info(f"[OK] RemoveHydrogenPocketDataset (2) access successful")
        except Exception as e:
            logger.error(f"[FAIL] RemoveHydrogenPocketDataset (2) failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤6: NormalizeDataset
        logger.info("\n" + "=" * 60)
        logger.info("Step 6: Adding NormalizeDataset wrappers...")
        logger.info("=" * 60)
        
        from unimol.data import NormalizeDataset
        try:
            apo_dataset = NormalizeDataset(dataset, 'coordinates')
            apo_dataset = NormalizeDataset(apo_dataset, 'pocket_coordinates')
            logger.info(f"[OK] NormalizeDataset wrappers created")
            
            item = apo_dataset[0]
            logger.info(f"[OK] NormalizeDataset access successful")
        except Exception as e:
            logger.error(f"[FAIL] NormalizeDataset failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤7: ReAlignLigandDataset
        logger.info("\n" + "=" * 60)
        logger.info("Step 7: Adding ReAlignLigandDataset wrapper...")
        logger.info("=" * 60)
        
        from unimol.data import ReAlignLigandDataset
        try:
            apo_dataset = ReAlignLigandDataset(
                dataset, 'coordinates', 'pocket_coordinates'
            )
            logger.info(f"[OK] ReAlignLigandDataset created")
            
            item = apo_dataset[0]
            logger.info(f"[OK] ReAlignLigandDataset access successful")
        except Exception as e:
            logger.error(f"[FAIL] ReAlignLigandDataset failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        # 步骤8: 测试完整的数据集加载（使用task.load_dataset）
        logger.info("\n" + "=" * 60)
        logger.info("Step 8: Testing complete dataset loading via task.load_dataset...")
        logger.info("=" * 60)
        
        try:
            task.load_dataset("ligand_predict", combine=False, epoch=1)
            final_dataset = task.dataset("ligand_predict")
            logger.info(f"[OK] Complete dataset loaded, size: {len(final_dataset) if hasattr(final_dataset, '__len__') else 'unknown'}")
            
            # 尝试访问第一个项目
            item = final_dataset[0]
            logger.info(f"[OK] Final dataset access successful, type: {type(item)}")
            if isinstance(item, dict):
                logger.info(f"  Keys: {list(item.keys())[:10]}...")
        except Exception as e:
            logger.error(f"[FAIL] Complete dataset loading failed: {e}")
            logger.error(f"Traceback:\n{traceback.format_exc()}")
            raise
        
        logger.info("\n" + "=" * 60)
        logger.info("All dataset chain tests passed!")
        logger.info("=" * 60)
        
    except Exception as e:
        logger.error(f"\nDataset chain test failed: {e}")
        logger.error(f"Traceback:\n{traceback.format_exc()}")
        raise

def main():
    logger.info("=" * 60)
    logger.info("Dataset Chain Step-by-Step Test")
    logger.info("=" * 60)
    logger.info(f"Platform: {platform.system()} {platform.release()}")
    logger.info(f"Python version: {sys.version}")
    
    try:
        test_dataset_chain_step_by_step()
    except Exception as e:
        logger.error(f"\nTest suite failed: {e}")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

