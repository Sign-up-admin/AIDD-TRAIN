'''
Configuration settings for the data processing and model training script.
'''

CONFIG = {
    # --- Data Source ---
    # Force the reprocessing of the entire dataset, ignoring any cached files.
    # This is useful after changing data processing logic in `src/data_processing.py`.
    # The script will automatically set this to False after a successful run.
    'force_data_reprocessing': True,

    # --- Data Processing ---
    # 'strict': Skips any PDB file with unresolvable issues (e.g., atom overlaps). This is the safest mode.
    # 'permissive': Attempts to apply smart heuristics to repair problematic PDB files. This may increase data yield
    #               but introduces a small risk of altering molecular structures.
    'data_processing_mode': 'strict',
    # Directory to save detailed reports and file snapshots for PDBs that fail processing.
    'failed_cases_dir': r'failed_cases',

    # --- Path Settings ---
    'index_file': r'index/INDEX_general_PL.2020R1.lst',
    'dataset_path': r'PDBbind-2025.8.4/P-L/',
    'processed_data_dir': r'processed_data',

    # --- Model & Training Hyperparameters ---
    'epochs': 20,
    # The actual batch size used by the hardware. The effective batch size will be this value multiplied by gradient_accumulation_steps.
    'batch_size': 1, 
    'learning_rate': 0.0001,
    'train_split': 0.8,
    'dropout_rate': 0.5,
    # Gradient Accumulation: Simulate a larger batch size to stabilize training.
    # Effective Batch Size = batch_size * gradient_accumulation_steps
    'gradient_accumulation_steps': 16, # Results in an effective batch size of 1 * 16 = 16

    # --- ViSNet Model Specific ---
    # These parameters control the size and complexity of the ViSNet model.
    'visnet_hidden_channels': 128, # 增加模型宽度以获得更多特征容量
    'visnet_num_layers': 4,      # 增加模型深度以学习更复杂的相互作用
    'max_num_neighbors': 32,     # 允许模型看到更丰富的原子邻域环境 (ViSNet 默认值)
    'visnet_cutoff': 8.0,
    'visnet_num_rbf': 64,      # 增加距离表示的分辨率

    # --- Hardware & Performance ---
    'processing_num_workers': 28,
    'loader_num_workers': 0,
    # --- Performance Profiling ---
    # Set to True to run the profiler on the first epoch to identify bottlenecks.
    'profile': False,
    # --- Manual Save Trigger ---
    # Create a file with this name in the project root to trigger an immediate checkpoint save.
    'manual_save_trigger_file': 'SAVE_NOW.flag',
    # --- Periodic Checkpointing ---
    # Save a checkpoint every N batches to prevent losing progress on long runs.
    # Set to 0 to disable and only save at the end of each epoch.
    # A value of 1000 is a reasonable starting point.
    'save_every_n_batches': 1000,
    # --- Reproducibility ---
    # Set a random seed for consistent data splitting and model initialization.
    'seed': 42,

    # --- Debugging ---
    # Set to True to save the batch that causes a NaN/Inf and stop training.
    'debug_mode': True,
}
