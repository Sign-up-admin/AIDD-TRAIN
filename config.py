'''
Configuration settings for the data processing and model training script.
'''

CONFIG = {
    # --- Path Settings ---
    'index_file': r'E:/Qinchaojun/AIDD-TRAIN/index/INDEX_general_PL.2020R1.lst',
    'dataset_path': r'E:/Qinchaojun/PDBbind-2025.8.4/P-L/',
    'processed_data_dir': r'E:/Qinchaojun/AIDD-TRAIN/processed_data',

    # --- Model & Training Hyperparameters ---
    'epochs': 20,
    'batch_size': 128,  # Increased for better GPU utilization (was 32)
    'learning_rate': 0.001,
    'train_split': 0.8,
    'dropout_rate': 0.5,

    # --- Hardware & Performance ---
    # Set the number of CPU cores for parallel data preprocessing.
    # Recommended: Set to the number of physical cores for optimal performance.
    'num_workers': 0, # Increased for faster data loading (was 4)
}
