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
    'batch_size': 128,
    'learning_rate': 0.001,
    'train_split': 0.8,
    'dropout_rate': 0.5,

    # --- Hardware & Performance ---
    # Number of CPU cores for the one-time, initial data preprocessing.
    # A high value caused memory crashes. A safer, more conservative value is set.
    # You can try tuning this value up (e.g., 12, 16) if you have sufficient RAM.
    'processing_num_workers': 8,

    # Number of CPU cores for the DataLoader during training.
    # This is an I/O-bound task. A smaller number is recommended.
    'loader_num_workers': 4,
}
