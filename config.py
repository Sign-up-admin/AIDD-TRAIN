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
    'batch_size': 5, # Final adjustment to find the sweet spot
    'learning_rate': 0.0001,
    'train_split': 0.8,
    'dropout_rate': 0.5,

    # --- ViSNet Model Specific ---
    # These parameters control the size and complexity of the ViSNet model.
    # Reducing them can help with CUDA out-of-memory errors.
    'visnet_hidden_channels': 128,  # Reduced from 256
    'visnet_num_layers': 4,         # Reduced from 6
    'visnet_max_neighbors': 16,
    'visnet_cutoff': 8.0, # Cutoff distance for neighbors
    'visnet_num_rbf': 64,           # Reduced from 128

    # --- Hardware & Performance ---
    # Number of CPU cores for the one-time, initial data preprocessing.
    # WARNING: High values can lead to I/O errors on some systems during file saving.
    # If you encounter RuntimeErrors, consider reducing this number.
    'processing_num_workers': 28,

    # Number of CPU cores for the DataLoader during training.
    # This is an I/O-bound task. A smaller number is recommended.
    'loader_num_workers': 12, # Kept at 12 to speed up data loading
}
