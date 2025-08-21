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
    # The actual batch size used by the hardware. The effective batch size will be this value multiplied by gradient_accumulation_steps.
    'batch_size': 4, 
    'learning_rate': 0.0001,
    'train_split': 0.8,
    'dropout_rate': 0.5,
    # Gradient Accumulation: Simulate a larger batch size to stabilize training.
    # Effective Batch Size = batch_size * gradient_accumulation_steps
    'gradient_accumulation_steps': 4, # Results in an effective batch size of 4 * 4 = 16

    # --- ViSNet Model Specific ---
    # These parameters control the size and complexity of the ViSNet model.
    'visnet_hidden_channels': 128,
    'visnet_num_layers': 4,
    'visnet_max_neighbors': 16,
    'visnet_cutoff': 8.0,
    'visnet_num_rbf': 64,

    # --- Hardware & Performance ---
    'processing_num_workers': 28,
    'loader_num_workers': 0,
}
