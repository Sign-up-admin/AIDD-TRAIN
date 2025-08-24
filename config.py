'''
Configuration settings for the data processing and model training script.
'''

import os

# --- Development Mode Switch ---
# Select the operating mode. This single switch controls all major hyperparameters.
# - 'smoke_test': A minimal configuration to ensure the entire pipeline runs without errors.
#                 Use this for initial setup and debugging code logic. Runs in minutes.
# - 'prototyping': A lightweight configuration for rapid experimentation and idea validation.
#                  Aims for a fast feedback loop (minutes per epoch) on smaller datasets.
# - 'validation': A medium-sized configuration to validate model performance before full-scale training.
#                 Optimized for GPUs with 6-8GB of VRAM.
# - 'production': The full-scale configuration for generating final results. Tuned for GPUs
#                 with >8GB of VRAM.
DEVELOPMENT_MODE = 'prototyping'  # Options: 'smoke_test', 'prototyping', 'validation', 'production'


# --- Mode-Specific Hyperparameter Sets ---
# Defines the core settings for each development mode.
MODES = {
    'smoke_test': {
        'epochs': 2,
        'batch_size': 1,
        'gradient_accumulation_steps': 2,
        'loader_num_workers': 2,
        'visnet_hidden_channels': 8,
        'visnet_num_layers': 1,
        'visnet_num_rbf': 8,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 1000,
    },
    # Recommended for rapid idea testing. Aims for epochs in minutes.
    'prototyping': {
        'epochs': 10,
        'batch_size': 2,
        'gradient_accumulation_steps': 4,
        'loader_num_workers': 4,
        'visnet_hidden_channels': 32,
        'visnet_num_layers': 2,
        'visnet_num_rbf': 32,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 5000,
    },
    # Recommended for GPUs with ~6-8GB VRAM (e.g., RTX 3060).
    'validation': {
        'epochs': 20,
        'batch_size': 4,
        'gradient_accumulation_steps': 8,
        'loader_num_workers': 4,
        'visnet_hidden_channels': 48,
        'visnet_num_layers': 3,
        'visnet_num_rbf': 48,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 10000,
    },
    # Tuned for GPUs with >8GB VRAM. May run on 6GB but can be slow or unstable.
    'production': {
        'epochs': 100,
        'batch_size': 8,
        'gradient_accumulation_steps': 8,
        'loader_num_workers': 4,
        'visnet_hidden_channels': 64,
        'visnet_num_layers': 4,
        'visnet_num_rbf': 64,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 10000,
    }
}

# --- Base Configuration ---
# Contains settings that are common across all modes.
CONFIG = {
    # --- Development Mode ---
    'development_mode': DEVELOPMENT_MODE,

    # --- Data Source ---
    'force_data_reprocessing': False,

    # --- Data Processing ---
    'data_processing_mode': 'strict',
    'max_atoms': 10000,
    'failed_cases_dir': r'failed_cases',

    # --- Path Settings ---
    'index_file': r'index/INDEX_general_PL.2020R1.lst',
    'dataset_path': r'PDBbind-2025.8.4/P-L/',
    'processed_data_dir': r'processed_data',

    # --- Model & Training Hyperparameters ---
    'learning_rate': 0.0001,
    'train_split': 0.8,
    'dropout_rate': 0.5,
    'gradient_clip_val': 1.0,

    # --- ViSNet Model Specific ---
    'max_num_neighbors': 32,
    'visnet_cutoff': 5.0,
    'visnet_lmax': 2,
    'visnet_vecnorm_type': 'max_min',

    # --- Hardware & Performance ---
    'processing_num_workers': 16,
    'profile': False,
    'manual_save_trigger_file': 'SAVE_NOW.flag',
    'save_every_n_batches': 1000,

    # --- Reproducibility ---
    'seed': 42,

    # --- Debugging ---
    'debug_mode': True,
}

# --- Merge Mode-Specific Settings ---
mode_settings = MODES.get(DEVELOPMENT_MODE)
if mode_settings is None:
    raise ValueError(f"Invalid DEVELOPMENT_MODE: '{DEVELOPMENT_MODE}'. Please choose from {list(MODES.keys())}.")

CONFIG.update(mode_settings)

# --- Post-Processing and Validation ---
CONFIG['run_name'] = f"visnet_{DEVELOPMENT_MODE}_lr{CONFIG['learning_rate']}_bs{CONFIG['batch_size']}x{CONFIG['gradient_accumulation_steps']}"

if not os.path.exists(CONFIG['failed_cases_dir']):
    os.makedirs(CONFIG['failed_cases_dir'])
