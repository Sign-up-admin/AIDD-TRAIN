"""
 Configuration settings for the data processing and model training script.
 """
 
import os
 
from compass.optimizer.utils import load_optimized_settings
 
# --- Execution Mode ---
# Selects the operational mode, controlling all major hyperparameters.
# - 'smoke_test': A minimal configuration to ensure the entire pipeline runs without errors.
#                 Ideal for initial setup and debugging code logic, runs in minutes.
# - 'prototyping': A lightweight setup for rapid experimentation and idea validation.
#                  Aims for a fast feedback loop (a few minutes per epoch) on a smaller dataset.
# - 'validation': A medium-scale configuration for validating model performance before full-scale training.
#                 Optimized for GPUs with 6-8GB of VRAM.
# - 'production': A comprehensive configuration for generating final results.
#                 Tuned for GPUs with more than 8GB of VRAM.
#
# Options: 'smoke_test', 'prototyping', 'validation', 'production'
EXECUTION_MODE = 'validation'
 
 
# --- Hyperparameter Sets for Each Mode ---
# Defines the core settings for each execution mode. Parameters can be adjusted
# based on hardware and experimental needs.
MODES = {
    # 'smoke_test': For quickly verifying that the code logic and pipeline are functional.
    'smoke_test': {
        'epochs': 2,  # Total number of training epochs.
        'batch_size': 1,  # Batch size per GPU.
 
        'profile': False,  # Whether to enable performance profiling.
        'max_atoms': 1000,  # Maximum number of atoms supported in this mode.
    },
    # 'prototyping': Recommended for rapid idea validation, targeting a few minutes per epoch.
    'prototyping': {
        'epochs': 50, # Total number of training epochs.
        'batch_size': 2, # Batch size per GPU. Adjust based on VRAM.
 
        'profile': False, # Whether to enable performance profiling.
        'max_atoms': 5000, # Maximum number of atoms supported in this mode.
    },
    # 'validation': Recommended for GPUs with ~6-8GB of VRAM (e.g., RTX 3060).
    'validation': {
        'epochs': 20,
        'batch_size': 4,
 
        'profile': False,
        'max_atoms': 10000,
    },
    # 'production': Tuned for GPUs with >8GB of VRAM. May be slow or unstable on 6GB cards.
    'production': {
        'epochs': 100,
        'batch_size': 8,
 
    }
}
 
def get_config(mode_name=EXECUTION_MODE):
    """
    Generates the final configuration dictionary by merging base settings
    with mode-specific overrides and optimized hardware profiles.
    """
    # --- Base Configuration ---
    # Contains settings common to all modes. These are the defaults and can be
    # overridden by the mode-specific settings.
    config = {
        # --- Execution Mode ---
        'execution_mode': mode_name,
 
        # --- Data Source & Processing ---
        'force_data_reprocessing': False,
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
        'gradient_accumulation_steps': 1, # Number of steps to accumulate gradients before updating weights.

        # --- ViSNet Model Specifics ---
        # Note: ViSNet's attention mechanism requires `visnet_hidden_channels` to be
        # divisible by `visnet_num_heads` (default is 8).
        'max_num_neighbors': 32,
        'visnet_cutoff': 5.0,
        'visnet_lmax': 2,
        'visnet_vecnorm_type': 'max_min',
 
        # --- Hardware & Performance ---
        'processing_num_workers': os.cpu_count() or 16,
        'profile': False,
        'manual_save_trigger_file': 'SAVE_NOW.flag',
        'save_every_n_batches': 100,
 
        # --- Reproducibility ---
        'seed': 42,
 
        # --- Debugging ---
        'debug_mode': True,
    }
 
    # --- Merge Mode-Specific Settings ---
    mode_settings = MODES.get(mode_name)
    if mode_settings is None:
        raise ValueError(f"Invalid EXECUTION_MODE: '{mode_name}'. Please choose from {list(MODES.keys())}.")
    config.update(mode_settings)
 
    # --- Apply Hardware-Optimized Settings ---
    # This will override any conflicting settings with values from the optimizer.
    load_optimized_settings(config)
 
    # --- Post-processing and Validation ---
    # Create a unique run name for logging and checkpointing.
    config['run_name'] = f"visnet_{mode_name}_lr{config['learning_rate']}_bs{config['batch_size']}x{config['gradient_accumulation_steps']}"
 
    # Create dynamic directories based on the run name.
    run_name = config['run_name']
    config['log_dir'] = os.path.join('logs', run_name)
    config['checkpoint_dir'] = os.path.join('checkpoints', run_name)
 
    # Ensure the failed cases directory exists.
    if not os.path.exists(config['failed_cases_dir']):
        os.makedirs(config['failed_cases_dir'])
 
    return config
 
# --- Global Config ---
# The main configuration object used throughout the application.
CONFIG = get_config()