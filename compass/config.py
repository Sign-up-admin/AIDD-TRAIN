"""
 Configuration settings for the data processing and model training script.
 """
 
import os

os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'
 
from compass.optimizer.utils import load_optimized_settings
 
# --- Execution Mode ---
# Selects the operational mode, controlling all major hyperparameters.
# - 'validation_tuned': A tuned version for better GPU utilization on 6-8GB cards.
# - 'validation': A safe, conservative mode for basic validation.
# - 'prototyping': A lightweight setup for rapid experimentation.
# ... and others
#
# Options: 'validation_tuned', 'validation', 'prototyping', 'smoke_test', 'production'
EXECUTION_MODE = 'validation_tuned'
 
 
# --- Hyperparameter Sets for Each Mode ---
# Defines the core settings for each execution mode. Parameters can be adjusted
# based on hardware and experimental needs.
MODES = {
    'validation_tuned': {
        'epochs': 200,
        'batch_size': 2,
        'gradient_accumulation_steps': 8, # Increased for larger effective batch size
        'visnet_hidden_channels': 64,   # Increased model width
        'visnet_num_layers': 4,         # Increased model depth
        'visnet_num_rbf': 80,           # Increased RBF complexity
        'visnet_lmax': 2,               # Increased angular complexity
        'profile': False,
        'max_atoms': 4000,              # Increased data complexity
    },
    'validation': {
        'epochs': 200,
        'batch_size': 2,
        'gradient_accumulation_steps': 4,
        'visnet_hidden_channels': 48,
        'visnet_num_layers': 3,
        'visnet_num_rbf': 64,
        'profile': False,
        'max_atoms': 2000,
    },
    'prototyping': {
        'epochs': 50,
        'batch_size': 2,
        'gradient_accumulation_steps': 2,
        'profile': False,
        'max_atoms': 1500,
    },
    'smoke_test': {
        'epochs': 2,
        'batch_size': 1,
        'gradient_accumulation_steps': 1,
        'profile': False,
        'max_atoms': 1000,
    },
    'production': {
        'epochs': 100,
        'batch_size': 8,
        'gradient_accumulation_steps': 8,
    }
}
 
def get_config(mode_name=EXECUTION_MODE):
    """
    Generates the final configuration dictionary by merging base settings
    with mode-specific overrides and optimized hardware profiles.

    Args:
        mode_name (str): The execution mode, which determines the set of
                         hyperparameters to use. Defaults to `EXECUTION_MODE`.

    Returns:
        dict: A dictionary containing the complete configuration settings.

    Raises:
        ValueError: If an invalid `mode_name` is provided.
    """
    # --- Base Configuration ---
    config = {
        'execution_mode': mode_name,
        'force_data_reprocessing': True, # Must be true when max_atoms changes
        'data_processing_mode': 'strict',
        'max_atoms': 4000, # Match the tuned mode
        'failed_cases_dir': r'failed_cases',
        'index_file': r'index/INDEX_general_PL.2020R1.lst',
        'dataset_path': r'PDBbind-2025.8.4/P-L/',
        'processed_data_dir': r'processed_data',
        'learning_rate': 0.0001,
        'train_split': 0.8,
        'warmup_epochs': 3, # 增加学习率预热参数，3-5个epoch是很好的起点
        'weight_decay': 1e-5, # 将权重衰减也放入配置中，方便调整
        'dropout_rate': 0.5,
        'gradient_clip_val': 1.0,
        'max_num_neighbors': 32, # Restored to a more standard value
        'visnet_cutoff': 5.0,    # Restored
        'visnet_vecnorm_type': 'max_min',
        'processing_num_workers': os.cpu_count() or 16,
        'profile': False,
        'manual_save_trigger_file': 'SAVE_NOW.flag',
        'save_every_n_batches': 100,
        'seed': 42,
        'debug_mode': False,
        'diffusion': {
            'use_two_stage_diffusion': True,
            'stage1': {
                'enabled': True,
                'noise_level': 0.02,
            },
            'stage2': {
                'enabled': True,
                'noise_level': 0.1,
            }
        }
    }
 
    # --- Merge Mode-Specific Settings ---
    mode_settings = MODES.get(mode_name)
    if mode_settings is None:
        raise ValueError(f"Invalid EXECUTION_MODE: '{mode_name}'. Please choose from {list(MODES.keys())}.")
    config.update(mode_settings)
 
    # --- Apply Hardware-Optimized Settings ---
    # The following line is disabled to prevent overriding manual settings with old JSON files.
    # load_optimized_settings(config)
 
    # --- Post-processing and Validation ---
    run_name = f"visnet_{mode_name}_lr{config['learning_rate']}_bs{config['batch_size']}x{config['gradient_accumulation_steps']}"
    config['run_name'] = run_name
    config['log_dir'] = os.path.join('logs', run_name)
    config['checkpoint_dir'] = os.path.join('checkpoints', run_name)
 
    if not os.path.exists(config['failed_cases_dir']):
        os.makedirs(config['failed_cases_dir'])
 
    return config
 
# --- Global Config ---
CONFIG = get_config()
