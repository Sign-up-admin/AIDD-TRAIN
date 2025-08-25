"""
Hardware Optimizer Configuration File

This file centralizes all the tunable parameters for the hardware optimizer.
Modifying these values will directly influence the search strategy and performance
of the find_optimal_configs script.
"""

# ==============================================================================
# 1. CORE SEARCH SPACE DEFINITIONS
# ==============================================================================
# Philosophy:
# The optimizer tailors its search based on the GPU's VRAM and the desired
# development mode. We define a hierarchy of architectures from the most complex
# ('production') to the simplest ('smoke_test'). This ensures that the optimizer
# searches within a reasonable space for each use case.
#
# - 'production':  Aims for the highest quality model. Architectures are complex,
#                    intended for final model training and deployment.
# - 'validation':  Seeks a balance between performance and speed. Used for
#                    hyperparameter tuning and validation loops where training time
#                    is a consideration.
# - 'prototyping': For rapid development and debugging. Architectures are simple
#                    and train very quickly.
# - 'smoke_test':  A minimal configuration to ensure the entire pipeline runs
#                    without errors. It's the fastest and simplest check.
# ==============================================================================
ARCH_DEFINITIONS = {
    'high_vram': {
        'vram_threshold': 12,
        'description': "Standard architecture for high VRAM GPUs (>12GB)",
        'architectures': {
            'production':  {'layers': [8, 7, 6, 5, 4, 3], 'channels': [256, 192, 128, 96, 64]},
            'validation':  {'layers': [6, 5, 4, 3], 'channels': [128, 96, 64]},
            'prototyping': {'layers': [4, 3, 2],   'channels': [64, 48, 32]},
            'smoke_test':  {'layers': [2, 1],     'channels': [32, 24, 16]}
        }
    },
    'medium_vram': {
        'vram_threshold': 8,
        'description': "Conservative architecture for medium VRAM GPUs (8-12GB)",
        'architectures': {
            'production':  {'layers': [6, 5, 4], 'channels': [96, 64, 48]},
            'validation':  {'layers': [4, 3, 2], 'channels': [64, 48, 32]},
            'prototyping': {'layers': [3, 2, 1], 'channels': [32, 24, 16]},
            'smoke_test':  {'layers': [1],     'channels': [16, 8]}
        }
    },
    'low_vram': {
        'vram_threshold': 6,
        'description': "More conservative architecture for low VRAM GPUs (6-8GB)",
        'architectures': {
            'production':  {'layers': [4, 3, 2], 'channels': [48, 32, 24]},
            'validation':  {'layers': [3, 2],   'channels': [32, 24, 16]},
            'prototyping': {'layers': [2, 1],   'channels': [16, 8]},
            'smoke_test':  {'layers': [1],     'channels': [8]}
        }
    },
    'ultra_low_vram': {
        'vram_threshold': 0,
        'description': "Ultra-conservative architecture for very low VRAM GPUs (<=6GB)",
        'architectures': {
            'production':  {'layers': [3, 2], 'channels': [32, 24]},
            'validation':  {'layers': [2, 1], 'channels': [24, 16]},
            'prototyping': {'layers': [1],   'channels': [16, 8]},
            'smoke_test':  {'layers': [1],   'channels': [8]}
        }
    }
}

# ==============================================================================
# 2. VRAM-BASED SCALING FACTORS
# ==============================================================================
# Philosophy:
# To make the most of available VRAM, we apply a scaling factor to the number of
# layers and channels. GPUs with more VRAM can handle larger models, so this
# allows the optimizer to be more ambitious on higher-end hardware.
# ==============================================================================
VRAM_SCALING_FACTORS = {
    24: 1.5,
    16: 1.2,
    12: 1.1,
    0:  1.0
}

# ==============================================================================
# 3. MODE-SPECIFIC PARAMETERS
# ==============================================================================
# Philosophy:
# Each mode has different requirements for batch size (bs) and stress testing.
# - 'bs': A starting batch size for the search. Larger models ('production') need
#         smaller batch sizes to fit in memory.
# - 'stress': The number of iterations to run when measuring cycle time. More
#             iterations lead to more accurate time estimates, which is crucial
#             for the time-sensitive 'validation' mode.
# ==============================================================================
MODE_PARAMS = {
    'production':  {'bs': 16, 'stress': 20},
    'validation':  {'bs': 32, 'stress': 29},
    'prototyping': {'bs': 64, 'stress': 20},
    'smoke_test':  {'bs': 128, 'stress': 1}
}

# ==============================================================================
# 4. TIME-GUIDED OPTIMIZATION PARAMETERS
# ==============================================================================
# Philosophy:
# For 'validation' and 'prototyping', the goal is not just performance, but
# developer productivity. We define a "sweet spot" time range (in minutes) for
# a standard training cycle. The optimizer uses Bayesian methods to find the most
# efficient model within this time frame. 'production' mode does not have a time
# limit, as the priority is finding the absolute best model.
# ==============================================================================
TIME_RANGES = {
    'prototyping': (10, 40),
    'validation':  (60, 240) # Increased upper bound to find better models
}
CYCLE_BATCHES = 450

# ==============================================================================
# 5. DATA-AWARE PARAMETER CAPPING
# ==============================================================================
# Philosophy:
# To prevent overfitting, we cap the maximum number of model parameters based on
# the dataset size. A larger dataset can support a more complex model without
# memorizing the training data. This acts as a crucial guardrail, ensuring the
# optimizer doesn't select a model that is too powerful for the data it's being
# trained on. A 4MB model file roughly corresponds to 1 million parameters.
# ==============================================================================
PARAMETER_CAPS = {
    1000: 250_000,          # For tiny datasets, keep models simple (~1MB file).
    10000: 750_000,        # For small datasets (~3MB file).
    25000: 1_000_000,       # For medium datasets like ~19.5k, cap at 1M params (~4MB file).
    float('inf'): 2_000_000 # For large datasets, allow complex models (~8MB file).
}

# ==============================================================================
# 6. GENERAL EXECUTION SETTINGS
# ==============================================================================
OPTIMIZATION_HIERARCHY = ['production', 'validation', 'prototyping', 'smoke_test']

# ==============================================================================
# 7. TEST SAMPLE CONFIGURATION
# ==============================================================================
# Defines the real-world data sample used for all hardware stress tests.
# Using a sample that is representative of your dataset (e.g., average or
# 95th-percentile complexity) is crucial for accurate results.
#
# - 'pdb_code': A descriptive name for the sample.
# - 'protein_path': Relative path from the project root to the protein PDB file.
# - 'ligand_path': Relative path from the project root to the ligand SDF file.
# - 'binding_data': A string representing the binding data (optional).
# - 'max_atoms': The 'max_atoms' setting required to process this specific
#   sample. This should be set slightly higher than the actual number of
#   atoms in the sample to ensure it can be processed correctly.
# ==============================================================================
TEST_SAMPLE_CONFIG = {
    'pdb_code': '1jmf',
    'protein_path': '1jmf/1jmf_protein.pdb',
    'ligand_path': '1jmf/1jmf_ligand.sdf',
    'binding_data': 'Kd=1.0nM',
    'max_atoms': 10200
}
