"""
Hardware Optimizer Configuration File

This file centralizes all the tunable parameters for the hardware optimizer.
Modifying these values will directly influence the search strategy and performance
of the find_optimal_configs script.
"""

# ==============================================================================
# 1. CORE SEARCH SPACE DEFINITIONS
# ==============================================================================
ARCH_DEFINITIONS = {
    'high_vram': {
        'vram_threshold': 12,
        'description': "Standard architecture for high VRAM GPUs (>12GB)",
        'architectures': {
            'production':  {'layers': [8, 7, 6, 5], 'channels': [256, 192, 128]},
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
VRAM_SCALING_FACTORS = {
    24: 1.5,
    16: 1.2,
    12: 1.1,
    0:  1.0
}

# ==============================================================================
# 3. MODE-SPECIFIC PARAMETERS
# ==============================================================================
MODE_PARAMS = {
    'production':  {'bs': 16, 'stress': 20},
    'validation':  {'bs': 32, 'stress': 20},
    'prototyping': {'bs': 64, 'stress': 20},
    'smoke_test':  {'bs': 128, 'stress': 1}
}

# ==============================================================================
# 4. TIME-GUIDED OPTIMIZATION PARAMETERS
# ==============================================================================
TIME_RANGES = {
    'prototyping': (10, 40),
    'validation':  (60, 120)
}
CYCLE_BATCHES = 450

# ==============================================================================
# 5. DATA-AWARE PARAMETER CAPPING
# ==============================================================================
PARAMETER_CAPS = {
    1000:  50_000,
    10000: 100_000,
    float('inf'): 250_000
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
    'protein_path': 'compass/optimizer/1jmf/1jmf_pocket.pdb',
    'ligand_path': 'compass/optimizer/1jmf/1jmf_ligand.sdf',
    'binding_data': 'Kd=1.0nM',
    'max_atoms': 10200
}
