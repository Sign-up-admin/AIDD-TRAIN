"""
Hardware Optimizer Configuration File

This file centralizes all the tunable parameters for the hardware optimizer.
Modifying these values will directly influence the search strategy and performance
of the find_optimal_configs script.
"""

# ==============================================================================
# 1. CORE SEARCH SPACE DEFINITIONS
# ==============================================================================
# This is the most critical part of the configuration. It defines the search
# space for different GPU VRAM tiers. The optimizer will select one of these
# based on the detected GPU memory.
#
# Structure:
#   - Keys: A descriptive name for the VRAM tier (e.g., 'low_vram', 'high_vram').
#   - Values: A dictionary containing:
#     - 'vram_threshold': The minimum VRAM in GB for this tier to be selected.
#     - 'architectures': A dictionary defining the model architectures for each
#       optimization mode ('production', 'validation', etc.).
#       - 'layers': A list of `visnet_num_layers` to test.
#       - 'channels': A list of `visnet_hidden_channels` to test.
#
# The tiers are evaluated in descending order of 'vram_threshold'.
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
# These factors are applied to the lists defined in ARCH_DEFINITIONS to subtly
# scale the search space based on more granular VRAM availability.
# This allows for a more continuous adaptation than the discrete tiers above.
#
# Structure:
#   - Keys: Minimum VRAM in GB for the factor to be applied.
#   - Values: A multiplier for `layers`, `channels`, and `start_batch_size`.
#
# The factors are evaluated in descending order.
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
# Defines the starting batch size and stress test intensity for each mode.
#
# Structure:
#   - 'bs': The initial batch size to start the search from.
#   - 'stress': The number of iterations for the binary search stability check.
# ==============================================================================
MODE_PARAMS = {
    'production':  {'bs': 16, 'stress': 20},
    'validation':  {'bs': 32, 'stress': 20},
    'prototyping': {'bs': 64, 'stress': 20},
    'smoke_test':  {'bs': 128, 'stress': 1} # Minimal stress for a quick check.
}

# ==============================================================================
# 4. TIME-GUIDED OPTIMIZATION PARAMETERS (for Prototyping/Validation)
# ==============================================================================
# Defines the "sweet spot" for training cycle times. The optimizer's goal for
# prototyping and validation is to find the most efficient model that fits
# within this time window.
#
# Structure:
#   - Keys: The optimization mode.
#   - Values: A tuple of (minimum_time, maximum_time) in minutes.
# ==============================================================================
TIME_RANGES = {
    'prototyping': (10, 40),
    'validation':  (60, 120)
}

# The number of batches used to estimate the total cycle time.
# This should represent a typical training epoch or a significant portion thereof.
CYCLE_BATCHES = 450

# ==============================================================================
# 5. DATA-AWARE PARAMETER CAPPING
# ==============================================================================
# To prevent overfitting, the optimizer caps the maximum number of model
# parameters based on the dataset size.
#
# Structure:
#   - Keys: The maximum dataset size for this cap to apply.
#   - Values: The maximum number of trainable parameters allowed.
#
# The caps are evaluated in ascending order.
# ==============================================================================
PARAMETER_CAPS = {
    1000:  50_000,   # For small datasets (<1k samples)
    10000: 100_000,  # For medium datasets (<10k samples)
    float('inf'): 250_000 # For large datasets
}

# ==============================================================================
# 6. GENERAL EXECUTION SETTINGS
# ==============================================================================
# The hierarchical order in which to run the optimizations if multiple modes
# are requested. This ensures that higher-priority modes can inform or
# influence lower-priority ones if needed in future implementations.
OPTIMIZATION_HIERARCHY = ['production', 'validation', 'prototyping', 'smoke_test']
