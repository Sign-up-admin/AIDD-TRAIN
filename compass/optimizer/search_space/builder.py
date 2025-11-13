import logging
from ..config import ARCH_DEFINITIONS, VRAM_SCALING_FACTORS, MODE_PARAMS, PARAMETER_CAPS

logger = logging.getLogger("HardwareOptimizer")


def get_search_space(vram_gb, mode_to_optimize):
    """
    Generates the search space for a given mode and VRAM size.

    Args:
        vram_gb (float): The amount of VRAM in GB.
        mode_to_optimize (str): The optimization mode.

    Returns:
        tuple: A tuple containing the starting batch size, a list of hidden channels,
               a list of number of layers, and the number of stress iterations.
    """
    logger.info(f"Generating dedicated search space for '{mode_to_optimize}' mode...")
    selected_tier = next(
        (
            ARCH_DEFINITIONS[t]
            for t in sorted(
                ARCH_DEFINITIONS, key=lambda k: ARCH_DEFINITIONS[k]["vram_threshold"], reverse=True
            )
            if vram_gb >= ARCH_DEFINITIONS[t]["vram_threshold"]
        ),
        None,
    )
    if not selected_tier:
        logger.critical("Could not find a suitable architecture definition.")
        exit()
    logger.info(f"--- Using '{selected_tier['description']}' ---")
    arch_def = selected_tier["architectures"][mode_to_optimize]
    vram_factor = next(
        (
            VRAM_SCALING_FACTORS[vt]
            for vt in sorted(VRAM_SCALING_FACTORS.keys(), reverse=True)
            if vram_gb >= vt
        ),
        1.0,
    )
    logger.info(
        f"--- Applying VRAM scaling factor of {vram_factor} (based on {vram_gb:.2f} GB) ---"
    )
    num_attention_heads = 8
    hidden_channels_list = sorted(
        [
            c
            for c in {
                ((int(c * vram_factor) // num_attention_heads) * num_attention_heads)
                for c in arch_def["channels"]
            }
            if c > 0
        ],
        reverse=True,
    )
    num_layers_list = sorted(
        [
            layer
            for layer in {int(layer * vram_factor) for layer in arch_def["layers"]}
            if layer > 0
        ],
        reverse=True,
    )
    start_batch_size = int(MODE_PARAMS[mode_to_optimize]["bs"] * vram_factor)
    stress_iterations = MODE_PARAMS[mode_to_optimize]["stress"]
    if not num_layers_list or not hidden_channels_list:
        logger.critical("Search space is empty.")
        exit()
    logger.info(f"Search space: Layers={num_layers_list}, Channels={hidden_channels_list}")
    logger.info(f"Starting BS: {start_batch_size}, Stress Iterations: {stress_iterations}")
    return start_batch_size, hidden_channels_list, num_layers_list, stress_iterations


def get_parameter_cap(current_dataset_size):
    """
    Gets the parameter cap for a given dataset size.

    Args:
        current_dataset_size (int): The size of the current dataset.

    Returns:
        int: The parameter cap.
    """
    level, cap = "large", PARAMETER_CAPS[float("inf")]
    for size_thresh in sorted(PARAMETER_CAPS.keys()):
        if current_dataset_size < size_thresh:
            cap, level = PARAMETER_CAPS[size_thresh], f"small (<{size_thresh})"
            break
    logger.info(
        f"--- Dataset size: {current_dataset_size} ({level}). Capping model parameters at ~{cap/1e3:.0f}k to prevent overfitting. ---"
    )
    return cap
