import os
import platform
import torch

def get_hardware_recommendations(config, logger):
    """
    Analyzes system hardware and the current configuration, providing actionable
    recommendations and warnings directly to the user and log files.
    """
    logger.log("--- Analyzing Hardware and Configuration ---")
    recommendations = []
    warnings = []
    num_cpu_cores = os.cpu_count() or 1

    # --- Worker Recommendations ---
    if platform.system() == "Windows":
        if config.get('loader_num_workers', 0) > 0:
            warnings.append(
                f"On Windows, it's highly recommended to set 'loader_num_workers' to 0 for stability. "
                f"Your current value is {config.get('loader_num_workers')}."
            )
    else:
        if config.get('loader_num_workers') != num_cpu_cores:
            recommendations.append(
                f"Your system has {num_cpu_cores} CPU cores. For optimal data loading, "
                f"consider setting 'loader_num_workers' to {num_cpu_cores} in config.py."
            )

    # --- GPU and Development Mode Recommendations ---
    if torch.cuda.is_available():
        gpu_name = torch.cuda.get_device_name(0)
        total_mem_gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        logger.log(f"GPU Detected: {gpu_name} with {total_mem_gb:.2f} GB VRAM.")

        current_mode = config.get('development_mode', 'N/A')
        
        # Suggest the best mode based on VRAM
        suggested_mode = 'prototyping'
        if total_mem_gb < 8: # Covers RTX 3060 6GB
            suggested_mode = 'prototyping'
            if current_mode == 'production':
                warnings.append(
                    f"Your VRAM is {total_mem_gb:.2f} GB, which may be insufficient for the 'production' mode. "
                    f"This can lead to CUDA out-of-memory errors. It is highly recommended to switch "
                    f"DEVELOPMENT_MODE to '{suggested_mode}' in config.py for stable training."
                )
            recommendations.append(
                f"For a GPU with {total_mem_gb:.2f} GB VRAM, the '{suggested_mode}' mode offers the best balance of speed and stability."
            )
        elif 8 <= total_mem_gb < 16:
            suggested_mode = 'production' # Can handle production
            recommendations.append(
                "Your VRAM is sufficient for all modes. 'prototyping' is great for speed, "
                "'production' for best results."
            )
        else: # >= 16GB
            suggested_mode = 'production'
            recommendations.append(
                "Your high VRAM allows you to run 'production' mode comfortably. You could even "
                "increase batch_size or model size within the 'production' settings in config.py."
            )

    else:
        logger.log("No CUDA GPU detected. Training will be performed on the CPU.")
        warnings.append(
            "No CUDA GPU detected. Training on CPU will be extremely slow. "
            "The 'smoke_test' mode is recommended to verify functionality."
        )
        if config.get('development_mode') != 'smoke_test':
            warnings.append(
                f"You are currently in '{config.get('development_mode')}' mode without a GPU. "
                f"Consider switching DEVELOPMENT_MODE to 'smoke_test' in config.py."
            )

    # --- Final Output ---
    if recommendations:
        logger.log("--- Recommendations ---")
        for rec in recommendations:
            logger.log(f"- {rec}")
    
    if warnings:
        logger.log_warning("--- Configuration Warnings ---")
        for warn in warnings:
            logger.log_warning(f"- {warn}")
    
    logger.log("--------------------------------------")
