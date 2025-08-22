# AIDD-TRAIN

## Project Overview

This repository, `AIDD-TRAIN`, documents a deep learning project focused on training models for drug discovery, specifically addressing common and complex challenges encountered during development. A significant part of this repository highlights the debugging journey and architectural improvements made to the training pipeline, transforming it into a robust and reliable system.

## Key Learnings and Architectural Improvements

This project served as a practical case study for building and debugging deep learning training systems. Key insights and improvements include:

*   **Data Pipeline Unification**: A critical architectural flaw was identified and fixed where the data processing for **training** and **prediction** were inconsistent. The training pipeline was modified to generate explicit protein, ligand, and interaction graphs, matching the structure required for prediction. This ensures the model learns relevant physical interactions and that its predictions are valid.

*   **Robust Checkpointing & Recovery**: The system was hardened to save the complete training state (model, optimizer, scheduler, scaler) and to intelligently load the most recent checkpoint, whether from a normal save or a graceful shutdown.

*   **Graceful Shutdown**: The training loop now handles `KeyboardInterrupt` (Ctrl+C) to perform a final, safe save of the training state, preventing loss of progress on manual termination.

*   **Reproducible Data Splits**: The train/validation data split was made deterministic to ensure that validation loss is a consistent and comparable metric across different training runs.

*   **Systematic Debugging**: The project demonstrates a clear methodology for debugging, from surface-level errors (like CUDA OOM) to deep, logical bugs in the data pipeline.

For a detailed, step-by-step account of the debugging process and agent collaboration, please refer to the [Agent Collaboration Log](AGENT.md).

## Setup and Usage

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/AIDD-TRAIN.git
    cd AIDD-TRAIN
    ```
2.  **Create a virtual environment and install dependencies:**
    ```bash
    # It is recommended to use a conda environment for managing PyTorch and RDKit
    conda create -n aidd-train python=3.9
    conda activate aidd-train
    # Install PyTorch, PyG, RDKit, etc. from requirements.txt
    pip install -r requirements.txt
    ```
3.  **Prepare Data:**
    *   Place your dataset (e.g., PDBbind) in the directory specified in `config.py`.
    *   **Important**: The data processing logic is versioned. If you make any changes to the functions in `src/data_processing.py`, you **must** update the `data_processing_version` string in `main.py`. This will trigger a warning and require you to delete the old `processed_data/` directory to force regeneration of the dataset with the correct logic.

4.  **Run Training:**
    ```bash
    python main.py
    ```

5.  **Make Predictions:**
    ```bash
    python predict.py --protein_path /path/to/protein.pdb --ligand_path /path/to/ligand.sdf --checkpoint /path/to/model_best.pth.tar
    ```

## Hardware Configuration Guide

The model's memory footprint, especially for ViSNet, is highly dependent on your GPU's VRAM. The settings in `config.py` need to be adjusted to prevent `CUDA out of memory` errors. Here are some practical guidelines based on different hardware.

### Key Parameters

*   `batch_size`: The number of data samples processed by the GPU at once. **This is the first parameter you should lower when facing memory errors.**
*   `gradient_accumulation_steps`: Simulates a larger batch size. The "effective batch size" is `batch_size * gradient_accumulation_steps`. If you decrease `batch_size`, you should increase this value proportionally to keep the effective batch size consistent, which helps stabilize training.
*   `visnet_hidden_channels`: Controls the size of the model's hidden layers. If lowering `batch_size` to 1 is still not enough, **reducing this value (e.g., from 128 to 64) is the next step.** Note that this reduces the model's capacity and may impact its final accuracy.

### Example Configurations

#### Case 1: Low VRAM (e.g., NVIDIA RTX 3060 6GB)

With only 6GB of VRAM, the ViSNet model can easily cause out-of-memory errors even with a single sample. The following configuration was found to work on this hardware:

*   **Problem**: `CUDA out of memory` errors occur even with `batch_size: 2`.
*   **Solution**:
    1.  Set `batch_size` to the absolute minimum: `1`.
    2.  Increase `gradient_accumulation_steps` to `16` to maintain an effective batch size of `1 * 16 = 16`.
    3.  If memory errors persist, reduce the model's complexity by setting `visnet_hidden_channels` to `64`.

```python
# config.py for a 6GB VRAM GPU
CONFIG = {
    # ...
    'batch_size': 1,
    'gradient_accumulation_steps': 16,
    'visnet_hidden_channels': 64,
    # ...
}
```

#### Case 2: High VRAM (e.g., NVIDIA RTX 4060 Ti 16GB)

With more VRAM, you can afford a larger batch size and a more complex model, which can speed up training and potentially improve results.

*   **Advantage**: Can handle larger batches and a full-sized model.
*   **Configuration**:
    1.  Start with a larger `batch_size`, for example, `4` or `8`.
    2.  Adjust `gradient_accumulation_steps` accordingly. For an effective batch size of 16, you could use `batch_size: 4` and `gradient_accumulation_steps: 4`.
    3.  You can use the full model complexity with `visnet_hidden_channels: 128`.

```python
# config.py for a 16GB VRAM GPU
CONFIG = {
    # ...
    'batch_size': 4,
    'gradient_accumulation_steps': 4,
    'visnet_hidden_channels': 128,
    # ...
}
```

### Dealing with Numerical Instability (NaN errors)

During training, you might encounter NaN (Not a Number) or Inf (Infinity) values in the loss, especially when using mixed-precision (float16). This is a sign of numerical instability. While the training script is designed to skip these problematic batches, the root cause should be addressed.

1.  **Increase Optimizer Epsilon**: The Adam optimizer uses a small value `eps` to prevent division by zero. The default (`1e-8`) can be too small for `float16`'s limited range. Increasing it can significantly improve stability. In `main.py`, the optimizer is already initialized with a more stable value:
    ```python
    optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], eps=1e-7)
    ```

2.  **Identify Problematic Data**: The training script will now log the PDB codes of batches that produce `NaN` values. Look for these warnings in the output. If you see the same PDB codes appearing repeatedly, it indicates a problem with those specific data files (e.g., unusual geometry). You may need to inspect or exclude them from your dataset.

3.  **Disable Mixed Precision (for Debugging)**: If `NaN` issues persist, you can temporarily disable Automatic Mixed Precision (AMP) to confirm if `float16` is the cause. This involves commenting out the `autocast` and `GradScaler` code in `src/training.py`. This will increase VRAM usage and slow down training but is a reliable way to diagnose the problem.

## Debugging and Development Logs

Comprehensive logs detailing the debugging process and the collaboration with the AI agent are available in the following files:

*   **Agent Collaboration Log**: [AGENT.md](AGENT.md) - A high-level summary of the key problems identified by the user and the solutions implemented by the agent.
*   **Detailed Debugging Log**: [GEMINI.md](GEMINI.md) - A verbose, step-by-step log of the entire debugging and development conversation.

---
*This README was last updated by Gemini.*