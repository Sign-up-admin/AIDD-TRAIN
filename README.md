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

## Debugging and Development Logs

Comprehensive logs detailing the debugging process and the collaboration with the AI agent are available in the following files:

*   **Agent Collaboration Log**: [AGENT.md](AGENT.md) - A high-level summary of the key problems identified by the user and the solutions implemented by the agent.
*   **Detailed Debugging Log**: [GEMINI.md](GEMINI.md) - A verbose, step-by-step log of the entire debugging and development conversation.

---
*This README was last updated by Gemini.*