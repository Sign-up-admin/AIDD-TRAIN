# AIDD-TRAIN

## Project Overview

This repository, `AIDD-TRAIN`, documents a deep learning project focused on training models, specifically addressing common challenges encountered during the development process. A significant part of this repository highlights the debugging journey to resolve complex issues like CUDA out of memory errors and subsequent logical bugs in a deep learning training pipeline.

## Key Learnings and Debugging Insights

This project served as a practical case study for debugging deep learning training processes. Key insights gained include:

*   **Systematic Debugging**: Starting from common issues (e.g., batch size for OOM) and progressively delving into more complex root causes (model architecture, data characteristics).
*   **Handling Extreme Cases**: Understanding that data outliers or "edge cases" can often be the primary cause of resource exhaustion (like GPU memory) even when average loads are manageable.
*   **Code Mechanism Awareness**: Recognizing and utilizing built-in mechanisms (e.g., `DATA_PROCESSING_VERSION` for forcing data re-processing) to ensure changes are effectively applied.
*   **Precise Error Interpretation**: The importance of carefully reading and understanding specific error messages (e.g., `AttributeError: 'tuple' object has no attribute 'size'`) to pinpoint the exact problem.

For a detailed, step-by-step account of the debugging process, including initial diagnostics, deep dives into model and data issues, and resolution of subsequent bugs, please refer to the [Gemini Debugging Log](#gemini-debugging-log) below.

## Setup and Usage

*(This section needs to be filled with instructions on how to set up the environment, install dependencies, and run the training scripts. For example:)*

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/AIDD-TRAIN.git
    cd AIDD-TRAIN
    ```
2.  **Create a virtual environment and install dependencies:**
    ```bash
    conda create -n aidd-train python=3.x
    conda activate aidd-train
    pip install -r requirements.txt # Assuming a requirements.txt exists
    ```
3.  **Prepare Data:**
    *(Instructions on how to download or preprocess data, if applicable. Mention the `DATA_PROCESSING_VERSION` if relevant for users.)*
    *Note: If you modify data processing logic, remember to update `DATA_PROCESSING_VERSION` in `main.py` and delete the `processed_data` folder to force reprocessing.*

4.  **Run Training:**
    ```bash
    python main.py # Or specific command to start training
    ```

## Gemini Debugging Log

A comprehensive log detailing the debugging process for resolving CUDA out of memory errors and other related issues is available here:

*   [GEMINI.md](GEMINI.md)

---
*This README was last updated by Gemini.*