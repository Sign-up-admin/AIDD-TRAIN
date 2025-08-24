# COMPASS: Navigating the Complexities of Protein-Ligand Binding

[![Python Version](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-orange.svg)](https://pytorch.org/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

**COMPASS** is a deep learning project dedicated to accurately predicting protein-ligand binding affinities. It leverages a state-of-the-art Graph Neural Network (GNN), ViSNet, to learn from the intricate 3D geometry of molecular complexes, aiming to accelerate the process of drug discovery.

This project is built not just on a powerful model, but on a core philosophy of extreme data robustness and a highly efficient, mode-driven development workflow. Its name, COMPASS, reflects its ability to navigate the often-problematic landscape of real-world structural biology data, ensuring reliable and reproducible results.

---

## Key Features

- **State-of-the-Art Model**: Implements the ViSNet architecture for high-precision, geometry-aware predictions.
- **Robust Data Pipeline**: The cornerstone of COMPASS. The data processing pipeline is meticulously designed to handle common and obscure issues found in PDB data, including the automatic detection and resolution of duplicate atom coordinates between protein-ligand pairs.
- **Mode-Driven Workflow**: A four-stage development process (`smoke_test`, `prototyping`, `validation`, `production`) that allows you to seamlessly switch between sanity checks, rapid experimentation, and full-scale training with a single line change.
- **High-Performance Training**: Utilizes Automatic Mixed Precision (AMP) via `torch.amp` for significant speed-ups on compatible hardware.
- **Resilient & Manageable**: Features robust checkpointing, graceful exit handling, automated log/checkpoint organization, and a clear configuration system.

---

## The Four-Stage Development Workflow

To solve the conflict between rapid iteration and time-consuming training, COMPASS implements a phased workflow. You can switch between modes by changing a single variable in `config.py`, and all relevant hyperparameters will adjust automatically.

1.  **`smoke_test` (Smoke Test Mode)**
    -   **Goal**: *"Does the code run?"*
    -   **Use For**: A minimal configuration to perform a quick health check on the entire pipeline after code changes. It runs in minutes and is not intended for meaningful training.

2.  **`prototyping` (Prototyping Mode)**
    -   **Goal**: *"Is my idea promising?"*
    -   **Use For**: The scientist's lab. A lightweight configuration for agile experimentation. It's designed for a fast feedback loop (minutes per epoch) to quickly validate new ideas and observe training trends.

3.  **`validation` (Validation Mode)**
    -   **Goal**: *"How does my idea perform under realistic conditions?"*
    -   **Use For**: A medium-sized configuration that serves as a bridge between prototyping and production. Use it to validate a promising idea on a more representative dataset before committing to a full-scale run.

4.  **`production` (Production Mode)**
    -   **Goal**: *"What are the final, best-effort results?"*
    -   **Use For**: The final build. This mode uses the full-scale model and data to generate the final, reliable results for publication or deployment. It is the most time and resource-intensive mode.

---

## Setup and Installation

1.  **Clone the repository:**
    ```sh
    git clone <your-repo-url>
    cd AIDD-TRAIN
    ```

2.  **Create a virtual environment (recommended):**
    ```sh
    python -m venv venv
    source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
    ```

3.  **Install dependencies:**
    *This project relies on PyTorch and PyTorch Geometric. Please follow their official installation instructions for your specific CUDA version first.*
    ```sh
    # Example for CUDA 11.8
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
    pip install torch_geometric
    pip install rdkit-pypi biopython tqdm scipy
    ```

4.  **Data:**
    - Download the PDBbind dataset (v.2020 or other).
    - Update the `dataset_path` and `index_file` variables in `config.py` to point to your dataset location.

## Configuration and Usage

The entire project is controlled from the `config.py` file.

1.  **Select Your Mode**: Open `config.py` and set the `DEVELOPMENT_MODE` variable to one of the four modes: `'smoke_test'`, `'prototyping'`, `'validation'`, or `'production'`.

2.  **Run Training**: Execute the main script from your terminal.
    ```sh
    python main.py
    ```

The script will automatically select the correct hyperparameters based on your chosen mode, run the training, and save all logs and model checkpoints into a uniquely named directory (e.g., `checkpoints/visnet_prototyping_.../`).

---

## The COMPASS Philosophy: A Case Study in Data Robustness

This project was forged through a deep-dive debugging session to solve the sudden appearance of `NaN` (Not a Number) values during training. Instead of simply skipping problematic data, we developed a strategy of **"Pause and Autopsy"**.

By implementing a `debug_mode` and a dedicated analysis toolkit (`src/analyze_problem_batch.py`), we discovered the root cause: a collision of atom coordinates between protein and ligand files upon their merger. This led to the development of a robust data processing pipeline that intelligently validates and cleans the data before it ever reaches the model.

This journey underscores the COMPASS philosophy: true progress in scientific machine learning comes not just from powerful architectures, but from a relentless commitment to understanding and purifying the data that fuels them.

## License

This project is licensed under the GNU AGPLv3 License. See the `LICENSE` file for details.

## Acknowledgments

This project utilizes the PDBbind dataset. We gratefully acknowledge the creators and maintainers of this valuable resource.
