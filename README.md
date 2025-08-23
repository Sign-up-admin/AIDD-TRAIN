# COMPASS: Navigating the Complexities of Protein-Ligand Binding

[![Python Version](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-orange.svg)](https://pytorch.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**COMPASS** is a deep learning project dedicated to accurately predicting protein-ligand binding affinities. It leverages a state-of-the-art Graph Neural Network (GNN), ViSNet, to learn from the intricate 3D geometry of molecular complexes, aiming to accelerate the process of drug discovery.

This project is built not just on a powerful model, but on a core philosophy of extreme data robustness. Its name, COMPASS, reflects its ability to navigate the often-problematic landscape of real-world structural biology data, ensuring reliable and reproducible results.

---

## Key Features

- **State-of-the-Art Model**: Implements the ViSNet architecture for high-precision, geometry-aware predictions.
- **Robust Data Pipeline**: The cornerstone of COMPASS. The data processing pipeline is meticulously designed to handle common and obscure issues found in PDB data, including:
    - Filtering of excessively large or malformed structures.
    - **Automatic detection and resolution of duplicate atom coordinates**, both within a single molecule and, crucially, between protein-ligand pairs upon merging. This prevents numerical instability (`NaN` values) during training.
- **High-Performance Training**: Utilizes Automatic Mixed Precision (AMP) via `torch.amp` for significant speed-ups on compatible hardware.
- **Resilient & Manageable**: Features robust checkpointing, graceful exit handling for interruptions, and a clear configuration system.

## The COMPASS Philosophy: A Case Study in Data Robustness

This project was forged through a deep-dive debugging session aimed at solving a common yet elusive problem in deep learning: the sudden appearance of `NaN` (Not a Number) values during training.

#### The Challenge

Initial training runs were plagued by `NaN` losses, which would appear unpredictably for certain data points (e.g., PDB ID `1p06`). A common but ineffective approach is to simply skip these problematic batches. We chose a different path.

#### The Investigation

We adopted a strategy of **"Pause and Autopsy"** over "Skip and Ignore".

1.  **`debug_mode`**: A debug flag was implemented in the training script. When activated, the script would automatically save the exact data batch that caused the `NaN` error and halt execution.
2.  **Analysis Toolkit**: A dedicated analysis script, `src/analyze_problem_batch.py`, was created to load these saved "crime scenes" and perform a deep inspection of every tensor.

#### The Discovery

Initial analysis confirmed the input data files themselves were clean. The `NaN`s were being generated *during* the model's forward pass. This led to the hypothesis of numerical instability caused by a **division-by-zero** error, likely from two atoms sharing the same coordinates.

Our specialized analysis script confirmed this, detecting duplicate coordinates in the problematic batch. However, the final piece of the puzzle was uncovered through manual inspection: the coordinate collision was not within a single file, but **between a protein atom and a ligand atom**. When the two molecules were merged into a single graph, these two atoms occupied the exact same point in 3D space.

#### The Solution

The data processing pipeline (`src/data_processing.py`) was upgraded to act as a true "compass" for our data. It now intelligently checks for and removes these overlapping atoms at the source, ensuring that every single data point entering the model is geometrically sound and numerically stable.

This journey underscores the COMPASS philosophy: true progress in scientific machine learning comes not just from powerful architectures, but from a relentless commitment to understanding and purifying the data that fuels them.

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

## Configuration

All project settings, including file paths, model hyperparameters, and training parameters, are managed in the `config.py` file. Key settings include:

- `dataset_path`: Path to the root of the PDBbind dataset.
- `processed_data_dir`: Directory to save the processed `.pt` graph files.
- `batch_size`, `epochs`, `learning_rate`: Standard training hyperparameters.
- `debug_mode`: Set to `True` to enable the "Pause and Autopsy" feature for debugging new data issues.

## Usage

To start the training process, simply run:

```sh
python main.py
```

The script will first process the raw PDB/SDF files into graph data and save them to the `processed_data_dir`. Subsequent runs will load these processed files directly, unless the data processing script (`src/data_processing.py`) is modified.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgments

This project utilizes the PDBbind dataset. We gratefully acknowledge the creators and maintainers of this valuable resource.
