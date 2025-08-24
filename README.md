# COMPASS: Navigating the Complexities of Protein-Ligand Binding

[![Python Version](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0%2B-orange.svg)](https://pytorch.org/)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

**COMPASS** is a deep learning project dedicated to accurately predicting protein-ligand binding affinities. It leverages a state-of-the-art Graph Neural Network (GNN), ViSNet, to learn from the intricate 3D geometry of molecular complexes, aiming to accelerate the process of drug discovery.

This project is built not just on a powerful model, but on a core philosophy of extreme data robustness and a highly efficient, mode-driven development workflow. Its name, COMPASS, reflects its ability to navigate the often-problematic landscape of real-world structural biology data, ensuring reliable and reproducible results.

---

## Key Features

- **State-of-the-Art Model**: Implements the ViSNet architecture for high-precision, geometry-aware predictions.
- **Robust Data Pipeline**: The cornerstone of COMPASS. The data processing pipeline is meticulously designed to handle common and obscure issues found in PDB data.
- **Automated Hardware Optimization**: A built-in tool to automatically find the best-performing configuration for your specific hardware, eliminating memory errors.
- **Mode-Driven Workflow**: A four-stage development process (`smoke_test`, `prototyping`, `validation`, `production`) that allows for seamless switching between quick checks, rapid experimentation, and full-scale training.
- **High-Performance Training**: Utilizes Automatic Mixed Precision (AMP) for significant speed-ups.
- **Resilient & Manageable**: Features robust checkpointing, graceful exit handling, and automated log/checkpoint organization.

---

## Step 1: Setup and Installation

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

4.  **Prepare Data:**
    - Download the PDBbind dataset (v.2020 or other).
    - Open `config.py` and update the `dataset_path` and `index_file` variables to point to your dataset location.

---

## Step 2: Automated Hardware Optimization (Recommended)

Before running a full training, it is highly recommended to calibrate the model for your specific hardware. This step eliminates `CUDA out of memory` errors and ensures you are always using the most performant configuration.

This process only needs to be run once for each mode you intend to use (e.g., `production`, `validation`).

**Optimize for Production Mode:** To find the best parameters for generating final results, run:
```sh
python src/hardware_optimizer.py --mode production
```

**Optimize for Validation Mode:** To find the best parameters for the validation stage, run:
```sh
python src/hardware_optimizer.py --mode validation
```

The script will test different configurations and save the optimal results to a `hardware_profile.json` file.

**Note:** Running the optimizer for a new mode will **add to** the existing `hardware_profile.json` file, not overwrite it. This allows you to store optimized settings for all your different workflow modes (e.g., `production`, `validation`, etc.) in a single file. The main training script will then automatically use these settings based on your selected mode.

---

## Step 3: Configuration and Daily Usage

Once the one-time setup and optimization are complete, your daily workflow is very simple.

1.  **Select Your Mode**: Open `config.py` and set the `DEVELOPMENT_MODE` variable to one of the four modes: `'smoke_test'`, `'prototyping'`, `'validation'`, or `'production'`.

2.  **Run Training**: Execute the main script from your terminal.
    ```sh
    python main.py
    ```

The script will automatically use the best settings for your chosen mode—either the optimized parameters from your `hardware_profile.json` or the default settings if no optimization was run for that mode.

All logs and model checkpoints will be saved into a uniquely named directory (e.g., `checkpoints/visnet_prototyping_.../`).

---

## Understanding the Four-Stage Workflow

COMPASS implements a phased workflow to balance speed and rigor. You can switch between modes by changing a single variable in `config.py`.

1.  **`smoke_test`**: *"Does the code run?"* A minimal check that runs in minutes.
2.  **`prototyping`**: *"Is my idea promising?"* A lightweight configuration for rapid experimentation.
3.  **`validation`**: *"How does my idea perform under realistic conditions?"* A medium-sized configuration for pre-production validation.
4.  **`production`**: *"What are the final, best-effort results?"* The full-scale configuration for generating final results.

---

## The COMPASS Philosophy: A Case Study in Data Robustness

This project was forged through a deep-dive debugging session to solve the sudden appearance of `NaN` (Not a Number) values during training. Instead of simply skipping problematic data, we developed a strategy of **"Pause and Autopsy"**.

This journey underscores the COMPASS philosophy: true progress in scientific machine learning comes not just from powerful architectures, but from a relentless commitment to understanding and purifying the data that fuels them.

## License

This project is licensed under the GNU AGPLv3 License. See the `LICENSE` file for details.

## Acknowledgments

This project utilizes the PDBbind dataset. We gratefully acknowledge the creators and maintainers of this valuable resource.
