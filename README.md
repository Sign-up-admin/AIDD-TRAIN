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
- **Important ViSNet Constraint**: The ViSNet model requires that the number of hidden channels be divisible by the number of attention heads. This is a key consideration when configuring the model.

---

## Step 1: Setup and Installation

1.  **Clone the repository:**
    ```sh
    git clone https://github.com/Sign-up-admin/AIDD-TRAIN.git
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

## Step 2: Automated Hardware Optimization (Highly Recommended)

Before the first training run, you should calibrate the model for your specific hardware. This crucial one-time step runs a series of benchmarks to find the most performant and stable configuration (batch size, model complexity, etc.) your GPU can handle, effectively eliminating `CUDA out of memory` errors.

The optimizer is now fully automated and uses a hierarchical strategy. It starts by finding the best configuration for the most demanding `production` mode and then uses those results to intelligently and rapidly find the optimal settings for the less demanding modes.

### Running the Optimizer

1.  Open your terminal.
2.  Navigate to the project's root directory (`AIDD-TRAIN`).
3.  Run the following single command:

    ```sh
    python src/hardware_optimizer.py
    ```

This command will optimize for all four modes (`production`, `validation`, `prototyping`, and `smoke_test`) in the correct order. The process may take some time, especially the initial `production` mode test, which runs for a high number of iterations (e.g., 500+) to ensure the configuration is truly stable under sustained load.

The script will create or update a `hardware_profile.json` file. The main training script (`main.py`) will automatically load the appropriate settings from this file based on the `DEVELOPMENT_MODE` you select in `config.py`.

### (Optional) Optimizing for Specific Modes

If you wish to re-run the optimization for only specific modes, you can use the `--modes` argument:

```sh
# Example: Optimize only for production and validation
python src/hardware_optimizer.py --modes production validation
```

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
