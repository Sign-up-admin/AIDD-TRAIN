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

This project contains an intelligent hardware optimizer, `hardware_optimizer.py`, designed to find the perfect model configuration for different stages of the development lifecycle.

### The Core Philosophy

This optimizer was built upon a clear, hierarchical development philosophy defined by its architect, in collaboration with the Gemini agent. It recognizes that the "best" configuration is not a single setting, but a set of trade-offs tailored to the specific goal of each development phase. The optimizer's intelligence lies in its ability to weigh these trade-offs, using real-time performance estimation.

The workflow is defined by three core stages, each with a unique optimization target:

1.  **`prototyping` (Soft Target: ~20 min/cycle)**
    *   **Philosophy**: **Time is the ruler.** This stage is for rapid trial and error. The configuration must be fast enough to allow developers to quickly test ideas. The optimizer targets a ~20 minute cycle time (for a fixed 450-batch run) but allows for a **20-minute flexibility window**. 
    *   **Strategy**: It searches its dedicated **small model space** to find the configuration with the **highest throughput (max batch size)** that fits within this flexible time budget. This ensures the fastest possible iteration speed without prematurely discarding a slightly slower but much more powerful configuration.

2.  **`validation` (Soft Target: ~90 min/cycle)**
    *   **Philosophy**: **Balance is the key.** This stage acts as the crucial bridge between a promising prototype and a full-scale production run. It must be close enough to production quality to give meaningful results, but fast enough to not halt the development flow. It serves to seriously validate the findings from the `prototyping` stage.
    *   **Strategy**: It targets a ~90 minute cycle time, also with a **20-minute flexibility window**. It searches its dedicated **large model space** for the configuration with the **highest throughput**, striking the perfect balance between speed and quality.

3.  **`production` (Goal: Time-unlimited)**
    *   **Philosophy**: **Quality is the ultimate goal.** Time is no longer the primary constraint. This stage is for building the best possible model that the hardware and data can support, ready for deployment.
    *   **Strategy**: It employs a **two-stage optimization**: first, it finds the highest-quality (largest) model that respects both data and hardware limits. Second, it squeezes all remaining performance out of the hardware by finding the maximum possible batch size for that single best model.

This structured approach ensures that from the earliest idea to the final deployment, there is a perfectly optimized configuration to support the task at hand.

### Running the Optimizer

1.  Open your terminal.
2.  Navigate to the project's root directory (`AIDD-TRAIN`).
3.  Run the following single command:

    ```sh
    python src/hardware_optimizer.py
    ```

This command will optimize for all four modes (`production`, `validation`, `prototyping`, and `smoke_test`) in the correct order. The process may take some time.

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
