# AIDD-TRAIN Hardware Optimizer (`hardware_optimizer.py`)

This document provides a comprehensive overview of the hardware optimization script, `hardware_optimizer.py`. Its primary purpose is to automatically determine the optimal model configurations for your specific GPU hardware, ensuring the best possible performance without encountering Out-of-Memory (OOM) errors.

## 1. Purpose

The script intelligently tests various model architectures and batch sizes to find the sweet spot for your GPU. The final output is a `hardware_profile.json` file, which can be used by other parts of the project to train and run models with settings perfectly tailored to your machine.

## 2. The Intelligent Search Strategy

Instead of a brute-force search, the optimizer employs a multi-faceted, intelligent strategy to find the best configuration efficiently.

### a. Dynamic VRAM Scaling

The script's behavior is not one-size-fits-all. It begins by detecting your GPU's available VRAM and uses this to define a sensible search space.

- **Low VRAM (e.g., 6GB Laptop GPU):** The script uses a conservative set of base parameters (e.g., fewer layers, smaller channel sizes) to avoid immediate failures.
- **Mid-Range VRAM (e.g., 16GB Gaming GPU):** It applies a `vram_factor` to scale up the base parameters, starting its search with more ambitious and powerful model configurations.
- **High-End VRAM (e.g., 80GB A100):** The scaling factor is higher, allowing the script to immediately begin testing large, complex models, saving time by skipping configurations that would under-utilize the hardware.

### b. Hierarchical Optimization

The script optimizes for different modes in a specific order: `production` -> `validation` -> `prototyping` -> `smoke_test`.

This hierarchy is leveraged to save time. The script performs a full search for the best model architecture (`layers` and `channels`) **only once**, for the most demanding `production` mode. For all subsequent modes, it **inherits this optimal architecture** and simply finds the largest possible `batch_size` that the mode can support. This avoids redundant searches.

### c. Efficient Batch Size Search

Once a viable model architecture is found (i.e., it runs with `batch_size=1`), the script does not linearly test every batch size. Instead, it uses a two-phase approach:

1.  **Find OOM Ceiling:** It exponentially increases the batch size (1, 2, 4, 8, 16...) to quickly find a point of failure.
2.  **Binary Search:** It then performs a binary search between the last successful size and the failure point to precisely and rapidly pinpoint the maximum stable batch size.

## 3. Parameter Selection: The What and The Why

A key design decision is choosing which parameters to include in the automatic search. The goal is to optimize hardware-dependent variables while keeping model-defining variables fixed.

### Optimized Parameters (The "What")

These parameters are searched dynamically as they have the most significant impact on VRAM and computational load:

-   `visnet_num_layers`: The depth of the model.
-   `visnet_hidden_channels`: The width of the model.
-   `batch_size`: The number of samples processed in one iteration.

### Fixed Parameters (The "Why")

Several important parameters are intentionally kept fixed during optimization. Changing them is treated as a scientific modeling decision, not a hardware tuning one.

-   **`lmax=2`**: This controls the complexity of the spherical harmonic basis functions, which are crucial for representing 3D geometry. `lmax=2` (including s, p, and d orbitals) is a widely accepted standard that provides an excellent balance between model accuracy and computational cost. It defines the *capability* of the model, which we assume is fixed.

-   **`cutoff`**: This defines the maximum distance for atomic interactions. It is a physical assumption about the system being modeled. Changing it affects what the model can "see" and should be guided by domain expertise, not hardware specs.

-   **`num_heads=8`**: The number of attention heads is tightly coupled with the `hidden_channels` (which must be divisible by `num_heads`). Fixing it to a proven default value (8) dramatically simplifies the search space and stabilizes the optimization process.

## 4. How to Run

Simply execute the script from the project root:

```bash
python src/hardware_optimizer.py
```

The script will print its progress and, upon completion, will save the results to `hardware_profile.json`.

## 5. Practical Adjustments

Based on practical experience, some of the default stress test iteration counts have been adjusted to balance reliability with speed. For example, the `production` mode test and the final stability check now use **35 iterations**—a number found to be robust enough for validation while significantly speeding up the optimization process.
