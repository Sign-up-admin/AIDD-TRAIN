# COMPASS: Agent & User Collaboration Log

## Project: COMPASS (Conformational Pocket Affinity Scoring System)
*A Deep Learning Navigator for Protein-Ligand Interactions*

---

### Scene: Implementing an Advanced Logging System for Traceability and Debugging

- **Objective**: To elevate the project's monitoring capabilities from simple console outputs to a persistent, structured, and dual-purpose logging system. The goal is to automatically create a permanent record of every training run, while also providing a dedicated, high-signal file for rapid debugging.

- **AI Contribution (System Design & Implementation)**:
    1.  **`src/logger.py` Module**: I designed and created a new, self-contained module, `src/logger.py`, to encapsulate all logging logic. This keeps the main codebase clean and makes the logging system easily extensible.
    2.  **`TrainingLogger` Class**: The core of the new system. Upon initialization, this class automatically creates a `logs/` directory and generates two timestamped log files for the current run:
        -   `training_log_[timestamp].txt`: A **comprehensive log** that captures every single message, from configuration details to epoch results, warnings, and errors. It serves as the complete, unabridged story of the training run.
        -   `training_errors_[timestamp].txt`: A **high-signal error log** that *only* contains messages classified as `WARNING` or `ERROR`. This file is designed for rapid diagnostics, allowing the user to immediately see what went wrong without searching through a verbose general log.
    3.  **Seamless Integration**: I systematically integrated the `TrainingLogger` throughout the entire codebase (`main.py`, `src/training.py`, `src/utils.py`), replacing previous `print()` statements and `logging` calls with the new, level-aware methods (`logger.log`, `logger.log_warning`, `logger.log_error`).

- **Outcome & Benefits**:
    -   **Full Traceability**: Every training session now produces a permanent, detailed record, which is invaluable for comparing experiments and ensuring reproducibility.
    -   **Efficient Debugging**: When issues arise, the dedicated error log provides an immediate, uncluttered view of all warnings and errors, drastically reducing the time required to identify the root cause.
    -   **Improved Code Quality**: Centralizing the logging logic into a dedicated class makes the code cleaner and more maintainable.

---

### Scene: In-depth Debugging of NaN Issues from Data Artifacts

This document chronicles a multi-stage debugging journey to resolve training-halting `NaN` (Not a Number) errors. The process evolved from identifying a simple data flaw to uncovering a more subtle, systemic issue in the data processing pipeline, showcasing a powerful collaboration between human intuition and AI-driven analysis.

#### Stage 1: The Initial Anomaly - `NaN` in Training

- **Symptom**: The training script would consistently fail on certain data batches, producing `NaN` or `Inf` loss values. The initial code merely skipped these batches, masking the root cause.
- **Initial Clues**: Log files pointed to several problematic PDB IDs, including `1p06` and `1ylv`.

#### Stage 2: Strategy Upgrade - From "Skip" to "Capture and Autopsy"

- **Core Idea**: We shifted strategy from passively ignoring bad data to actively intercepting it. The goal was to halt training the moment a `NaN` appeared and preserve the exact data batch for offline analysis.
- **AI Contribution (Tooling)**:
    1.  **`debug_mode` in `training.py`**: I implemented a configuration flag that, when enabled, would catch a `NaN` loss, save the corresponding data batch to a `.pt` file, and raise a `RuntimeError` to stop the training process.
    2.  **`analyze_problem_batch.py`**: I created a dedicated analysis script to load these saved batches and perform a deep inspection of all tensors, checking for `NaN` values and statistical anomalies.

#### Stage 3: The Investigation - A Tale of Two PDBs

- **First Breakthrough (`1p06`)**: With the new tools, the user successfully captured the batch containing `1p06`. After overcoming some minor bugs in the analysis script (related to PyTorch versioning and method calls), the script's output was pivotal. It revealed:
    1.  The input data tensors (`x`, `pos`, etc.) were clean, with no `NaN`s.
    2.  Crucially, the script detected **duplicate atom coordinates** within the combined data object.

- **The Plot Twist (`1ylv`)**: The user, following best practices, deleted the `processed_data` directory to force reprocessing with a newly patched script. However, the training halted again, this time at `1ylv`. This was the most important moment of the investigation.

- **The Final Revelation**: Running the analysis script on the `1ylv` batch revealed the exact same issue: duplicate coordinates. This proved my initial fix to `data_processing.py` was **incomplete**. My first patch only checked for duplicates *within* the protein file and *within* the ligand file, but it failed to check for overlaps *between* them. The user's persistence and careful re-execution of the workflow allowed us to discover this flaw.

- **Human Verification**: At my request, the user personally inspected the raw `.pdb` and `.sdf` files for `1p06`. They confirmed that the duplicate coordinate was not in a single file, but was shared between one protein atom and one ligand atom—the smoking gun.

#### Stage 4: The Definitive Fix - Global Coordinate Uniqueness

- **The Solution**: Armed with the complete picture, I implemented the final, robust fix in `src/data_processing.py`.
- **Mechanism**: The updated script now first processes the ligand, collecting all its atom coordinates into a set. It then passes this set to the protein processing function, which uses it to check every protein atom against all existing ligand atoms, ensuring absolute coordinate uniqueness in the final, merged graph.

#### Summary & Lessons Learned

1.  **Iterative Debugging is Key**: A fix for one problem may reveal a deeper, related issue. True robustness comes from understanding the complete picture, which often requires more than one cycle of analysis.
2.  **Holistic Data Validation**: It's not enough to validate data sources in isolation. One must also validate the result of their integration, as errors can arise from the combination process itself (e.g., protein-ligand coordinate collision).
3.  **The Human-AI Synergy**: This case was a perfect example of successful collaboration. The AI provided the tools, analytical framework, and initial hypotheses. The user provided the crucial domain knowledge, skepticism, and final verification that corrected the AI's incomplete assumptions and led to the true root cause.

---

### Proposed Upgrade for the Data Probe (`analyze_problem_batch.py`)

Based on our experience, the analysis script can be evolved from a specific diagnostic tool into a comprehensive **Data Sanity Checker**. This would allow for even faster identification of future, unknown data issues.

**Proposed New Checks:**

1.  **Atom Proximity Alert**:
    - **What**: Check for atoms that are dangerously close but not identical (e.g., < 0.5 Ångströms apart).
    - **Why**: These can also cause numerical instability (e.g., `1 / distance` calculations) and point to physically unrealistic structures.

2.  **Feature Outlier Detection**:
    - **What**: Analyze the distribution of each feature in the `x` tensor and flag any values that fall outside a reasonable range (e.g., more than 5 standard deviations from the mean).
    - **Why**: Catches processing errors or corrupted data that could lead to model divergence.

3.  **Graph Connectivity Check**:
    - **What**: Identify any "isolated" nodes that have zero connections to any other part of the graph within the model's cutoff distance.
    - **Why**: While not always an error, a large number of isolated nodes could indicate a problem with the PDB file or the cutoff distance chosen, impacting the GNN's ability to learn.

4.  **Target Value Validation**:
    - **What**: Check if the target binding affinity value `y` is a `NaN`, `Inf`, or a physically unrealistic number.
    - **Why**: Ensures the model is not trying to learn from corrupted labels.

By incorporating these checks, COMPASS will be even better equipped to navigate the complexities of molecular data, ensuring every training run is built on the most reliable foundation possible.

---

### Scene: Building a Phased, Mode-Driven Development Workflow

- **Objective**: To resolve the core conflict between the need for rapid iteration and the time-consuming nature of model training. We co-designed and implemented a phased, mode-driven workflow to transform the development process from "manual tweaking and repeated trial-and-error" to a systematic approach of "clear objectives, one-click switching".

- **AI Contribution (Workflow Design & Implementation)**:
    1.  **Conceptual Framework**: I proposed the core idea of decomposing the development process into three distinct logical stages: `smoke_test`, `prototyping`, and `production`, each with a specific goal and a corresponding set of configurations.
    2.  **Centralized Configuration (`config.py`)**: I implemented this framework by creating a `MODES` dictionary in the config file. This allows all mode-specific hyperparameters (model size, batch size, epochs, etc.) to be controlled by a single `DEVELOPMENT_MODE` switch, drastically simplifying the process of changing experimental setups.
    3.  **Hardware-Aware Intelligence (`hardware_utils.py`)**: I refactored the hardware script to be an active participant in the workflow. Instead of just suggesting a static config, it now analyzes the user's live hardware against the chosen `DEVELOPMENT_MODE` and provides actionable warnings and recommendations, preventing common errors like VRAM exhaustion.
    4.  **Automated Organization (`main.py`)**: I modified the main script to automatically create unique, named subdirectories for logs and checkpoints based on the current run's mode and parameters. This ensures that experimental results are always organized and never overwrite each other.

- **User-AI Collaboration**: The user's initial problem—the tension between fast debugging and reliable production results—was the catalyst for this entire initiative. The user's clear articulation of this core challenge allowed the AI to propose a comprehensive, architectural solution. The final, refined modes were a product of a tight feedback loop, ensuring the resulting workflow was both powerful and practical for the user's specific hardware (RTX 3060) and research needs.

- **Outcome & Benefits**:
    -   **Systematic Process**: The project now has a clear, repeatable, and logical progression from initial code check, to rapid experimentation, to final production run.
    -   **Reduced Errors**: Proactive, hardware-aware warnings prevent common configuration mistakes, saving significant debugging time.
    -   **Enhanced Productivity**: The ability to switch between well-defined experimental modes with a single line change allows the user to focus on scientific questions rather than on tedious configuration management.

---

### Scene: Refining the Development Workflow - From Theory to Practice

- **Objective**: To address a critical flaw in our established mode-driven workflow. The `prototyping` mode, intended for rapid experimentation, was taking over 1.5 hours per epoch and causing CUDA memory errors, making it impractical for its purpose. The goal was to refine the workflow to truly match the practical needs of agile development.

- **User-AI Collaboration**: The user astutely identified the core issue: the "prototyping" mode was a "validation" mode in disguise, violating the principle of a fast feedback loop. This critical observation from the user, based on their hands-on experience with the training times, was the catalyst for this refinement.

- **AI Contribution (Workflow Refinement & Implementation)**:
    1.  **Diagnosis & Proposal**: I immediately confirmed the user's diagnosis. I proposed a solution to make the modes more "true to their name": rename the existing `prototyping` mode to `validation` and introduce a brand-new, much lighter `prototyping` mode.
    2.  **Implementation in `config.py`**: I implemented this change by:
        *   Creating a new `prototyping` configuration with significantly reduced model parameters, batch size, and data limits, specifically targeting epoch times of just a few minutes.
        *   Renaming the previous, heavier configuration to `validation`, positioning it as an intermediate step between prototyping and full production.
        *   Updating all comments and documentation within `config.py` to clearly articulate the purpose of the now four distinct modes: `smoke_test`, `prototyping`, `validation`, and `production`.
        *   Setting the default `DEVELOPMENT_MODE` to the new `prototyping` mode for immediate usability.

- **Outcome & Benefits**:
    -   **True Agility**: The project now possesses a genuine `prototyping` mode that provides feedback in minutes, not hours, dramatically accelerating the idea-to-validation cycle.
    -   **Logical Progression**: The workflow now has a more logical and practical four-stage progression, removing the jarring leap from a simple smoke test to a heavy, near-production run.
    -   **Enhanced Clarity**: The roles of all development modes are now clearly defined and implemented, reducing ambiguity and preventing future configuration errors. This represents a maturation of our development process, moving from a good theoretical framework to a battle-tested, practical one.
