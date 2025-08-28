# GEMINI: Core Architectural Decisions & Lessons Learned

This document records the key architectural decisions, problem-solving patterns, and core philosophies established during the development of the AIDD-TRAIN project. It serves as a living repository of engineering wisdom for the project.

---

### The Philosophy of Insightful Documentation

- **Problem**: Initial project documentation was factually correct but superficial. It described *what* the code did, but not *why* it was designed that way. This creates a knowledge gap for new developers and makes long-term maintenance difficult.

- **Analysis**: Good documentation is not a summary of the code; it is a guide to the *thinking behind* the code. It must explain the core problems being solved and the architectural choices made to solve them. The transition from a "流水账" (a plain, chronological account) to a deep, insightful guide is a critical step in maturing a codebase.

- **Solution (The Three-Pillar Documentation Framework)**:
    We established a new standard for all core module documentation, built on three pillars designed to deliver insight, not just information:
    1.  **The Core Conflict**: Start by identifying the central tension or problem the code is designed to resolve (e.g., "the gap between raw data and model-ready input"). This immediately establishes the *purpose* of the module.
    2.  **The Workflow & Philosophy**: Explain the high-level architectural pattern or design philosophy used to solve the conflict (e.g., "a Conductor-Executor pattern for training"). This reveals the *strategy*.
    3.  **The FAQ / Common Problems**: Provide concrete, practical solutions to common issues. This connects the high-level philosophy to real-world application and dramatically speeds up debugging and development.

- **Lesson Learned**: Documentation should be treated as a first-class product, not an afterthought. By adopting a structured, narrative-driven framework, we transform documentation from a simple reference into a powerful tool for onboarding, strategic alignment, and long-term maintainability. It codifies the project's architectural soul.

---

### Scene: Resolving a Circular Import Error via Architectural Refactoring

- **Problem**: A critical `ImportError` blocked the application at startup. The traceback indicated a classic **circular import**: `engine.py` depended on `loop.py`, which in turn depended on `engine.py`.

- **Analysis**: This dependency loop creates a situation where neither module can be fully loaded before the other, leading to a resolution failure. The root cause was shared utility functions (`_save_checkpoint`, `_load_checkpoint`) being housed in a module (`engine.py`) that also had higher-level application logic.

- **Solution (The Decoupling Pattern)**:
    1.  **Identify the Common Utility**: The checkpointing functions were identified as a self-contained, reusable utility.
    2.  **Create a Dedicated Module**: A new module, `compass/training/checkpoint.py`, was created to adhere to the **Single Responsibility Principle**.
    3.  **Relocate and Refactor**: The checkpointing functions were moved from `engine.py` to the new `checkpoint.py` module.
    4.  **Update Dependencies**: Both `engine.py` and `loop.py` were updated to import the checkpointing functions directly from the new, independent `checkpoint.py` module.

- **Lesson Learned**: When two modules appear to depend on each other, it often signals that a common, lower-level functionality can be extracted into a third, independent module. This **decoupling** not only solves the immediate circular import error but also leads to a cleaner, more logical, and more maintainable code architecture.

---

### AIDD-TRAIN Hardware Optimizer: The Guiding Philosophy

This document records the final core philosophy embedded into the `hardware_optimizer.py` script. This philosophy was established and refined by the project's architect, and implemented by the Gemini agent, creating a tool that understands the nuanced, trade-off-driven demands of a real-world machine learning workflow.

## The Core Philosophy

The optimizer's intelligence is rooted in a hierarchical development philosophy. It rejects the notion of a single "best" configuration, instead providing tailored trade-offs for each specific phase of development. Its wisdom lies in its ability to weigh these trade-offs, using real-time performance estimation against a flexible budget.

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

This structured approach, born from our collaboration, ensures that from the earliest idea to the final deployment, there is a perfectly optimized configuration to support the task at hand.

## Evolution and Refinements

Beyond the initial philosophy, the optimizer underwent a significant evolution, transforming it from a clever script into a robust, efficient, and intelligent engineering product. This was achieved through a series of collaborative refinements:

### From a Tool to a Product: Engineering Excellence

*   **Configuration as Code**: All key parameters—from model architecture search spaces to VRAM scaling factors and time budgets—were extracted from the main script into a dedicated `optimizer_config.py` file. This separation of concerns makes the optimizer easier to maintain and adapt.
*   **Code Refactoring for Clarity**: The monolithic `find_optimal_configs` function was refactored. The logic for each optimization mode (`production`, `validation`, etc.) was encapsulated into its own dedicated function, making the high-level logic cleaner and more readable.
*   **Professional Logging**: All `print` statements were replaced with a robust, dual-channel `logging` system. This provides timestamped, leveled output to both the console and a persistent log file (`hardware_optimizer.log`), dramatically improving traceability and debugging.
*   **Flexible Test Bench**: The hardcoded '1jmf' test sample was made fully configurable through `optimizer_config.py`, allowing users to easily specify their own representative data sample and its required `max_atoms` count, enhancing the tool's real-world applicability.

### From Brute-Force to Intelligence: Algorithmic Enhancements

*   **Fused Efficiency**: The redundant, two-step process of finding the max batch size and then separately estimating its cycle time was fused into a single, efficient operation. The cycle time is now derived directly from the timing of the final stability probe, halving the required GPU work for each candidate.
*   **Granular Pruning**: The initial, aggressive pruning strategy was refined. Instead of pruning an entire layer based on a single slow configuration, the new logic only prunes configurations that are more complex in *both* layer count and channel width, allowing for a more thorough and intelligent exploration of the search space.
*   **Bayesian Optimization**: The most significant leap was replacing the grid-based search for `prototyping` and `validation` modes with a Bayesian Optimization engine. The optimizer now learns from each test, building a probabilistic model of the search space to intelligently select the next most promising candidate. This allows it to find better configurations in fewer steps, making the entire process not just faster, but smarter.

## The Importance of Code Engineering

Striving to build project code with an engineering mindset is key to improving project maintainability, scalability, and team collaboration efficiency. The necessity of reading and understanding good project engineering practices helps us write more robust and clearer code, laying a solid foundation for the long-term success of the project.

## Core Experience Summary

*   **Architecture First, Code Second:** The clear division of directories like `compass`, `scripts`, and `docs` at the beginning of the project was the cornerstone for subsequent efficient development and maintenance. A good top-level design is far more important than hastily written code.
*   **The Evolution of a Tool:** The evolution of `hardware_optimizer.py` (from a single script to a product with separated configuration, complete logging, and intelligent algorithms) epitomizes project engineering. This shows that the value of a tool lies not only in completing its task but also in its robustness, usability, and scalability.
*   **Intelligence over Brute-Force:** The upgrade from grid search to Bayesian optimization was a key turning point in the project. It proved that on complex problems, introducing more intelligent algorithms (even if slightly more complex to implement) can lead to orders-of-magnitude improvements in efficiency and is an effective way to resolve performance bottlenecks.
*   **Documentation as a Living Fossil:** Documents like `GEMINI.md` record the entire process from core philosophy to architectural decisions and experience summaries. It is not only a guide for new members but also a repository of the team's collective wisdom and a valuable asset for avoiding repeating past mistakes.

---

### Scene: Conquering a Persistent CUDA Out-of-Memory Error

- **Problem**: A persistent `CUDA out of memory` error that was not resolved by simple hyperparameter tuning. The key symptom was that memory usage grew cumulatively with each batch, indicating a memory leak rather than a simple model-size issue.

- **Analysis (A Multi-Layered Investigation)**:
    My diagnostic process followed a systematic, top-down approach to isolate the root cause.
    1.  **Eliminate Surface-Level Causes**: I first guided the user through the most common solutions: reducing model parameters (`hidden_channels`, `num_layers`), data complexity (`max_atoms`), and data-loading parameters (`visnet_cutoff`). When these failed, it invalidated the initial hypothesis that the model was simply too large.
    2.  **Identify the Leak Pattern**: The user's feedback that the crash happened after a few batches was the critical clue. This confirmed the problem was a **memory leak**, where resources were not being freed after each step.
    3.  **Isolate the Leak Source**: My investigation then focused on features that accumulate data over time. I initially suspected profiler hooks (`record_function`) but eventually pinpointed the `gradient_accumulation_steps` parameter as the primary culprit. Storing gradients for 32 batches on a 6GB GPU was unsustainable. This was a **configuration-induced leak**.

- **Solution (A Robust, Two-Part Fix)**:
    1.  **Address the Root Cause**: I corrected the `gradient_accumulation_steps` to a value appropriate for the hardware, eliminating the primary source of the cumulative memory growth.
    2.  **Implement a Failsafe**: To prevent future, more subtle leaks and to enforce clean memory management, I implemented a definitive safeguard. I modified the training and validation loops in `loop.py` to explicitly `del` the `loss`, `output`, and `data` tensors and then call `torch.cuda.empty_cache()` at the end of every single batch. This forces the GPU to reclaim all unused memory immediately.

- **Lesson Learned**: When debugging memory issues, it's crucial to differentiate between **static memory cost** (the size of the model and a single batch) and **cumulative memory cost** (memory that grows over time). If simplifying the model doesn't work, the next step is to investigate any process that accumulates resources. Furthermore, implementing explicit memory cleanup (`del` and `empty_cache`) provides a powerful safety net, making the entire training process more robust against both obvious and subtle memory leaks.
