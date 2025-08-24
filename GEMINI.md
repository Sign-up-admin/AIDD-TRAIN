# AIDD-TRAIN Hardware Optimizer: The Guiding Philosophy

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
