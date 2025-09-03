"""
Two-Stage Diffusion Training Recipe

A recipe that divides the training into two distinct stages.
"""

# A user-friendly name for the recipe, used in the selection menu.
recipe_name = "Two-Stage Diffusion"

def run(trainer):
    """
    Executes the training recipe.

    Args:
        trainer (Trainer): The trainer instance (the "oven") to run epochs.
    """
    trainer.logger.log(f"--- Using {recipe_name} Recipe ---")
    
    total_epochs = trainer.config['epochs']
    stage1_epochs = total_epochs // 2

    # --- Stage 1 --- #
    if trainer.config.get('diffusion', {}).get('stage1', {}).get('enabled', False):
        trainer.logger.log("--- Starting Diffusion Stage 1 --- ")
        for epoch in range(trainer.start_epoch, stage1_epochs + 1):
            trainer.run_epoch(epoch, current_stage=1)
        trainer.start_epoch = stage1_epochs + 1 # Set for stage 2

    # --- Stage 2 --- #
    if trainer.config.get('diffusion', {}).get('stage2', {}).get('enabled', False):
        trainer.logger.log("--- Starting Diffusion Stage 2 --- ")
        for epoch in range(trainer.start_epoch, total_epochs + 1):
            trainer.run_epoch(epoch, current_stage=2)
