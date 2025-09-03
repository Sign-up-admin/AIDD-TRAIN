"""
Standard Training Recipe

A simple, single-stage training recipe that runs for the total number of epochs.
"""

# A user-friendly name for the recipe, used in the selection menu.
recipe_name = "Standard Training"

def run(trainer):
    """
    Executes the training recipe.

    Args:
        trainer (Trainer): The trainer instance (the "oven") to run epochs.
    """
    trainer.logger.log(f"--- Using {recipe_name} Recipe ---")

    # The entire training process is defined here.
    # The recipe has full control over the trainer.
    for epoch in range(trainer.start_epoch, trainer.config['epochs'] + 1):
        trainer.run_epoch(epoch, current_stage=0) # 0 for standard
