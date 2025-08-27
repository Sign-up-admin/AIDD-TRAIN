import signal
import sys

import torch
from torch.amp import GradScaler
from torch.optim.lr_scheduler import CosineAnnealingLR
from torch.profiler import profile, ProfilerActivity

from .checkpoint import save_checkpoint, create_checkpoint_data, resume_from_checkpoint
from .loop import train_epoch, validate_epoch


class Trainer:
    """
    Main class for handling the model training and validation process.

    This class encapsulates the entire training loop, including learning rate scheduling,
    checkpointing, profiling, and graceful shutdown handling.
    """
    def __init__(self, config, model, train_loader, val_loader, device, logger):
        """
        Initializes the Trainer object.

        Args:
            config (dict): Configuration dictionary with training parameters.
            model (torch.nn.Module): The model to be trained.
            train_loader (torch.utils.data.DataLoader): DataLoader for the training set.
            val_loader (torch.utils.data.DataLoader): DataLoader for the validation set.
            device (torch.device): The device to run training on (e.g., 'cuda' or 'cpu').
            logger (Logger): The logger for recording training progress.
        """
        self.config = config
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.device = device
        self.logger = logger

        self.base_lr = config['learning_rate']
        self.warmup_epochs = self.config.get('warmup_epochs', 0)
        self.optimizer = torch.optim.Adam(model.parameters(), lr=self.base_lr, eps=1e-6, weight_decay=config.get('weight_decay', 1e-5))
        self.scaler = GradScaler()
        self.scheduler = CosineAnnealingLR(self.optimizer, T_max=self.config['epochs'] - self.warmup_epochs, eta_min=1e-7)
        
        self.start_epoch = 1
        self.start_batch = 0
        self.best_val_loss = float('inf')
        self.training_state = {'epoch': 0, 'batch_idx': 0}

    def _save_interrupt_checkpoint(self, epoch, batch_idx):
        """
        Saves the model state when training is interrupted.

        Args:
            epoch (int): The current epoch number.
            batch_idx (int): The current batch index.
        """
        if epoch > 0:
            interrupt_data = create_checkpoint_data(
                epoch, batch_idx, self.model, self.optimizer, self.scaler, 
                self.scheduler, self.best_val_loss, interrupted=True
            )
            save_checkpoint(interrupt_data, self.config['checkpoint_dir'], 'INTERRUPTED.pth.tar', self.logger)
            self.logger.log(f"--- Final state for epoch {epoch}, batch {batch_idx} saved. Exiting. ---")
        else:
            self.logger.log_warning("--- Interruption occurred before training. No state to save. Exiting. ---")

    def _setup_signal_handlers(self):
        """Handles graceful shutdown on SIGTERM or KeyboardInterrupt."""
        def graceful_exit_handler(sig, _frame):
            log_msg = "\n\n--- "
            log_msg += "SIGTERM received" if sig == signal.SIGTERM else "SIGINT (Ctrl+C) received"
            log_msg += ". Saving final state... ---"
            self.logger.log_warning(log_msg)
            self._save_interrupt_checkpoint(self.training_state.get('epoch', 0), self.training_state.get('batch_idx', 0))
            sys.exit(0)

        signal.signal(signal.SIGTERM, graceful_exit_handler)
        signal.signal(signal.SIGINT, graceful_exit_handler)

    def _handle_lr_warmup(self, epoch):
        """
        Adjusts learning rate based on warmup schedule.

        Args:
            epoch (int): The current epoch number.
        """
        if self.warmup_epochs > 0 and epoch <= self.warmup_epochs:
            lr = self.base_lr * (epoch / self.warmup_epochs)
            for param_group in self.optimizer.param_groups:
                param_group['lr'] = lr
            self.logger.log(f"Warmup Epoch {epoch}/{self.warmup_epochs}: Set LR to {lr:.6f}")
        elif epoch == self.warmup_epochs + 1:
            for param_group in self.optimizer.param_groups:
                param_group['lr'] = self.base_lr
            self.logger.log(f"Warmup finished. LR set to base value: {self.base_lr}")

    def _profile_epoch(self):
        """
        Runs one training epoch with the profiler enabled.

        Returns:
            float: The training loss for the profiled epoch.
        """
        self.logger.log("--- Profiling enabled for the first epoch ---")
        with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA], record_shapes=True) as prof:
            train_loss = train_epoch(self.model, self.train_loader, self.optimizer, self.device, self.scaler, self.config, self, self.logger)
        self.logger.log("--- Profiler Results ---")
        self.logger.log(prof.key_averages().table(sort_by="cuda_time_total", row_limit=15))
        return train_loss

    def _save_regular_checkpoint(self, epoch, val_loss):
        """
        Saves a regular end-of-epoch checkpoint and the best model if applicable.

        Args:
            epoch (int): The current epoch number.
            val_loss (float): The validation loss for the current epoch.
        """
        is_best = val_loss < self.best_val_loss
        if is_best:
            self.best_val_loss = val_loss

        checkpoint_data = create_checkpoint_data(
            epoch, 0, self.model, self.optimizer, self.scaler, 
            self.scheduler, val_loss, interrupted=False
        )

        save_checkpoint(checkpoint_data, self.config['checkpoint_dir'], 'checkpoint.pth.tar', self.logger)
        if is_best:
            save_checkpoint(checkpoint_data, self.config['checkpoint_dir'], 'model_best.pth.tar', self.logger)
            self.logger.log(f"-> New best model saved with validation loss: {val_loss:.4f}")

    def _run_epoch(self, epoch):
        """
        Runs a single training and validation epoch.

        Args:
            epoch (int): The current epoch number.
        """
        self.training_state['epoch'] = epoch
        self.training_state['batch_idx'] = 0

        self._handle_lr_warmup(epoch)

        is_profiling_epoch = (epoch == 1 and self.start_batch == 0 and self.config.get('profile', False))
        if is_profiling_epoch:
            train_loss = self._profile_epoch()
        else:
            train_loss = train_epoch(self.model, self.train_loader, self.optimizer, self.device, self.scaler, self.config, self, self.logger)

        self.start_batch = 0  # Reset after any potential resume

        val_loss = validate_epoch(self.model, self.val_loader, self.device, self.logger)

        if epoch > self.warmup_epochs:
            self.scheduler.step()

        current_lr = self.optimizer.param_groups[0]['lr']
        self.logger.log(f"Epoch {epoch:02d}/{self.config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f} | LR: {current_lr:.6f}")

        self._save_regular_checkpoint(epoch, val_loss)

    def run(self):
        """
        Starts the main training loop.

        The loop runs for the number of epochs specified in the config,
        handling signal interrupts, checkpoint resumption, and epoch execution.
        """
        self._setup_signal_handlers()
        resume_from_checkpoint(self)

        self.logger.log("Step 5: Starting training...")
        self.logger.log(f"-> Model initialized with {sum(p.numel() for p in self.model.parameters()):,} parameters.")

        for epoch in range(self.start_epoch, self.config['epochs'] + 1):
            self._run_epoch(epoch)

        self.logger.log("--- Training Finished ---")
