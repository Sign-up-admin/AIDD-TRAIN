import os
import sys
import signal
import torch
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.amp import GradScaler
from torch.profiler import profile, ProfilerActivity

from .loop import train_epoch, validate_epoch


def _save_checkpoint(state, directory, filename, logger=None):
    """Saves checkpoint to file."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    torch.save(state, filepath)
    if logger:
        logger.log(f"Checkpoint saved to {filepath}")


def _load_checkpoint(directory, device, logger=None):
    """Loads the most recent checkpoint from a directory."""
    if not os.path.isdir(directory):
        if logger: logger.log(f"=> No checkpoint directory found at '{directory}'")
        return None

    checkpoint_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.pth.tar')]
    if not checkpoint_files:
        if logger: logger.log(f"=> No checkpoint files (.pth.tar) found in '{directory}'")
        return None

    latest_checkpoint_path = max(checkpoint_files, key=os.path.getmtime)
    if not latest_checkpoint_path:
        if logger: logger.log(f"=> No valid checkpoints found.")
        return None

    if logger: logger.log(f"=> Attempting to load latest checkpoint: '{os.path.basename(latest_checkpoint_path)}'")
    try:
        checkpoint = torch.load(latest_checkpoint_path, map_location=device)
        if logger: logger.log(f"=> Successfully loaded checkpoint from epoch {checkpoint.get('epoch', 'N/A')}.")
        return checkpoint
    except Exception as e:
        if logger:
            logger.log_error(f"Error loading checkpoint '{os.path.basename(latest_checkpoint_path)}': {e}")
        return None


class Trainer:
    def __init__(self, config, model, train_loader, val_loader, device, logger):
        self.config = config
        self.model = model
        self.train_loader = train_loader
        self.val_loader = val_loader
        self.device = device
        self.logger = logger

        self.optimizer = torch.optim.Adam(model.parameters(), lr=config['learning_rate'], eps=1e-7, weight_decay=1e-5)
        self.scaler = GradScaler()
        self.scheduler = ReduceLROnPlateau(self.optimizer, 'min', patience=3, factor=0.5)

        self.start_epoch = 1
        self.start_batch = 0
        self.best_val_loss = float('inf')
        self.training_state = {'epoch': 0, 'batch_idx': 0}

    def _setup_signal_handlers(self):
        """Handles graceful shutdown on SIGTERM or KeyboardInterrupt."""
        def graceful_exit_handler(sig, _frame):
            log_msg = "\n\n--- "
            if sig == signal.SIGTERM:
                log_msg += "SIGTERM received"
            else: # SIGINT
                log_msg += "SIGINT (Ctrl+C) received"
            log_msg += ". Saving final state... ---"
            self.logger.log_warning(log_msg)

            epoch_to_save = self.training_state.get('epoch', 0)
            batch_idx_to_save = self.training_state.get('batch_idx', 0)

            if epoch_to_save > 0:
                interrupt_checkpoint_data = {
                    'epoch': epoch_to_save, 'batch_idx': batch_idx_to_save,
                    'model_state_dict': self.model.state_dict(),
                    'optimizer_state_dict': self.optimizer.state_dict(),
                    'scaler_state_dict': self.scaler.state_dict(),
                    'scheduler_state_dict': self.scheduler.state_dict(),
                    'val_loss': self.best_val_loss, 'interrupted': True
                }
                _save_checkpoint(interrupt_checkpoint_data, self.config['checkpoint_dir'], 'INTERRUPTED.pth.tar', self.logger)
                self.logger.log(f"--- Final state for epoch {epoch_to_save}, batch {batch_idx_to_save} saved. Exiting. ---")
            else:
                self.logger.log_warning("--- Interruption occurred before training. No state to save. Exiting. ---")

            sys.exit(0)

        signal.signal(signal.SIGTERM, graceful_exit_handler)
        signal.signal(signal.SIGINT, graceful_exit_handler)

    def _resume_from_checkpoint(self):
        checkpoint = _load_checkpoint(self.config['checkpoint_dir'], self.device, self.logger)
        if not checkpoint:
            self.logger.log("--> No checkpoint found, starting from scratch.")
            return

        try:
            self.model.load_state_dict(checkpoint['model_state_dict'])
            self.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

            if checkpoint.get('interrupted', False):
                self.start_epoch = checkpoint['epoch']
                self.start_batch = checkpoint.get('batch_idx', 0)
                self.logger.log_warning(f"--> Resuming from an interrupted run. Restarting epoch {self.start_epoch} from batch {self.start_batch}.")
            else:
                self.start_epoch = checkpoint['epoch'] + 1
                self.start_batch = 0

            self.best_val_loss = checkpoint.get('val_loss', float('inf'))
            if 'scaler_state_dict' in checkpoint:
                self.scaler.load_state_dict(checkpoint['scaler_state_dict'])
                self.logger.log("--> Resumed GradScaler state.")
            if 'scheduler_state_dict' in checkpoint:
                self.scheduler.load_state_dict(checkpoint['scheduler_state_dict'])
                self.logger.log("--> Resumed LR Scheduler state.")

            self.logger.log(f"--> Checkpoint loaded. Starting at epoch {self.start_epoch}. Best val loss: {self.best_val_loss:.4f}")

        except (KeyError, RuntimeError) as e:
            self.logger.log_error(f"Could not load checkpoint due to an error: {e}. Starting from scratch.")
            self.start_epoch = 1
            self.start_batch = 0
            self.best_val_loss = float('inf')

    def run(self):
        self._setup_signal_handlers()
        self._resume_from_checkpoint()

        self.logger.log("Step 5: Starting training...")
        self.logger.log(f"-> Model initialized with {sum(p.numel() for p in self.model.parameters()):,} parameters.")

        for epoch in range(self.start_epoch, self.config['epochs'] + 1):
            self.training_state['epoch'] = epoch
            self.training_state['batch_idx'] = 0

            is_profiling_epoch = (epoch == 1 and self.start_batch == 0 and self.config.get('profile', False))

            if is_profiling_epoch:
                self.logger.log("--- Profiling enabled for the first epoch ---")
                with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA], record_shapes=True) as prof:
                    train_loss = train_epoch(self.model, self.train_loader, self.optimizer, self.device, self.scaler, self.config, self, self.logger)
                self.logger.log("--- Profiler Results ---")
                self.logger.log(prof.key_averages().table(sort_by="cuda_time_total", row_limit=15))
            else:
                train_loss = train_epoch(self.model, self.train_loader, self.optimizer, self.device, self.scaler, self.config, self, self.logger)

            # Reset start_batch for subsequent epochs after a resume
            self.start_batch = 0

            val_loss = validate_epoch(self.model, self.val_loader, self.device, self.logger)

            self.scheduler.step(val_loss)

            self.logger.log(f"Epoch {epoch:02d}/{self.config['epochs']} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

            is_best = val_loss < self.best_val_loss
            if is_best:
                self.best_val_loss = val_loss

            checkpoint_data = {
                'epoch': epoch,
                'model_state_dict': self.model.state_dict(),
                'optimizer_state_dict': self.optimizer.state_dict(),
                'scaler_state_dict': self.scaler.state_dict(),
                'scheduler_state_dict': self.scheduler.state_dict(),
                'val_loss': val_loss,
            }

            _save_checkpoint(checkpoint_data, self.config['checkpoint_dir'], 'checkpoint.pth.tar', self.logger)
            if is_best:
                _save_checkpoint(checkpoint_data, self.config['checkpoint_dir'], 'model_best.pth.tar', self.logger)
                self.logger.log(f"-> New best model saved with validation loss: {val_loss:.4f}")

        self.logger.log("--- Training Finished ---")