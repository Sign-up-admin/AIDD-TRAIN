import torch
import gc
import time
import logging

from torch_geometric.data import Batch # type: ignore

from ...training.model import ViSNetPDB

logger = logging.getLogger("HardwareOptimizer")
progress_logger = logging.getLogger("ProgressBar")


def _setup_probe_environment(config, processed_test_data):
    """
    Sets up the environment for a single probe, including model, data, and optimizer.

    Args:
        config (dict): The configuration for the probe, including model and batch size.
        processed_test_data (Data): The pre-processed data sample to use for the probe.

    Returns:
        tuple: A tuple containing the initialized model, a batch of data, and the optimizer.
    """
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ViSNetPDB(
        hidden_channels=config['visnet_hidden_channels'],
        num_layers=config['visnet_num_layers'],
        lmax=2, vecnorm_type='max_min'
    ).to(device)
    data_list = [processed_test_data] * config['batch_size']
    batch = Batch.from_data_list(data_list).to(device)  # type: ignore
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
    return model, batch, optimizer


def probe_config(config, processed_test_data, stress_iterations=20, prefix=""):
    """
    Probes a given configuration for stability and performance.

    Args:
        config (dict): The configuration to probe.
        processed_test_data (Data): The data sample to use for the probe.
        stress_iterations (int): The number of iterations to run the probe for.
        prefix (str): A prefix for the log message.

    Returns:
        tuple: A tuple containing a boolean indicating success and the average time per iteration.
    """
    base_text = f"{prefix}--- Probing: batch_size={config['batch_size']:<3}, layers={config['visnet_num_layers']}, channels={config['visnet_hidden_channels']:<3} for {stress_iterations} iterations..."
    bar_length = 20
    model, batch, optimizer = None, None, None
    try:
        model, batch, optimizer = _setup_probe_environment(config, processed_test_data)
        torch.cuda.synchronize()
        start_time = time.time()
        for i in range(stress_iterations):
            optimizer.zero_grad(); output = model(batch); loss = output.sum() * 100; loss.backward(); optimizer.step()
            progress = (i + 1) / stress_iterations
            filled_length = int(bar_length * progress)
            bar = '█' * filled_length + '-' * (bar_length - filled_length)
            progress_logger.info(f'\r{base_text} [{bar}] {progress: >4.0%}')
        torch.cuda.synchronize()
        end_time = time.time()
        avg_time_per_iter = (end_time - start_time) / stress_iterations if stress_iterations > 0 else float('inf')
        logger.info(f'\r{base_text} SUCCESS' + ' ' * (bar_length + 10))
        return True, avg_time_per_iter
    except torch.cuda.OutOfMemoryError:
        logger.warning(f'\r{base_text} FAILED (OOM)' + ' ' * (bar_length + 10))
        return False, float('inf')
    finally:
        if model is not None: del model, batch, optimizer; gc.collect(); torch.cuda.empty_cache()


def find_max_batch_size_by_stressing(base_config, processed_test_data, start_batch_size, stress_iters):
    """
    Finds the maximum stable batch size for a given base configuration.

    Args:
        base_config (dict): The base configuration (model architecture).
        processed_test_data (Data): The data sample to use for stress testing.
        start_batch_size (int): The initial batch size to start the search from.
        stress_iters (int): The number of iterations for the stress test.

    Returns:
        tuple: A tuple containing the maximum stable batch size and the time per iteration at that size.
    """
    logger.debug("--- Strategy: Finding max batch size and associated cycle time ---")
    logger.debug("[1/2] Finding OOM ceiling...")
    probe_bs = start_batch_size
    while True:
        config = {**base_config, 'batch_size': probe_bs}
        success, _ = probe_config(config, processed_test_data, stress_iterations=5, prefix="	")
        if success:
            logger.debug(f"> OK. Trying batch_size={probe_bs * 2}")
            probe_bs = probe_bs * 2 if probe_bs > 1 else 2
        else:
            oom_ceiling = probe_bs
            logger.debug(f"> OOM at batch_size={oom_ceiling}. Ceiling found."); break
    
    logger.debug(f"[2/2] Binary searching from ceiling ({oom_ceiling}) to find stable edge...")
    bs_candidate, time_at_candidate = 0, float('inf')
    low, high = 1, oom_ceiling - 1
    while low <= high:
        mid = (low + high) // 2
        if mid == 0: break
        config = {**base_config, 'batch_size': mid}
        success, avg_time = probe_config(config, processed_test_data, stress_iterations=stress_iters, prefix="	")
        if success:
            bs_candidate, time_at_candidate = mid, avg_time
            low = mid + 1
            logger.debug(f"> Stable at batch_size={mid}. Avg time/iter: {avg_time:.4f}s. Trying higher.")
        else:
            high = mid - 1
            logger.debug(f"> OOM at batch_size={mid}. Trying lower.")

    if bs_candidate == 0:
        logger.warning("Could not find any stable batch size for this architecture.")
        return None, float('inf')
    
    logger.debug(f"> Max stable batch size found: {bs_candidate}.")
    return bs_candidate, time_at_candidate
