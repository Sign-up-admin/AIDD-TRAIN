'''
Configuration settings for the data processing and model training script.
'''

import os

# --- 开发模式切换 ---
# 选择操作模式。此开关控制所有主要超参数。
# - 'smoke_test': 最小化配置，确保整个流程能无错运行。
#                 用于初始设置和调试代码逻辑，运行时间仅需几分钟。
# - 'prototyping': 轻量级配置，用于快速实验和验证想法。
#                  目标是在较小的数据集上实现快速的反馈循环（每轮数分钟）。
# - 'validation': 中等规模配置，用于在全面训练前验证模型性能。
#                 为显存为 6-8GB 的 GPU 优化。
# - 'production': 用于生成最终结果的全面配置。
#                 为显存大于 8GB 的 GPU 调优。
#
# 可选选项: 'smoke_test', 'prototyping', 'validation', 'production'
DEVELOPMENT_MODE = 'prototyping'


# --- 各模式的超参数设置 ---
# 定义了每种开发模式的核心设置。不同模式下的参数可以根据硬件和实验需求进行调整。
MODES = {
    # 'smoke_test': 用于快速验证代码逻辑和流程是否通畅。
    'smoke_test': {
        'epochs': 2,  # 训练的总轮数。对于 smoke_test，设置为一个很小的值。
        'batch_size': 1,  # 每个 GPU 的批处理大小。对于 smoke_test，设置为 1。
        'gradient_accumulation_steps': 2,  # 梯度累积步数。模拟更大的批次大小，例如 batch_size * steps = 1 * 2 = 2。
        'loader_num_workers': 0,  # 数据加载器的工作进程数。设置为0以避免多进程问题。
        'visnet_hidden_channels': 8,  # ViSNet 模型的隐藏层通道数。必须是注意力头数（默认为8）的倍数。
        'visnet_num_layers': 1,  # ViSNet 模型的层数。值越小，模型越简单。
        'visnet_num_rbf': 8,  # 径向基函数的数量。用于表示原子间的距离。
        'force_data_reprocessing': False,  # 是否强制重新处理数据。在 smoke_test 中通常为 False。
        'profile': False,  # 是否启用性能分析。
        'max_atoms': 1000,  # 此模式下支持的最大原子数。
    },
    # 'prototyping': 推荐用于快速验证想法，目标是每轮训练耗时几分钟。
    'prototyping': {
        'epochs': 50, # 训练的总轮数。
        'batch_size': 2, # 每个 GPU 的批处理大小。根据显存大小调整。
        'gradient_accumulation_steps': 4, # 梯度累积步数。
        'loader_num_workers': 0, # 数据加载器的工作进程数。设置为0以避免多进程问题。
        'visnet_hidden_channels': 32, # ViSNet 模型的隐藏层通道数。必须是注意力头数（默认为8）的倍数。
        'visnet_num_layers': 2, # ViSNet 模型的层数。
        'visnet_num_rbf': 32, # 径向基函数的数量。
        'force_data_reprocessing': False, # 是否强制重新处理数据。
        'profile': False, # 是否启用性能分析。
        'max_atoms': 5000, # 此模式下支持的最大原子数。
    },
    # 'validation': 推荐用于显存约 6-8GB 的 GPU (例如, RTX 3060)。
    'validation': {
        'epochs': 20,
        'batch_size': 4,
        'gradient_accumulation_steps': 8,
        'loader_num_workers': 0,
        'visnet_hidden_channels': 48, # ViSNet 模型的隐藏层通道数。必须是注意力头数（默认为8）的倍数。
        'visnet_num_layers': 3,
        'visnet_num_rbf': 48,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 10000,
    },
    # 'production': 为显存大于 8GB 的 GPU 调优。在 6GB 显存上可能运行缓慢或不稳定。
    'production': {
        'epochs': 100,
        'batch_size': 8,
        'gradient_accumulation_steps': 8,
        'loader_num_workers': 0,
        'visnet_hidden_channels': 64, # ViSNet 模型的隐藏层通道数。必须是注意力头数（默认为8）的倍数。
        'visnet_num_layers': 4,
        'visnet_num_rbf': 64,
        'force_data_reprocessing': False,
        'profile': False,
        'max_atoms': 10000,
    }
}

# --- 基础配置 ---
# 包含所有模式通用的设置。这些是基础参数，会被上面选择的模式中的同名参数覆盖。
CONFIG = {
    # --- 开发模式 ---
    # 当前的开发模式，由 DEVELOPMENT_MODE 变量设置。
    'development_mode': DEVELOPMENT_MODE,

    # --- 数据源 ---
    # 是否强制重新处理数据。如果为 True，将忽略已有的缓存文件，重新生成所有数据。
    # 在更改数据处理逻辑或切换数据集时建议设为 True。
    # 可选值: True, False
    'force_data_reprocessing': False,

    # --- 数据处理 ---
    # 数据处理模式。'strict' 模式下，遇到任何处理错误都会中断程序。未来可能支持 'relaxed' 模式。
    # 可选值: 'strict'
    'data_processing_mode': 'strict',
    # 结构中允许的最大原子数。超过此数量的结构将被忽略，以防止内存溢出。
    'max_atoms': 10000,
    # 保存处理失败案例信息的目录路径。
    'failed_cases_dir': r'failed_cases',

    # --- 路径设置 ---
    # 数据集索引文件的路径。此文件列出了所有需要处理的 PDB ID。
    'index_file': r'index/INDEX_general_PL.2020R1.lst',
    # PDBbind 数据集根目录的路径。
    'dataset_path': r'PDBbind-2025.8.4/P-L/',
    # 存放预处理后数据的目录路径。
    'processed_data_dir': r'processed_data',

    # --- 模型与训练超参数 ---
    # 模型训练的学习率。Adam 优化器通常建议使用 1e-3 到 1e-4 范围内的值。
    'learning_rate': 0.0001,
    # 训练集在整个数据集中的比例。剩余部分 (1.0 - train_split) 将用作验证集。
    # 可选范围: 0.0 - 1.0
    'train_split': 0.8,
    # Dropout 比率，用于模型正则化，防止过拟合。在全连接层后使用。
    # 可选范围: 0.0 - 1.0 (通常为 0.5)
    'dropout_rate': 0.5,
    # 梯度裁剪值。将梯度的范数限制在此值内，用于防止梯度爆炸问题。
    'gradient_clip_val': 1.0,

    # --- ViSNet 模型特定参数 ---
    # 注意：ViSNet 的注意力机制要求 `visnet_hidden_channels` 必须能被 `visnet_num_heads` 整除。
    # `visnet_num_heads` 在模型内部定义，默认为 8。请确保 `visnet_hidden_channels` 是 8 的倍数。
    # 在构建分子图时，每个原子考虑的最大邻居数。值越大，计算成本越高，但可能捕捉更复杂的相互作用。
    'max_num_neighbors': 32,
    # ViSNet 模型中邻居搜索的截断半径（单位：埃）。只有在此半径内的原子才被视为邻居。
    'visnet_cutoff': 5.0,
    # 球面谐波的最大阶数 (l)。l=0: s轨道, l=1: p轨道, l=2: d轨道。
    # lmax=2 表示使用零阶、一阶和二阶球面谐波，可以捕捉更复杂的角度信息。
    # 可选值: 0, 1, 2, ...
    'visnet_lmax': 2,
    # 向量归一化的类型。'max_min' 通过最大最小值进行归一化。'none' 表示不进行归一化。
    # 可选值: 'max_min', 'none'
    'visnet_vecnorm_type': 'max_min',

    # --- 硬件与性能 ---
    # 数据预处理时使用的并行工作进程数。建议设置为 CPU 核心数或稍小的值。
    'processing_num_workers': os.cpu_count() or 16,
    # 是否启用性能分析。启用后会记录时间和内存使用情况，但会稍微降低性能。
    # 可选值: True, False
    'profile': False,
    # 手动触发模型保存的标志文件名。如果此文件存在，将在下一个检查点强制保存模型。
    'manual_save_trigger_file': 'SAVE_NOW.flag',
    # 每处理 N 个批次后保存一次模型。用于在长时间训练中定期备份。
    'save_every_n_batches': 100,

    # --- 可复现性 ---
    # 随机种子。设置一个固定的种子可以确保数据划分、模型初始化等随机过程的结果一致，便于复现实验。
    'seed': 42,

    # --- 调试 ---
    # 是否启用调试模式。在调试模式下，可能会打印更多日志信息。
    # 可选值: True, False
    'debug_mode': True,
}

# --- 合并模式特定设置 ---
# 从 MODES 中获取当前开发模式的设置
mode_settings = MODES.get(DEVELOPMENT_MODE)
if mode_settings is None:
    # 如果设置了无效的模式，则抛出错误
    raise ValueError(f"Invalid DEVELOPMENT_MODE: '{DEVELOPMENT_MODE}'. Please choose from {list(MODES.keys())}.")

# 使用模式特定的设置更新基础配置。
# 如果基础配置和模式设置中有相同的键，模式设置中的值会覆盖基础配置中的值。
CONFIG.update(mode_settings)

# --- 硬件自适应配置加载 ---
# 检查是否存在由 hardware_optimizer.py 生成的优化配置文件。
# 如果存在，并且配置文件中有当前 DEVELOPMENT_MODE 的设置，则加载这些优化后的参数。
try:
    import json
    OPTIMIZED_CONFIG_PATH = 'hardware_profile.json'
    if os.path.exists(OPTIMIZED_CONFIG_PATH):
        print(f"--- Found hardware profile: '{OPTIMIZED_CONFIG_PATH}' ---")
        with open(OPTIMIZED_CONFIG_PATH, 'r') as f:
            optimized_params = json.load()

        # 如果当前模式在优化配置中存在
        if DEVELOPMENT_MODE in optimized_params:
            print(f"--- Applying optimized settings for '{DEVELOPMENT_MODE}' mode. ---")
            # 使用优化后的参数更新当前配置
            CONFIG.update(optimized_params[DEVELOPMENT_MODE])
            print("--- Optimized settings applied. ---")
        else:
            print(f"--- No optimized settings found for '{DEVELOPMENT_MODE}' mode in profile. Using default settings. ---")

except Exception as e:
    print(f"[Warning] Could not load or apply hardware profile: {e}")


# --- 后处理与验证 ---
# 根据配置生成一个唯一的运行名称，用于日志和模型保存。
CONFIG['run_name'] = f"visnet_{DEVELOPMENT_MODE}_lr{CONFIG['learning_rate']}_bs{CONFIG['batch_size']}x{CONFIG['gradient_accumulation_steps']}"

# 确保用于存放处理失败案例的目录存在。
if not os.path.exists(CONFIG['failed_cases_dir']):
    os.makedirs(CONFIG['failed_cases_dir'])
