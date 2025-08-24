import torch
import os
import json
import gc
from itertools import product
import argparse

# 这是一个简化的模型定义，仅用于内存测试，无需加载整个项目代码
class MockViSNet(torch.nn.Module):
    def __init__(self, hidden_channels, num_layers, lmax, vecnorm_type):
        super().__init__()
        # 模拟模型参数量，这是影响显存占用的关键
        self.layers = torch.nn.ModuleList([
            torch.nn.Linear(hidden_channels, hidden_channels) for _ in range(num_layers)
        ])
        self.embedding = torch.nn.Embedding(100, hidden_channels)
        print(f"  - Mock model initialized: {num_layers} layers, {hidden_channels} channels.")

    def forward(self, x, pos, batch):
        # 模拟前向传播
        h = self.embedding(x)
        for layer in self.layers:
            h = layer(h)
        return torch.sum(h) # 返回一个标量以模拟损失计算

def probe_config(config):
    """
    使用给定的配置尝试运行一个训练步骤，测试是否会导致显存溢出。
    """
    print(f"\n--- Probing config: batch_size={config['batch_size']}, hidden_channels={config['visnet_hidden_channels']}, layers={config['visnet_num_layers']} ---")
    try:
        # 1. 创建模型并移至GPU
        model = MockViSNet(
            hidden_channels=config['visnet_hidden_channels'],
            num_layers=config['visnet_num_layers'],
            lmax=2,
            vecnorm_type='max_min'
        ).cuda()

        # 2. 创建模拟数据并移至GPU
        # 模拟一个批次中包含的大分子图
        num_atoms = config['max_atoms']
        batch_size = config['batch_size']
        
        # 模拟原子特征 (整数)
        x = torch.randint(0, 100, (num_atoms * batch_size,)).cuda()
        # 模拟原子坐标 (浮点数)
        pos = torch.randn(num_atoms * batch_size, 3).cuda()
        # 模拟批次索引
        batch = torch.repeat_interleave(torch.arange(batch_size), num_atoms).cuda()
        
        print(f"  - Simulating data: {batch_size} graphs, {num_atoms} atoms each.")

        # 3. 模拟训练步骤 (前向传播 + 反向传播)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-4)
        optimizer.zero_grad()
        
        output = model(x, pos, batch)
        output.backward() # 反向传播是显存占用的主要部分
        optimizer.step()

        print("  - SUCCESS: Configuration is stable.")
        return True

    except torch.cuda.OutOfMemoryError:
        print("  - FAILED: CUDA out of memory. This configuration is too large.")
        return False
    finally:
        # 清理显存，为下一次测试做准备
        gc.collect()
        torch.cuda.empty_cache()

def find_optimal_config(mode_to_optimize):
    """
    主函数，用于搜索给定模式的最佳配置。
    """
    if not torch.cuda.is_available():
        print("CUDA is not available. Aborting hardware optimization.")
        return

    # 获取GPU信息
    gpu_properties = torch.cuda.get_device_properties(0)
    vram_gb = gpu_properties.total_memory / (1024**3)
    print(f"Detected GPU: {gpu_properties.name} with {vram_gb:.2f} GB VRAM.")
    print(f"--- Starting optimization for '{mode_to_optimize}' mode ---")


    # 定义搜索空间 (从大到小，优先保证模型复杂度)
    # 对于验证模式，我们可能希望使用更大的批处理大小
    if mode_to_optimize == 'validation':
        batch_sizes = [32, 24, 16, 12, 8, 4, 2, 1]
        max_atoms_target = 5000 # 验证可能使用较小的分子
    else: # production 或其他模式
        batch_sizes = [16, 12, 8, 6, 4, 2, 1]
        max_atoms_target = 10000

    search_space = {
        'visnet_hidden_channels': [128, 96, 64, 48],
        'visnet_num_layers': [6, 5, 4, 3],
        'batch_size': batch_sizes,
        'max_atoms': [max_atoms_target]
    }

    best_config = None

    # 从最大的模型开始尝试
    param_combinations = product(
        search_space['visnet_num_layers'],
        search_space['visnet_hidden_channels'],
        search_space['batch_size']
    )

    for layers, channels, bs in param_combinations:
        current_config = {
            'visnet_hidden_channels': channels,
            'visnet_num_layers': layers,
            'batch_size': bs,
            'max_atoms': search_space['max_atoms'][0]
        }
        
        if probe_config(current_config):
            best_config = current_config
            print(f"\n>>> Found optimal configuration for '{mode_to_optimize}' mode!")
            break
    
    if best_config:
        # 加载现有的配置文件（如果存在）
        output_path = 'hardware_profile.json'
        if os.path.exists(output_path):
            try:
                with open(output_path, 'r') as f:
                    final_profile = json.load(f)
            except json.JSONDecodeError:
                final_profile = {} # 处理空文件或损坏的文件
        else:
            final_profile = {}

        # 更新或添加当前模式的配置
        final_profile[mode_to_optimize] = best_config

        print(f"\n--- Optimization Complete for '{mode_to_optimize}' ---")
        print(f"Optimal configuration for your hardware:")
        print(json.dumps(best_config, indent=4))
        print(f"Updating settings in '{output_path}'...")

        with open(output_path, 'w') as f:
            json.dump(final_profile, f, indent=4)
        
        print(f"\nRun 'python src/hardware_optimizer.py --mode <another_mode>' to optimize for other modes.")
    else:
        print(f"\n--- Optimization Failed for '{mode_to_optimize}' ---")
        print("Could not find any stable configuration even with the smallest settings.")
        print("Please consider reducing 'max_atoms' or using a smaller model.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find the optimal hardware configuration for a given mode.")
    parser.add_argument(
        '--mode',
        type=str,
        default='production',
        help="The development mode to optimize for (e.g., 'production', 'validation')."
    )
    args = parser.parse_args()
    find_optimal_config(args.mode)