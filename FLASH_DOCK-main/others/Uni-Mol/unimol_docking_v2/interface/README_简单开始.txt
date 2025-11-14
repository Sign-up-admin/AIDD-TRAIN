========================================
最简单的开始方法
========================================

方法1: 双击运行（最简单）
----------------------------------------
1. 双击 "进入WSL2.bat"
2. 在 WSL2 中输入以下命令：

   cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface
   bash setup_wsl2_env.sh

3. 等待安装完成（可能需要几分钟）

4. 然后运行：
   conda activate flash_dock
   python3 demo.py --mode single --input-protein 你的文件路径 --input-ligand 你的文件路径 --input-docking-grid 你的文件路径 --output-ligand-name result --output-ligand-dir 输出路径


方法2: 一条命令（在 PowerShell 中）
----------------------------------------
wsl -d Ubuntu-24.04 -e bash -c "cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface && bash setup_wsl2_env.sh"


========================================
就这么简单！
========================================



