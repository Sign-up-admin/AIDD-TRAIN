@echo off
echo ========================================
echo 启动 WSL2 并运行 Uni-Mol 对接预测
echo ========================================
echo.

wsl -d Ubuntu-24.04 bash -c "cd /mnt/e/Qinchaojun/AIDD-TRAIN/FLASH_DOCK-main/others/Uni-Mol/unimol_docking_v2/interface && bash setup_wsl2_env.sh && echo. && echo ======================================== && echo 环境设置完成！ && echo ======================================== && echo. && echo 现在可以运行: && echo conda activate flash_dock && echo python3 demo.py [你的参数]"

pause



