# Stop all COMPASS services
# 停止 COMPASS 项目所有服务

$ErrorActionPreference = "SilentlyContinue"
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "停止 COMPASS 项目所有服务" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

function Stop-PortService {
    param(
        [int]$Port,
        [string]$ServiceName
    )
    
    Write-Host "[步骤] 停止 $ServiceName (端口 $Port)..." -ForegroundColor Yellow
    
    # Get processes using the port
    $netstat = netstat -ano | Select-String ":$Port" | Select-String "LISTENING"
    
    if ($netstat) {
        foreach ($line in $netstat) {
            $parts = $line -split '\s+'
            $pid = $parts[-1]
            
            if ($pid -match '^\d+$') {
                Write-Host "  停止进程 $pid" -ForegroundColor Gray
                $process = Get-Process -Id $pid -ErrorAction SilentlyContinue
                if ($process) {
                    Stop-Process -Id $pid -Force -ErrorAction SilentlyContinue
                    Write-Host "  [OK] 成功停止进程 $pid" -ForegroundColor Green
                } else {
                    Write-Host "  [INFO] 进程 $pid 已停止或不存在" -ForegroundColor Gray
                }
            }
        }
    } else {
        Write-Host "  [INFO] 端口 $Port 上没有运行的服务" -ForegroundColor Gray
    }
}

# Stop Service Registry (port 8500)
Stop-PortService -Port 8500 -ServiceName "服务注册中心"

# Stop COMPASS Service (port 8080)
Stop-PortService -Port 8080 -ServiceName "COMPASS 服务"

# Stop FLASH-DOCK Frontend (port 8501)
Stop-PortService -Port 8501 -ServiceName "FLASH-DOCK 前端"

Write-Host ""
Write-Host "等待所有进程完全停止..." -ForegroundColor Yellow
Start-Sleep -Seconds 2

Write-Host ""
Write-Host "检查剩余进程..." -ForegroundColor Yellow

$remaining = netstat -ano | Select-String ":8500|:8080|:8501" | Select-String "LISTENING"

if ($remaining) {
    Write-Host "警告：仍有服务在运行：" -ForegroundColor Red
    $remaining | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
} else {
    Write-Host "========================================" -ForegroundColor Green
    Write-Host "所有服务已成功停止！" -ForegroundColor Green
    Write-Host "========================================" -ForegroundColor Green
}

Write-Host ""
Read-Host "按 Enter 键退出"








