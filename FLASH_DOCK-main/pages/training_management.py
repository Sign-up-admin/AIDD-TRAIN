"""
Training management page for FLASH-DOCK.
"""

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import sys
from pathlib import Path
import time
from datetime import datetime
import json
import threading
import requests

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
# Add FLASH_DOCK-main/services to path
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from compass_client import CompassClient, CompassError


def _format_error_message(e: Exception) -> str:
    """
    Format error message with detailed information if it's a CompassError.

    Args:
        e: Exception to format

    Returns:
        Formatted error message string
    """
    if isinstance(e, CompassError):
        parts = [f"**错误**: {e.message}"]
        if e.error_code:
            parts.append(f"**错误代码**: {e.error_code}")
        if e.status_code:
            parts.append(f"**HTTP状态码**: {e.status_code}")
        if e.detail:
            detail_str = json.dumps(e.detail, indent=2, ensure_ascii=False)
            parts.append(f"**详细信息**:\n```json\n{detail_str}\n```")
        return "\n\n".join(parts)
    else:
        return str(e)


def _handle_stop_task(client: CompassClient, task_id: str, task_status: str):
    """
    Handle stopping a training task with comprehensive error handling and status checking.

    Args:
        client: CompassClient instance
        task_id: Task ID to stop
        task_status: Current task status
    """
    if task_status not in ["running", "initializing"]:
        st.warning(f"任务状态为 {task_status}，无法停止。只有运行中或初始化中的任务可以停止。")
        return

    with st.spinner("正在停止任务..."):
        try:
            # 直接发送停止请求
            result = client.stop_training_task(task_id)

            # 验证响应格式
            if not result or not isinstance(result, dict):
                st.error("**停止任务失败**")
                st.warning("⚠️ 后端返回的响应格式异常，无法确认停止状态")
                # 尝试获取任务状态以确认
                try:
                    latest_task = client.get_training_task(task_id)
                    latest_status = latest_task.get("status")
                    st.info(f"📊 当前任务状态: {latest_status}")
                    if latest_status in ["cancelled", "completed", "failed"]:
                        st.success("✓ 任务实际上已经停止")
                        time.sleep(1)
                        st.rerun()
                    elif latest_status == "cancelling":
                        st.info("⏳ 任务正在取消中，请稍候...")
                except Exception:
                    pass
            else:
                # 检查响应消息，确认后端是否真的接受了停止请求
                response_message = result.get("message", "")
                if (
                    "stopped" not in response_message.lower()
                    and "stop" not in response_message.lower()
                ):
                    st.warning("⚠️ 后端响应消息异常，无法确认停止请求是否被接受")

                # 立即检查任务状态，确认是否真的停止成功
                status_checked = False
                try:
                    updated_task = client.get_training_task(task_id)
                    current_status = updated_task.get("status")

                    # 如果任务状态已经改变，说明停止成功
                    if current_status in ["cancelled", "completed", "failed"]:
                        status_messages = {
                            "cancelled": "已取消",
                            "completed": "已完成",
                            "failed": "已失败",
                        }
                        st.success(
                            f"✓ 任务已成功停止，当前状态: {status_messages.get(current_status, current_status)}"
                        )
                        status_checked = True
                        time.sleep(1)
                        st.rerun()
                    elif current_status == "cancelling":
                        # 任务正在取消中，显示提示并继续轮询
                        st.info("⏳ 任务正在取消中，请稍候...")
                        status_checked = False  # 继续轮询
                    elif current_status not in ["running", "initializing"]:
                        st.warning(f"⚠️ 任务状态已改变为: {current_status}，但非预期的停止状态")
                        status_checked = True
                        time.sleep(1)
                        st.rerun()
                except Exception:
                    # 无法立即检查状态，继续轮询
                    pass

                # 如果状态还未改变，进行轮询检查
                if not status_checked:
                    # 显示提示信息，说明正在验证停止状态
                    status_placeholder = st.empty()
                    status_placeholder.info("⏳ 停止请求已发送，正在验证任务状态...")

                    # 轮询检查任务状态
                    max_poll_time = 30.0
                    poll_interval = 0.5  # 减少轮询间隔，提高响应速度
                    poll_elapsed = 0.0
                    final_status = None
                    last_status = None

                    # 创建进度条显示轮询进度
                    progress_bar = st.progress(0)

                    while poll_elapsed < max_poll_time:
                        try:
                            updated_task = client.get_training_task(task_id)
                            final_status = updated_task.get("status")

                            # 更新进度条和状态显示
                            progress = min(poll_elapsed / max_poll_time, 1.0)
                            progress_bar.progress(progress)

                            # 显示状态信息，包括cancelling状态
                            status_display = {
                                "cancelling": "正在取消中",
                                "running": "运行中",
                                "initializing": "初始化中",
                                "cancelled": "已取消",
                                "completed": "已完成",
                                "failed": "已失败",
                            }.get(final_status, final_status)

                            status_placeholder.info(
                                f"正在检查任务状态... ({poll_elapsed:.1f}s / {max_poll_time:.0f}s) - 当前状态: {status_display}"
                            )

                            if final_status in ["cancelled", "completed", "failed"]:
                                status_messages = {
                                    "cancelled": "已取消",
                                    "completed": "已完成",
                                    "failed": "已失败",
                                }
                                progress_bar.progress(1.0)
                                status_placeholder.success(
                                    f"✓ 任务已成功停止，当前状态: {status_messages.get(final_status, final_status)}"
                                )
                                break
                            elif final_status == "cancelling":
                                # 任务正在取消中，继续轮询
                                if last_status != "cancelling":
                                    status_placeholder.info(
                                        f"⏳ 任务正在取消中... ({poll_elapsed:.1f}s / {max_poll_time:.0f}s)"
                                    )
                                last_status = final_status
                            elif final_status not in ["running", "initializing", "cancelling"]:
                                progress_bar.progress(1.0)
                                status_placeholder.warning(
                                    f"⚠️ 任务状态已改变为: {final_status}，但非预期的停止状态"
                                )
                                break
                        except Exception as poll_e:
                            # 轮询时出错，继续轮询
                            status_placeholder.warning(
                                f"⚠️ 轮询时出错: {str(poll_e)[:50]}... (继续轮询)"
                            )

                        time.sleep(poll_interval)
                        poll_elapsed += poll_interval

                    # 清理进度条
                    progress_bar.empty()

                    # 根据最终状态显示结果
                    if final_status and final_status in ["cancelled", "completed", "failed"]:
                        # 停止成功，已在上面显示成功消息
                        time.sleep(1)
                        st.rerun()
                    elif final_status == "cancelling":
                        status_placeholder.warning(
                            f"⚠️ 任务仍在取消中（已等待 {poll_elapsed:.1f} 秒）。"
                            f"任务可能正在清理资源，请稍后刷新页面查看最新状态。"
                        )
                    elif final_status and final_status in ["running", "initializing"]:
                        status_placeholder.warning(
                            "⚠️ 停止请求已发送，但任务仍在运行。任务可能正在停止中，请稍后刷新页面查看最新状态。"
                        )
                    else:
                        status_placeholder.warning(
                            f"⚠️ 无法确认任务停止状态。当前状态: {final_status if final_status else '未知'}"
                        )

                    time.sleep(1)
                    st.rerun()

        except CompassError as e:
            # 根据错误类型显示不同的消息
            st.error("**停止任务失败**")
            st.markdown(_format_error_message(e))

            # 显示可能的解决建议
            if e.status_code == 400:
                # 400错误通常表示任务状态不允许停止
                error_detail = str(e.detail) if e.detail else str(e.message)
                if "status" in error_detail.lower():
                    st.info("💡 提示: 任务状态可能已改变，请刷新页面查看最新状态")
                else:
                    st.info("💡 提示: 任务可能已经停止或状态已改变，请刷新页面查看最新状态")
                # 尝试获取最新状态
                try:
                    latest_task = client.get_training_task(task_id)
                    latest_status = latest_task.get("status")
                    st.info(f"📊 当前任务状态: {latest_status}")
                except Exception:
                    pass
            elif e.status_code == 404:
                st.info("💡 提示: 任务不存在，请刷新任务列表")
            elif e.status_code == 500:
                st.error("💡 提示: 服务器内部错误，请查看后端日志或联系管理员")
            else:
                st.info("💡 提示: 请检查任务状态或刷新页面")

        except requests.exceptions.Timeout as e:
            # 处理超时错误
            st.warning("**停止请求超时**")
            st.info(f"💡 提示: 请求超时（超时时间: {getattr(e, 'timeout', '未知')}秒）")
            st.info("💡 提示: 停止请求可能仍在处理中，请稍后刷新页面查看任务状态")
            st.info("💡 如果任务仍在运行，可以再次尝试停止")

            # 尝试获取任务状态以确认
            try:
                latest_task = client.get_training_task(task_id)
                latest_status = latest_task.get("status")
                st.info(f"📊 当前任务状态: {latest_status}")
                if latest_status in ["cancelled", "completed", "failed"]:
                    st.success("✓ 任务实际上已经停止")
                    time.sleep(1)
                    st.rerun()
            except Exception:
                st.warning("无法获取最新任务状态，请刷新页面查看")

        except ConnectionError as e:
            # 处理连接错误
            st.error("**无法连接到服务**")
            st.info("💡 提示: 请检查COMPASS服务是否正在运行")
            st.info(f"💡 错误详情: {str(e)}")

        except requests.exceptions.RequestException as e:
            # 处理其他请求异常
            error_str = str(e).lower()
            if "connection" in error_str or "无法连接" in error_str:
                st.error("**无法连接到服务**")
                st.info("💡 提示: 请检查COMPASS服务是否正在运行")
            else:
                st.error(f"**停止失败: {type(e).__name__}**")
                st.exception(e)

        except Exception as e:
            # 处理其他未知异常
            st.error(f"**停止失败: {type(e).__name__}**")
            st.exception(e)
            st.info("💡 提示: 发生未知错误，请查看错误详情或联系管理员")


def _create_terminal_html(terminal_key: str, task_id: str, ws_url: str) -> str:
    """Create HTML with xterm.js terminal and direct WebSocket connection."""
    # Escape WebSocket URL for JavaScript
    ws_url_escaped = json.dumps(ws_url)
    task_id_escaped = json.dumps(task_id)

    # Create unique error reporting key for this terminal instance
    # Note: This key is used in JavaScript only. Python code uses ws_errors_ prefix.
    # The JavaScript will send errors via postMessage to be collected by the parent.

    return f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/xterm@5.3.0/css/xterm.css" />
        <script src="https://cdn.jsdelivr.net/npm/xterm@5.3.0/lib/xterm.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/xterm-addon-fit@0.8.0/lib/xterm-addon-fit.js"></script>
        <style>
            body {{
                margin: 0;
                padding: 10px;
                background-color: #1e1e1e;
                font-family: monospace;
            }}
            #terminal {{
                width: 100%;
                height: 480px;
            }}
            #resources {{
                display: grid;
                grid-template-columns: repeat(4, 1fr);
                gap: 10px;
                padding: 10px;
                background-color: #2d2d2d;
                border-bottom: 1px solid #3e3e3e;
            }}
            .resource-item {{
                background-color: #1e1e1e;
                padding: 10px;
                border-radius: 4px;
                border: 1px solid #3e3e3e;
            }}
            .resource-label {{
                color: #d4d4d4;
                font-size: 12px;
                margin-bottom: 5px;
            }}
            .resource-value {{
                color: #0dbc79;
                font-size: 18px;
                font-weight: bold;
                margin-bottom: 5px;
            }}
            .resource-detail {{
                color: #888;
                font-size: 11px;
            }}
            .progress-bar {{
                width: 100%;
                height: 6px;
                background-color: #3e3e3e;
                border-radius: 3px;
                overflow: hidden;
                margin-top: 5px;
            }}
            .progress-fill {{
                height: 100%;
                background-color: #0dbc79;
                transition: width 0.3s ease;
            }}
            #status {{
                color: #d4d4d4;
                font-size: 12px;
                padding: 5px;
                background-color: #2d2d2d;
                border-bottom: 1px solid #3e3e3e;
            }}
            .status-connected {{
                color: #0dbc79;
            }}
            .status-connecting {{
                color: #e5e510;
            }}
            .status-disconnected {{
                color: #cd3131;
            }}
        </style>
    </head>
    <body>
        <div id="status">状态: <span id="status-text" class="status-connecting">连接中...</span></div>
        <div id="resources">
            <div class="resource-item">
                <div class="resource-label">CPU</div>
                <div class="resource-value" id="cpu-value">0.0%</div>
                <div class="progress-bar"><div class="progress-fill" id="cpu-progress" style="width: 0%"></div></div>
            </div>
            <div class="resource-item">
                <div class="resource-label">内存</div>
                <div class="resource-value" id="memory-value">0.0%</div>
                <div class="resource-detail" id="memory-detail">0.0 GB / 0.0 GB</div>
                <div class="progress-bar"><div class="progress-fill" id="memory-progress" style="width: 0%"></div></div>
            </div>
            <div class="resource-item">
                <div class="resource-label">GPU</div>
                <div class="resource-value" id="gpu-value">不可用</div>
                <div class="resource-detail" id="gpu-detail"></div>
                <div class="progress-bar"><div class="progress-fill" id="gpu-progress" style="width: 0%"></div></div>
            </div>
            <div class="resource-item">
                <div class="resource-label">存储</div>
                <div class="resource-value" id="storage-value">0.0%</div>
                <div class="resource-detail" id="storage-detail">0.0 GB / 0.0 GB</div>
                <div class="progress-bar"><div class="progress-fill" id="storage-progress" style="width: 0%"></div></div>
            </div>
        </div>
        <div id="terminal"></div>
        <script>
            const wsUrl = {ws_url_escaped};
            const taskId = {task_id_escaped};

            // Initialize terminal
            const term = new Terminal({{
                theme: {{
                    background: '#1e1e1e',
                    foreground: '#d4d4d4',
                    cursor: '#aeafad',
                    selection: '#3e3e3e',
                    black: '#1e1e1e',
                    red: '#cd3131',
                    green: '#0dbc79',
                    yellow: '#e5e510',
                    blue: '#2472c8',
                    magenta: '#bc3fbc',
                    cyan: '#11a8cd',
                    white: '#e5e5e5',
                    brightBlack: '#666666',
                    brightRed: '#f14c4c',
                    brightGreen: '#23d18b',
                    brightYellow: '#f5f543',
                    brightBlue: '#3b8eea',
                    brightMagenta: '#d670d6',
                    brightCyan: '#29b8db',
                    brightWhite: '#e5e5e5'
                }},
                fontSize: 12,
                fontFamily: 'Consolas, "Courier New", monospace',
                cursorBlink: true,
                cursorStyle: 'block',
                scrollback: 10000,
                convertEol: true,
                disableStdin: false  // Allow input for future command support
            }});

            const fitAddon = new FitAddon.FitAddon();
            term.loadAddon(fitAddon);

            term.open(document.getElementById('terminal'));
            fitAddon.fit();

            // WebSocket connection
            let ws = null;
            let reconnectDelay = 1000;
            let maxReconnectDelay = 60000;
            let reconnectTimer = null;
            let isReconnecting = false;
            let lastPongTime = Date.now();
            let reconnectAttempts = 0;
            let maxReconnectAttempts = 10;
            const pingInterval = 30000; // 30 seconds
            const pongTimeout = 60000; // 60 seconds

            // Status element
            const statusText = document.getElementById('status-text');

            // Function to report errors to parent window (Streamlit)
            function reportError(errorType, errorMessage, errorDetails) {{
                try {{
                    // Send error to parent window
                    if (window.parent && window.parent !== window) {{
                        window.parent.postMessage({{
                            type: 'websocket_error',
                            taskId: taskId,
                            errorType: errorType,
                            errorMessage: errorMessage,
                            errorDetails: errorDetails,
                            timestamp: new Date().toISOString()
                        }}, '*');
                    }}
                    // Also log to console
                    console.error(`[${{errorType}}] ${{errorMessage}}`, errorDetails);
                }} catch (e) {{
                    console.error('Failed to report error:', e);
                }}
            }}

            function updateStatus(text, className) {{
                statusText.textContent = text;
                statusText.className = className;
            }}

            function updateResources(resources) {{
                // Update CPU
                const cpuPercent = resources.cpu_percent || 0;
                document.getElementById('cpu-value').textContent = cpuPercent.toFixed(1) + '%';
                document.getElementById('cpu-progress').style.width = cpuPercent + '%';

                // Update Memory
                const memory = resources.memory || {{}};
                const memPercent = memory.percent || 0;
                const memUsed = memory.used_gb || 0;
                const memTotal = memory.total_gb || 0;
                document.getElementById('memory-value').textContent = memPercent.toFixed(1) + '%';
                document.getElementById('memory-detail').textContent = memUsed.toFixed(1) + ' GB / ' + memTotal.toFixed(1) + ' GB';
                document.getElementById('memory-progress').style.width = memPercent + '%';

                // Update GPU
                const gpu = resources.gpu || {{}};
                if (gpu.available && gpu.memory) {{
                    const gpuPercent = gpu.memory.allocated_percent || 0;
                    document.getElementById('gpu-value').textContent = gpuPercent.toFixed(1) + '%';
                    document.getElementById('gpu-detail').textContent = gpu.device_name || 'N/A';
                    document.getElementById('gpu-progress').style.width = gpuPercent + '%';
                }} else {{
                    document.getElementById('gpu-value').textContent = '不可用';
                    document.getElementById('gpu-detail').textContent = '';
                    document.getElementById('gpu-progress').style.width = '0%';
                }}

                // Update Storage
                const storage = resources.storage || {{}};
                if (storage.available || storage.total_gb > 0) {{
                    const storagePercent = storage.percent || 0;
                    const storageUsed = storage.used_gb || 0;
                    const storageTotal = storage.total_gb || 0;
                    document.getElementById('storage-value').textContent = storagePercent.toFixed(1) + '%';
                    document.getElementById('storage-detail').textContent = storageUsed.toFixed(1) + ' GB / ' + storageTotal.toFixed(1) + ' GB';
                    document.getElementById('storage-progress').style.width = storagePercent + '%';
                }} else {{
                    document.getElementById('storage-value').textContent = '不可用';
                    document.getElementById('storage-detail').textContent = '';
                    document.getElementById('storage-progress').style.width = '0%';
                }}
            }}

            function connect() {{
                if (isReconnecting) {{
                    return;
                }}

                try {{
                    updateStatus('连接中...', 'status-connecting');
                    term.write('\\r\\n\\x1b[33m[WebSocket] 正在连接到 ' + wsUrl + '...\\x1b[0m\\r\\n');

                    // Chrome-specific: Ensure WebSocket URL is valid
                    if (!wsUrl || wsUrl === 'undefined' || wsUrl === 'null') {{
                        const errorMsg = 'WebSocket URL无效: ' + wsUrl;
                        term.write('\\r\\n\\x1b[31m[错误] ' + errorMsg + '\\x1b[0m\\r\\n');
                        updateStatus('连接失败: URL无效', 'status-disconnected');
                        reportError('invalid_url', errorMsg, {{ url: wsUrl, taskId: taskId }});
                        return;
                    }}

                    // Create WebSocket connection
                    // Chrome may have issues with WebSocket in iframes, so we add error handling
                    try {{
                        ws = new WebSocket(wsUrl);
                    }} catch (wsError) {{
                        const errorMsg = '无法创建WebSocket连接: ' + wsError.message;
                        term.write('\\r\\n\\x1b[31m[错误] ' + errorMsg + '\\x1b[0m\\r\\n');
                        updateStatus('连接失败', 'status-disconnected');
                        reportError('websocket_creation_failed', errorMsg, {{
                            error: wsError.toString(),
                            url: wsUrl,
                            taskId: taskId,
                            browser: navigator.userAgent
                        }});
                        return;
                    }}

                    ws.onopen = function(event) {{
                        updateStatus('已连接', 'status-connected');
                        reconnectDelay = 1000; // Reset reconnect delay on successful connection
                        reconnectAttempts = 0; // Reset reconnect attempts on successful connection
                        isReconnecting = false;
                        lastPongTime = Date.now();
                        term.write('\\r\\n\\x1b[32m[WebSocket] 已连接到训练任务流\\x1b[0m\\r\\n');
                        term.write('\\x1b[32m[WebSocket] 等待任务输出...\\x1b[0m\\r\\n');
                        reportError('connection_success', 'WebSocket连接成功', {{ url: wsUrl, taskId: taskId }});
                    }};

                    ws.onmessage = function(event) {{
                        try {{
                            const message = JSON.parse(event.data);
                            const msgType = message.type;
                            const msgData = message.data;

                            if (msgType === 'log') {{
                                // Write log data directly to terminal (preserves ANSI escape codes)
                                term.write(msgData);
                            }} else if (msgType === 'resources') {{
                                // Update resources display directly
                                updateResources(msgData);
                            }} else if (msgType === 'connected') {{
                                term.write('\\r\\n\\x1b[32m[WebSocket] ' + msgData + '\\x1b[0m\\r\\n');
                            }} else if (msgType === 'ping') {{
                                // Respond to ping
                                ws.send(JSON.stringify({{ type: 'pong', data: 'pong' }}));
                            }} else if (msgType === 'pong') {{
                                lastPongTime = Date.now();
                            }} else if (msgType === 'error') {{
                                term.write('\\r\\n\\x1b[31m[错误] ' + msgData + '\\x1b[0m\\r\\n');
                                updateStatus('错误: ' + msgData, 'status-disconnected');
                                reportError('server_error', '服务器错误', {{
                                    message: msgData,
                                    taskId: taskId,
                                    url: wsUrl
                                }});
                            }}
                        }} catch (e) {{
                            console.error('Error processing message:', e);
                            // If not JSON, treat as plain text log
                            term.write(event.data);
                        }}
                    }};

                    ws.onerror = function(error) {{
                        console.error('WebSocket error:', error);
                        // Chrome may not provide error.message, so we check various properties
                        const errorMsg = error.message || error.reason || 'WebSocket连接错误';
                        const errorDetails = {{
                            message: errorMsg,
                            error: error.toString(),
                            taskId: taskId,
                            url: wsUrl,
                            readyState: ws ? ws.readyState : 'unknown',
                            browser: navigator.userAgent,
                            timestamp: new Date().toISOString()
                        }};

                        updateStatus('连接错误', 'status-disconnected');
                        term.write('\\r\\n\\x1b[31m[WebSocket] 连接错误: ' + errorMsg + '\\x1b[0m\\r\\n');

                        // Provide helpful error messages for common issues
                        if (wsUrl.includes('localhost') || wsUrl.includes('127.0.0.1')) {{
                            term.write('\\x1b[33m[提示] 如果是本地连接，请确保COMPASS服务正在运行\\x1b[0m\\r\\n');
                        }}
                        if (ws && ws.readyState === WebSocket.CONNECTING) {{
                            term.write('\\x1b[33m[提示] 连接正在进行中，如果持续失败，请检查任务状态\\x1b[0m\\r\\n');
                        }}

                        reportError('connection_error', 'WebSocket连接错误', errorDetails);
                    }};

                    ws.onclose = function(event) {{
                        const closeReason = event.reason || '未知原因';
                        const closeCode = event.code || '未知';
                        updateStatus('已断开', 'status-disconnected');
                        term.write('\\r\\n\\x1b[33m[WebSocket] 连接已断开 (代码: ' + closeCode + ', 原因: ' + closeReason + ')\\x1b[0m\\r\\n');

                        // Provide helpful messages for specific close codes
                        if (closeCode === 1006) {{
                            term.write('\\x1b[33m[提示] 连接异常关闭，可能的原因：\\x1b[0m\\r\\n');
                            term.write('\\x1b[33m  1. 任务未启动或已停止\\x1b[0m\\r\\n');
                            term.write('\\x1b[33m  2. COMPASS服务未运行\\x1b[0m\\r\\n');
                            term.write('\\x1b[33m  3. 网络连接问题\\x1b[0m\\r\\n');
                        }} else if (closeCode === 1008) {{
                            term.write('\\x1b[31m[错误] 任务不存在或无法访问\\x1b[0m\\r\\n');
                        }}

                        ws = null;

                        // Report close event
                        reportError('connection_closed', 'WebSocket连接已关闭', {{
                            code: closeCode,
                            reason: closeReason,
                            wasClean: event.wasClean,
                            taskId: taskId,
                            url: wsUrl,
                            timestamp: new Date().toISOString()
                        }});

                        // Only attempt to reconnect if it was an unexpected close
                        // Don't reconnect if the close was clean or if task might not be running
                        const shouldReconnect = !event.wasClean && reconnectAttempts < maxReconnectAttempts;

                        if (shouldReconnect && !isReconnecting) {{
                            reconnectAttempts++;
                            isReconnecting = true;
                            term.write('\\r\\n\\x1b[33m[WebSocket] 正在重连... (尝试 ' + reconnectAttempts + '/' + maxReconnectAttempts + ')\\x1b[0m\\r\\n');
                            reconnectTimer = setTimeout(function() {{
                                reconnectDelay = Math.min(reconnectDelay * 2, maxReconnectDelay);
                                isReconnecting = false;
                                connect();
                            }}, reconnectDelay);
                        }} else if (reconnectAttempts >= maxReconnectAttempts) {{
                            const maxAttemptsMsg = '已达到最大重连次数 (' + maxReconnectAttempts + ')，停止重连';
                            term.write('\\r\\n\\x1b[31m[WebSocket] ' + maxAttemptsMsg + '\\x1b[0m\\r\\n');
                            term.write('\\x1b[33m[提示] 请检查任务状态，确保任务已启动\\x1b[0m\\r\\n');
                            updateStatus('连接失败', 'status-disconnected');
                            reportError('max_reconnect_exceeded', maxAttemptsMsg, {{
                                attempts: reconnectAttempts,
                                taskId: taskId,
                                url: wsUrl
                            }});
                        }} else if (event.wasClean) {{
                            term.write('\\x1b[33m[提示] 连接正常关闭，如果任务正在运行，请刷新页面重新连接\\x1b[0m\\r\\n');
                        }}
                    }};
                }} catch (error) {{
                    console.error('Failed to create WebSocket:', error);
                    const errorMsg = error.message || '无法创建WebSocket连接';
                    updateStatus('连接失败', 'status-disconnected');
                    term.write('\\r\\n\\x1b[31m[错误] ' + errorMsg + '\\x1b[0m\\r\\n');

                    // Report creation error
                    reportError('connection_failed', '无法创建WebSocket连接', {{
                        message: errorMsg,
                        error: error.toString(),
                        stack: error.stack,
                        taskId: taskId,
                        url: wsUrl
                    }});

                    // Retry connection if not exceeded max attempts
                    if (!isReconnecting && reconnectAttempts < maxReconnectAttempts) {{
                        reconnectAttempts++;
                        isReconnecting = true;
                        term.write('\\r\\n\\x1b[33m[WebSocket] 正在重试连接... (尝试 ' + reconnectAttempts + '/' + maxReconnectAttempts + ')\\x1b[0m\\r\\n');
                        reconnectTimer = setTimeout(function() {{
                            reconnectDelay = Math.min(reconnectDelay * 2, maxReconnectDelay);
                            connect();
                        }}, reconnectDelay);
                    }} else if (reconnectAttempts >= maxReconnectAttempts) {{
                        const maxAttemptsMsg = '已达到最大重连次数 (' + maxReconnectAttempts + ')，停止重连';
                        term.write('\\r\\n\\x1b[31m[WebSocket] ' + maxAttemptsMsg + '\\x1b[0m\\r\\n');
                        updateStatus('连接失败', 'status-disconnected');
                        reportError('max_reconnect_exceeded', maxAttemptsMsg, {{
                            attempts: reconnectAttempts,
                            taskId: taskId,
                            url: wsUrl
                        }});
                    }}
                }}
            }}

            // Chrome compatibility: Ensure DOM is fully loaded before connecting
            if (document.readyState === 'loading') {{
                document.addEventListener('DOMContentLoaded', function() {{
                    // Small delay to ensure terminal is fully initialized
                    setTimeout(connect, 100);
                }});
            }} else {{
                // DOM already loaded, connect immediately
                // Small delay to ensure terminal is fully initialized
                setTimeout(connect, 100);
            }}

            // Send periodic ping to keep connection alive
            const pingIntervalId = setInterval(function() {{
                if (ws && ws.readyState === WebSocket.OPEN) {{
                    const now = Date.now();
                    // Check if we haven't received a pong in too long
                    if (now - lastPongTime > pongTimeout) {{
                        console.warn('Pong timeout, closing connection');
                        term.write('\\r\\n\\x1b[33m[WebSocket] 心跳超时，关闭连接\\x1b[0m\\r\\n');
                        ws.close();
                        return;
                    }}
                    // Send ping
                    try {{
                        ws.send(JSON.stringify({{ type: 'ping', data: 'ping' }}));
                    }} catch (e) {{
                        console.error('Failed to send ping:', e);
                        clearInterval(pingIntervalId);
                    }}
                }}
            }}, pingInterval);

            // Clean up ping interval on page unload
            window.addEventListener('beforeunload', function() {{
                if (pingIntervalId) {{
                    clearInterval(pingIntervalId);
                }}
            }});

            // Auto-scroll to bottom when new content is written
            const originalWrite = term.write.bind(term);
            term.write = function(data) {{
                originalWrite(data);
                term.scrollToBottom();
            }};

            // Handle window resize - Chrome may need debouncing
            let resizeTimer = null;
            window.addEventListener('resize', function() {{
                if (resizeTimer) {{
                    clearTimeout(resizeTimer);
                }}
                resizeTimer = setTimeout(function() {{
                    try {{
                        fitAddon.fit();
                    }} catch (e) {{
                        console.error('Error fitting terminal:', e);
                    }}
                }}, 100);
            }});

            // Handle terminal input (for future command support)
            term.onData(function(data) {{
                // Echo input to terminal
                term.write(data);
                // Send command to server (if WebSocket is open)
                if (ws && ws.readyState === WebSocket.OPEN) {{
                    // For now, just log it. Future: send commands to server
                    console.log('Terminal input:', data);
                }}
            }});

            // Cleanup on page unload
            window.addEventListener('beforeunload', function() {{
                if (reconnectTimer) {{
                    clearTimeout(reconnectTimer);
                }}
                if (resizeTimer) {{
                    clearTimeout(resizeTimer);
                }}
                if (pingIntervalId) {{
                    clearInterval(pingIntervalId);
                }}
                if (ws) {{
                    try {{
                        ws.close();
                    }} catch (e) {{
                        console.error('Error closing WebSocket:', e);
                    }}
                }}
            }});

            // Scroll to bottom initially and on content updates
            setTimeout(function() {{
                try {{
                    term.scrollToBottom();
                }} catch (e) {{
                    console.error('Error scrolling terminal:', e);
                }}
            }}, 100);

            // Chrome compatibility: Force terminal to render
            setTimeout(function() {{
                try {{
                    fitAddon.fit();
                    term.focus();
                }} catch (e) {{
                    console.error('Error initializing terminal:', e);
                }}
            }}, 200);
        </script>
    </body>
    </html>
    """


st.title("训练管理")
st.write("管理COMPASS训练任务")

# Initialize client
try:
    client = CompassClient()
    st.success("已连接到COMPASS服务")
except Exception as e:
    st.error(f"无法连接到COMPASS服务: {_format_error_message(e)}")
    st.info("请确保COMPASS服务已启动并注册到服务注册中心")
    st.stop()

# Tabs
tab1, tab2, tab3 = st.tabs(["创建训练任务", "任务列表", "任务详情"])

# Tab 1: Create Training Task
with tab1:
    st.subheader("创建新的训练任务")

    with st.form("create_task_form"):
        # Execution mode
        execution_mode = st.selectbox(
            "执行模式",
            ["validation_tuned", "validation", "prototyping", "smoke_test", "production"],
            help="validation_tuned: 验证调优模式\nvalidation: 验证模式\nprototyping: 原型模式\nsmoke_test: 快速测试\nproduction: 生产模式",
        )

        col1, col2 = st.columns(2)

        with col1:
            # Epochs
            epochs = st.number_input(
                "训练轮数 (Epochs)",
                min_value=1,
                max_value=10000,
                value=100,
                help="训练的总轮数 (1-10000)",
            )

            # Batch size
            batch_size = st.number_input(
                "批次大小 (Batch Size)",
                min_value=1,
                max_value=128,
                value=32,
                help="每个批次的样本数 (1-128)",
            )

        with col2:
            # Learning rate
            learning_rate = st.number_input(
                "学习率 (Learning Rate)",
                min_value=0.0001,
                max_value=1.0,
                value=0.001,
                step=0.0001,
                format="%.4f",
                help="学习率 (0-1.0)",
            )

            # Optimizer
            optimizer = st.selectbox(
                "优化器 (Optimizer)", ["adam", "adamw", "sgd", "rmsprop"], help="选择优化算法"
            )

        # Dataset selection
        try:
            datasets = client.list_datasets()
            dataset_options = ["无（使用默认数据集）"] + [
                f"{ds['dataset_id']} - {ds['name']}"
                for ds in datasets
                if ds.get("status") == "ready"
            ]
            selected_dataset = st.selectbox("选择数据集（可选）", dataset_options)

            if selected_dataset and selected_dataset != "无（使用默认数据集）":
                dataset_id = selected_dataset.split(" - ")[0]
            else:
                dataset_id = None
        except Exception as e:
            st.warning(f"无法加载数据集列表: {_format_error_message(e)}")
            dataset_id = None

        # Description
        description = st.text_area("任务描述（可选）", help="描述此训练任务的用途或目标")

        # Submit button
        submitted = st.form_submit_button("创建训练任务")

        if submitted:
            # Build config
            config = {
                "execution_mode": execution_mode,
                "epochs": int(epochs),
                "batch_size": int(batch_size),
                "learning_rate": float(learning_rate),
                "optimizer": optimizer,
            }

            with st.spinner("正在创建训练任务..."):
                try:
                    task = client.create_training_task(
                        config=config,
                        dataset_id=dataset_id,
                        description=description if description else None,
                    )
                    # Store task info in session_state for display outside form
                    st.session_state["last_created_task"] = task
                    st.session_state["last_created_task_id"] = task["task_id"]
                    st.session_state["selected_task_id"] = task["task_id"]
                    st.session_state["show_task_actions"] = True
                    st.rerun()

                except CompassError as e:
                    # Detailed error display for CompassError
                    st.error("**创建训练任务失败**")
                    with st.expander("查看详细错误信息", expanded=True):
                        st.markdown(_format_error_message(e))
                        if e.original_exception:
                            st.code(str(e.original_exception), language="text")
                except Exception as e:
                    # Generic error display
                    st.error("**创建训练任务失败**")
                    with st.expander("查看详细错误信息", expanded=True):
                        st.write(f"**错误类型**: {type(e).__name__}")
                        st.write(f"**错误消息**: {str(e)}")
                        import traceback

                        st.code(traceback.format_exc(), language="text")

    # Display success message and quick actions outside the form
    if st.session_state.get("show_task_actions", False) and st.session_state.get(
        "last_created_task"
    ):
        task = st.session_state["last_created_task"]
        st.session_state["show_task_actions"] = False  # Reset flag

        st.success("训练任务创建成功！")
        st.info(f"任务ID: {task['task_id']}\n状态: {task['status']}")

        # Quick actions after task creation
        st.subheader("快速操作")
        col1, col2, col3 = st.columns(3)

        with col1:
            if st.button("查看任务详情", key="view_after_create"):
                st.session_state["navigate_to_details"] = True
                st.rerun()

        with col2:
            # Auto-start option for pending tasks
            if task["status"] == "pending":
                if st.button("启动任务", key="start_after_create"):
                    try:
                        client.start_training_task(task["task_id"])
                        st.success("任务已启动！")
                        st.session_state[f"realtime_{task['task_id']}"] = True
                        st.session_state["navigate_to_details"] = True
                        time.sleep(1)
                        st.rerun()
                    except CompassError as e:
                        st.error("**启动任务失败**")
                        with st.expander("查看详细错误信息", expanded=False):
                            st.markdown(_format_error_message(e))
                    except Exception as e:
                        st.error(f"启动失败: {_format_error_message(e)}")
                        st.exception(e)
            else:
                st.info(f"任务状态: {task['status']}")

        with col3:
            if st.button("刷新任务列表", key="refresh_after_create"):
                st.rerun()

        # Auto-navigate to details if requested
        if st.session_state.get("navigate_to_details", False):
            st.session_state["navigate_to_details"] = False
            st.info("💡 提示：任务详情已自动填充，请切换到'任务详情'标签页查看。")

# Tab 2: Task List
with tab2:
    st.subheader("训练任务列表")

    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("刷新列表"):
            st.rerun()

    try:
        tasks = client.list_training_tasks()

        if tasks:
            # Status color mapping
            status_colors = {
                "creating": "🟡",
                "initializing": "🟡",
                "pending": "⚪",
                "running": "🟢",
                "paused": "🟠",
                "completed": "✅",
                "failed": "❌",
                "cancelled": "⛔",
            }

            # Prepare data for display
            display_tasks = []
            for task in tasks:
                status_icon = status_colors.get(task["status"], "❓")
                progress_info = task.get("progress", {})
                current_epoch = progress_info.get("current_epoch", 0)
                total_epochs = progress_info.get(
                    "total_epochs", task.get("config", {}).get("epochs", 0)
                )

                display_tasks.append(
                    {
                        "状态": f"{status_icon} {task['status']}",
                        "任务ID": task["task_id"][:8] + "...",  # Short ID for display
                        "完整ID": task["task_id"],
                        "描述": task.get("description", "无"),
                        "执行模式": task.get("config", {}).get("execution_mode", "N/A"),
                        "轮数": (
                            f"{current_epoch}/{total_epochs}"
                            if total_epochs > 0
                            else f"{task.get('config', {}).get('epochs', 'N/A')}"
                        ),
                        "创建时间": task.get("created_at", "N/A"),
                        "开始时间": task.get("started_at", "N/A"),
                        "完成时间": task.get("completed_at", "N/A"),
                    }
                )

            df = pd.DataFrame(display_tasks)

            # Display table
            st.dataframe(
                df[["状态", "任务ID", "描述", "执行模式", "轮数", "创建时间"]],
                width="stretch",
                hide_index=True,
            )

            # Display action buttons for each task
            st.subheader("任务操作")
            st.markdown("---")

            # Create a task dictionary for quick lookup
            task_dict = {task["task_id"]: task for task in tasks}

            # Display action buttons for each task in a more compact way
            for idx, display_task in enumerate(display_tasks):
                task_id = display_task["完整ID"]
                task = task_dict.get(task_id)
                if not task:
                    continue

                task_status = task["status"]
                short_id = display_task["任务ID"]
                description = display_task["描述"]

                # Create a container for each task's actions
                with st.container():
                    # Use columns to display task info and actions side by side
                    col_info, col_actions = st.columns([3, 2])

                    with col_info:
                        status_icon = status_colors.get(task_status, "❓")
                        st.markdown(
                            f"**{status_icon} {short_id}** | {description} | 状态: `{task_status}`"
                        )

                    with col_actions:
                        action_cols = st.columns(4)

                        # Start button
                        with action_cols[0]:
                            if task_status in ["pending", "paused"]:
                                if st.button(
                                    "▶️ 启动",
                                    key=f"start_task_{task_id}_{idx}",
                                    use_container_width=True,
                                ):
                                    try:
                                        client.start_training_task(task_id)
                                        st.success("任务已启动")
                                        st.session_state[f"realtime_{task_id}"] = True
                                        time.sleep(1)
                                        st.rerun()
                                    except CompassError as e:
                                        st.error("**启动任务失败**")
                                        with st.expander("查看详细错误信息", expanded=False):
                                            st.markdown(_format_error_message(e))
                                    except Exception as e:
                                        st.error(f"启动失败: {_format_error_message(e)}")
                                        st.exception(e)

                        # Stop button
                        with action_cols[1]:
                            if task_status in ["running", "initializing"]:
                                if st.button(
                                    "⏹️ 停止",
                                    key=f"stop_task_{task_id}_{idx}",
                                    use_container_width=True,
                                ):
                                    _handle_stop_task(client, task_id, task_status)

                        # Pause button
                        with action_cols[2]:
                            if task_status == "running":
                                if st.button(
                                    "⏸️ 暂停",
                                    key=f"pause_task_{task_id}_{idx}",
                                    use_container_width=True,
                                ):
                                    try:
                                        client.pause_training_task(task_id)
                                        st.success("任务已暂停")
                                        time.sleep(1)
                                        st.rerun()
                                    except CompassError as e:
                                        st.error("**暂停任务失败**")
                                        with st.expander("查看详细错误信息", expanded=False):
                                            st.markdown(_format_error_message(e))
                                    except Exception as e:
                                        st.error(f"暂停失败: {_format_error_message(e)}")
                                        st.exception(e)

                        # View details button
                        with action_cols[3]:
                            if st.button(
                                "📋 详情",
                                key=f"view_task_{task_id}_{idx}",
                                use_container_width=True,
                            ):
                                st.session_state["selected_task_id"] = task_id
                                st.session_state["active_tab"] = "tab3"
                                st.rerun()

                    st.markdown("---")

            # Task selection for details (keep for backward compatibility)
            task_ids = [task["完整ID"] for task in display_tasks]
            selected_task_id = st.selectbox(
                "选择任务查看详情（或使用上方操作按钮）", task_ids, key="task_list_select"
            )

            if selected_task_id:
                st.session_state["selected_task_id"] = selected_task_id

                # Quick actions
                st.subheader("快速操作")
                col1, col2, col3, col4 = st.columns(4)

                # Find task status
                selected_task = next((t for t in tasks if t["task_id"] == selected_task_id), None)
                if selected_task:
                    task_status = selected_task["status"]

                    with col1:
                        if task_status in ["pending", "paused"]:
                            if st.button("启动任务", key="start_quick"):
                                try:
                                    client.start_training_task(selected_task_id)
                                    st.success("任务已启动")
                                    # Auto-enable real-time terminal for running tasks
                                    st.session_state[f"realtime_{selected_task_id}"] = True
                                    time.sleep(1)
                                    st.rerun()
                                except CompassError as e:
                                    st.error("**启动任务失败**")
                                    with st.expander("查看详细错误信息", expanded=False):
                                        st.markdown(_format_error_message(e))
                                except Exception as e:
                                    st.error(f"启动失败: {_format_error_message(e)}")
                                    st.exception(e)

                    with col2:
                        if task_status == "running" or task_status == "initializing":
                            if st.button("停止任务", key="stop_quick"):
                                _handle_stop_task(client, selected_task_id, task_status)

                    with col3:
                        if task_status == "running":
                            if st.button("暂停任务", key="pause_quick"):
                                try:
                                    client.pause_training_task(selected_task_id)
                                    st.success("任务已暂停")
                                    time.sleep(1)
                                    st.rerun()
                                except CompassError as e:
                                    st.error("**暂停任务失败**")
                                    with st.expander("查看详细错误信息", expanded=False):
                                        st.markdown(_format_error_message(e))
                                except Exception as e:
                                    st.error(f"暂停失败: {_format_error_message(e)}")
                                    st.exception(e)

                    with col4:
                        if st.button("查看详情", key="view_details_quick"):
                            st.session_state["active_tab"] = "tab3"
                            st.rerun()
        else:
            st.info("暂无训练任务")
    except Exception as e:
        st.error(f"获取任务列表失败: {_format_error_message(e)}")

# Tab 3: Task Details
with tab3:
    st.subheader("任务详情")

    # Get task ID from session state or input
    if "selected_task_id" in st.session_state:
        default_task_id = st.session_state["selected_task_id"]
    elif "last_created_task_id" in st.session_state:
        default_task_id = st.session_state["last_created_task_id"]
    else:
        default_task_id = ""

    task_id_input = st.text_input("输入任务ID", value=default_task_id, key="task_detail_input")

    if task_id_input:
        try:
            task = client.get_training_task(task_id_input)

            # Basic info
            col1, col2 = st.columns(2)
            with col1:
                st.metric("任务ID", task["task_id"][:16] + "...")
                st.metric("状态", task["status"])
            with col2:
                if task.get("created_at"):
                    st.metric(
                        "创建时间",
                        (
                            task["created_at"][:19]
                            if len(task["created_at"]) > 19
                            else task["created_at"]
                        ),
                    )
                if task.get("started_at"):
                    st.metric(
                        "开始时间",
                        (
                            task["started_at"][:19]
                            if len(task["started_at"]) > 19
                            else task["started_at"]
                        ),
                    )

            # Description
            if task.get("description"):
                st.info(f"**描述**: {task['description']}")

            # Configuration
            st.subheader("训练配置")
            config = task.get("config", {})
            config_df = pd.DataFrame(
                [
                    {"参数": "执行模式", "值": config.get("execution_mode", "N/A")},
                    {"参数": "训练轮数", "值": config.get("epochs", "N/A")},
                    {"参数": "批次大小", "值": config.get("batch_size", "N/A")},
                    {"参数": "学习率", "值": config.get("learning_rate", "N/A")},
                    {"参数": "优化器", "值": config.get("optimizer", "N/A")},
                ]
            )
            # Convert all values to string to avoid Arrow serialization issues
            config_df["值"] = config_df["值"].astype(str)
            st.dataframe(config_df, width="stretch", hide_index=True)

            # Progress
            st.subheader("训练进度")
            progress = task.get("progress", {})

            if progress:
                # Progress metrics
                current_epoch = progress.get("current_epoch", 0)
                total_epochs = progress.get("total_epochs", config.get("epochs", 0))

                if total_epochs > 0:
                    progress_percent = (current_epoch / total_epochs) * 100
                    st.progress(progress_percent / 100)
                    st.write(
                        f"当前轮数: {current_epoch} / {total_epochs} ({progress_percent:.1f}%)"
                    )

                # Display all progress info
                if len(progress) > 0:
                    st.json(progress)
            else:
                st.info("暂无进度信息")

            # Error info - display task errors prominently
            if task.get("error"):
                st.subheader("❌ 任务错误")
                with st.expander("查看任务错误详情", expanded=True):
                    error_msg = task.get("error", "未知错误")
                    st.error(f"**错误**: {error_msg}")

                    # Try to get more error details if available
                    if isinstance(error_msg, dict):
                        st.json(error_msg)
                    elif isinstance(error_msg, str) and error_msg.startswith("{"):
                        try:
                            error_dict = json.loads(error_msg)
                            st.json(error_dict)
                        except Exception:
                            pass

            # Task control
            st.subheader("任务控制")
            task_status = task["status"]

            col1, col2, col3, col4 = st.columns(4)

            with col1:
                if task_status in ["pending", "paused"]:
                    if st.button("启动任务", key="start_detail"):
                        try:
                            client.start_training_task(task_id_input)
                            st.success("任务已启动")
                            # Auto-enable real-time terminal for running tasks
                            st.session_state[f"realtime_{task_id_input}"] = True
                            time.sleep(1)
                            st.rerun()
                        except CompassError as e:
                            st.error("**启动任务失败**")
                            with st.expander("查看详细错误信息", expanded=False):
                                st.markdown(_format_error_message(e))
                        except Exception as e:
                            st.error(f"启动失败: {_format_error_message(e)}")
                            st.exception(e)

            with col2:
                if task_status == "running" or task_status == "initializing":
                    if st.button("停止任务", key="stop_detail"):
                        _handle_stop_task(client, task_id_input, task_status)

            with col3:
                if task_status == "running":
                    if st.button("暂停任务", key="pause_detail"):
                        try:
                            client.pause_training_task(task_id_input)
                            st.success("任务已暂停")
                            time.sleep(1)
                            st.rerun()
                        except CompassError as e:
                            st.error("**暂停任务失败**")
                            with st.expander("查看详细错误信息", expanded=False):
                                st.markdown(_format_error_message(e))
                        except Exception as e:
                            st.error(f"暂停失败: {_format_error_message(e)}")
                            st.exception(e)

            with col4:
                if st.button("删除任务", key="delete_detail"):
                    try:
                        client.delete_training_task(task_id_input)
                        st.success("任务已删除")
                        time.sleep(1)
                        st.rerun()
                    except CompassError as e:
                        st.error("**删除任务失败**")
                        with st.expander("查看详细错误信息", expanded=False):
                            st.markdown(_format_error_message(e))
                    except Exception as e:
                        st.error(f"删除失败: {_format_error_message(e)}")
                        st.exception(e)

            # Error display section - show errors from WebSocket and other operations
            # Use consistent error key name (ws_errors_)
            error_key = f"ws_errors_{task_id_input}"
            if error_key not in st.session_state:
                st.session_state[error_key] = []

            # Add JavaScript to listen for postMessage from iframe
            # This will capture WebSocket errors from the terminal component
            error_listener_script = f"""
            <script>
            (function() {{
                // Listen for messages from iframe (terminal component)
                window.addEventListener('message', function(event) {{
                    // Check if this is a WebSocket error message from our terminal
                    if (event.data && event.data.type === 'websocket_error') {{
                        const errorData = {{
                            errorType: event.data.errorType || '未知错误',
                            errorMessage: event.data.errorMessage || '无错误消息',
                            errorDetails: event.data.errorDetails || {{}},
                            timestamp: event.data.timestamp || new Date().toISOString()
                        }};

                        // Store error in Streamlit session state via URL parameters
                        // Note: This is a workaround since we can't directly access Streamlit session state from JavaScript
                        // The errors will be collected and displayed when the page reruns
                        console.log('WebSocket error received:', errorData);

                        // Try to notify parent Streamlit app
                        if (window.parent && window.parent !== window) {{
                            window.parent.postMessage({{
                                type: 'streamlit_websocket_error',
                                taskId: '{task_id_input}',
                                errorData: errorData
                            }}, '*');
                        }}
                    }}
                }});
            }})();
            </script>
            """
            components.html(error_listener_script, height=0)

            # Display accumulated errors
            if st.session_state[error_key]:
                st.subheader("⚠️ 错误信息")
                with st.expander("查看错误详情", expanded=True):
                    for idx, error_info in enumerate(st.session_state[error_key]):
                        error_type = error_info.get("errorType", "未知错误")
                        error_message = error_info.get("errorMessage", "无错误消息")
                        error_details = error_info.get("errorDetails", {})
                        timestamp = error_info.get("timestamp", "未知时间")

                        st.error(f"**{error_type}** ({timestamp})")
                        st.write(f"**错误消息**: {error_message}")
                        if error_details:
                            st.json(error_details)
                        st.divider()

                if st.button("清除错误日志", key=f"clear_errors_{task_id_input}"):
                    st.session_state[error_key] = []
                    st.rerun()

            # Real-time terminal and resource monitoring
            st.subheader("实时终端和资源监控")

            # Auto-enable real-time terminal for running tasks
            task_status_for_display = task.get("status", "")
            if task_status_for_display == "running" and not st.session_state.get(
                f"realtime_{task_id_input}", False
            ):
                st.session_state[f"realtime_{task_id_input}"] = True

            # Toggle between traditional logs and real-time terminal
            # Allow terminal to be enabled even for pending tasks (will show connection status)
            use_realtime = st.checkbox(
                "启用实时终端转播",
                value=st.session_state.get(
                    f"realtime_{task_id_input}",
                    task_status_for_display in ["running", "initializing"],
                ),
                key=f"realtime_toggle_{task_id_input}",
                help="使用 xterm.js 终端模拟器实时显示训练输出，支持 ANSI 转义码。即使任务未运行，也可以启用以查看连接状态。",
            )

            if use_realtime:
                st.session_state[f"realtime_{task_id_input}"] = True

                # Validate task status and show appropriate messages
                valid_statuses = ["running", "initializing"]
                if task_status_for_display not in valid_statuses:
                    # Show terminal even for non-running tasks, but with a warning
                    st.warning(
                        f"⚠️ 任务状态为 `{task_status_for_display}`，实时终端将在任务启动后自动连接。"
                    )
                    st.info(
                        "💡 提示：\n"
                        "1. 终端已准备就绪，等待任务启动\n"
                        '2. 点击"启动任务"按钮后，终端将自动连接\n'
                        "3. 连接状态会显示在终端顶部"
                    )

                    # Still show the terminal component so user can see connection attempts
                    # Get WebSocket URL from client
                    try:
                        service_url = client._get_service_url()
                        from urllib.parse import urlparse

                        parsed = urlparse(service_url)
                        ws_scheme = "wss" if parsed.scheme == "https" else "ws"
                        ws_url = f"{ws_scheme}://{parsed.netloc}/api/v1/training/tasks/{task_id_input}/stream"
                    except Exception as e:
                        st.error(f"无法获取服务 URL: {_format_error_message(e)}")
                        st.info("请确保COMPASS服务已启动并注册到服务注册中心。")
                        ws_url = None

                    if ws_url:
                        terminal_key = f"terminal_{task_id_input}"
                        with st.expander("🔧 调试信息", expanded=False):
                            st.code(
                                f"WebSocket URL: {ws_url}\n任务ID: {task_id_input}\n任务状态: {task_status_for_display}"
                            )
                            st.info(
                                "💡 终端将等待任务启动后自动连接。\n"
                                "当前状态: 等待任务状态变为 running 或 initializing"
                            )

                        # Show terminal with waiting message
                        terminal_html = _create_terminal_html(terminal_key, task_id_input, ws_url)
                        components.html(
                            terminal_html,
                            height=600,
                            scrolling=False,
                        )

                        # Add refresh button for terminal
                        col_refresh1, col_refresh2 = st.columns([4, 1])
                        with col_refresh2:
                            if st.button(
                                "🔄 刷新终端",
                                key=f"refresh_terminal_pending_{task_id_input}",
                                help="刷新终端连接",
                            ):
                                st.rerun()
                else:
                    # Task is running or initializing - show terminal normally
                    # Get WebSocket URL from client
                    try:
                        service_url = client._get_service_url()
                        # Convert HTTP URL to WebSocket URL
                        from urllib.parse import urlparse

                        parsed = urlparse(service_url)
                        # Use wss:// for HTTPS, ws:// for HTTP
                        ws_scheme = "wss" if parsed.scheme == "https" else "ws"
                        ws_url = f"{ws_scheme}://{parsed.netloc}/api/v1/training/tasks/{task_id_input}/stream"
                    except Exception as e:
                        st.error(f"无法获取服务 URL: {_format_error_message(e)}")
                        st.info("请确保COMPASS服务已启动并注册到服务注册中心。")
                        ws_url = None

                    if ws_url:
                        # Real-time terminal using xterm.js with direct WebSocket connection
                        # Includes both terminal and resource monitoring in one component
                        terminal_key = f"terminal_{task_id_input}"

                        # Show WebSocket URL for debugging (can be removed in production)
                        with st.expander("🔧 调试信息", expanded=False):
                            st.code(
                                f"WebSocket URL: {ws_url}\n任务ID: {task_id_input}\n任务状态: {task_status_for_display}"
                            )
                            st.info(
                                "💡 如果连接失败，请检查：\n1. COMPASS服务是否正在运行\n2. 任务是否已启动（状态为running或initializing）\n3. 网络连接是否正常\n4. 浏览器控制台是否有错误信息"
                            )

                        # Terminal HTML component with JavaScript WebSocket client
                        terminal_html = _create_terminal_html(terminal_key, task_id_input, ws_url)
                        components.html(
                            terminal_html,
                            height=600,
                            scrolling=False,
                        )

                        # Add refresh button for terminal
                        col_refresh1, col_refresh2 = st.columns([4, 1])
                        with col_refresh2:
                            if st.button(
                                "🔄 刷新终端",
                                key=f"refresh_terminal_running_{task_id_input}",
                                help="刷新终端连接",
                            ):
                                st.rerun()

                        # Note: WebSocket errors are reported via postMessage and will be displayed
                        # in the error section above. The JavaScript code in the HTML component
                        # sends error messages to the parent window.
                    else:
                        st.error("无法建立WebSocket连接。请检查服务状态。")
                        st.info(
                            "**可能的原因：**\n1. COMPASS服务未启动或未注册到服务注册中心\n2. 服务注册中心不可访问\n3. 网络连接问题"
                        )
            else:
                st.session_state[f"realtime_{task_id_input}"] = False
                # Clean up any terminal-related state
                terminal_key = f"terminal_{task_id_input}"
                if terminal_key in st.session_state:
                    del st.session_state[terminal_key]

                # Traditional logs view
                st.subheader("任务日志")

                col1, col2 = st.columns([3, 1])
                with col1:
                    auto_refresh = st.checkbox("自动刷新", value=False, key="auto_refresh_logs")
                with col2:
                    log_limit = st.number_input(
                        "日志行数",
                        min_value=10,
                        max_value=1000,
                        value=100,
                        step=10,
                        key="log_limit",
                    )

                try:
                    logs = client.get_task_logs(task_id_input, limit=log_limit)

                    if logs:
                        log_text = "\n".join(logs)
                        st.text_area(
                            "日志内容", log_text, height=400, key="log_display", disabled=True
                        )
                    else:
                        st.info("暂无日志")

                    if auto_refresh:
                        time.sleep(2)
                        st.rerun()
                except CompassError as e:
                    st.error("**获取日志失败**")
                    with st.expander("查看详细错误信息", expanded=False):
                        st.markdown(_format_error_message(e))
                except Exception as e:
                    st.error(f"获取日志失败: {_format_error_message(e)}")
                    st.exception(e)

            # Metrics
            st.subheader("训练指标")
            try:
                metrics = client.get_task_progress(task_id_input)
                if metrics:
                    st.json(metrics)
                else:
                    st.info("暂无指标数据")
            except CompassError as e:
                st.warning("**获取指标失败**")
                with st.expander("查看详细错误信息", expanded=False):
                    st.markdown(_format_error_message(e))
            except Exception as e:
                st.warning(f"获取指标失败: {_format_error_message(e)}")
                st.exception(e)
        except CompassError as e:
            st.error("**获取任务信息失败**")
            with st.expander("查看详细错误信息", expanded=True):
                st.markdown(_format_error_message(e))
        except Exception as e:
            st.error(f"获取任务信息失败: {_format_error_message(e)}")
            st.exception(e)
