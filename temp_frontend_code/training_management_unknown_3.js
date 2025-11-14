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