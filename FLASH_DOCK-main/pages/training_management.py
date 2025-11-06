"""
Training management page for FLASH-DOCK with real-time progress and comprehensive monitoring.
"""
import streamlit as st
import pandas as pd
import time
from datetime import datetime
import sys
from pathlib import Path
import traceback

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
# Add FLASH_DOCK-main/services to path
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

# Initialize debug monitor
try:
    from debug_monitor import ServiceMonitor, format_diagnostic
    monitor = ServiceMonitor()
except Exception as e:
    st.error(f"Failed to import debug monitor: {e}")
    monitor = None

st.title("训练管理")
st.write("管理和监控COMPASS训练任务")

# Service status monitoring
with st.expander("🔍 服务连接状态", expanded=False):
    if monitor:
        if st.button("检查服务状态"):
            with st.spinner("正在检查服务状态..."):
                diagnostic = monitor.full_diagnostic()
                st.code(format_diagnostic(diagnostic), language="text")
                
                # Show summary
                registry_ok = diagnostic["registry"]["available"]
                compass_available = diagnostic["compass_services"]["count"] > 0
                
                col1, col2 = st.columns(2)
                with col1:
                    if registry_ok:
                        st.success("[OK] 服务注册中心: 正常")
                    else:
                        st.error(f"[FAIL] 服务注册中心: {diagnostic['registry'].get('error', '未知错误')}")
                
                with col2:
                    if compass_available:
                        st.success(f"[OK] COMPASS服务: {diagnostic['compass_services']['count']} 个可用")
                    else:
                        st.error("[FAIL] COMPASS服务: 不可用")
    else:
        st.warning("调试监控功能不可用")

# Initialize client with detailed error handling
client = None
connection_error = None
registry_available = False
service_available = False

try:
    from compass_client import CompassClient
    client = CompassClient()
    
    # Check service availability
    try:
        registry_available = client.service_manager.is_registry_available()
        if registry_available:
            service_available = client.is_service_available()
            if service_available:
                # Test connection by listing tasks
                tasks = client.list_tasks()
                st.success("[OK] 已连接到COMPASS服务")
            else:
                st.warning("[WARNING] 服务注册中心可用，但未找到COMPASS服务")
                st.info("请确保COMPASS服务已启动并注册到服务注册中心")
        else:
            st.error("[FAIL] 服务注册中心不可用")
            st.warning("请检查服务注册中心是否正在运行（端口8500）")
    except ConnectionError as e:
        connection_error = str(e)
        st.error(f"[FAIL] 连接错误: {e}")
        if "注册中心" in str(e) or "registry" in str(e).lower():
            st.warning("💡 请确保服务注册中心正在运行")
            st.info("运行服务注册中心: `python services/registry/server.py --port 8500`")
        else:
            st.warning("请检查服务状态并重试")
        if st.checkbox("显示详细错误", key="conn_error_detail"):
            st.exception(e)
    except TimeoutError as e:
        connection_error = str(e)
        st.error(f"[FAIL] 连接超时: {e}")
        st.warning("服务可能响应缓慢或不可用，请稍后重试")
        if st.checkbox("显示详细错误", key="timeout_error_detail"):
            st.exception(e)
    except Exception as e:
        connection_error = str(e)
        st.error(f"[FAIL] 无法连接到COMPASS服务: {e}")
        st.warning("请检查服务状态并重试")
        if st.checkbox("显示详细错误", key="conn_error_detail"):
            st.exception(e)
except Exception as e:
    connection_error = str(e)
    st.error(f"[FAIL] 初始化COMPASS客户端失败: {e}")
    if st.checkbox("显示详细错误", key="init_error_detail"):
        st.exception(e)

# Show detailed error if connection failed
if connection_error and st.checkbox("显示详细错误信息", key="detailed_error"):
    st.code(traceback.format_exc(), language="python")

# Dataset format info
with st.expander("📋 数据集格式说明", expanded=False):
    st.markdown("""
    ### PDBbind数据集格式要求
    
    COMPASS使用PDBbind数据集格式。如果未上传自定义数据集，系统将自动使用默认的PDBbind-2025.8.4数据集。
    
    **数据集结构要求：**
    ```
    dataset/
    ├── index/
    │   └── INDEX_general_PL.2020R1.lst  (索引文件)
    └── P-L/  (蛋白质-配体对目录)
        ├── 1981-2000/
        │   ├── 1a30/
        │   │   ├── 1a30_protein.pdb
        │   │   └── 1a30_ligand.mol2
        │   └── ...
        └── ...
    ```
    
    **索引文件格式：**
    每行包含：PDB代码、分辨率、年份、结合亲和力等信息
    
    **上传数据集：**
    - 将数据集打包为ZIP或TAR格式
    - 确保包含索引文件和P-L目录结构
    - 上传后系统会自动处理数据
    """)

# Stop if client not available
if not client:
    st.stop()

# Tabs
tab1, tab2, tab3 = st.tabs(["创建任务", "任务列表与监控", "任务详情"])

with tab1:
    st.subheader("创建新的训练任务")
    
    # Dataset selection - Show default PDBbind dataset in dropdown
    st.info("💡 提示：选择要使用的数据集")
    
    dataset_id = None
    use_default_dataset = False
    
    try:
        datasets = client.list_datasets()
        # Add default dataset as first option
        dataset_options = ["默认PDBbind-2025.8.4数据集"] + [f"{ds['name']} ({ds['dataset_id']})" for ds in datasets] if datasets else ["默认PDBbind-2025.8.4数据集"]
        selected = st.selectbox("选择数据集", dataset_options, index=0)
        
        if selected == "默认PDBbind-2025.8.4数据集":
            use_default_dataset = True
            dataset_id = None
        elif datasets:
            # Find the selected dataset
            for ds in datasets:
                if f"{ds['name']} ({ds['dataset_id']})" == selected:
                    dataset_id = ds['dataset_id']
                    use_default_dataset = False
                    break
    except Exception as e:
        st.warning(f"获取数据集列表失败: {e}，将使用默认数据集")
        use_default_dataset = True
        dataset_id = None
    
    with st.form("create_task_form"):
        execution_mode = st.selectbox(
            "执行模式",
            ["validation_tuned", "validation", "prototyping", "smoke_test", "production"],
            help="validation_tuned: 优化的验证模式 | validation: 标准验证模式 | prototyping: 快速原型 | smoke_test: 快速测试 | production: 生产模式"
        )
        
        epochs = st.number_input("训练轮数", min_value=1, max_value=1000, value=200)
        batch_size = st.number_input("批次大小", min_value=1, max_value=32, value=2)
        learning_rate = st.number_input("学习率", min_value=1e-6, max_value=1e-2, value=0.0001, format="%.6f")
        
        description = st.text_area("任务描述（可选）")
        
        submitted = st.form_submit_button("创建任务")
        
        if submitted:
            # Check service availability before creating task
            if not registry_available:
                st.error("❌ 无法创建任务：服务注册中心不可用")
                st.warning("请先启动服务注册中心，然后重试")
                st.info("💡 运行服务注册中心: `python services/registry/server.py --port 8500`")
                st.stop()
            
            if not service_available:
                st.error("❌ 无法创建任务：COMPASS服务不可用")
                st.warning("请确保COMPASS服务已启动并注册到服务注册中心")
                st.info("💡 运行COMPASS服务: `python compass/service_main.py --port 8080`")
                st.stop()
            
            config = {
                "execution_mode": execution_mode,
                "epochs": epochs,
                "batch_size": batch_size,
                "learning_rate": learning_rate
            }
            
            if not use_default_dataset and dataset_id:
                config["dataset_id"] = dataset_id
            
            try:
                # Initialize session state for auto-navigation
                if 'new_task_id' not in st.session_state:
                    st.session_state.new_task_id = None
                if 'switch_to_tab2' not in st.session_state:
                    st.session_state.switch_to_tab2 = False
                
                # Create progress container for better UX
                progress_container = st.empty()
                
                # Create task with timeout and progress indication
                try:
                    with progress_container.container():
                        st.info("🔄 正在验证配置和检查系统资源...")
                        st.info("⏳ 这可能需要几秒钟，请稍候...")
                    
                    # Use longer timeout for task creation (60 seconds)
                    task_id = client.create_training_task(
                        config, 
                        dataset_id=dataset_id if not use_default_dataset else None, 
                        description=description,
                        timeout=60  # 60 second timeout
                    )
                    
                    # Clear progress indicators
                    progress_container.empty()
                    
                    # Store task ID and set flag for auto-start and navigation
                    st.session_state.new_task_id = task_id
                    st.session_state.switch_to_tab2 = True
                    
                    # Auto-start training with fast verification
                    try:
                        with st.spinner("正在启动训练（快速验证模式）..."):
                            # Check connection before starting
                            conn_status = client.check_connection()
                            if not conn_status['compass_service_available']:
                                raise ConnectionError(
                                    "无法启动训练：COMPASS服务连接失败\n" +
                                    (conn_status.get('error_message', '未知错误') or '')
                                )
                            
                            # Start training with fast verification (优化后的验证逻辑，最多2秒)
                            start_result = client.start_training(task_id, verify_start=True, max_retries=3)
                            
                            # 客户端已经做了快速验证，这里直接显示成功消息
                            # 如果验证失败，客户端会记录警告但不阻塞
                            st.success(f"✅ 任务创建成功并已开始训练！")
                            st.info(f"📌 任务ID: `{task_id}`")
                            st.warning("⚠️ 请切换到 **'任务列表与监控'** 标签页查看实时进度和监控信息")
                        
                        # Clear the form by rerunning (立即返回，不等待)
                        st.rerun()
                    except TimeoutError as e:
                        st.error(f"⏱️ 启动训练超时：{str(e)}")
                        st.info(f"📌 任务ID: `{task_id}` - 已创建但未启动")
                        st.markdown("""
                        **可能原因：**
                        - COMPASS服务响应缓慢
                        - 任务初始化时间过长
                        - 网络连接问题
                        
                        **建议：**
                        - 等待片刻后检查任务状态
                        - 切换到'任务列表与监控'标签页手动启动任务
                        - 检查COMPASS服务日志
                        """)
                        if st.checkbox("显示详细错误", key="start_timeout_error_detail"):
                            st.exception(e)
                    except ConnectionError as e:
                        st.error(f"🔌 连接错误：{str(e)}")
                        st.info(f"📌 任务ID: `{task_id}` - 已创建但未启动")
                        
                        # Show connection diagnostics
                        with st.expander("📊 连接诊断信息", expanded=False):
                            try:
                                conn_status = client.check_connection()
                                if conn_status['registry_available']:
                                    st.success("✓ 服务注册中心：可用")
                                else:
                                    st.error(f"✗ 服务注册中心：不可用 - {conn_status.get('error_message', '')}")
                                
                                if conn_status['compass_service_available']:
                                    st.success("✓ COMPASS服务：已发现")
                                    if conn_status.get('compass_service_info'):
                                        st.json(conn_status['compass_service_info'])
                                else:
                                    st.error(f"✗ COMPASS服务：未找到 - {conn_status.get('error_message', '')}")
                                
                                if conn_status.get('connection_test'):
                                    st.success("✓ 连接测试：成功")
                                else:
                                    st.warning("⚠ 连接测试：失败或未执行")
                            except Exception as diag_e:
                                st.error(f"无法获取诊断信息: {diag_e}")
                        
                        st.markdown("""
                        **解决方案：**
                        1. 检查COMPASS服务是否正在运行（端口8080）
                        2. 检查服务注册中心是否可用（端口8500）
                        3. 查看服务状态检查器（上方展开"服务连接状态"）
                        4. 检查服务日志以获取更多信息
                        """)
                        if st.checkbox("显示详细错误", key="start_connection_error_detail"):
                            st.exception(e)
                    except ValueError as e:
                        error_msg = str(e)
                        st.error(f"❌ 启动失败：{error_msg}")
                        st.info(f"📌 任务ID: `{task_id}` - 已创建但未启动")
                        
                        # Check if it's a status verification issue
                        if "无法验证" in error_msg or "状态未在预期时间内变化" in error_msg:
                            st.warning("""
                            **说明：** 启动命令已发送到COMPASS服务，但无法确认任务是否成功启动。
                            
                            **建议：**
                            - 切换到'任务列表与监控'标签页查看任务状态
                            - 如果任务状态仍为pending，尝试手动启动
                            - 检查COMPASS服务日志
                            """)
                        else:
                            st.markdown("""
                            **可能原因：**
                            - 任务状态不正确（不能启动）
                            - 任务已失败或不存在
                            - 服务器返回错误
                            
                            **建议：**
                            - 切换到'任务列表与监控'标签页查看任务详细信息
                            - 检查任务状态和错误信息
                            """)
                        if st.checkbox("显示详细错误", key="start_value_error_detail"):
                            st.exception(e)
                    except Exception as e:
                        st.error(f"❌ 启动训练时发生未知错误: {str(e)}")
                        st.info(f"📌 任务ID: `{task_id}` - 已创建但未启动")
                        st.markdown("""
                        **建议：**
                        - 切换到'任务列表与监控'标签页查看任务状态
                        - 如果任务状态仍为pending，尝试手动启动
                        - 检查服务日志以获取更多信息
                        """)
                        if st.checkbox("显示详细错误", key="start_error_detail"):
                            st.exception(e)
                            st.code(traceback.format_exc(), language="python")
                except TimeoutError as e:
                    progress_container.empty()
                    st.error(f"⏱️ 超时错误：{str(e)}")
                    st.warning("💡 建议：")
                    st.markdown("""
                    - 检查系统资源使用情况（CPU、内存、GPU）
                    - 等待一段时间后重试
                    - 减少并发任务数量
                    - 检查服务日志以获取更多信息
                    """)
                    if st.checkbox("显示详细错误", key="timeout_error_detail"):
                        st.exception(e)
                except ConnectionError as e:
                    progress_container.empty()
                    error_msg = str(e)
                    st.error(f"🔌 连接错误：{error_msg}")
                    
                    # Provide specific guidance based on error message
                    if "注册中心" in error_msg or "registry" in error_msg.lower():
                        st.warning("💡 问题：服务注册中心不可用")
                        st.markdown("""
                        **解决方案：**
                        1. 检查服务注册中心是否正在运行
                        2. 检查端口8500是否被占用
                        3. 运行服务注册中心: `python services/registry/server.py --port 8500`
                        """)
                    elif "COMPASS服务" in error_msg or "compass" in error_msg.lower():
                        st.warning("💡 问题：COMPASS服务不可用")
                        st.markdown("""
                        **解决方案：**
                        1. 检查COMPASS服务是否正在运行
                        2. 检查服务是否已在注册中心注册
                        3. 运行COMPASS服务: `python compass/service_main.py --port 8080`
                        4. 查看服务状态检查器（上方展开"服务连接状态"）
                        """)
                    else:
                        st.warning("💡 建议：")
                        st.markdown("""
                        - 检查COMPASS服务是否正在运行
                        - 检查服务注册中心是否可用
                        - 查看服务状态检查器（上方展开"服务连接状态"）
                        """)
                    if st.checkbox("显示详细错误", key="connection_error_detail"):
                        st.exception(e)
                except ValueError as e:
                    progress_container.empty()
                    error_msg = str(e)
                    # Check if it's a resource error
                    if "资源不足" in error_msg or "资源" in error_msg:
                        st.error(f"⚠️ 系统资源不足：")
                        st.markdown(error_msg.replace("\n", "\n\n"))
                        st.warning("💡 建议：")
                        st.markdown("""
                        - 等待当前运行的任务完成
                        - 停止一些正在运行的任务
                        - 减少新任务的资源配置（如批次大小）
                        - 检查系统资源使用情况
                        """)
                    else:
                        st.error(f"❌ 配置错误：{error_msg}")
                    if st.checkbox("显示详细错误", key="value_error_detail"):
                        st.exception(e)
            except Exception as e:
                st.error(f"❌ 创建任务失败：{str(e)}")
                st.warning("💡 如果问题持续，请检查服务日志或联系管理员")
                if st.checkbox("显示详细错误", key="create_error_detail"):
                    st.exception(e)
                    st.code(traceback.format_exc(), language="python")

with tab2:
    st.subheader("训练任务列表与实时监控")
    
    # Check if we should auto-switch to this tab and highlight a new task
    if st.session_state.get('switch_to_tab2', False):
        new_task_id = st.session_state.get('new_task_id')
        if new_task_id:
            st.info(f"✨ 新任务已创建并启动: {new_task_id}")
            # Clear the flag after showing the message
            st.session_state.switch_to_tab2 = False
    
    # Auto-refresh for running tasks
    auto_refresh = st.checkbox("自动刷新（运行中的任务）", value=True)
    refresh_interval = st.slider("刷新间隔（秒）", min_value=1, max_value=10, value=3) if auto_refresh else None
    
    if st.button("手动刷新列表"):
        st.rerun()
    
    try:
        with st.spinner("正在获取任务列表..."):
            try:
                tasks = client.list_tasks()
            except TimeoutError as e:
                st.error(f"⏱️ 获取任务列表超时：{str(e)}")
                st.warning("💡 请稍后重试或刷新页面")
                tasks = []
            except Exception as e:
                st.error(f"❌ 获取任务列表失败：{str(e)}")
                if st.checkbox("显示详细错误", key="list_tasks_error"):
                    st.exception(e)
                tasks = []
        running_tasks = []  # Initialize outside if block
        
        # Get the new task ID if auto-switching
        new_task_id = st.session_state.get('new_task_id')
        
        if tasks:
            # Filter and display running tasks prominently
            running_tasks = [t for t in tasks if t['status'] == 'running']
            other_tasks = [t for t in tasks if t['status'] != 'running']
            
            if running_tasks:
                st.subheader("🔄 运行中的任务")
                for task in running_tasks:
                    # Highlight new task if it's running
                    is_new_task = new_task_id and task['task_id'] == new_task_id
                    if is_new_task:
                        st.success(f"📌 新任务: {task['task_id']}")
                    
                    with st.container():
                        col1, col2 = st.columns([3, 1])
                        with col1:
                            st.write(f"**任务ID:** {task['task_id']}")
                            if task.get('description'):
                                st.caption(task['description'])
                        
                        with col2:
                            status_color = {"running": "🟢", "pending": "🟡", "completed": "🟢", "failed": "🔴"}
                            st.write(f"{status_color.get(task['status'], '⚪')} {task['status'].upper()}")
                        
                        # Show progress if available
                        progress_info = task.get('progress', {})
                        if progress_info:
                            stage = progress_info.get('stage', 'unknown')
                            progress = progress_info.get('progress', 0.0)
                            message = progress_info.get('message', '')
                            
                            # Progress bar
                            st.progress(progress, text=f"{stage}: {message}")
                            
                            # Resource monitoring
                            try:
                                resources = client.get_task_resources(task['task_id'])
                                if resources:
                                    col1, col2, col3 = st.columns(3)
                                    with col1:
                                        st.metric("CPU使用率", f"{resources.get('cpu_percent', 0):.1f}%")
                                    with col2:
                                        mem_info = resources.get('memory', {})
                                        st.metric("内存使用", f"{mem_info.get('used_gb', 0):.2f} GB / {mem_info.get('total_gb', 0):.2f} GB")
                                    with col3:
                                        gpu_info = resources.get('gpu', {})
                                        if gpu_info.get('available'):
                                            gpu_mem = gpu_info.get('memory', {})
                                            st.metric("GPU内存", f"{gpu_mem.get('allocated_mb', 0):.0f} MB / {gpu_mem.get('total_mb', 0):.0f} MB")
                                        else:
                                            st.metric("GPU", "不可用")
                            except Exception as e:
                                # Silently fail if resource monitoring is not available
                                pass
                            
                            # Detailed progress with console output and resources
                            expander_key = f"task_details_{task['task_id']}"
                            expanded = is_new_task  # Auto-expand new tasks
                            with st.expander("📊 详细监控", expanded=expanded):
                                # Progress details
                                if stage == 'data_processing':
                                    dp_info = progress_info.get('data_processing', {})
                                    completed = dp_info.get('completed', 0)
                                    total = dp_info.get('total', 0)
                                    percentage = dp_info.get('percentage', 0)
                                    st.write(f"数据处理: {completed}/{total} ({percentage:.1f}%)")
                                
                                elif stage == 'training':
                                    train_info = progress_info.get('training', {})
                                    epoch = train_info.get('current_epoch', 0)
                                    total_epochs = train_info.get('total_epochs', 0)
                                    batch = train_info.get('current_batch', 0)
                                    total_batches = train_info.get('total_batches', 0)
                                    train_loss = train_info.get('train_loss', 0.0)
                                    val_loss = train_info.get('val_loss', 0.0)
                                    
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        st.metric("当前轮次", f"{epoch}/{total_epochs}")
                                        st.metric("当前批次", f"{batch}/{total_batches}")
                                    with col2:
                                        if train_loss > 0:
                                            st.metric("训练损失", f"{train_loss:.4f}")
                                        if val_loss > 0:
                                            st.metric("验证损失", f"{val_loss:.4f}")
                                
                                # Elapsed time
                                elapsed = progress_info.get('elapsed_time', 0)
                                if elapsed:
                                    hours = int(elapsed // 3600)
                                    minutes = int((elapsed % 3600) // 60)
                                    seconds = int(elapsed % 60)
                                    st.caption(f"⏱️ 已运行时间: {hours:02d}:{minutes:02d}:{seconds:02d}")
                                
                                # Console output (real-time logs)
                                st.subheader("📝 控制台输出")
                                try:
                                    logs = client.get_task_logs(task['task_id'], limit=500)
                                    if logs:
                                        # Display logs in a code block with auto-scroll
                                        log_text = "\n".join(logs)
                                        st.code(log_text, language=None)
                                    else:
                                        st.info("暂无控制台输出")
                                except Exception as e:
                                    st.warning(f"获取日志失败: {e}")
                                
                                # Resource details
                                st.subheader("💻 资源使用详情")
                                try:
                                    resources = client.get_task_resources(task['task_id'])
                                    if resources:
                                        col1, col2 = st.columns(2)
                                        with col1:
                                            st.write("**CPU和内存**")
                                            st.write(f"- CPU使用率: {resources.get('cpu_percent', 0):.1f}%")
                                            mem_info = resources.get('memory', {})
                                            st.write(f"- 内存使用: {mem_info.get('used_gb', 0):.2f} GB / {mem_info.get('total_gb', 0):.2f} GB ({mem_info.get('percent', 0):.1f}%)")
                                        
                                        with col2:
                                            st.write("**GPU**")
                                            gpu_info = resources.get('gpu', {})
                                            if gpu_info.get('available'):
                                                st.write(f"- 设备: {gpu_info.get('device_name', 'N/A')}")
                                                gpu_mem = gpu_info.get('memory', {})
                                                st.write(f"- 已分配: {gpu_mem.get('allocated_mb', 0):.0f} MB ({gpu_mem.get('allocated_percent', 0):.1f}%)")
                                                st.write(f"- 已保留: {gpu_mem.get('reserved_mb', 0):.0f} MB ({gpu_mem.get('reserved_percent', 0):.1f}%)")
                                                st.write(f"- 总计: {gpu_mem.get('total_mb', 0):.0f} MB")
                                            else:
                                                st.write("- GPU不可用")
                                except Exception as e:
                                    st.warning(f"获取资源信息失败: {e}")
                        
                        st.divider()
            
            # All tasks table
            st.subheader("所有任务")
            df = pd.DataFrame([
                {
                    "任务ID": task['task_id'],
                    "状态": task['status'],
                    "创建时间": task['created_at'],
                    "描述": task.get('description', '')
                }
                for task in tasks
            ])
            st.dataframe(df, width='stretch')
            
            # Task actions
            selected_task_id = st.selectbox("选择任务进行操作", [task['task_id'] for task in tasks])
            
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                if st.button("开始", key="start"):
                    try:
                        # Check connection before starting
                        conn_status = client.check_connection()
                        if not conn_status['compass_service_available']:
                            st.error(f"❌ 无法启动：COMPASS服务不可用")
                            st.warning(f"错误信息：{conn_status.get('error_message', '未知错误')}")
                            with st.expander("📊 连接诊断", expanded=False):
                                st.json(conn_status)
                            st.stop()
                        
                        with st.spinner("正在启动任务（快速验证模式）..."):
                            # Start with fast verification (优化后的验证逻辑，最多2秒)
                            start_result = client.start_training(selected_task_id, timeout=30, verify_start=True, max_retries=3)
                            
                            # 客户端已经做了快速验证，直接显示成功消息
                            st.success(f"✅ 任务 {selected_task_id} 启动命令已发送")
                            st.info("请刷新任务列表查看最新状态")
                        
                        st.rerun()
                    except TimeoutError as e:
                        st.error(f"⏱️ 启动超时：{str(e)}")
                        st.markdown("""
                        **可能原因：**
                        - COMPASS服务响应缓慢
                        - 任务初始化时间过长
                        - 网络连接问题
                        
                        **建议：**
                        - 等待片刻后刷新任务列表
                        - 检查任务状态
                        - 查看COMPASS服务日志
                        """)
                    except ConnectionError as e:
                        st.error(f"🔌 连接错误：{str(e)}")
                        with st.expander("📊 连接诊断", expanded=False):
                            try:
                                conn_status = client.check_connection()
                                st.json(conn_status)
                            except Exception as diag_e:
                                st.error(f"无法获取诊断信息: {diag_e}")
                    except ValueError as e:
                        error_msg = str(e)
                        st.error(f"❌ 启动失败：{error_msg}")
                        
                        # Try to get task status for diagnostics
                        try:
                            task_status = client.get_task_status(selected_task_id)
                            st.info(f"**当前任务状态：** {task_status.get('status')}")
                            if task_status.get('error'):
                                st.warning(f"**错误信息：** {task_status.get('error')}")
                        except Exception:
                            pass
                    except Exception as e:
                        st.error(f"❌ 启动失败: {str(e)}")
                        if st.checkbox("显示详细错误", key="start_button_error"):
                            st.exception(e)
                            st.code(traceback.format_exc(), language="python")
            
            with col2:
                if st.button("停止", key="stop"):
                    try:
                        client.stop_training(selected_task_id)
                        st.success(f"任务 {selected_task_id} 已停止")
                        st.rerun()
                    except Exception as e:
                        st.error(f"停止失败: {e}")
                        st.exception(e)
            
            with col3:
                if st.button("暂停", key="pause"):
                    try:
                        client.pause_training(selected_task_id)
                        st.success(f"任务 {selected_task_id} 已暂停")
                        st.rerun()
                    except Exception as e:
                        st.error(f"暂停失败: {e}")
                        st.exception(e)
            
            with col4:
                if st.button("删除", key="delete"):
                    try:
                        client.delete_task(selected_task_id)
                        st.success(f"任务 {selected_task_id} 已删除")
                        st.rerun()
                    except Exception as e:
                        st.error(f"删除失败: {e}")
                        st.exception(e)
        else:
            st.info("暂无训练任务")
    except Exception as e:
        st.error(f"[FAIL] 获取任务列表失败: {e}")
        if st.checkbox("显示详细错误信息", key="error_details"):
            st.exception(e)
            st.code(traceback.format_exc(), language="python")
    
    # Auto-refresh using Streamlit's built-in mechanism
    if auto_refresh and refresh_interval and running_tasks:
        time.sleep(refresh_interval)
        st.rerun()

with tab3:
    st.subheader("任务详情与日志")
    
    task_id = st.text_input("输入任务ID")
    
    if task_id:
        try:
            with st.spinner("正在获取任务详情..."):
                try:
                    task_status = client.get_task_status(task_id)
                except TimeoutError as e:
                    st.error(f"⏱️ 获取任务详情超时：{str(e)}")
                    st.stop()
                except Exception as e:
                    st.error(f"❌ 获取任务详情失败：{str(e)}")
                    if st.checkbox("显示详细错误", key="get_task_status_error"):
                        st.exception(e)
                    st.stop()
            
            # Task info
            col1, col2 = st.columns(2)
            with col1:
                st.metric("任务ID", task_status['task_id'])
                st.metric("状态", task_status['status'])
            with col2:
                if task_status.get('started_at'):
                    st.metric("开始时间", task_status['started_at'])
                if task_status.get('completed_at'):
                    st.metric("完成时间", task_status['completed_at'])
            
            # Progress section
            if task_status.get('progress'):
                st.subheader("训练进度")
                progress_info = task_status['progress']
                stage = progress_info.get('stage', 'unknown')
                progress = progress_info.get('progress', 0.0)
                message = progress_info.get('message', '')
                
                st.progress(progress, text=f"{stage}: {message}")
                st.json(progress_info)
            
            # Configuration
            with st.expander("训练配置"):
                st.json(task_status.get('config', {}))
            
            # Logs section
            st.subheader("训练日志")
            log_limit = st.slider("显示日志行数", min_value=10, max_value=500, value=100)
            
            if st.button("刷新日志"):
                st.rerun()
            
            try:
                with st.spinner("正在获取日志..."):
                    try:
                        logs = client.get_task_logs(task_id, limit=log_limit)
                    except TimeoutError as e:
                        st.error(f"⏱️ 获取日志超时：{str(e)}")
                        st.warning("💡 日志可能较大，请稍后重试或减少日志行数")
                        logs = []
                    except Exception as e:
                        st.error(f"❌ 获取日志失败: {str(e)}")
                        if st.checkbox("显示详细错误", key="get_logs_error"):
                            st.exception(e)
                        logs = []
                
                if logs:
                    log_text = "\n".join(logs)
                    st.text_area("日志", log_text, height=400, key=f"logs_{task_id}")
                else:
                    st.info("暂无日志")
            except Exception as e:
                st.error(f"获取日志失败: {e}")
                if st.checkbox("显示详细错误", key="logs_general_error"):
                    st.exception(e)
        except Exception as e:
            st.error(f"[FAIL] 获取任务详情失败: {e}")
            if st.checkbox("显示详细错误", key="get_task_error_detail"):
                st.exception(e)
