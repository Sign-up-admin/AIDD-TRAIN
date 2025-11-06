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

try:
    from compass_client import CompassClient
    client = CompassClient()
    
    # Test connection
    try:
        tasks = client.list_tasks()
        st.success("[OK] 已连接到COMPASS服务")
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
if connection_error and st.checkbox("显示详细错误信息"):
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
    
    # Dataset selection
    st.info("💡 提示：如果不上传数据集，系统将使用默认的PDBbind-2025.8.4数据集")
    
    use_default_dataset = st.checkbox("使用默认PDBbind-2025.8.4数据集", value=True)
    dataset_id = None
    
    if not use_default_dataset:
        try:
            datasets = client.list_datasets()
            if datasets:
                dataset_options = ["无"] + [f"{ds['name']} ({ds['dataset_id']})" for ds in datasets]
                selected = st.selectbox("选择数据集", dataset_options)
                if selected != "无":
                    dataset_id = datasets[dataset_options.index(selected) - 1]['dataset_id']
            else:
                st.warning("暂无可用数据集，将使用默认数据集")
                use_default_dataset = True
        except Exception as e:
            st.error(f"获取数据集列表失败: {e}")
            use_default_dataset = True
    
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
            config = {
                "execution_mode": execution_mode,
                "epochs": epochs,
                "batch_size": batch_size,
                "learning_rate": learning_rate
            }
            
            if not use_default_dataset and dataset_id:
                config["dataset_id"] = dataset_id
            
            try:
                with st.spinner("正在创建任务..."):
                    task_id = client.create_training_task(
                        config, 
                        dataset_id=dataset_id if not use_default_dataset else None, 
                        description=description
                    )
                st.success(f"[OK] 任务创建成功！任务ID: {task_id}")
                
                # Auto-start option
                auto_start = st.checkbox("立即开始训练", value=True, key="auto_start")
                if auto_start:
                    try:
                        with st.spinner("正在启动训练..."):
                            client.start_training(task_id)
                        st.success(f"[OK] 任务 {task_id} 已开始训练")
                        st.info("请切换到'任务列表与监控'标签页查看实时进度")
                        time.sleep(1)
                        st.rerun()
                    except Exception as e:
                        st.error(f"启动训练失败: {e}")
                        st.exception(e)
            except Exception as e:
                st.error(f"[FAIL] 创建任务失败: {e}")
                if st.checkbox("显示详细错误", key="create_error_detail"):
                    st.exception(e)

with tab2:
    st.subheader("训练任务列表与实时监控")
    
    # Auto-refresh for running tasks
    auto_refresh = st.checkbox("自动刷新（运行中的任务）", value=True)
    refresh_interval = st.slider("刷新间隔（秒）", min_value=1, max_value=10, value=3) if auto_refresh else None
    
    if st.button("手动刷新列表"):
        st.rerun()
    
    try:
        with st.spinner("正在获取任务列表..."):
            tasks = client.list_tasks()
        running_tasks = []  # Initialize outside if block
        
        if tasks:
            # Filter and display running tasks prominently
            running_tasks = [t for t in tasks if t['status'] == 'running']
            other_tasks = [t for t in tasks if t['status'] != 'running']
            
            if running_tasks:
                st.subheader("🔄 运行中的任务")
                for task in running_tasks:
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
                            
                            # Detailed progress
                            with st.expander("详细进度", expanded=True):
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
                        client.start_training(selected_task_id)
                        st.success(f"任务 {selected_task_id} 已开始")
                        st.rerun()
                    except Exception as e:
                        st.error(f"启动失败: {e}")
                        st.exception(e)
            
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
                task_status = client.get_task_status(task_id)
            
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
                    logs = client.get_task_logs(task_id, limit=log_limit)
                
                if logs:
                    log_text = "\n".join(logs)
                    st.text_area("日志", log_text, height=400, key=f"logs_{task_id}")
                else:
                    st.info("暂无日志")
            except Exception as e:
                st.error(f"获取日志失败: {e}")
                st.exception(e)
        except Exception as e:
            st.error(f"[FAIL] 获取任务详情失败: {e}")
            if st.checkbox("显示详细错误", key="get_task_error_detail"):
                st.exception(e)
