"""
Training management page for FLASH-DOCK.
"""

import streamlit as st
import pandas as pd
import sys
from pathlib import Path
import time
from datetime import datetime

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
# Add FLASH_DOCK-main/services to path
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from compass_client import CompassClient

st.title("训练管理")
st.write("管理COMPASS训练任务")

# Initialize client
try:
    client = CompassClient()
    st.success("已连接到COMPASS服务")
except Exception as e:
    st.error(f"无法连接到COMPASS服务: {e}")
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
            st.warning(f"无法加载数据集列表: {e}")
            dataset_id = None

        # Description
        description = st.text_area("任务描述（可选）", help="描述此训练任务的用途或目标")

        # Submit button
        submitted = st.form_submit_button("创建训练任务", use_container_width=True)

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
                    st.success(f"训练任务创建成功！")
                    st.info(f"任务ID: {task['task_id']}\n状态: {task['status']}")
                    st.session_state["last_created_task_id"] = task["task_id"]
                except Exception as e:
                    st.error(f"创建训练任务失败: {e}")

# Tab 2: Task List
with tab2:
    st.subheader("训练任务列表")

    col1, col2 = st.columns([1, 4])
    with col1:
        if st.button("刷新列表", use_container_width=True):
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
                use_container_width=True,
                hide_index=True,
            )

            # Task selection for details
            task_ids = [task["完整ID"] for task in display_tasks]
            selected_task_id = st.selectbox("选择任务查看详情", task_ids, key="task_list_select")

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
                            if st.button("启动任务", use_container_width=True, key="start_quick"):
                                try:
                                    client.start_training_task(selected_task_id)
                                    st.success("任务已启动")
                                    time.sleep(1)
                                    st.rerun()
                                except Exception as e:
                                    st.error(f"启动失败: {e}")

                    with col2:
                        if task_status == "running":
                            if st.button("停止任务", use_container_width=True, key="stop_quick"):
                                try:
                                    client.stop_training_task(selected_task_id)
                                    st.success("任务已停止")
                                    time.sleep(1)
                                    st.rerun()
                                except Exception as e:
                                    st.error(f"停止失败: {e}")

                    with col3:
                        if task_status == "running":
                            if st.button("暂停任务", use_container_width=True, key="pause_quick"):
                                try:
                                    client.pause_training_task(selected_task_id)
                                    st.success("任务已暂停")
                                    time.sleep(1)
                                    st.rerun()
                                except Exception as e:
                                    st.error(f"暂停失败: {e}")

                    with col4:
                        if st.button(
                            "查看详情", use_container_width=True, key="view_details_quick"
                        ):
                            st.session_state["active_tab"] = "tab3"
                            st.rerun()
        else:
            st.info("暂无训练任务")
    except Exception as e:
        st.error(f"获取任务列表失败: {e}")

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
            st.dataframe(config_df, use_container_width=True, hide_index=True)

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

            # Error info
            if task.get("error"):
                st.error(f"**错误信息**: {task['error']}")

            # Task control
            st.subheader("任务控制")
            task_status = task["status"]

            col1, col2, col3, col4 = st.columns(4)

            with col1:
                if task_status in ["pending", "paused"]:
                    if st.button("启动任务", use_container_width=True, key="start_detail"):
                        try:
                            client.start_training_task(task_id_input)
                            st.success("任务已启动")
                            time.sleep(1)
                            st.rerun()
                        except Exception as e:
                            st.error(f"启动失败: {e}")

            with col2:
                if task_status == "running":
                    if st.button("停止任务", use_container_width=True, key="stop_detail"):
                        try:
                            client.stop_training_task(task_id_input)
                            st.success("任务已停止")
                            time.sleep(1)
                            st.rerun()
                        except Exception as e:
                            st.error(f"停止失败: {e}")

            with col3:
                if task_status == "running":
                    if st.button("暂停任务", use_container_width=True, key="pause_detail"):
                        try:
                            client.pause_training_task(task_id_input)
                            st.success("任务已暂停")
                            time.sleep(1)
                            st.rerun()
                        except Exception as e:
                            st.error(f"暂停失败: {e}")

            with col4:
                if st.button("删除任务", use_container_width=True, key="delete_detail"):
                    try:
                        client.delete_training_task(task_id_input)
                        st.success("任务已删除")
                        time.sleep(1)
                        st.rerun()
                    except Exception as e:
                        st.error(f"删除失败: {e}")

            # Logs
            st.subheader("任务日志")

            col1, col2 = st.columns([3, 1])
            with col1:
                auto_refresh = st.checkbox("自动刷新", value=False, key="auto_refresh_logs")
            with col2:
                log_limit = st.number_input(
                    "日志行数", min_value=10, max_value=1000, value=100, step=10, key="log_limit"
                )

            try:
                logs = client.get_task_logs(task_id_input, limit=log_limit)

                if logs:
                    log_text = "\n".join(logs)
                    st.text_area("日志内容", log_text, height=400, key="log_display", disabled=True)
                else:
                    st.info("暂无日志")

                if auto_refresh:
                    time.sleep(2)
                    st.rerun()
            except Exception as e:
                st.error(f"获取日志失败: {e}")

            # Metrics
            st.subheader("训练指标")
            try:
                metrics = client.get_task_progress(task_id_input)
                if metrics:
                    st.json(metrics)
                else:
                    st.info("暂无指标数据")
            except Exception as e:
                st.warning(f"获取指标失败: {e}")

        except Exception as e:
            st.error(f"获取任务详情失败: {e}")
            st.info("请检查任务ID是否正确，或任务是否已被删除")
    else:
        st.info("请输入任务ID以查看详情，或从任务列表中选择一个任务")
