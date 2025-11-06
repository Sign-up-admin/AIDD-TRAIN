"""
Data management page for FLASH-DOCK.
"""

import streamlit as st
import pandas as pd
import sys
from pathlib import Path

# Add parent directory to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
# Add FLASH_DOCK-main/services to path
flashdock_services = Path(__file__).parent.parent / "services"
sys.path.insert(0, str(flashdock_services))

from compass_client import CompassClient

st.title("数据管理")
st.write("管理COMPASS训练数据集")

# Initialize client
try:
    client = CompassClient()
    st.success("已连接到COMPASS服务")
except Exception as e:
    st.error(f"无法连接到COMPASS服务: {e}")
    st.stop()

# Tabs
tab1, tab2 = st.tabs(["上传数据集", "数据集列表"])

with tab1:
    st.subheader("上传数据集")

    uploaded_file = st.file_uploader("选择数据集文件", type=["zip", "tar", "tar.gz"])

    if uploaded_file:
        dataset_name = st.text_input("数据集名称", value=uploaded_file.name)
        dataset_description = st.text_area("数据集描述（可选）")

        if st.button("上传"):
            with st.spinner("正在上传数据集..."):
                try:
                    # Save to temp file
                    import tempfile
                    import os

                    with tempfile.NamedTemporaryFile(
                        delete=False, suffix=os.path.splitext(uploaded_file.name)[1]
                    ) as tmp_file:
                        tmp_file.write(uploaded_file.getbuffer())
                        tmp_path = tmp_file.name

                    dataset_id = client.upload_dataset(
                        tmp_path, name=dataset_name, description=dataset_description
                    )
                    os.remove(tmp_path)

                    st.success(f"数据集上传成功！数据集ID: {dataset_id}")
                except Exception as e:
                    st.error(f"上传失败: {e}")

with tab2:
    st.subheader("数据集列表")

    if st.button("刷新列表"):
        st.rerun()

    try:
        datasets = client.list_datasets()

        if datasets:
            df = pd.DataFrame(
                [
                    {
                        "数据集ID": ds["dataset_id"],
                        "名称": ds["name"],
                        "大小 (MB)": f"{ds['size'] / (1024*1024):.2f}",
                        "文件数": ds["file_count"],
                        "状态": ds["status"],
                        "创建时间": ds["created_at"],
                    }
                    for ds in datasets
                ]
            )
            st.dataframe(df, width="stretch")

            # Delete dataset
            selected_dataset_id = st.selectbox(
                "选择要删除的数据集", [ds["dataset_id"] for ds in datasets]
            )

            if st.button("删除数据集"):
                try:
                    client.delete_dataset(selected_dataset_id)
                    st.success(f"数据集 {selected_dataset_id} 已删除")
                    st.rerun()
                except Exception as e:
                    st.error(f"删除失败: {e}")
        else:
            st.info("暂无数据集")
    except Exception as e:
        st.error(f"获取数据集列表失败: {e}")
