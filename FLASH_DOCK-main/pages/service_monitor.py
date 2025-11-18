"""
Service monitor page for FLASH-DOCK.
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

from service_manager import ServiceManager
from compass_client import CompassClient
from registry_url_helper import get_registry_url

st.title("服务监控")
st.write("监控COMPASS服务状态")

# 自动检测注册中心 URL（在 WSL 中时使用 Windows 主机 IP）
registry_url = get_registry_url()

# Initialize service manager
try:
    service_manager = ServiceManager(registry_url=registry_url)
    client = CompassClient(registry_url=registry_url)
    st.success(f"已连接到服务注册中心 ({registry_url})")
except Exception as e:
    st.error(f"无法连接到服务注册中心 ({registry_url}): {e}")
    st.info("提示: 如果在 WSL 中运行，请确保可以访问 Windows 主机的注册中心")
    st.stop()

# Refresh services
if st.button("刷新服务列表"):
    service_manager.refresh_services()
    st.rerun()

# Service status
st.subheader("COMPASS服务状态")

try:
    services = service_manager.registry_client.discover_compass_services(healthy_only=False)

    if services:
        df = pd.DataFrame(
            [
                {
                    "服务ID": s.service_id,
                    "地址": f"{s.host}:{s.port}",
                    "状态": s.status.value,
                    "版本": s.version,
                    "最后心跳": s.last_heartbeat.isoformat() if s.last_heartbeat else "N/A",
                }
                for s in services
            ]
        )
        st.dataframe(df, width="stretch")

        # Health status summary
        healthy_count = sum(1 for s in services if s.status.value == "healthy")
        total_count = len(services)

        col1, col2 = st.columns(2)
        with col1:
            st.metric("健康服务数", healthy_count)
        with col2:
            st.metric("总服务数", total_count)
    else:
        st.warning("未发现COMPASS服务")
except Exception as e:
    st.error(f"获取服务状态失败: {e}")

# Inference status
st.subheader("推理服务状态")

try:
    inference_status = client.get_inference_status()
    st.json(inference_status)
except Exception as e:
    st.error(f"获取推理服务状态失败: {e}")

# Model list
st.subheader("可用模型")

try:
    models = client.list_models()

    if models:
        df = pd.DataFrame(
            [
                {
                    "模型ID": m["model_id"],
                    "名称": m["name"],
                    "版本": m["version"],
                    "大小 (MB)": f"{m['file_size'] / (1024*1024):.2f}",
                    "创建时间": m["created_at"],
                }
                for m in models
            ]
        )
        st.dataframe(df, width="stretch")
    else:
        st.info("暂无可用模型")
except Exception as e:
    st.error(f"获取模型列表失败: {e}")
