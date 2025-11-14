# Windows兼容性：修复sh模块问题
import sys
import os
import shutil
import subprocess

if sys.platform == "win32":
    # Windows上sh模块不可用，创建兼容的sh模块
    class ShWrapper:
        @staticmethod
        def rm(*args, **kwargs):
            """删除文件或目录"""
            if "-r" in args:
                path = args[args.index("-r") + 1] if args.index("-r") + 1 < len(args) else None
                if path and os.path.exists(path):
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
            else:
                for arg in args:
                    if os.path.exists(arg):
                        if os.path.isdir(arg):
                            shutil.rmtree(arg)
                        else:
                            os.remove(arg)

        @staticmethod
        def cp(src, dst):
            """复制文件或目录"""
            if os.path.isdir(src):
                if os.path.exists(dst):
                    shutil.rmtree(dst)
                shutil.copytree(src, dst)
            else:
                shutil.copy2(src, dst)

        class Command:
            """命令包装器"""

            def __init__(self, cmd):
                self.cmd = cmd

            def __call__(self, *args, **kwargs):
                cmd_list = [self.cmd] + list(args)
                cwd = kwargs.get("_cwd", None)
                fg = kwargs.get("_fg", False)
                if fg:
                    return subprocess.run(cmd_list, cwd=cwd, check=False)
                else:
                    return subprocess.Popen(cmd_list, cwd=cwd)

    # 在导入streamlit_molstar之前，先注册sh模块
    sys.modules["sh"] = ShWrapper()

import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from streamlit_molstar import st_molstar, st_molstar_rcsb, st_molstar_remote

# 如需使用口袋预测相关函数
from streamlit_molstar.pocket import (
    select_pocket_from_local_protein,
    # 如果你的项目需要也可以 import select_pocket_from_upload_protein
)

# docking 模块
from streamlit_molstar.docking import st_molstar_docking

import json
import tempfile  # 用于创建临时文件
import re
import tqdm

# 如果没有在 session_state 中记录 page，就初始化一个默认值
if "page" not in st.session_state:
    st.session_state["page"] = "主页"

# 在侧边栏使用按钮来切换页面
st.sidebar.title("Navigation")
if st.sidebar.button("主页"):
    st.session_state["page"] = "主页"

# 新增“准备配体”按钮（插在“主页”和“口袋预测”之间）
if st.sidebar.button("准备配体"):
    st.session_state["page"] = "准备配体"

if st.sidebar.button("口袋预测"):
    st.session_state["page"] = "口袋预测"
if st.sidebar.button("分子对接"):
    st.session_state["page"] = "分子对接"
if st.sidebar.button("批量口袋预测与对接"):
    st.session_state["page"] = "批量口袋预测与对接"
    # 新增“预测亲和力”按钮
if st.sidebar.button("预测亲和力"):
    st.session_state["page"] = "预测亲和力"

# 添加数据管理、服务监控、训练管理按钮
st.sidebar.markdown("---")
st.sidebar.subheader("服务管理")
if st.sidebar.button("数据管理"):
    st.session_state["page"] = "数据管理"
if st.sidebar.button("服务监控"):
    st.session_state["page"] = "服务监控"
if st.sidebar.button("训练管理"):
    st.session_state["page"] = "训练管理"

# 获取当前页面
page = st.session_state["page"]

# ------------------------------------------------------------------------------
# 主页
# ------------------------------------------------------------------------------
if page == "主页":
    # 使用 HTML 和 Markdown 居中标题
    st.markdown(
        "<h1 style='text-align: center;'>⚡️欢迎使用⚡️</h1>",
        unsafe_allow_html=True,
    )
    st.markdown("<br>", unsafe_allow_html=True)

    # 显示字符画
    try:
        with open("./others/logo.txt", "r", encoding="utf-8") as file:
            ascii_art = file.read()
            styled_ascii_art = ascii_art.replace(" ", "&nbsp;").replace("\n", "<br>")
            html_code = f"""
            <div style='text-align: center; font-family: monospace; font-size: 14px; line-height: 1;'>
                {styled_ascii_art}
            </div>
            """
            st.markdown(html_code, unsafe_allow_html=True)

    except FileNotFoundError:
        st.error("logo.txt 文件未找到，请确保它与脚本位于同一目录下。")
    except UnicodeDecodeError:
        st.error("无法解码 logo.txt 文件，请确认文件编码格式是否为 UTF-8。")

    # 在字符画和图片之间插入若干空行
    st.markdown("<br><br><br><br><br>", unsafe_allow_html=True)

    # 显示 logo.png
    if os.path.exists("./others/logo.png"):
        st.image("./others/logo.png", use_container_width=True)
    else:
        st.error("logo.png 文件未找到，请确保它与脚本位于同一目录下。")

# ------------------------------------------------------------------------------
# 准备配体
# ------------------------------------------------------------------------------
elif page == "准备配体":
    st.title("准备配体")

    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    from rdkit.Chem.Draw import rdMolDraw2D
    import py3Dmol
    from stmol import showmol

    # 尝试导入 streamlit_ketcher
    ketcher_available = False
    st_ketcher = None
    try:
        from streamlit_ketcher import st_ketcher

        ketcher_available = True
    except ImportError:
        # 导入失败，库未安装
        ketcher_available = False
        st_ketcher = None
    except Exception:
        # 其他导入错误
        ketcher_available = False
        st_ketcher = None

    # 1. 允许用户上传一个 SDF 文件
    st.info("请上传一个 SDF 文件，或在画布中绘制分子结构/粘贴SMILES")
    st.markdown("**上传**分子文件（SDF 格式）：")
    sdf_file = st.file_uploader("上传 SDF 文件", type=["sdf"], label_visibility="hidden")

    # 2. 允许用户使用 Ketcher 绘制或输入 SMILES
    st.markdown("**或者** 在下方绘制分子结构/粘贴SMILES：")
    smiles_input = None

    if ketcher_available and st_ketcher is not None:
        # 注入 JavaScript 补丁来修复 eventBus 初始化问题
        # 这个补丁会在组件加载前确保 eventBus 正确初始化
        ketcher_fix_script = """
        <script>
        (function() {
            // 立即创建 eventBus，不等待任何事件
            function createEventBus() {
                return {
                    listeners: {},
                    on: function(event, callback) {
                        if (!this.listeners[event]) {
                            this.listeners[event] = [];
                        }
                        this.listeners[event].push(callback);
                    },
                    emit: function(event, data) {
                        if (this.listeners[event]) {
                            this.listeners[event].forEach(callback => {
                                try {
                                    callback(data);
                                } catch (e) {
                                    console.error('EventBus callback error:', e);
                                }
                            });
                        }
                    },
                    off: function(event, callback) {
                        if (this.listeners[event]) {
                            this.listeners[event] = this.listeners[event].filter(cb => cb !== callback);
                        }
                    }
                };
            }
            
            // 创建全局 eventBus 实例
            var globalEventBus = createEventBus();
            
            // 在多个位置设置 eventBus，确保组件能找到它
            if (typeof window !== 'undefined') {
                // 设置全局 eventBus（多个可能的名称）
                window.__ketcher_eventBus = globalEventBus;
                window.eventBus = globalEventBus;
                window.ketcherEventBus = globalEventBus;
                
                // 确保在 document 上也有引用
                if (typeof document !== 'undefined') {
                    document.__ketcher_eventBus = globalEventBus;
                    document.eventBus = globalEventBus;
                }
                
                // 使用 Object.defineProperty 确保 eventBus 始终可用
                try {
                    Object.defineProperty(window, 'eventBus', {
                        get: function() {
                            return globalEventBus;
                        },
                        set: function(value) {
                            globalEventBus = value || createEventBus();
                        },
                        configurable: true
                    });
                } catch (e) {
                    console.warn('Could not define eventBus property:', e);
                }
                
                // 使用 Proxy 拦截对 undefined 对象的 eventBus 访问
                try {
                    var originalGet = Object.prototype.__lookupGetter__;
                    var eventBusProxy = new Proxy({}, {
                        get: function(target, prop) {
                            if (prop === 'eventBus') {
                                return globalEventBus;
                            }
                            return undefined;
                        }
                    });
                    
                    // 将 eventBus 设置为全局可访问
                    window.__streamlit_ketcher_eventBus = globalEventBus;
                } catch (e) {
                    console.warn('Could not create eventBus proxy:', e);
                }
                
                // 监听 iframe 加载，确保 iframe 内部也能访问 eventBus
                function setupIframeEventBus(iframe) {
                    try {
                        if (iframe.contentWindow) {
                            var iframeWindow = iframe.contentWindow;
                            iframeWindow.__ketcher_eventBus = globalEventBus;
                            iframeWindow.eventBus = globalEventBus;
                            iframeWindow.ketcherEventBus = globalEventBus;
                            
                            // 也在 iframe 的 document 上设置
                            try {
                                if (iframeWindow.document) {
                                    iframeWindow.document.__ketcher_eventBus = globalEventBus;
                                    iframeWindow.document.eventBus = globalEventBus;
                                }
                            } catch (e) {
                                // 可能无法访问 iframe document
                            }
                        }
                    } catch (e) {
                        // 跨域限制，无法访问 iframe 内容
                    }
                }
                
                // 立即检查现有的 iframe
                if (typeof document !== 'undefined') {
                    var existingIframes = document.querySelectorAll('iframe');
                    existingIframes.forEach(setupIframeEventBus);
                    
                    // 使用 MutationObserver 监听新添加的 iframe
                    var observer = new MutationObserver(function(mutations) {
                        mutations.forEach(function(mutation) {
                            mutation.addedNodes.forEach(function(node) {
                                if (node.nodeType === 1) { // Element node
                                    if (node.tagName === 'IFRAME') {
                                        setupIframeEventBus(node);
                                        // 等待 iframe 加载完成后再设置
                                        node.addEventListener('load', function() {
                                            setupIframeEventBus(node);
                                        });
                                    }
                                    var iframes = node.querySelectorAll ? node.querySelectorAll('iframe') : [];
                                    iframes.forEach(function(iframe) {
                                        setupIframeEventBus(iframe);
                                        iframe.addEventListener('load', function() {
                                            setupIframeEventBus(iframe);
                                        });
                                    });
                                }
                            });
                        });
                    });
                    
                    if (document.body) {
                        observer.observe(document.body, {
                            childList: true,
                            subtree: true
                        });
                    } else {
                        // 如果 body 还没加载，等待 DOMContentLoaded
                        document.addEventListener('DOMContentLoaded', function() {
                            observer.observe(document.body, {
                                childList: true,
                                subtree: true
                            });
                            // 再次检查所有 iframe
                            var allIframes = document.querySelectorAll('iframe');
                            allIframes.forEach(setupIframeEventBus);
                        });
                    }
                }
            }
        })();
        </script>
        """

        # 注入修复脚本（必须在组件加载前执行）
        components.html(ketcher_fix_script, height=0)

        try:
            # 尝试使用 Ketcher 组件
            smiles_input = st_ketcher()
        except Exception as e:
            # 如果仍然失败，显示错误信息
            st.error(
                f"⚠️ **分子编辑器组件初始化失败**: {str(e)}\n\n"
                "请尝试刷新页面，或使用下方的文本输入框输入 SMILES 字符串。"
            )
            smiles_input = None

    # 如果 Ketcher 不可用或失败，提供文本输入作为备选方案
    if not smiles_input:
        smiles_input = st.text_input(
            "输入 SMILES 字符串、分子式或结构式（如果上方编辑器不可用）",
            placeholder="例如: CCO (乙醇), CC(=O)O (乙酸), C9H8O4 (分子式), CH3COOC6H4COOH (结构式) 等",
            help="支持格式：\n1. SMILES字符串（如：CCO, CC(=O)O）\n2. 分子式（如：C9H8O4）\n3. 结构式（如：CH3COOC6H4COOH）",
        )

    def convert_formula_to_smiles(formula: str) -> str:
        """
        尝试将分子式或结构式转换为SMILES字符串
        
        支持的格式：
        1. 标准分子式：C9H8O4
        2. 结构式：CH3COOC6H4COOH
        
        Args:
            formula: 分子式或结构式字符串
            
        Returns:
            SMILES字符串，如果转换失败则返回None
        """
        # 首先尝试直接作为SMILES解析
        mol = Chem.MolFromSmiles(formula)
        if mol is not None:
            return formula
        
        # 清理输入，移除空格
        formula_clean = formula.strip().replace(" ", "")
        
        # 尝试一些常见的结构式转换规则（内置字典，无需网络）
        # 例如：CH3COOC6H4COOH 可能是阿司匹林的结构式
        common_structures = {
            "CH3COOC6H4COOH": "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林
            "C9H8O4": "CC(=O)OC1=CC=CC=C1C(=O)O",  # 阿司匹林（阿司匹林的分子式）
        }
        
        if formula_clean.upper() in common_structures:
            return common_structures[formula_clean.upper()]
        
        # 尝试通过PubChem API查找分子式对应的SMILES（需要网络连接和requests库）
        try:
            import requests
        except ImportError:
            # requests库未安装，跳过API调用
            return None
        
        try:
            # 检查是否是标准分子式格式（如C9H8O4）
            # 标准分子式通常以元素符号开头，后跟数字
            if re.match(r'^[A-Z][a-z]?\d*([A-Z][a-z]?\d*)*$', formula_clean):
                # 使用PubChem API通过分子式查找
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/formula/{formula_clean}/property/CanonicalSMILES/JSON"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                        properties = data['PropertyTable']['Properties']
                        if properties and len(properties) > 0:
                            smiles = properties[0].get('CanonicalSMILES')
                            if smiles:
                                # 验证SMILES是否有效
                                test_mol = Chem.MolFromSmiles(smiles)
                                if test_mol is not None:
                                    return smiles
            else:
                # 可能是结构式（如CH3COOC6H4COOH），尝试通过名称查找
                # 这里我们尝试通过PubChem的名称搜索
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{formula_clean}/property/CanonicalSMILES/JSON"
                response = requests.get(url, timeout=10)
                if response.status_code == 200:
                    data = response.json()
                    if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                        properties = data['PropertyTable']['Properties']
                        if properties and len(properties) > 0:
                            smiles = properties[0].get('CanonicalSMILES')
                            if smiles:
                                test_mol = Chem.MolFromSmiles(smiles)
                                if test_mol is not None:
                                    return smiles
        except Exception:
            # API调用失败（网络问题、超时等），继续返回None
            pass
        
        return None

    def process_and_show_mol(
        mol: Chem.Mol, uploaded_sdf_name: str = None, user_defined_filename: str = None
    ):
        """
        对分子进行加氢、3D 嵌入、MMFF 优化并展示 2D/3D 结构；
        根据不同来源决定最终保存的文件名：
        - 如果有 uploaded_sdf_name，则用 "原文件名去除.sdf + '_prepared.sdf'"
        - 如果没有 uploaded_sdf_name，但用户给了自定义文件名，则用 "用户自定义文件名 + '.sdf'"
        """
        if not mol:
            return

        # 2D 可视化
        st.subheader("2D 分子结构")
        try:
            # 方法1: 尝试使用 PIL 图像显示（更可靠）
            try:
                img = Draw.MolToImage(mol, size=(400, 400))
                if img:
                    st.image(img, caption="2D 分子结构", use_container_width=False)
                else:
                    raise ValueError("图像生成失败")
            except Exception as img_error:
                # 方法2: 如果图像失败，尝试使用 SVG
                st.warning(f"图像显示失败，尝试使用 SVG 方式: {str(img_error)}")
                try:
                    svg = rdMolDraw2D.MolToSVG(mol, width=400, height=400)
                    if svg and len(svg) > 0:
                        # 确保 SVG 格式正确，添加居中样式
                        if svg.strip().startswith('<svg'):
                            # 完整的 SVG，添加容器
                            svg_html = f'<div style="text-align: center; margin: 10px 0;">{svg}</div>'
                        else:
                            # 不完整的 SVG，尝试包装
                            svg_html = f'<div style="text-align: center; margin: 10px 0;">{svg}</div>'
                        st.markdown(svg_html, unsafe_allow_html=True)
                    else:
                        raise ValueError("SVG 生成失败：返回空内容")
                except Exception as svg_error:
                    # 方法3: 如果都失败，显示错误信息
                    st.error(f"2D 分子结构显示失败: {str(svg_error)}")
                    st.info("请检查分子结构是否有效，或尝试使用其他可视化工具。")
                    # 尝试显示分子的基本信息
                    try:
                        st.write(f"分子原子数: {mol.GetNumAtoms()}")
                        st.write(f"分子键数: {mol.GetNumBonds()}")
                    except:
                        pass
        except Exception as e:
            st.error(f"2D 可视化过程中出现错误: {str(e)}")
            st.info("请检查分子结构是否有效。")

        # 生成 3D 构象并能量优化
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol_3d)

        # 3D 可视化
        st.subheader("3D 分子结构")
        mol_block = Chem.MolToMolBlock(mol_3d)
        xyzview = py3Dmol.view(width=500, height=400)
        xyzview.addModel(mol_block, "mol")
        xyzview.setStyle({"stick": {}})
        xyzview.zoomTo()
        showmol(xyzview, height=400, width=500)

        # 提供保存按钮，将 3D 结构写出为 SDF 文件
        if st.button("保存 3D 结构为 SDF"):
            os.makedirs("./Result/Prepare_Ligand", exist_ok=True)

            if uploaded_sdf_name:
                # 如果用户上传了 SDF，就使用该 SDF 名（去 .sdf）并加上 _prepared
                base_name = os.path.splitext(uploaded_sdf_name)[0]
                out_filename = base_name + "_prepared.sdf"
            else:
                # 如果没有上传的 SDF，就使用用户输入的文件名（不含 .sdf 后缀），再加上 .sdf
                if user_defined_filename:
                    out_filename = user_defined_filename.strip() + ".sdf"
                else:
                    # 如果用户也没有输入任何自定义文件名，可给一个默认值
                    out_filename = "ligand_3d.sdf"

            sdf_path = os.path.join("./Result/Prepare_Ligand", out_filename)
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol_3d)
            writer.close()
            st.success(f"已将 3D 结构保存到 {sdf_path}")

    # 首先解析用户上传的 SDF
    mol_from_sdf = None
    uploaded_sdf_name = None

    def fix_mol_dimension(mol):
        """修复分子的维度标记，避免 RDKit 警告"""
        if mol is None:
            return mol
        # 检查分子是否有 3D 坐标（非零 Z 坐标）
        conf = mol.GetConformer()
        if conf is not None:
            has_3d = False
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                if abs(pos.z) > 0.001:  # 如果 Z 坐标不为零
                    has_3d = True
                    break
            # 如果有 3D 坐标，确保分子被标记为 3D
            if has_3d:
                conf.Set3D(True)
        return mol

    if sdf_file is not None:
        uploaded_sdf_name = sdf_file.name  # 记录用户上传的文件名
        try:
            # 使用 removeHs=False 和 sanitize=True 来正确读取 SDF 文件
            sdf_supplier = Chem.ForwardSDMolSupplier(sdf_file, removeHs=False, sanitize=True)
            mols = []
            for mol in sdf_supplier:
                if mol is not None:
                    # 立即修复维度标记，避免警告
                    mol = fix_mol_dimension(mol)
                    mols.append(mol)
            if len(mols) > 0:
                mol_from_sdf = mols[0]
            else:
                st.error("无法从 SDF 文件中解析出分子，请检查文件格式或内容。")
        except Exception as e:
            st.error(f"读取 SDF 文件出现错误: {e}")

    if mol_from_sdf:
        # 如果成功解析出上传的 SDF，则展示并保存
        process_and_show_mol(mol_from_sdf, uploaded_sdf_name=uploaded_sdf_name)
    else:
        # 如果用户没有上传 SDF 或上传的 SDF 解析失败，则查看 Ketcher 中有没有输入 SMILES
        if smiles_input:
            # 首先尝试直接作为SMILES解析
            mol_from_smiles = Chem.MolFromSmiles(smiles_input)
            
            # 如果直接解析失败，尝试作为分子式或结构式转换
            if mol_from_smiles is None:
                with st.spinner("正在尝试将分子式/结构式转换为SMILES..."):
                    converted_smiles = convert_formula_to_smiles(smiles_input)
                    if converted_smiles:
                        st.info(f"已成功转换：{smiles_input} → {converted_smiles}")
                        mol_from_smiles = Chem.MolFromSmiles(converted_smiles)
                    else:
                        st.warning(f"无法将 '{smiles_input}' 转换为有效的SMILES。尝试作为SMILES直接解析...")
            
            if mol_from_smiles:
                user_defined_filename = st.text_input(
                    "请输入保存时的 SDF 文件名（不含 .sdf）", value="my_mol"
                )
                process_and_show_mol(
                    mol_from_smiles,
                    uploaded_sdf_name=None,
                    user_defined_filename=user_defined_filename,
                )
            else:
                # 检查是否安装了requests库
                try:
                    import requests
                    requests_available = True
                except ImportError:
                    requests_available = False
                
                error_msg = (
                    f"无法解析输入 '{smiles_input}'。\n\n"
                    "请确保输入格式正确：\n"
                    "1. SMILES字符串（如：CCO, CC(=O)O）\n"
                    "2. 标准分子式（如：C9H8O4）\n"
                    "3. 结构式（如：CH3COOC6H4COOH）\n\n"
                )
                
                if not requests_available:
                    error_msg += (
                        "⚠️ 注意：检测到未安装 `requests` 库，无法使用PubChem API查询分子式。\n"
                        "如需支持分子式转换，请安装：`pip install requests`\n\n"
                        "目前仅支持内置的常见分子式（如：C9H8O4, CH3COOC6H4COOH）。"
                    )
                else:
                    error_msg += (
                        "注意：分子式和结构式需要通过PubChem API查询，需要网络连接。\n"
                        "如果网络不可用，请使用SMILES字符串或上传SDF文件。"
                    )
                
                st.error(error_msg)

# ------------------------------------------------------------------------------
# 口袋预测
# ------------------------------------------------------------------------------
elif page == "口袋预测":
    st.title("口袋预测")

    # 让用户选择如何加载蛋白质
    option = st.radio("Select how to load the protein:", ("上传蛋白质", "加载示例文件"))

    # 用于保存用户上传的蛋白文件名称（用于替换 Pocket Name）
    uploaded_pdb_filename = None

    if option == "上传蛋白质":
        try:
            # 用户上传蛋白质（只出现一次，不会再弹二次上传）
            pdb_file = st.file_uploader("请上传蛋白质文件 (.pdb)", type=["pdb"])

            if pdb_file is not None:
                # 记下上传的名称
                uploaded_pdb_filename = pdb_file.name

                # 使用临时文件的方式进行口袋预测
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                    tmp.write(pdb_file.getvalue())
                    tmp.flush()
                    file_path = tmp.name

                # 调用 p2rank (或其他函数) ，读取该临时文件进行预测
                selected = select_pocket_from_local_protein(
                    file_path, p2rank_home="./others/p2rank_2.5/"
                )
                # 预测完成后删除该临时文件
                os.remove(file_path)

                if selected:
                    pocket = selected
                    st.write("预测到的口袋信息: ", pocket)

                    # 如果 rank=1 的口袋
                    if pocket["rank"] == "1":
                        # 如果上传了文件名，则用之，否则用 pocket['name']
                        final_name = (
                            uploaded_pdb_filename if uploaded_pdb_filename else pocket["name"]
                        )
                        data = {
                            "Pocket Name": [final_name],
                            "Center": [pocket["center"]],
                        }
                        df = pd.DataFrame(data)

                        st.write("最优口袋信息预览：")
                        st.dataframe(df)

                        # 用户点击按钮后，才将CSV保存到指定文件夹
                        if st.button("保存 best_pocket.csv"):
                            os.makedirs("./Result/Predict_Pocket", exist_ok=True)
                            csv_path = "./Result/Predict_Pocket/best_pocket.csv"
                            df.to_csv(csv_path, index=False)
                            st.success(f"best_pocket.csv 已保存到 {csv_path}")

        except Exception as e:
            st.warning(f"处理上传蛋白时发生错误: {e}")

    elif option == "加载示例文件":
        try:
            # 用示例文件名
            uploaded_pdb_filename = "protein_example.pdb"
            # 调用 p2rank 做预测
            selected = select_pocket_from_local_protein(
                "examples/pocket/protein.pdb", p2rank_home="./others/p2rank_2.5/"
            )
            if selected:
                pocket = selected
                st.write("预测到的口袋信息: ", pocket)

                if pocket["rank"] == "1":
                    data = {
                        "Pocket Name": [uploaded_pdb_filename],
                        "Center": [pocket["center"]],
                    }
                    df = pd.DataFrame(data)

                    st.write("最优口袋信息预览：")
                    st.dataframe(df)

                    if st.button("保存 best_pocket.csv"):
                        os.makedirs("./Result/Predict_Pocket", exist_ok=True)
                        csv_path = "./Result/Predict_Pocket/best_pocket.csv"
                        df.to_csv(csv_path, index=False)
                        st.success(f"best_pocket.csv 已保存到 {csv_path}")

        except Exception as e:
            st.warning(f"加载示例文件时发生错误: {e}")

# ------------------------------------------------------------------------------
# 分子对接
# ------------------------------------------------------------------------------
elif page == "分子对接":
    st.title("分子对接")
    st.write("请上传蛋白质 (PDB 格式) 和配体 (SDF 格式)，并设置对接参数。")

    # 让用户上传蛋白质和配体文件
    protein_file = st.file_uploader("上传蛋白质文件 (.pdb)", type=["pdb"])
    ligand_file = st.file_uploader("上传配体文件 (.sdf)", type=["sdf"])

    # 让用户上传口袋预测结果文件（可选）
    st.write("可选：上传口袋预测结果 CSV 文件，将自动填充对接网格参数。")
    pocket_csv_file = st.file_uploader("上传口袋预测结果文件 (CSV)", type=["csv"])

    # 默认网格参数
    center_x = 0.0
    center_y = 0.0
    center_z = 0.0

    if pocket_csv_file is not None:
        try:
            # 读取 CSV 文件并获取中心坐标
            pocket_df = pd.read_csv(pocket_csv_file)
            if "Center" in pocket_df.columns:
                center_coords = pocket_df.loc[0, "Center"]  # 获取第一个口袋的中心坐标
                if isinstance(center_coords, str):
                    coords = [float(c) for c in re.findall(r"[-+]?[0-9]*\.?[0-9]+", center_coords)]
                    if len(coords) == 3:
                        center_x, center_y, center_z = coords
                    else:
                        st.warning("CSV 文件中的 Center 格式不正确，无法自动填充网格参数。")
                else:
                    st.warning("CSV 文件中的 Center 格式不正确，无法自动填充网格参数。")
            else:
                st.warning("CSV 文件中未找到 Center 列，无法自动填充网格参数。")
        except Exception as e:
            st.error(f"读取 CSV 文件时出现错误: {e}")

    # 显示网格参数输入框，无论是否上传 CSV 文件
    st.subheader("设置对接口袋参数")
    center_x = st.number_input("Center X", value=center_x)
    center_y = st.number_input("Center Y", value=center_y)
    center_z = st.number_input("Center Z", value=center_z)

    size_x = st.number_input("Size X", value=100.0)
    size_y = st.number_input("Size Y", value=100.0)
    size_z = st.number_input("Size Z", value=100.0)

    # 当用户点击“开始分子对接”时，生成 docking_grid.json 文件并调用对接命令
    if st.button("开始分子对接"):
        # 如果没有上传蛋白质或配体，提示错误
        if not protein_file or not ligand_file:
            st.error("请先上传蛋白质 (pdb) 和配体 (sdf) 文件。")
        else:
            try:
                # 创建临时文件夹保存缓存文件
                with tempfile.TemporaryDirectory() as temp_dir:
                    docking_grid = {
                        "center_x": center_x,
                        "center_y": center_y,
                        "center_z": center_z,
                        "size_x": size_x,
                        "size_y": size_y,
                        "size_z": size_z,
                    }
                    docking_grid_path = os.path.join(temp_dir, "docking_grid.json")

                    with open(docking_grid_path, "w") as f:
                        json.dump(docking_grid, f, indent=4)

                    # 保存蛋白质和配体文件到临时目录
                    protein_path = os.path.join(temp_dir, "protein.pdb")
                    ligand_path = os.path.join(temp_dir, "ligand.sdf")

                    with open(protein_path, "wb") as f:
                        f.write(protein_file.getvalue())

                    with open(ligand_path, "wb") as f:
                        f.write(ligand_file.getvalue())

                    # 设置结果保存目录
                    result_dir = "./Result/Docking_Result"
                    os.makedirs(result_dir, exist_ok=True)

                    # 构造命令（使用列表形式避免shell注入风险）
                    command = [
                        "python",
                        "./others/Uni-Mol/unimol_docking_v2/interface/demo.py",
                        "--mode",
                        "single",
                        "--conf-size",
                        "10",
                        "--cluster",
                        "--input-protein",
                        str(protein_path),
                        "--input-ligand",
                        str(ligand_path),
                        "--input-docking-grid",
                        str(docking_grid_path),
                        "--output-ligand-name",
                        "ligand_predict",
                        "--output-ligand-dir",
                        str(result_dir),
                        "--steric-clash-fix",
                        "--model-dir",
                        "./others/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt",
                    ]

                    # 执行命令（使用列表形式，避免shell注入风险）
                    result = subprocess.run(command, capture_output=True, text=True, timeout=300)

                    # 根据命令返回值判断是否执行成功
                    if result.returncode == 0:
                        st.success("分子对接完成！")
                        st.text_area("对接输出日志", value=result.stdout, height=150)

                        # 分子对接完成后，处理结果文件
                        try:
                            ligand_output_path = os.path.join(result_dir, "ligand_predict.sdf")

                            # 删除结果目录中除 ligand_predict.sdf 外的所有文件
                            for file_name in os.listdir(result_dir):
                                file_path = os.path.join(result_dir, file_name)
                                if file_name != "ligand_predict.sdf" and os.path.isfile(file_path):
                                    os.remove(file_path)

                            # 重命名 ligand_predict.sdf
                            output_name = f"{os.path.splitext(ligand_file.name)[0]}_{os.path.splitext(protein_file.name)[0]}_docked.sdf"
                            renamed_path = os.path.join(result_dir, output_name)
                            os.rename(ligand_output_path, renamed_path)

                            # 提示用户结果保存位置
                            st.success(f"对接结果保存为 {renamed_path}")

                            # 可视化对接结果
                            st_molstar_docking(protein_path, renamed_path, key="5", height=600)
                        except Exception:
                            st.error("处理结果文件时出错，请检查路径或权限。")

                    else:
                        st.error("分子对接失败！")
                        st.text_area("错误信息", value=result.stderr, height=150)

            except Exception as e:
                st.error(f"对接过程出现错误: {e}")

# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
elif page == "批量口袋预测与对接":
    import os
    import pandas as pd
    import subprocess
    import tempfile
    import json
    from pathlib import Path
    import streamlit as st
    from streamlit_molstar.pocket import select_pocket_from_local_protein

    st.title("批量口袋预测与分子对接")

    # 定义固定路径
    batch_docking_dir = Path("./Batch_Docking")
    result_dir = Path("./Batch_Docking")
    result_dir.mkdir(parents=True, exist_ok=True)

    # 检查 Batch_Docking 目录是否存在
    if not batch_docking_dir.exists():
        st.error(f"目录 {batch_docking_dir} 不存在。请创建该目录并添加 PDB 和 SDF 文件。")
    else:
        # 自动生成任务 CSV 文件
        def generate_task_csv():
            pdb_files = list(batch_docking_dir.glob("*.pdb"))
            sdf_files = list(batch_docking_dir.glob("*.sdf"))

            if not pdb_files:
                st.error("在 ./Batch_Docking 文件夹中未找到 PDB 文件。请添加至少一个 PDB 文件。")
                return None
            if not sdf_files:
                st.error("在 ./Batch_Docking 文件夹中未找到 SDF 文件。请添加至少一个 SDF 文件。")
                return None

            tasks = []
            for pdb_file in pdb_files:
                for sdf_file in sdf_files:
                    tasks.append(
                        {
                            "Protein": pdb_file.name,
                            "Ligand": sdf_file.name,
                            "Run": "Yes",  # 默认所有任务为 "Yes"
                        }
                    )

            task_df = pd.DataFrame(tasks)
            return task_df

        # 生成任务 DataFrame
        task_df = generate_task_csv()

        if task_df is not None:
            # 提供下载任务 CSV 的按钮
            csv = task_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                label="下载任务 CSV 文件", data=csv, file_name="docking_tasks.csv", mime="text/csv"
            )

            st.markdown("---")
            st.info(
                """
                1. 下载上方的任务 CSV 文件。
                2. 在本地编辑 CSV 文件，修改 `Run` 列为 `Yes` 的任务将被执行，`No` 列的任务将被跳过。
                3. 修改完成后，上传修改后的 CSV 文件并点击“开始批量预测和对接”按钮。
            """
            )

            # 上传修改后的任务 CSV 文件
            uploaded_csv = st.file_uploader(
                "上传修改后的任务 CSV 文件", type=["csv"], key="upload_task_csv"
            )

            if uploaded_csv is not None:
                try:
                    uploaded_tasks_df = pd.read_csv(uploaded_csv)

                    # 检查必要的列是否存在
                    required_columns = {"Protein", "Ligand", "Run"}
                    if not required_columns.issubset(uploaded_tasks_df.columns):
                        st.error(
                            f"上传的任务文件缺少必要的列：{required_columns - set(uploaded_tasks_df.columns)}"
                        )
                    else:
                        # 过滤需要运行的任务
                        tasks_to_run = uploaded_tasks_df[
                            uploaded_tasks_df["Run"].str.lower() == "yes"
                        ]

                        if tasks_to_run.empty:
                            st.warning(
                                "没有任务需要运行，请确保至少有一项任务的 `Run` 列为 `Yes`。"
                            )
                        else:
                            st.write(f"发现 {len(tasks_to_run)} 个任务需要运行。")

                            # 显示需要运行的任务表格
                            st.subheader("待运行的任务列表")
                            st.dataframe(tasks_to_run[["Protein", "Ligand"]].reset_index(drop=True))

                            # 开始批量预测和对接按钮
                            if st.button("开始批量预测和对接", key="start_batch_processing"):
                                log_messages = []
                                progress_bar = st.progress(0)
                                status_text = st.empty()

                                for i, task in tasks_to_run.iterrows():
                                    protein_path = batch_docking_dir / task["Protein"]
                                    ligand_path = batch_docking_dir / task["Ligand"]

                                    # 口袋预测前更新状态
                                    status_text.text(
                                        f"任务 {i + 1}/{len(tasks_to_run)}: 正在为 {task['Protein']} 预测口袋..."
                                    )

                                    # 口袋预测
                                    try:
                                        # 为每个任务传递唯一的 key
                                        pocket_result = select_pocket_from_local_protein(
                                            str(protein_path),
                                            p2rank_home="./others/p2rank_2.5/",
                                            key=f"select_pocket_{i}",  # 添加唯一 key
                                        )
                                    except Exception as e:
                                        log_messages.append(
                                            f"任务 {task['Protein']} 和 {task['Ligand']} 口袋预测失败：{e}"
                                        )
                                        progress_bar.progress((i + 1) / len(tasks_to_run))
                                        continue

                                    if pocket_result:
                                        center_coords = [
                                            float(coord) for coord in pocket_result["center"]
                                        ]
                                        docking_grid = {
                                            "center_x": center_coords[0],
                                            "center_y": center_coords[1],
                                            "center_z": center_coords[2],
                                            "size_x": 100.0,
                                            "size_y": 100.0,
                                            "size_z": 100.0,
                                        }

                                        # 创建临时目录存储对接网格
                                        with tempfile.TemporaryDirectory() as temp_dir:
                                            docking_grid_path = Path(temp_dir) / "docking_grid.json"
                                            with open(docking_grid_path, "w") as f:
                                                json.dump(docking_grid, f, indent=4)

                                            # 更新状态：开始对接
                                            status_text.text(
                                                f"任务 {i + 1}/{len(tasks_to_run)}: 正在对接 {task['Protein']} 和 {task['Ligand']}..."
                                            )

                                            # 构造对接命令（使用列表形式避免shell注入风险）
                                            command = [
                                                "python",
                                                "./others/Uni-Mol/unimol_docking_v2/interface/demo.py",
                                                "--mode",
                                                "single",
                                                "--conf-size",
                                                "10",
                                                "--cluster",
                                                "--input-protein",
                                                str(protein_path),
                                                "--input-ligand",
                                                str(ligand_path),
                                                "--input-docking-grid",
                                                str(docking_grid_path),
                                                "--output-ligand-name",
                                                "ligand_predict",
                                                "--output-ligand-dir",
                                                str(result_dir),
                                                "--steric-clash-fix",
                                                "--model-dir",
                                                "./others/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt",
                                            ]

                                            # 执行对接命令（使用列表形式，避免shell注入风险）
                                            result = subprocess.run(
                                                command,
                                                capture_output=True,
                                                text=True,
                                                timeout=300,
                                            )

                                            if result.returncode == 0:
                                                ligand_output_path = (
                                                    result_dir / "ligand_predict.sdf"
                                                )
                                                output_name = f"{protein_path.stem}_{ligand_path.stem}_docked.sdf"
                                                renamed_path = result_dir / output_name

                                                try:
                                                    os.rename(ligand_output_path, renamed_path)
                                                    log_messages.append(
                                                        f"任务 {task['Protein']} 和 {task['Ligand']} 对接完成。结果保存为 {renamed_path}"
                                                    )
                                                except Exception as e:
                                                    log_messages.append(
                                                        f"任务 {task['Protein']} 和 {task['Ligand']} 结果保存失败：{e}"
                                                    )
                                            else:
                                                log_messages.append(
                                                    f"任务 {task['Protein']} 和 {task['Ligand']} 对接失败。错误信息：{result.stderr}"
                                                )

                                    else:
                                        log_messages.append(
                                            f"任务 {task['Protein']} 的口袋信息未找到。"
                                        )

                                    # 更新进度条
                                    progress_bar.progress((i + 1) / len(tasks_to_run))

                                # 所有任务完成后更新状态
                                status_text.text("所有任务已完成。")

                                # 显示日志
                                st.success("所有任务已完成。")
                                st.text_area("任务日志", value="\n".join(log_messages), height=300)
                except pd.errors.EmptyDataError:
                    st.error("上传的 CSV 文件为空，请检查文件内容。")
                except pd.errors.ParserError:
                    st.error("上传的 CSV 文件格式错误，请确保文件为有效的 CSV 格式。")
                except Exception as e:
                    st.error(f"读取任务文件时出错：{e}")

# ------------------------------------------------------------------------------
# 预测亲和力
# ------------------------------------------------------------------------------

elif page == "预测亲和力":
    import os
    import tempfile
    import subprocess
    import pandas as pd
    import streamlit as st
    import time
    import matplotlib.pyplot as plt
    import seaborn as sns

    if page == "预测亲和力":
        st.title("预测亲和力")
        st.write("在此页面，你可以进行小分子与蛋白质的结合亲和力预测。选择单个预测或批量预测模式。")

        # 模式选择
        mode = st.radio("选择模式", ("单个预测", "批量预测"))

        if mode == "单个预测":
            st.subheader("单个蛋白与小分子的亲和力预测")

            # 用户上传蛋白质 PDB 文件
            protein_file = st.file_uploader("上传蛋白质 PDB 文件", type=["pdb"])

            # 用户上传小分子 SDF 文件
            ligand_file = st.file_uploader("上传小分子 SDF 文件", type=["sdf"])

            # 按钮触发预测
            if st.button("开始预测"):
                if protein_file is None:
                    st.error("请上传蛋白质 PDB 文件。")
                elif ligand_file is None:
                    st.error("请上传小分子 SDF 文件。")
                else:
                    with st.spinner("正在进行亲和力预测，请稍候..."):
                        try:
                            # 创建临时目录
                            with tempfile.TemporaryDirectory() as tmpdir:
                                # 保存上传的蛋白质文件
                                protein_path = os.path.join(tmpdir, protein_file.name)
                                with open(protein_path, "wb") as f:
                                    f.write(protein_file.getbuffer())

                                # 保存上传的小分子文件
                                ligand_path = os.path.join(tmpdir, ligand_file.name)
                                with open(ligand_path, "wb") as f:
                                    f.write(ligand_file.getbuffer())

                                # 输出 CSV 文件路径
                                output_csv_path = os.path.join(tmpdir, "single_prediction.csv")

                                # 调用预测脚本
                                pred_dir = "./others/PLANET"
                                pred_script = "pred.py"
                                pred_script_path = os.path.join(pred_dir, pred_script)

                                cmd = [
                                    "python",
                                    pred_script_path,
                                    "-p",
                                    protein_path,
                                    "-l",
                                    ligand_path,
                                    "-m",
                                    ligand_path,
                                    "-o",
                                    output_csv_path,
                                ]

                                result = subprocess.run(cmd, capture_output=True, text=True)

                                if result.returncode != 0:
                                    st.error(f"预测过程中发生错误:\n{result.stderr}")
                                else:
                                    if os.path.exists(output_csv_path):
                                        df = pd.read_csv(output_csv_path)
                                        st.success("预测完成！结果如下：")
                                        st.dataframe(df)
                                    else:
                                        st.error("预测完成但未找到输出 CSV 文件。")
                        except Exception as e:
                            st.error(f"发生异常: {e}")

        elif mode == "批量预测":
            st.subheader("批量蛋白与小分子亲和力预测")

            # 按钮触发预测
            if st.button("开始批量预测"):
                with st.spinner("正在进行批量亲和力预测，请稍候..."):
                    try:
                        batch_dir = "./Batch_Docking_Result"
                        if not os.path.exists(batch_dir):
                            st.error("批量预测目录不存在。")
                        else:
                            final_results = []

                            # 扫描文件夹中的 SDF 和 PDB 文件
                            sdf_files = [f for f in os.listdir(batch_dir) if f.endswith(".sdf")]
                            pdb_files = [f for f in os.listdir(batch_dir) if f.endswith(".pdb")]

                            st.write("发现以下蛋白质文件：")
                            st.write(pdb_files)
                            st.write("发现以下配体文件：")
                            st.write(sdf_files)

                            progress_bar = st.progress(0)
                            total_files = len(sdf_files)

                            for i, sdf_file in enumerate(sdf_files):
                                receptor_name = sdf_file.split("_")[0]
                                ligand_name = sdf_file.split("_")[1]
                                pdb_file = os.path.join(batch_dir, receptor_name + ".pdb")
                                sdf_file_path = os.path.join(batch_dir, sdf_file)

                                if os.path.exists(pdb_file):
                                    st.text(
                                        f"正在计算第 {i + 1}/{total_files} 对：蛋白 {pdb_file} 和 配体 {sdf_file} 的亲和力..."
                                    )
                                    with tempfile.TemporaryDirectory() as tmpdir:
                                        output_csv_path_tmp = os.path.join(
                                            tmpdir, "temp_result.csv"
                                        )

                                        cmd = [
                                            "python",
                                            "./others/PLANET/pred.py",
                                            "-p",
                                            pdb_file,
                                            "-l",
                                            sdf_file_path,
                                            "-m",
                                            sdf_file_path,
                                            "-o",
                                            output_csv_path_tmp,
                                        ]

                                        result = subprocess.run(cmd, capture_output=True, text=True)

                                        if result.returncode == 0 and os.path.exists(
                                            output_csv_path_tmp
                                        ):
                                            temp_df = pd.read_csv(output_csv_path_tmp)
                                            if "Binding_Affinity" in temp_df.columns:
                                                binding_affinity = temp_df["Binding_Affinity"].iloc[
                                                    0
                                                ]
                                                final_results.append(
                                                    {
                                                        "Protein_File": receptor_name,
                                                        "Ligand_File": ligand_name,
                                                        "Binding_Affinity": binding_affinity,
                                                    }
                                                )
                                        else:
                                            st.error(f"文件 {sdf_file} 处理失败。")

                                # 更新进度条
                                progress_bar.progress((i + 1) / total_files)
                                time.sleep(0.1)  # 模拟计算时间

                            if final_results:
                                results_df = pd.DataFrame(final_results)

                                # 保存结果到 Binding_Affinity 文件夹
                                binding_affinity_dir = "./Result/Binding_Affinity"
                                os.makedirs(binding_affinity_dir, exist_ok=True)

                                output_csv_path = os.path.join(
                                    binding_affinity_dir, "batch_prediction_results.csv"
                                )
                                results_df.to_csv(output_csv_path, index=False)

                                st.success("批量预测完成！结果已保存到以下目录：")
                                st.write(output_csv_path)

                                # 绘制热图
                                st.subheader("亲和力热图")
                                heatmap_data = results_df.pivot(
                                    index="Protein_File",
                                    columns="Ligand_File",
                                    values="Binding_Affinity",
                                )
                                plt.figure(figsize=(10, 8), dpi=600)
                                sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".2f")
                                plt.xlabel("Ligands")
                                plt.ylabel("Proteins")
                                plt.title("Binding Affinity Heatmap")
                                st.pyplot(plt)

                                heatmap_path = os.path.join(
                                    binding_affinity_dir, "binding_affinity_heatmap.png"
                                )
                                plt.savefig(heatmap_path, dpi=600)

                                st.write("热图已保存到以下目录：")
                                st.write(heatmap_path)

                            else:
                                st.error("未生成任何预测结果。")
                    except Exception as e:
                        st.error(f"发生异常: {e}")

# ------------------------------------------------------------------------------
# 数据管理
# ------------------------------------------------------------------------------
elif page == "数据管理":
    import sys
    from pathlib import Path

    # Add parent directory to path
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))
    # Add FLASH_DOCK-main/services to path
    flashdock_services = Path(__file__).parent / "services"
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

# ------------------------------------------------------------------------------
# 服务监控
# ------------------------------------------------------------------------------
elif page == "服务监控":
    import sys
    from pathlib import Path

    # Add parent directory to path
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))
    # Add FLASH_DOCK-main/services to path
    flashdock_services = Path(__file__).parent / "services"
    sys.path.insert(0, str(flashdock_services))

    from service_manager import ServiceManager
    from compass_client import CompassClient

    st.title("服务监控")
    st.write("监控COMPASS服务状态")

    # Initialize service manager
    try:
        service_manager = ServiceManager()
        client = CompassClient()
        st.success("已连接到服务注册中心")
    except Exception as e:
        st.error(f"无法连接到服务注册中心: {e}")
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

# ------------------------------------------------------------------------------
# 训练管理
# ------------------------------------------------------------------------------
elif page == "训练管理":
    import sys
    from pathlib import Path

    # Add parent directory to path
    project_root = Path(__file__).parent.parent
    sys.path.insert(0, str(project_root))
    # Add FLASH_DOCK-main/services to path
    flashdock_services = Path(__file__).parent / "services"
    sys.path.insert(0, str(flashdock_services))

    # Import and execute training management page
    import importlib.util

    training_management_path = Path(__file__).parent / "pages" / "training_management.py"
    spec = importlib.util.spec_from_file_location("training_management", training_management_path)
    training_management = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(training_management)
