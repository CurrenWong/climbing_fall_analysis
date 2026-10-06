---
title: 踝关节 FE 建模 —— 几何与网格（经验 + 交接）
date: 2026-10-03
type: 技术交接 / 经验文档
status: 几何网格阶段（3 部件已验证，软骨阻塞）
tags:
  - 踝关节
  - 有限元
  - gmsh
  - STL
  - 四面体网格
  - 吴恺2012
  - IITD
aliases:
  - 踝关节FE交接
  - ankle FE handoff
---

# 踝关节 FE 建模 —— 几何与网格（经验 + 交接）

> **本文档定位**：把「踝关节 FE 的几何/网格」这一段的知识、已验成果、阻塞点
> 与下一步**完整落到纸面**，供上下文压缩后无缝续接。
>
> **与上位文档的关系**：
> - [`OpenSim_FE复现方案.md`](OpenSim_FE复现方案.md) 是**上位方案**（OpenSim→FEBio 全链，
>   S0 已完成）。本文档只覆盖其中的**几何/网格**子问题，不重复其架构与阶段计划。
> - [`../ANKLE_FE_TODO.md`](../ANKLE_FE_TODO.md) 是同源的进度清单（本文档是其展开版）。
> - 论文出处：[`../paper/人体足踝部有限元模型的建立及有效性分析_吴恺.pdf`](../paper)
>   （中国骨与关节外科 2012, 5(4):352-356），转文本在
>   `../paper/_text/吴恺_足踝FE.txt`。

---

## 0. 一句话状态

**3 个部件四面体网格已完成且体积误差 0.000%（含 20 万面的胫骨）；两个关节软骨
因源 STL 几何缺陷不可用，需改为从骨面法向 offset 生成。**

---

## 1. 目标与对标基准

复现吴恺 2012 的踝关节建模方法。**注意论文用的是 1 名志愿者的 CT + Mimics +
Geomagic + Ansys**；我们**没有那套数据与商业软件**，改用：

| 论文要素 | 我们的替代 |
|---|---|
| CT + Mimics 重建骨骼 | **IITD 开源踝足模型的真实 STL/OBJ**（`paper/ankleModel-latest.zip`） |
| Geomagic 生成软骨 | ⚠️ 待做：从骨面法向 offset（见 §4） |
| Ansys Workbench | **pyfebio 0.3.0 + FEBio 4.13.0** |
| 论文表 1 材料参数 | ✅ 直接可用（模型本身没给骨/软骨材料） |
| 韧带"线弹性 + 只受拉" | ✅ IITD 给了**更精细的非线性本构** |

### 验证目标（论文自己的自我验证方式）

| 量 | 论文值 | 我们的对标方式 |
|---|---|---|
| 胫骨下关节面最大接触应力 | 3.74 MPa（文献域 2.7–3.97） | 看**量级**是否落域 |
| 接触面积 | 335.5 mm² | 看量级 |
| 四条韧带位移排序 | 胫跟 5.24 > 胫舟 5.04 > 跟腓 2.05 > 下胫腓前 1.92 mm | 看**定性排序** |

> **不要对 3.74 这个精确值**：那是他的具体几何算出来的。论文自己也只是
> "与既往研究具有相似性"。

---

## 2. ✅ 已验证成果

### 2.1 三个部件的四面体网格（体积误差 0.000%，无负体积）

| 部件 | 源 | 源面数 | 实际面数 | 节点 | 四面体 | 体积 | 误差 | 质量 min/中 |
|---|---|---|---|---|---|---|---|---|
| 距骨 talus | `talus_r.stl` | 3,670 | 3,670 | 2,628 | 9,529 | **27.5596 cm³** | **0.000%** | 0.087 / 0.402 |
| 跟骨 calcaneus | `calcaneus_r.stl` | 6,718 | 6,718 | 5,247 | 20,145 | **57.6893 cm³** | **0.000%** | 0.062 / 0.420 |
| 胫骨远端 tibia | `OK_tibia.obj` | 200,000 | **6,000**（降采样） | 5,576 | 21,287 | **76.9127 cm³** | **0.000%** | 0.003 / 0.364 |

- 产物：`temp/ankle_mesh/{talus_r_8, calcaneus_r_8, tibia_6000}.msh`（gmsh v2.2）
- 边长范围：距骨 0.404–6.43 mm / 跟骨 0.229–6.31 mm / 胫骨 0.276–7.28 mm
- **零负体积**（Jacobian 全正）

> **胫骨的意义**：20 万面 → 降采样 6,000 面 → 网格成功，证明**方案 A
> （降采样跑通全链路）在最难部件上成立**。后续加密只需调 `decimate` 比例。

### 2.2 坐标系一致性（装配前提，已验证）

三者同一坐标系，单位 **米**（存 .msh 后按 cm 打印）：

```
部件          x (cm)              y (cm)              z (cm)
胫骨远端   -2.46 ~  2.96     -7.99 ~ -2.67     -3.78 ~  3.12
距骨       -2.12 ~  3.84     -2.29 ~  0.93     -2.33 ~  1.93
跟骨       -5.54 ~  2.15     -3.33 ~  1.93     -2.37 ~  2.27
```

**关键几何关系**：胫骨远端最低面 `y=-2.67` 与距骨最高面 `y=-2.29` 之间
**有 0.38 cm 间隙** —— 这正是关节腔/软骨的位置，**几何关系正确**。
x/z 方向三者重叠良好 → **装配可行**。

### 2.3 网格质量实测（2026-10-03 复核，修正 §2.1 的口径）

⚠️ **口径更正**：§2.1 表格里的「质量 min/中」来自 `scripts/ankle_fe/msh_info.py`，
它算的是 **`12/√2 · V / Lmax³`**（用**最长边**做归一，对偏斜单元惩罚极重）。
下表是用标准 **gamma = `12(3V)^(2/3) / ΣL²`** 重算的结果（两者对正四面体都 = 1）：

| 网格 | tet | **gamma min** | p5 | **gamma 中位** | <0.05 | <0.10 | 长宽比中位/p95/max | 判定 |
|---|---|---|---|---|---|---|---|---|
| `calcaneus_r_8` | 20,145 | 0.3795 | 0.632 | **0.8410** | 0.00% | 0.00% | 1.63 / 2.57 / 4.50 | ✅ |
| `talus_r_8` | 9,529 | 0.4309 | 0.628 | **0.8254** | 0.00% | 0.00% | 1.68 / 2.58 / 4.52 | ✅ |
| `tibia_6000` | 21,287 | 0.0422 | 0.460 | **0.7925** | 0.01% | **0.03%** | 1.76 / 3.74 / **10.59** | ✅ 可用 |

**结论**：
- **网格质量不是当前瓶颈**。gamma 中位 0.79–0.84，sliver（<0.10）在胫骨上仅
  **6 个单元（0.03%）**，其余为 0。全部 **负体积 0**。
- 唯一需处理的是胫骨有 **1 个长宽比 10.59** 的偏斜单元 → 需**局部重划**，
  但**不阻塞**求解（`[LRN-20261003-013]`：不要拿偏严的自定义指标当结论）。
- ❌ **缺失的是「网格收敛性研究」**（加密一档 σ_vm 变化应 < 20%，见上位方案 §七 T6）
  —— 一次都没做过。



---

## 3. ❌ 阻塞：两个关节软骨 STL 是缺陷几何

### 3.1 现象

| 软骨 | 边界边 (OCC 判定) | 封口所需容差 | 封口后体积 |
|---|---|---|---|
| `tibial_talus_cartilage_r.stl` | 18（tol 1e-6）→ 10（1e-5） | **仅 1e-3（1 mm）** | **−0.1228 cm³（负！）** vs 表面 2.326 |
| `talus_cartilage_r.stl` | 93（1e-6）→ 25（1e-5） | 1e-4 | **0.0357 cm³** vs 表面 4.150（**−99%**） |

**封口成功 = 几何被毁。** 用这样的实体会得到比失败更糟的结果。

### 3.2 排除的假设（三轮实验，证据链完整）

| 假设 | 检验方法 | 结果 |
|---|---|---|
| 有开口 | 统计每条边被几个三角形共享 | **13,581 条全是 2 次，0 例外** → 无开口 |
| 非流形 | 统计 >2 次的边 | **0 条** |
| 绕向不一致 | 统计共享边是否被反向遍历 | **13,581 反向配对 / 0 冲突** → 完全一致 |
| 精度不够 | 扫合并容差 1e-8 → 1e-3 | 18→10→5→0，但**同时体积被毁** |

**结论**：这是**拓扑完全干净的闭壳**，OCC 却缝合不了。根因是几何层面
（细长片 19–243 μm 被合并后消失，留下 OCC 认得的洞），而拓扑统计
在合并前的点上做，看不出来。

### 3.3 解决的思路（且更忠于论文）

吴恺原文 §1.2 明确写：
> "根据各关节面的几何形状，利用软件**加厚功能在各骨面生成关节软骨**"

**软骨本来就是从骨面 offset 生成的，不是分割出来的。**
所以用胫骨/距骨的关节面沿法向 offset 出软骨层 ——
**既绕开缺陷，又贴合原方法**。这是下一步的首选（见 §4）。

### 3.4 备选方案

1. **放弃软骨实体**，直接让骨-骨接触（会改变接触应力物理，不推荐）
2. 用软骨 STL 仅作**接触面位置参考**，不建实体
3. 换一个更好的踝足模型（需另找数据）

---

## 4. 下一步（按优先级）

### 4.1 决策点（需要用户拍板）

| # | 决策 | 选项 |
|---|---|---|
| D1 | **腓骨缺失怎么处理** | IITD 模型**无独立腓骨文件**，`OK_tibia.obj` 只是 5.3 cm 远端胫骨。而吴恺加载是"胫、腓骨下端上截面"。→ (a) 接受无腓骨简化；(b) 另找带腓骨的几何 |
| D2 | 软骨生成方式 | (a) 骨面法向 offset（推荐）；(b) 放弃软骨实体 |

### 4.2 待做清单

- [ ] 由骨面法向 offset 生成软骨层（1–2 mm，对齐论文"软骨较厚、与骨有交叠"）
- [ ] 按 D1 处理腓骨
- [ ] pyfebio 写 `.feb`：材料 + 骨-软骨 `bonded` + 软骨面面接触（μ=0.01）
- [ ] 600 N 垂直加载（胫腓骨下端上截面）→ FEBio 4.13 求解
- [ ] 读回接触应力/接触面积 → 对标 2.7–3.97 MPa / 335.5 mm²
- [ ] 加四条韧带 → 位移排序对标
- [ ] **把已验证的网格路径固化进 `scripts/ankle_fe/geom.py`**（见 §5.4）

---

## 5. 完整可复现流程

### 5.1 环境（一次性）

```cmd
:: venv（已存在）
D:\Project\climbing_fall_analysis\.venv\Scripts\python.exe      :: Python 3.12.12

:: 依赖（venv 里没有 pip，必须用 uv 的 --python 指过去）
D:\Program\uv\uv.exe pip install pyfebio gmsh pyvista --python .venv\Scripts\python.exe
:: 已装：pyfebio 0.3.0 / gmsh 4.15.2(pip wheel) / pyvista 0.49.0 / meshio / h5py

:: ★ gmsh 官方 SDK（pip wheel 不够用，见 §6.1）
:: 下载 https://gmsh.info/bin/Windows/gmsh-4.15.2-Windows64-sdk.zip (42 MB)
:: 解压到 D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk\
copy /y D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk\lib\gmsh-4.15.dll ^
        D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk\bin\
```

**关键路径常量**：

| 名称 | 路径 |
|---|---|
| FEBio 求解器 | `D:\Program\FEBioStudio\bin\febio4.exe`（4.13.0，**不在 PATH**） |
| gmsh SDK | `D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk` |
| gmsh DLL 目录 | `<SDK>\lib`（**不在 bin**！） |
| gmsh.exe | `<SDK>\bin\gmsh.exe` |

### 5.2 ⭐ 已验证可用的 STL→四面体 流程

> **核心**：**逐三角形手工重建 BREP**。这是官方示例
> `share/doc/gmsh/examples/api/stl_to_brep.py` 的做法，
> 在**闭合干净的 STL** 上给出 **0.000% 体积误差**。

```python
import os, sys
from pathlib import Path
import numpy as np

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK/'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))
import gmsh

gmsh.initialize(str(LIB))                      # ★ 必须指向官方 DLL
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.Verbosity", 0)

gmsh.open(str(stl))                            # 读离散三角网
et, _, en = gmsh.model.mesh.getElements(2)
tri = next(np.array(n).reshape(-1, 3) for e, n in zip(et, en) if e == 2)
ntags, ncoord, _ = gmsh.model.mesh.getNodes()
ncoord = np.array(ncoord, float).reshape(-1, 3)
tag2row = {int(x): i for i, x in enumerate(ntags)}
rows = np.vectorize(tag2row.get)(tri)          # ★ 三角形存的是 tag, 不是行号

gmsh.model.add("m")
occ = gmsh.model.occ
pts = [occ.addPoint(*ncoord[i]) for i in range(len(ntags))]   # 每节点一个点

e2c, surfs = {}, []
for a, b, c in rows:
    cur = []
    for u, v in ((a, b), (b, c), (c, a)):
        k = (min(u, v), max(u, v))
        t = e2c.get(k)
        if t is None:                          # ★ 边去重, 否则环建不起来
            t = occ.addLine(pts[u], pts[v]); e2c[k] = t
        cur.append(t)
    surfs.append(occ.addPlaneSurface([occ.addCurveLoop(cur)]))

occ.synchronize()
bnd = gmsh.model.getBoundary([(2, t) for _, t in gmsh.model.getEntities(2)])
if len(bnd) == 0:                              # ★ 封闭才能成实体
    sl = occ.addSurfaceLoop([t for _, t in gmsh.model.getEntities(2)])
    occ.addVolume([sl]); occ.synchronize()

gmsh.option.setNumber("Mesh.MeshSizeMax", char_len)
gmsh.option.setNumber("Mesh.MeshSizeMin", char_len * 0.4)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)   # Delaunay
gmsh.model.mesh.generate(3)
gmsh.model.mesh.optimize("Netgen")
```

**降采样（大部件必用，方案 A）**：

```python
import pyvista as pv
m = pv.read("OK_tibia.obj")            # 100,002 点 / 200,000 面
dm = m.decimate(1.0 - 6000 / m.n_cells)          # 目标 6000 面
dm = dm.clean(point_merging=True, tolerance=1e-7, absolute=True)
# 再去掉退化三角形 + 压缩点索引, 然后送 §5.2 流程
```

### 5.3 已有脚本

| 脚本 | 作用 | 状态 |
|---|---|---|
| `scripts/ankle_fe/geom.py` | STL→tet（**当前用的是 `classifySurfaces` 路径**） | ⚠️ 未在骨头上验证；软骨报 `overlapping facets` |
| `scripts/ankle_fe/msh_info.py` | 解析 .msh，报节点/四面体/体积/质量/负体积 | ✅ 可用 |
| `temp/pyfebio_demo/scan_parts.py` | **已验证的** plane-surface 路径 + 多部件批量 | ✅ 产出了 3 个骨网格 |
| `temp/pyfebio_demo/tibia_test.py` | 降采样 + 网格化（胫骨 20 万面） | ✅ 验证通过 |
| `temp/pyfebio_demo/cart_close.py` | 软骨封口容差扫描（证明不可用） | ✅ 证据来源 |
| `temp/pyfebio_demo/zip_list.py` | 从 ankleModel zip 抽几何 + 解析 `.osim` 参数 | ✅ |

### 5.4 ⚠️ 技术债

`scripts/ankle_fe/geom.py` 目前实现的是**官方 `classifySurfaces` 路径**
（`merge → removeDuplicateNodes → classifySurfaces(pi/2) → createGeometry →
addSurfaceLoop → addVolume`），**该路径只在软骨上测过（失败）**，
没在骨头上验证过。**已验证的是 §5.2 的 plane-surface 路径**。

→ **待办**：把 §5.2 路径固化进 `geom.py`（保留 classifySurfaces 作为备选分支）。

---

## 6. 踩坑清单（重要，避免重复踩）

### 6.1 ★ pip 的 gmsh wheel ≠ 官方 SDK

**pip wheel 的 `occ` 绑定残缺**，实测：

| 调用 | 结果 |
|---|---|
| `occ.importShapes(stl, format="stl")` | `Unknown file type`（6 种 format 全试） |
| `occ.reclassify(1, 1)` | `AttributeError`（不存在） |
| `geo.addTriangle(...)` | `AttributeError`（geo 无此函数） |
| `mesh.importStl(path)` | `takes 0 positional arguments` |
| `gmsh.merge(stl)` → `healShapes` | `Unknown OpenCASCADE entity of dimension 2` |
| `gmsh.open(stl)` | 只得 "Discrete surface"，`getEntities(3)` 为空 |

但 `General.BuildInfo` 显示 **OCC 7.8.1 确实编译进去了** ——
所以**不是内核缺失，是绑定/SDK 不完整**。
pip wheel **只有 `gmsh.py`，没有 `gmsh.exe`**。

→ **必须用官方 Windows64 SDK。**

### 6.2 gmsh 其他坑

| 坑 | 正确做法 |
|---|---|
| DLL 在 `<SDK>\lib` **不在 bin** | `gmsh.initialize(str(LIB))`；`gmsh.exe` 需把 dll 拷到 `bin\`（否则 `3221225781` = DLL_NOT_FOUND） |
| `Geometry.ToleranceMergePoints` **不存在** | 正确的是 **`Geometry.Tolerance`** |
| `occ.addDiscreteSurface` / `occ.reclassify` 不存在 | 用 §5.2 的逐三角形重建 |
| 逐三角形 `addPlaneSurface` 的面**互相独立**，OCC 不认共面 | 这才是 `classifySurfaces` 的用途；但对已闭合 STL，前者反而更直接可靠 |
| `occ.healShapes` 认不了 "Discrete surface" | 先 `open` 拿网格自己建点，别指望 heal |

### 6.3 pyvista / VTK 坑

| 坑 | 正确做法 |
|---|---|
| `CellType.TETRA` 编号 | **pyvista 里 = 10**（VTK 原本 12，pyvista 重编过） |
| `delaunay_3d` 参数名 | 是 `alpha` / `tol` / `offset`，**没有** `tolerance` / `depth` |
| `delaunay_3d` 结果 | ❌ **会填满凹陷**：距骨 +32.5%、跟骨 +46.5%、软骨 +148.7% → **不可用** |
| `bounds` 是 tuple | 要 `np.asarray(bounds, float)` 才能相减 |
| h5py 读 xplt 结果 | `displacement` 是**纯 float32 (n,3)，列=xyz**，不是 recarray；`u[4,2]` 而非 `u[4]["z"]` |

### 6.4 cmd.exe 坑

- `python -c "多行 + 缩进"` → **`IndentationError`**（cmd 把换行折成空格）
  → **多行 Python 一律写 `.py` 再跑**（本轮踩了 4 次）
- `where` / `rm` / `tr` 不可用 → 用 `dir /b` / `del` / PowerShell

---

## 7. 已提取的参数（`temp/ankle_geom/osim_params.json`）

来源：`ankleModel-latest.zip` 内 `Model/IITD_Ankle_foot_model_NLSRclass.OSIM`。
单位 **m / N**。（该 `.osim` 是**纯运动学 + 力元件**，本身不含 FE 材料。）

### 7.1 骨/软骨材料（来自吴恺论文表 1，模型没给）

| 材料 | E (MPa) | ν | 建模方式 |
|---|---|---|---|
| 骨 | 7300 | 0.3 | 各向同性线弹性实体 |
| 软骨 | 260 | 0.4 | 各向同性线弹性实体 |
| 韧带 | 10 | 0.4 | 只受拉单轴梁，**每条只 1 个单元**（避免应力集中） |

其他：关节间摩擦系数 **0.01**；软骨间隙 **< 0.1 mm**；骨-软骨 **bonded**；
网格 automatic / coarse。

### 7.2 踝部韧带（IITD 非线性本构，比论文更精细）

共同参数：`eta_a=70`，`alpha=1.148`，`beta=2`，`k_a=1000`，`max_strain=0.45`

| IITD 名 | 对应吴恺的韧带 | L0 (m) | pcsa_force (N) |
|---|---|---|---|
| `pTTL_r` | **下胫腓联合前韧带** | 0.0200 | 39.1 |
| `aTTL_r` | 胫距前（aTTL） | 0.0250 | 19.6 |
| `aTFL_r` | 胫距前韧带 | 0.0290 | 47.5 |
| `CFL_r` | **跟腓韧带** | 0.0360 | 10.1 |
| `TNL_r` | **胫舟韧带** | 0.0400 | 450.0 |
| `TCL_r` | 胫距后 | 0.0310 | 19.0 |
| `pTFL_r` | 腓距后 | 0.0310 | 4.5 |
| `TPS_r` | 距腓后突 | 0.0181 | 1032.8 |
| `FPS_r` | 腓距后 | 0.0202 | 1032.8 |

> 每条韧带都带**精确附着点坐标**（在 `osim_params.json` 的 `attachment_points`
> 字段：`name` / `body` / `location_m`），可直接用作梁单元端点。

### 7.3 踝部肌肉（173 条中踝相关）

| 名 | Fmax (N) | 最优肌纤维长 (m) | 腱松弛长 (m) | 羽状角 (°) |
|---|---|---|---|---|
| `soleus_r` | 6194.8 | 0.0440 | 0.2768 | 0.38 |

### 7.4 关节运动学

| 关节 | 坐标名 | 范围 (rad) | 对应角度 |
|---|---|---|---|
| `Ankle_r` | `ankle_pfdf_r` | [-0.8727, 0.3491] | -50° ~ +20°（跖屈/背屈） |
| `Subtalar_r` | `subtalar_inev_r` | [-0.6109, 0.2094] | -35° ~ +12°（内翻/外翻） |
| `Chopart_r` | `chopart_obl_r` | [-0.1745, 0.0524] | -10° ~ +3° |

---

## 8. 几何资源清单

**抽出位置**：`temp/ankle_geom/`（用 `temp/pyfebio_demo/zip_list.py` 从
`paper/ankleModel-latest.zip` 内的 `Model.zip` 提取）

| 文件 | 内容 | 状态 |
|---|---|---|
| `talus_r.stl` | 距骨 | ✅ 已网格化 |
| `calcaneus_r.stl` | 跟骨 | ✅ 已网格化 |
| `OK_tibia.obj` | 胫骨远端（5.3 cm） | ✅ 已网格化（降采样） |
| `tibial_talus_cartilage_r.stl` | 胫侧关节软骨 | ❌ 缺陷几何 |
| `talus_cartilage_r.stl` | 距侧关节软骨 | ❌ 缺陷几何 |
| `midfoot_r.stl` | 中足 | 未测 |
| `metatarsal_*.stl` `toes_*.stl` | 跖骨/趾骨 | 未测 |
| `OK_femur.obj` | 股骨 | 未测 |
| `femur_cartilage_*.obj` `*_meniscus_*.obj` | 膝部软骨/半月板 | 未测 |
| `tibial_plane.obj` | 胫骨平面（98 顶点，102×10×40 mm） | 可能是加载面参考 |
| `IITD_Ankle_foot_model_NLSRclass.OSIM` | 运动学 + 力参数源 | ✅ 已解析 |
| `osim_params.json` | 上文 §7 的结构化参数 | ✅ 已生成 |

**⚠️ 缺腓骨**：整个 zip 里没有独立 `fibula` 几何。见 §4.1 D1。

---

## 9. 关键教训（已入 `.learnings/`）

| ID | 教训 |
|---|---|
| `LRN-20261003-010` | **「pip 装了同一个库」≠「你能用它的全部 API」**。同一版本号可以是两份不同二进制；判断依据是**实测 API 可调用性**，不是版本号。用户说"换方案"时要换**真的不同的东西**——在已确认失败的载体上换写法是伪装的重复尝试（判据：如果重试的还是同一个文件/进程，那就不是换方案）。 |
| `LRN-20261003-011` | **「拓扑干净」≠「几何可缝合」**。"边被恰好 2 个面共享"是拓扑条件；几何缝合还要求共享边处连续。**参数扫描要扫到能做/不能做的边界，并检查解的可用性**——只扫到"✅ 封闭"就收工会拿到**体积为负**的实体，比失败更糟。**看到成功要立刻质疑：这个成功的产物还有意义吗？** |
| `ERR-20261003-006` / `007` | pyfebio 0.3.0 只写不能读；FEBio 不在 PATH / xplt 布局坑 |

---

## 10. 快速自检清单（续接时先跑这个）

```cmd
:: 1. 环境
D:\Project\climbing_fall_analysis\.venv\Scripts\python.exe -c "import pyfebio,gmsh,pyvista;print('OK')"
dir D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk\lib\gmsh-4.15.dll
dir D:\Program\FEBioStudio\bin\febio4.exe

:: 2. 已有网格体检（应报 3 个部件, 体积误差 0, 负体积 0）
cd /d D:\Project\climbing_fall_analysis
.venv\Scripts\python.exe scripts\ankle_fe\msh_info.py ^
    temp\ankle_mesh\talus_r_8.msh temp\ankle_mesh\calcaneus_r_8.msh ^
    temp\ankle_mesh\tibia_6000.msh
```

**期望输出**：

```
talus_r_8.msh      节点   2,628 四面体    9,529 vol= 27.5596 cm³ | 负体积 0
calcaneus_r_8.msh  节点   5,247 四面体   20,145 vol= 57.6893 cm³ | 负体积 0
tibia_6000.msh     节点   5,576 四面体   21,287 vol= 76.9127 cm³ | 负体积 0
```

---

## 11. 模型选型评估：Maharaj2021 vs IITD（2026-10-03 实测）

> 起因：用户问 `model/FootwithinLowerlimbModel` 换掉 IITD 是否更好。
> **结论：看目的，一分为二 —— 对 FE 几何更差，对多体动力学更好。**

### 11.1 该目录是什么

`model/FootwithinLowerlimbModel/` = **Maharaj2021「Foot within lower limb」模型**
（source.txt → `https://simtk.org/projects/footankle_model`）。
含 3 个 `.osim` + 10 个足部 STL + **真实实验数据**（`static.c3d` / `walk.c3d` /
`walk.trc` / `walk_ik.mot`）+ Scale/IK 任务文件。**`paper/` 里没有它的 zip**（另行下载）。

### 11.2 实测对比表（OpenSim 4.6 真实加载）

| 模型 | bodies | joints | coords | muscles | mass(kg) | 能否加载 |
|---|---|---|---|---|---|---|
| Maharaj2021_BothLegs | 17 | 17 | 37 | 93 | 49.31 | ✅ |
| Maharaj2021_BothLegs_V2 | 17 | 17 | 37 | 93 | 74.44 | ✅ |
| Maharaj_SampleData_SCALED | 17 | 17 | 37 | 93 | 90.00 | ✅ |
| Rajagopal FullBodyModel-4.0 | 22 | 22 | 39 | 80 | 75.34 | ✅ |
| **IITD 踝足** | — | — | — | — | — | ❌ **加载失败** |

**⚠️ IITD 的 OSIM 在原生 OpenSim 4.6 无法加载**：

```
[error] Object::newInstanceOfType(): object type 'NonLinearLigament'
        is not a registered Object! It will be ignored.
→ RuntimeError: std::exception in 'OpenSim::Model::Model(std::string const&)'
```

根因就在文件名里：`IITD_Ankle_foot_model_NLSRclass.OSIM` —— **NLSR = NonLinear
Spring/Ligament，是自定义插件类**，stock OpenSim 不认识。

> **重大含义**：我们此前只是把该 `.osim` 当**文本**解析出参数（`osim_params.json`），
> **从未真正把它当模型加载过**。所以 **IITD 不能用于多体动力学**，
> 除非拿到并编译它的 NLSR 插件。这直接影响上位方案
> `OpenSim_FE复现方案.md` 的 S1/S3 —— 那条链需要的是**能加载的**模型。

**Maharaj 的 body 构成（17）**：pelvis + ``femur_r``, `tibia_r`, `rpatella`,
`talus_r`, `calcn_r`, `midfoot_r`, `forefoot_r`, `digits_r`（左腿同）。
**足部 5 关节**（`ankle_r` / `subtalar_r` / `midtarsal_r` / `tarsometatarsal_r` /
`digits_r`）比 IITD 的 3 个（ankle / subtalar / chopart）更细。

### 11.3 「是否包含胫骨腓骨」—— 分两层回答

| 层面 | 结论 |
|---|---|
| **运动学 body** | ✅ 有 `tibia_r` / `tibia_l`；**但无独立 fibula body**。腓骨是 tibia body 的**第 2 个显示网格**：`tibia_r_geom_2` → `<mesh_file>fibula_r.vtp</mesh_file>`（Rajagopal 惯例） |
| **可网格化几何** | ❌ **一份都没有**。21 个 `<geometry_file>` 引用里只有 **10 个**在本地（全是足部 STL），**8 个两处都缺**：`tibia_r.vtp`、`fibula_r.vtp`、`femur_r.vtp`、`femur_l.vtp`、`tibia_l.vtp`、`pat.vtp`、`l_pat.vtp`、`pelvis.vtp` |

`<vertices>` 计数 = **0** → `.osim` **不内嵌网格**，全是 `<mesh_file>` 外部引用
（`SampleData_SCALED.osim` 亦然）。**该目录没有 Geometry 文件夹**，
且 OpenSim 加载时如实报警：`Couldn't find file 'tibia_r.vtp'` ×11。

### 11.4 补上 Rajagopal 的 vtp 也不解决问题（关键）

`model/opensim/FullBodyModel-4.0/Geometry/` 里**确有**胫/腓骨 vtp，但**只是低精度可视化网格**：

| vtp | 三角面 | 闭合 | 体积 cm³ | ext cm |
|---|---|---|---|---|
| `r_tibia.vtp` | **390** | ❌ False (euler=13) | 211.5 | 5.75 × 37.60 × 6.89 |
| `r_fibula.vtp` | **240** | ✅ | 55.5 | 2.87 × 36.40 × 3.17 |
| `r_femur.vtp` | 908 | ✅ | 348.5 | 7.03 × 45.27 × 9.95 |
| `r_patella.vtp` | 140 | ✅ | 27.9 | 2.61 × 5.29 × 5.14 |
| `r_talus.vtp` | 194 | ✅ | 21.2 | 4.61 × 3.01 × 4.16 |

- 37.6 cm 的胫骨只有 **390 面** → 平均每面 ~2 cm；
  我们已网格化的 IITD 胫骨用了 **200,000 面**。**量级差 500×，做应力 FE 无意义。**
- `r_tibia.vtp` **还不闭合** → 连实体都建不成。
- ⚠️ **命名不一致**：Maharaj 引用 `tibia_r.vtp` / `fibula_r.vtp`（后缀式），
  Rajagopal 目录是 `r_tibia.vtp` / `r_fibula.vtp`（前缀式）→ 按名找不到。
  仅 `l_fibula.vtp` / `sacrum.vtp` / `l_pelvis.vtp` 三个刚好同名。

### 11.5 它自带的足部 STL 是降采样版（比我们手上的差）

| 骨 | Maharaj | IITD 原版 | 体积变化 |
|---|---|---|---|
| `FootModel_talus_r.stl` | 460 pts / **916 tris** | 1,837 / **3,670** | 20.461 vs 27.560 cm³（**−26%**） |
| `FootModel_calcn_r.stl` | 550 / **1,096** | 3,361 / **6,718** | 61.342 vs 57.689（+6%） |
| `FootModel_midfoot_r.stl` | 915 / **1,810** | 9,916 / **19,838** | 19.192 vs 35.809（**−46%**） |

体积大幅偏移 ⇒ 降采样**破坏了几何**。这些都是**为了可视化/多体驱动**而精简的版本，
**不能当 FE 几何源**。

### 11.6 净结论与建议

| 目的 | 判定 | 理由 |
|---|---|---|
| **FE 骨几何**（当前阻塞点） | ❌ **别换** | 胫腓骨几何缺失；能补的是 390/240 面低精度且不闭合；自带足部 STL 还是降采样版 |
| **OpenSim 多体动力学**（上位方案 S1/S3） | ✅ **收进来** | 能加载（IITD 加载失败）；足部 5 关节；含真实 c3d/trc/IK 数据 |

**对当前两个阻塞的净帮助 = 零**：不给更好的腓骨、不给软骨。

**但有一个有用副产品**：Maharaj2021 论文标题就是 *Foot within lower limb*，
说明「IITD 足 + Rajagopal 下肢」这套拼法是**有人正式发表过的** ——
若将来要拼腓骨，有对齐先例可查。

**可选的后续动作**（需用户确认，会在 `model/` 新建文件）：
把 Rajagopal 的 `r_tibia.vtp`→`tibia_r.vtp`、`r_fibula.vtp`→`fibula_r.vtp` 等
**改名拷贝**进 Maharaj 目录，让几何在 OpenSim 里能显示。
**注意**：这只解决"看得见"，**不改变 §11.4 的精度结论**。

---

## 12. ⭐「更好的有限元效果」的最大杠杆：THUMS V7.1 现成的 FE 骨网格

**发现日期**：2026-10-03。**结论：项目里已经躺着一套比 IITD 好得多的 FE 骨骼网格，
此前只被当成"参考"而没被当成几何源。**

### 12.1 在哪 / 是什么

`paper/AM50_V71_Occupant.zip` → `AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k`
= **219.4 MB 的 LS-DYNA 关键字文件，即 THUMS V7.1 全身有限元模型本体**。

> 此前 `OpenSim_FE复现方案.md` §2.1 把它标为「**仅作几何/材料/阈值参考**
> （无 LS-DYNA 许可，不求解）」——「不求解」是对的，但**"几何可提取"这点被漏掉了**。

### 12.2 实测：下肢骨是完整的、皮质/松质分离的 FE 部件

全文件扫描结果 —— `*NODE` × 1、`*ELEMENT_SOLID` × 1、`*ELEMENT_SHELL` × 1173、
`*PART` × 1834、`*SECTION_SOLID` × 1。**网格数据齐全，可提取。**

| 部件 ID | 名称 | 说明 |
|---|---|---|
| 81000700 | **`R_TIBIA_CORT`** | 胫骨**皮质壳**（shell） |
| 81000600 / 81000601 | `right_tibia_end_spon` / `_center_spon` | 胫骨**松质芯**（solid） |
| 81000900 | **`R_FIBULA_CORT`** | **腓骨皮质壳** ← **我们缺的东西** |
| 81000800 / 81000801 | `right_fibula_end_spon` / `_center_spon` | 腓骨松质芯 |
| 81001100 / 81001000 | `R_TALUS_CORT` / `R_TALUS_SPON` | 距骨 |
| 81001300 / 81001200 | `R_CALCANEUS_CORT` / `R_CALCANEUS_SPON` | 跟骨 |
| 81000100 / 81000000-01 | `R_FEMUR_CORT` / `right_femur_*_spon` | 股骨 |
| 81000200-81000301 | `R_PATELLA_SPON` / `R_PATELLA_CORT` | 髌骨 |
| 81001400-81003701 | 舟骨/骰骨/3 楔骨/5 跖骨/全部趾骨（各含 SPON+CORT） | 全足 |
| 81100801 | **`L_LIG_CALCANEOFIBULARE`** | **跟腓韧带** |
| 81100701 / 81100901 | `L_LIG_TALOFIBULARE_POSTERIUS` / `_ANTERIUS` | 距腓后/前韧带 |
| 81101101 | **`L_LIG_TIBIOFIBULARE`** | **下胫腓韧带** |
| 81100100-81100201 | `R_LIG_COLLATERAL_FIBULAR/TIBIAL_TISSUES` | 内外侧副韧带 |
| 81101301-09 | `R_LIG_FOOT 1..9` | 足部韧带组 |

（左右各一套；`*MAT_MUSCLE` × 808 为肌肉。）

### 12.3 为什么这是"更好的有限元效果"的最大杠杆

| 收益 | 说明 |
|---|---|
| ✅ **拿到腓骨** | 直接解决 §4.1 D1「IITD 无腓骨」的阻塞 |
| ✅ **皮质/松质分层** | THUMS 把 cortical / spongy 拆成**独立部件** → 白送**非均质骨**，而吴恺的单 E=7300 MPa 无法表达骨干 vs 骨端的刚度差 |
| ✅ **[Y25] 同一解剖** | [Y25] 用的就是 THUMS AM50 → 对标从"跨模型"变成"同模型"，**这是最大的一次有效性提升** |
| ✅ **已是 FE 网格** | 不用再 mesh STL，且元素尺寸是按冲击 FE 设计的（~1–3 mm 量级） |
| ✅ **韧带已建好** | 跟腓/距腓/下胫腓都在，可直接对标吴恺的 4 条韧带位移排序 |
| ✅ **材料卡齐全** | `mat_bone_AM50_V71_Occ_fracture.k` / `_no_fracture.k`（含骨折本构） |
| ✅ **文档** | `THUMS_AM50_V71_Documentation_202411.pdf`（11.7 MB）含材料与单元规格 |

### 12.4 代价与风险（必须先说清）

| 风险 | 说明 | 缓解 |
|---|---|---|
| **许可** | 文件头写明 THUMS USER POLICY：仅注册用户可"refer to, use and share"。**内部研究可用，但不得再分发/外传** | 抽出的网格只在本项目内用；发布论文时按 Toyota 要求引用 |
| **单位制** | LS-DYNA 惯用 mm–ms–kg–kN → 应力单位是 **GPa**，与 FEBio 侧 mm–N–MPa 差 **1000×** | 显式换算 + 单位自检（§七 T10，这是本类项目最典型的静默失效） |
| **格式** | 219 MB 单一 `*NODE` 块 + 单一 `*ELEMENT_SOLID` 块，需按 part ID 拆分 | 需要写一个流式 .k 解析器（复用 `temp/pyfebio_demo/scan_thums.py` 的流式读法） |
| **壳单元** | 皮质骨在 THUMS 里是 **shell**（有厚度属性），不是实体 | 两条路：FEBio 用 shell 单元；或把 shell 加厚+与 spon solid 合并重划 tet |
| **解剖年龄** | AM50 = 50 百分位成年男性，与 IITD 受试者不同 | 这是**换源**不是"混用"——全身同一模型，内部自洽 |

### 12.5 ✅ 提取验证已完成（2026-10-03 实测，结论：可用）

**"能扫到关键字" ≠ "能成功提取网格"** —— 这条敬畏是对的，但**验证已通过**。

#### 实测踩到的坑（差点静默全废）

同一个 `.k` 文件里**三种字段宽度并存**，按空白 `split()` 会让
**1,696,785 行（占 99.998%）只解析出 1 个 token** 而静默丢弃：

| 块 | 真实格式 |
|---|---|
| `*NODE` | I8 + 3×F16，**有空格** → 空白分隔可用 |
| `*ELEMENT_SOLID` | **8 字符定宽、完全无空格**（行=80 字符）`eid,pid,n1..n8` |
| `*ELEMENT_SHELL` | **8 字符定宽、完全无空格**（行=48 字符）`eid,pid,n1..n4` |
| `*PART` 数值卡 | **10 字符定宽**（行=90 字符） |

> 第一次探测时抓到的"带空格样本"只是极少数特例，**据此写解析器直接全废**。
> 教训：**探测格式要看分布（token 数直方图 / 行长度直方图），不能只看前几行。**
> （`[LRN-20261003-014]`）

#### 提取结果：26 个部件全部成功，0 解析失败、0 缺失节点

| pid | 类型 | 单元 | 名称 |
|---|---|---|---|
| **81000900** | **hex** | **842** | **`R_FIBULA_CORT`** ← 我们缺的腓骨 |
| 81000800 / 81000801 | tet | 762 / 2,216 | `right_fibula_end_spon` / `center_spon` |
| 81000700 | hex | 2,417 | `R_TIBIA_CORT` |
| 81000600 / 81000601 | tet | 6,651 / 14,482 | `right_tibia_end_spon` / `center_spon` |
| 81001100 / 81001000 | hex / tet | 394 / 1,766 | `R_TALUS_CORT` / `SPON` |
| 81001300 / 81001200 | hex / tet | 634 / 3,323 | `R_CALCANEUS_CORT` / `SPON` |

+ 踝部韧带：跟腓 12、距腓后 6、距腓前 18、**下胫腓 20**、距跟后 18、
三角韧带 4 组 45、足底长韧带 60、内外侧副韧带（实体）125

#### 五项验证全部通过

| # | 验证项 | 结果 |
|---|---|---|
| ① | 单元/节点数非零 | ✅ 保留 21,608 节点 |
| ② | **单位 = mm** | ✅ 股骨 45.20 / 胫骨 32.76 / 腓骨 30.43 / 距骨 5.73 cm，全符合解剖 |
| ③ | **体积可靠** | ✅ **3 种独立算法吻合 ≤0.15%**：自写 hex 分解 / 边界积分 / **VTK 第三方库**；单位立方体自检两法均 1.0000000000；均 < 凸包上界 |
| ④ | 单元形状 | ✅ 皮质 = **hex**、松质 = **tet**（混用已正确处理）；CORT/SPON **共享大量节点 = 共形贴合**（非两份独立网格，无重复计数） |
| ⑤ | **装配可行性** | ✅ **天然正确，无需任何对齐**（见下） |

#### ⭐ 最大的意外收获：装配关系天然成立

THUMS 是**全身模型**，各骨同处**一个全局坐标系**。实测部件间最近节点距：

| 一对 | 最近距 | 对应关节 |
|---|---|---|
| `R_TIBIA_CORT` ↔ `R_TALUS_CORT` | **1.41 mm** | 踝关节 |
| `R_FIBULA_CORT` ↔ `R_CALCANEUS_CORT` | **2.78 mm** | 跟腓 |
| `R_FIBULA_CORT` ↔ `R_TIBIA_CORT` | **1.41 mm** | 下胫腓 |
| `R_TALUS_CORT` ↔ `R_CALCANEUS_CORT` | **1.04 mm** | 距下 |

各骨相互贴合在 **1–3 mm**（正是关节软骨/韧带间隙的量级）→ **IITD 时代需要自己
拼装的装配工作，在 THUMS 里直接免掉**。

#### ⚠️ 一个必须声明的差异：THUMS 跗骨比 IITD 大 1.1–1.45×（线性）

| | THUMS | IITD | 线性比 | 体积比 |
|---|---|---|---|---|
| 距骨 PCA 主轴 (mm) | 68.8 / 52.1 / 43.3 | 59.0 / 40.7 / 30.1 | ~1.24 | **1.84** |
| 跟骨 PCA 主轴 (mm) | 89.9 / 54.6 / 48.6 | 81.2 / 51.0 / 41.1 | ~1.11 | **1.66** |

`1.24³ = 1.91 ≈ 观测 1.84` ⇒ **体积差完全由均匀线性尺度解释**，
是**真实的个体差异**（THUMS AM50 = 175 cm 美国 50 百分位男性 vs IITD 受试者），
**不是提取错误**（三条独立证据：3 法体积吻合、四肢长骨长度正确、体积比符合线性比）。

> **含义**：**两个模型的骨不可混用**（尺寸对不上）；跨模型比较应力**绝对值**时
> 必须声明解剖来源 —— 这也再次说明**不要对吴恺的 3.74 MPa 精确值**。

#### 产物

| 文件 | 内容 |
|---|---|
| `scripts/ankle_fe/thums_extract.py` | 流式提取器（改 `TARGETS` 即可抽其它部件） |
| `scripts/ankle_fe/thums_verify.py` | 验证器（自检 + 双算法体积 + 共形检查） |
| `temp/thums/thums_lowerlimb.npz` | 21,608 节点 + 27 部件单元连接 |
| `temp/thums/summary.json` / `verify.json` | 结构化结果 |

### 12.6 决策：要不要换几何源？

| 维度 | IITD（现状） | THUMS V7.1 |
|---|---|---|
| 胫骨 | 仅远端 5.3 cm，200k 面 | **全长 32.8 cm**，皮质+松质 |
| **腓骨** | ❌ **没有** | ✅ **有**（842 hex + 2,978 tet） |
| 软骨 | ❌ 源 STL 缺陷，需 offset 生成 | ❌ 同样没有独立软骨（需 offset） |
| 骨材料分层 | ❌ 单一 E=7300 MPa | ✅ **皮质/松质分离**（非均质骨） |
| 装配 | 需自己拼 | ✅ **天然对齐**（1–3 mm） |
| 韧带 | IITD 参数（但模型加载不了） | ✅ 8 组踝韧带实体/壳 |
| 解剖来源 | 特定受试者 | AM50（= [Y25] 同源） |
| 许可 | 开源 | ⚠️ THUMS USER POLICY，不可再分发 |
| 单位 | m | mm（且 LS-DYNA 应力是 **GPa**，进 FEBio 需 ×1000） |

**建议**：换。THUMS 一次性解掉「缺腓骨」「骨头均质」「装配」三个问题，
且与 [Y25] 同源使对标从"跨模型"变"同模型"。
**但软骨仍需自己 offset 生成**（两边都没有），且跗骨尺寸差异要在结论里声明。

---

## 13. 材料 / 软骨 / FEBio 管线（2026-10-03，见 `FE_PIPELINE_TODO.md`）

> **压缩上下文后接手 → 先读 `FE_PIPELINE_TODO.md` 的「🧭 接手入口」节**
> （状态 / 一键复现 / 下一步排序 / 负面结果清单 / 环境开关 / 别再踩的坑）。
> 本节 §13.1–13.7 是该入口的展开细节。

按用户要求按 **c → b → a** 执行。**C、B、A 全部完成** —— FEBio 求解已收敛并产出物理结果。

| 阶段 | 状态 | 产物 |
|---|---|---|
| **C** 材料参数 | ✅ | `scripts/ankle_fe/thums_materials.py` → `temp/thums/bone_materials.json` |
| **B** 软骨生成 | ✅（零穿透） | `scripts/ankle_fe/cartilage.py` → `temp/thums/cartilage{,.mesh}.json/npz` |
| **A** FEBio 求解 | ✅ **收敛** | `scripts/ankle_fe/build_feb.py` → `temp/thums/febio/ankle_tibiotalar.feb` |

### 13.1 C 的三个关键结论

1. **单位是 `mm-s-tonne-N-MPa`** ⇒ **E 已是 MPa，进 FEBio 零换算**
   （**推翻 §12 里"LS-DYNA 应力是 GPa"的提醒** —— 那只适用 mm-ms-kg-kN 制）。
   校验方式：ρ=2.0E-9 若为 tonne/mm³ 则 =2000 kg/m³（皮质骨 ✓）；E=18000 MPa = 18 GPa ✓。
2. **`*_center_spon` 不是松质骨，是骨髓腔**（MAT_VISCOELASTIC，E=12 MPa、ν=0.499、ρ=1000）。
   长骨的"松质"只在**两端**（`_end_spon`）。
3. **跗骨皮质建模方式与长骨不同**：长骨 SIGY=36 MPa+硬化曲线；跗骨 SIGY=224 MPa+零硬化。

### 13.2 B 的核心：软骨厚度必须由**实测间隙**约束

踝关节骨间最小距离仅 **1.406 mm**。若两侧各取 2 mm 会互相穿透 2.6 mm
（FE 里表现为初始过盈 → 虚假高接触应力）。
故取 `t = min(1.5 mm, 0.5 × 局部间隙)`，实测**零穿透**。

**THUMS 完全没有关节软骨部件**（只有肋骨软骨/椎间盘/半月板）⇒ 必须自己生成。

### 13.3 A 的收获：FEBio 4.13 的 `.feb` 语法已全部实测确定

踩了 8 轮 `unrecognized tag` / `invalid value` 才摸清（**详见 `FE_PIPELINE_TODO.md` §A
—— 那里有一份可直接复制的完整模板**）。三个最反直觉的点：

- `SolidDomain` 的 **`mat` 要填材料名**（不是数字 id）
- **`<Surface>` 的子标签按单元类型命名**（`<quad4>`/`<tri3>`，不是 `<elem>`）
- **段序**必须 `Loads → Boundary → Contact`，顺序错会被报成"标签不识别"

**方法论教训**：同类错误连续 3 次失败就该换信息源。
本次最后靠 **pyfebio 生成样例 XML** + **本机 `FEBio_User_Manual.pdf`** 一次全对
（`[LRN-20261003-015]`）。

### 13.4 A 的结果：已收敛，且与吴恺量级吻合

**运行**：`DRIVER=disp DISP_MM=0.1472` → `rc=0` / `NORMAL TERMINATION` / `converged at time : 1`（17 s）
**支反力**（沿加载轴求和）= **440.7 N**

| 域 | von Mises 中位 | p99 | max (MPa) |
|---|---|---|---|
| 距骨软骨 | 0.223 | **4.42** | 24.26 |
| 胫骨软骨 | 0.044 | **2.48** | 23.57 |
| 胫骨皮质骨 | 1.490 | 7.79 | 18.01 |
| 距骨皮质骨 | 6.270 | 31.89 | **38.68** |
| 胫骨/距骨松质 | 0.055 / 0.138 | 0.147 / 0.335 | 0.274 / 0.396 |
| 骨髓 | 0.0017 | 0.0062 | 0.0075 |

**与吴恺 2012 对照**（600 N 轴向，软骨/骨接触区 von Mises ≈ 3.74 MPa）：
440.7 N 下胫侧软骨 p99 = 2.48 MPa → 线性外推到 600 N ≈ **3.4 MPa**，
与吴恺的 **3.74 MPa 高度接近**。✅ 管线打通且量级可信。

**收敛包络**：稳定收敛区间 0.10–0.15 mm（407–441 N）；上限约 0.20–0.25 mm。
0.30 mm 细化到 60 载荷步可推到 96.5% 才停。成因：软骨大变形（≈18% 应变）
+ 部分近退化软骨薄片（det(J) 小至 0.008）+ penalty 接触刚度。

### 13.5 A 的五个硬骨头（全部已解，别再踩）

1. **`.feb` 语法**：`SolidDomain.mat` 填**材料名**；`<Surface>` 子标签**按单元类型命名**；
   段序必须是 `Mesh→MeshDomains→LoadData→Loads→Boundary→Contact→Output→Control`
   （顺序错会被报成 `unrecognized tag` 这种误导性错误）。**完整模板见 `FE_PIPELINE_TODO.md`**
2. **胫骨只靠接触约束 ⇒ 刚体模态 ⇒ 刚度矩阵奇异 ⇒ 2200 个单元翻转**
   → 改「顶面三自由度全约束（夹持端）+ 指定位移」驱动，**不再依赖接触来稳定**
3. **tet4 的翻转必须是奇置换**：4 个节点整体倒序 = 3 次交换 = 偶置换，**行列式符号不变**！
   → 122/462 + 183/618 个软骨单元被误剔除。正解：**交换两个节点**
4. **接触面必须与网格真实的面重合**：软骨已拆成 tet4，接触面若仍用原始 offset 四边形
   （三角化方式不一致）⇒ FEBio 接触搜索关联不上 ⇒ **完全不检测**
   （压缩 0.3 mm 仍零应力、两层软骨自由互穿）。正解：取「三节点全为 offset 节点」的 tet 面
5. **初始过盈会让第一步就爆**：让软骨过盈以"激活接触"→ t=0 就有接触力 → 压翻 4–5 个单元
   → 回到 FILL_FRAC=0.5，改用方案 2 解决刚体模态

**⚠️ 二分法的陷阱**：用「**无载荷**」变体做二分得到的是**假失败** ——
残差本已 8e-14（机器精度），但 FEBio 要求降到 `0.000000e+00`，永远满足不了
→ `Problem is diverging` → `Max nr of reformations reached`。
**二分时必须给非零载荷。**

**权威信息源**（比试错快 10 倍）：`D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf`
（4.7 MB，pdfplumber 直接提取）+ 让 `pyfebio` 生成样例 XML 读回。

### 13.6 读结果的坑
`pyfebio.xplt.to_hdf5` 对自定义 Surface/SolidDomain 集合会 `KeyError`
（`parse_state` 查 `mesh_dict[...][set_id]`）→ 需打容错补丁（见 `temp/pyfebio_demo/read_xplt3.py`）。
**且不要输出 `contact traction`**：它写出的 `PLT_FACE_DATA` 会让解析器崩掉。

### 13.7 第二轮（2026-10-04）：一个负面结果 + 一个新根因

**① "提高软骨网格质量"—— 假设是错的 ❌**
质量闸门（`CART_QUALITY=1`）把四面体 `|V| min` 从 0.0013 提到 **0.053 mm³（42 倍）**、
穿透仍 0，**但 FEBio 收敛更差**：0.30 mm/60 步从 0.965 掉到 0.642，工作点 17 s → 60 s。
**⇒ 收敛瓶颈是接触，不是薄片单元**。逐面过滤会在接触面上挖洞（见 ②）。
教训 `LRN-20261003-020`：**"看起来最像的原因"不等于瓶颈**；做优化前先做同指标 A/B。

**② 接 subtalar 关节 —— 网格建好了，接触面不合法（未完成）❌**
`FULL_ANKLE=1`：35,095 单元 / 10,324 节点，骨单元 0 剔除，两个接触界面都注册成功，
但**第一个载荷步就 `2 negative jacobians`**（载荷≈0，细化载荷步无效）。

**根因（关键工具：`temp/pyfebio_demo/diag_facets.py`）**：接触面**不是合法 2-流形**。

| 接触面 | 三角面 | 边使用次数 {1:开放边界, 2:正常, >2:非流形} |
|---|---|---|
| `surf_subtalar_calcaneus` | 787 | **{1: 333, 2: 592, 3: 163, 4: 76, 5: 9, 6: 1}** |
| `surf_subtalar_talus` | 633 | {1: 292, 2: 452, 3: 149, 4: 60, 5: 2, 6: 1} |

合法三角化曲面**每条边必须恰好被 2 个面用**。开放边界数 ≈ 被丢弃的 66 个软骨单元挖出的洞。
FEBio 也警告 `18 invalid facets`（**但单关节能收敛的模型同样有 2/3 个 invalid facet
却能跑完 ⇒ 警告本身不致命，致命的是非流形的规模**）。
教训 `LRN-20261003-019`。
→ **下一步**：先构造完整偏移曲面（合法 2-流形）**再**筛单元，或对接触面做边界修补。

**顺带一个已落地的有效修法**（`cartilage.py` 的 `CAP_SAFETY=0.97`）：
两侧软骨厚度各按自己那侧的 `0.5×gap` 算，而偏移方向各自指向对侧骨的最近点，
几何不对称时一侧会戳进另一侧。修法是后生成的一侧**厚度再被对侧已生成 offset 面的
距离限制**（留 3% 间隙）。实测压薄 10/89 与 35/253 个节点，体积几乎不变、穿透仍 0。

**③ 新增开关**：`CART_QUALITY` / `CAP_SAFETY` / `FULL_ANKLE` / `LAUGON` /
`NODE_RELOC` / `TWO_PASS` / `PENALTY` / `TIME_STEPS`（见 `FE_PIPELINE_TODO.md` 表格）。
试过 `AUGLAG`+`node_reloc`：第一个步就压翻 4 个，比 PENALTY 更糟。

**④ 单关节模型回归（改动后仍成立）✅**
`DRIVER=disp DISP_MM=0.1472 TIME_STEPS=20` → rc=0 / converged / 16.6 s，
支反力 **436.8 N**（改前 440.7 N），软骨 p99 **2.58 / 4.27 MPa**（改前 2.48 / 4.42），
差异 <1% ⇒ 第二轮所有改动**没有破坏已验证结果**。

**⑤ 仍然可做**
- 换回 force 驱动（600 N 面力）：需先解决刚体模态（弱弹簧 / 初始 tied 接触）
- 加 8 组踝韧带（`81100101~81101401_shell` 已在 npz 里）→ 复现吴恺的
  **定性排序**（胫跟＞胫舟＞跟腓＞下胫腓前），比 3.74 MPa 绝对值更有验证价值

---

**维护者**：Spark ⚡
**最后更新**：2026-10-03（踝关节几何/网格阶段：3 部件网格化完成，软骨阻塞待改 offset 方案；
**§2.3** 网格质量复核（gamma 中位 0.79–0.84，网格**不是**瓶颈）；
**§11** Maharaj2021 选型评估（对 FE 几何更差、对多体更好；IITD OSIM 在原生 OpenSim 4.6 加载失败）；
**§12 THUMS V7.1 现成 FE 骨网格** —— **提取验证已完成并全部通过**：
26 部件 0 失败、单位 mm、体积 3 法互验 ≤0.15%、皮质/松质共形、**装配天然对齐 1–3 mm**；
唯一差异是跗骨比 IITD 大 1.1–1.45×（经 PCA 验证为真实个体差异））
