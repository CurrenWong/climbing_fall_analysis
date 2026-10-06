---
title: OpenSim + Python 复现 [Y25] 有限元方法 —— 面向室内抱石软垫受力
date: 2026-10-03
type: 仿真方案
status: 主方案（对齐+场景+非垂直 S1–S4 已完成；执行状态见 OpenSim_FE交接.md）
tags:
  - OpenSim
  - 有限元
  - FEBio
  - 抱石
  - 坠落
  - 骨骼应力
  - 多体动力学
aliases:
  - OpenSim FE 复现方案
  - 抱石软垫有限元
---

# OpenSim + Python 复现 [Y25] 有限元方法

> **目标**：用 **OpenSim（多体动力学）+ Python/FEBio（有限元）** 复现
> [Y25]《成人直立位双足坠地跨高度骨骼损伤机制的有限元分析》的**方法链**，
> 然后把它从「刚性地面 · 1–50 m」**移植到「室内抱石软垫 · ≤4.5 m」**，
> 回答"摔在软垫上时各部位骨骼承受多少应力、会不会骨折"。

> **本方案与现有文档的关系**
> - 依赖 [`方案更新_v2.md`](方案更新_v2.md) 的 P0（分部位阈值）与 P1（场景库）成果；
> - 它实现的是 v2 §2.4 引述 [Y25] 自己点名的**「多体动力学-有限元耦合」**，
>   即 v2 里 **P3（转动/多体）+ P2（踝）+ P4（惯性分支）** 的工程化落地路径；
> - 注意：v2 §六 写过「❌ 不做全身有限元」。本方案**只做载荷路径链**（9 个部位，
>   非 THUMS 全身重复），增量价值是 **耦合 + 软垫**，与该结论不冲突。

---

# 一、锁定的范围（用户已确认）

| 维度 | 决定 |
|---|---|
| 架构 | **混合 OpenSim→FEBio**：OpenSim 出动载荷，FEBio 出 von Mises 应力 |
| 解剖范围 | **载荷路径链**：跟骨/足 → 胫骨 → 腓骨 → 股骨 → 骨盆 → 腰/胸/颈椎 → 颅骨（=[Y25] 报告的 9 个部位） |
| 验证顺序 | **先刚性地面标定 [Y25]**（1–50 m 直立双足），**再换软垫**转向抱石 |
| FE 几何 | **网格化 OpenSim VTP 骨面**（gmsh 表面→实体），不解析 THUMS `.k` |

## 1.1 ⚠️ 必须点明的关键前提：OpenSim 不是 FEM

**[Y25] 是纯有限元**（THUMS AM50 + LS-DYNA，直接算 von Mises 应力）。
**OpenSim 是多刚体动力学引擎，本身算不出 von Mises 应力。**
所以"用 OpenSim 复现 [Y25]"必须显式拆成两步：

```
OpenSim  →  关节反力/运动学（载荷）  →  FEBio  →  von Mises 应力  →  [Y25] 后处理
（多体）                              （有限元）
```

这不是妥协，而是 [Y25] 结论里自己写明的下一步：
> "……还可以引入**多体动力学-有限元耦合方法**，不仅可以研究瞬时冲击后的人体响应，
> 还可以得到高坠的整体动态过程，如坠落高度、初始速度、初始姿势等。"

换句话说，我们要做的正是 [Y25] 想做但没做的**耦合**部分。

---

# 二、资源盘点

## 2.1 仓库里已有的（不用重新找）

| 资源 | 位置 | 用途 |
|---|---|---|
| OpenSim 全身模型 | `paper/FullBodyModel-latest.zip` → `FullBodyModel-4.0/Rajagopal2015.osim` + `Geometry/*.vtp` | 多体动力学 + 骨几何 |
| THUMS V7.1（[Y25] 原模型） | `paper/AM50_V71_Occupant.zip` / `AF50_...zip`（LS-DYNA `.k`） | **仅作几何/材料/阈值参考**（无 LS-DYNA 许可，不求解） |
| IITD 踝足 OpenSim 模型 | `model/Model/*.OSIM` + STL/OBJ | 踝部高保真备用 |
| 有限元工具链 | gmsh 4.15.2 SDK、`meshio`、`pyvista`、`vtk`、`pyfebio` 0.3.0 | 网格 / 前处理 |
| FEBio 求解器 | `D:\Program\FEBioStudio\bin\febio4.exe`（4.13.0，**不在 PATH**） | FE 求解 |
| STL→tet 流水线 | `scripts/ankle_fe/geom.py`（体积误差 **0.00%**） | 骨网格化（复用） |
| 1D 冲击模型 | `src/climbing/pad.py`（双质点 + 非线性软垫 + 渐进接触） | **产生地面反力时程（GRF）** |
| 分部位阈值 | `src/climbing/bone.py`（已对标 [Y25]，命中 3/7） | 骨折判据 |
| 场景库 | `src/climbing/scenarios.py`（[B25] 真实频次） | 抱石加权 |
| 文献正文 | `paper/_text/*.txt`（[Y25]、踝 [T22]、跟骨 HIPC、[B25]…） | 阈值/对照 |

## 2.2 尚缺的

| 缺口 | 解决 |
|---|---|
| **OpenSim Python 未安装** | `pip install opensim`（见 §四） |
| 胫/腓/股/骨盆/脊柱/颅骨的 **FE 网格** | 用 `geom.py` 的 VTP→tet 流程批量生成 |
| **关节载荷 → FE 边界条件** 的映射代码 | 新写（§五 S3） |
| **[Y25] 后处理链**（骨折矩阵/Logistic/聚类） | 新写（§五 S4），可借鉴论文公式 |

---

# 三、总体架构

```
┌─ 输入 ────────────────────────────────────────────────────────────────────┐
│  高度 h（刚性 1–50 m / 软垫 0.5–4.5 m）· 质量 m · 姿势（直立双足）· 落点  │
└──────────────────────────────────────────────────────────────────────────┘
                                    │
              ┌─────────────────────┴─────────────────────┐
              ▼ ① 地面反力（用现有 1D 模型，不重造轮子）      │
┌─ src/climbing/pad.py ─────────────────────────────────┐ │
│  刚性地面  hard_surface()  → 短脉冲 GRF(t)             │ │
│  软垫  CrashPad(thickness=0.20) → 非线性长脉冲 GRF(t)   │ │
│  自检: ∫GRF dt = m·√(2gh)（冲量守恒，见 §七）          │ │
└───────────────────────────────────────────────────────┘ │
              │                                            │
              ▼ ② 多体动力学（OpenSim）                     │
┌─ OpenSim Rajagopal2015/2016 ──────────────────────────┐ │
│  站立姿态 + 骨盆 v0 = −√(2gh)                          │ │
│  GRF(t) 作为 ExternalForce 施加于 calcn/toes            │ │
│  Forward Dynamics 积分冲击窗口                          │ │
│  JointReaction → 踝/膝/髋/腰/胸/颈 关节反力（含力矩）    │ │
│  BodyKinematics → 各段加速度（惯性项）                  │ │
│  肌肉：基准用 "dead-drop"（见 §3.3）                    │ │
└───────────────────────────────────────────────────────┘ │
              │                                            │
              ▼ ③ 载荷映射（新模块）                        │
┌─ load_transfer.py ────────────────────────────────────┐ │
│  关节 wrench(F,M) → 关节面压力 p=F/A_cart + 弯矩力偶    │ │
│  并联构件载荷分配：胫骨 90% / 腓骨 10%（CSA 比）        │ │
│  坐标/符号统一（压缩为正），见 §七不变量                 │ │
└───────────────────────────────────────────────────────┘ │
              │                                            │
              ▼ ④ 有限元（FEBio）                          │
┌─ per-bone FE ─────────────────────────────────────────┐ │
│  gmsh: VTP 骨面 → 实体 tet 网格（复用 geom.py）         │ │
│  pyfebio/XML: 写 .feb（皮质/松质骨 E,ν）                │ │
│  febio4.exe 求解（static 峰值 / transient 波形）        │ │
│  .xplt → HDF5 → von Mises（Voigt 6 分量）               │ │
└───────────────────────────────────────────────────────┘ │
              │                                            │
              ▼ ⑤ [Y25] 后处理（新模块）                    │
┌─ y25_analysis.py ─────────────────────────────────────┐ │
│  各部位峰值 σ_vm vs 强度阈值 → 50×9 二值骨折矩阵         │
│  Logistic 回归（腓骨两端/骨盆/颅骨）                     │
│  层次聚类（Jaccard + 平均联结）→ 5 群集 / 3 阶段         │
│  双路径 PCA（[Y25] 的传导机制）                          │
└───────────────────────────────────────────────────────┘
              │
              ▼ ⑥ 对标验证 → ⑦ 换软垫转向抱石
```

---

# 四、关键技术决策

## 3.1 D1 · 载荷传递用 Forward Dynamics，不用 Inverse Dynamics

| 方案 | 判定 |
|---|---|
| **A. 1D 脉冲 → OpenSim 逆动力学** | ❌ 逆动力学需要**实测运动学**（`.mot`）；我们没有真实坠落动作捕捉，凭空给轨迹=伪造物理 |
| **B. OpenSim 正动力学 + ExternalForce（pad.py 提供 GRF）** | ✅ **采用**。初速度 `v0=−√(2gh)`，GRF(t) 作为外力，正向积分；关节反力由约束自然传导 |
| C. 规定运动学 + 分配惯性载荷 | ⚠️ 备选。会切断"高度→软垫→关节载荷"的因果链 |

**为什么 B 对**：载荷沿链向上的传导（正是 [Y25] 双路径的成因）由关节约束在正动力学里
自然产生；`pad.py` 提供**唯一诚实的接触物理边界**——OpenSim 不自己编接触本构。

> **载荷来源**：`JointReaction` 的**关节反力 wrench**（力 + 力矩）。
> 段惯性项 `m·a` 已包含在关节反力里（每个刚体满足自由体平衡）。
> **不用** `StaticOptimization`：肌肉反射潜伏期 10–50 ms ≫ 冲击时长 ~20 ms，SO 出的肌肉力不物理。

## 3.2 D2 · 关节载荷 → 单骨 FE 边界条件

一个关节反力**不能**直接整体压到一根骨头上（这是本项目 §3.14 已犯过的
"总力喂给单根构件"错误）。映射规则：

```
关节 wrench (F, M)  →  ① 关节面分布压力  p = F_n / A_cartilage
                       ② 弯矩力偶  ±F_c = M_bend / d（跨骨干质心两侧）
                       ③ 并联构件按 CSA 分配：胫骨 0.90 / 腓骨 0.10
```

- **关节面压力**优于刚性约束：两端都刚性钳制会**抑制弯曲模态、虚高骨干应力**。
- **力偶**表达弯矩，才能区分股骨颈 vs 股骨骨干的应力差异。
- **CSA 载荷分配**是仓库既有约定（`bone.py` 的 `load_fraction`）的延伸，
  防止"全力给胫骨、腓骨为零"导致腓骨两端 OR 排序失真。

**FEBio 里怎么把 wrench 施加得物理正确（已核实 pyfebio 0.3.0）**：

FEBio **没有**"节点力偶 / 面力矩"原语。惯用做法是**"幽灵刚体"**：

```
① add_simple_rigid_body(origin=关节中心)       # 在关节中心建一个参考刚体
② BCRigid(node_set="关节面节点", rb="ghost")    # 把骨关节面节点绑到该刚体
③ RigidForceLoad / RigidMomentLoad             # 施加 (F, M)，可带 loadcurve 时程
```

刚体求解器会把 F 与 M **正确分配到骨网格**，避免单节点加载的应力集中。
（备选：关节面分布压力 `p=F_n/A` + 弯矩力偶——防两端过约束，但力偶手搓易错。）
共享关节处仍按 CSA 分胫/腓。

**FEBio 分析类型（已核实）**：`Control analysis="DYNAMIC"`（隐式动力学）是
~ms 冲击的首选；收敛困难时退到显式 `ExplicitSolver`。载荷时程用
`load_controller type="loadcurve"`（`LINEAR`/`SMOOTH`，`dtmax` 要匹配最快特征 ~µs 级）。
骨材料用 `IsotropicElastic`（mm-N-MPa 下 E=7300–17000 MPa, ν=0.3）；
软垫用 `MooneyRivlinUC`/`OgdenUC`（pyfebio 无 hyperfoam，见附 A）。

## 3.3 D3 · 肌肉处理

| 方案 | 采用 |
|---|---|
| **dead-drop**：关闭激活动力学、激活=0（纯被动） | ✅ **基准**。冲击 ~20 ms 内肌肉来不及反应 |
| equilibrate + 默认肌肉 | ⚠️ 灵敏度对照 |
| Static Optimization / CMC | ❌ 无预定轨迹时无意义 |

## 3.4 D4 · 刚性 → 软垫只换 GRF

先验证阶段（S1–S4）用**刚性地面**：`ExternalForce` 用短脉冲。
抱石阶段（S5）**只把 GRF 换成 `pad.py` 软垫脉冲**，OpenSim/FEBio 管线不动。
理由：1D 软垫本构已是非线性泡沫物理的正确载体，脚-垫 FE 接触留到必要时再加。

> **[Y25] 的"阶跃 vs 渐进"对架构的含义**：
> - **阶跃部位**（足、胫骨两端、股骨颈、脊柱、颈椎）→ 应力由**外载荷**主导 → 单骨 FE 足够；
> - **渐进部位**（腓骨两端、骨盆、颅骨）→ 应力由**传导动力学**主导 → **必须有多体**（OpenSim）参与。
> 双路径本身也是**耦合**属性，不是单纯 FE 能给的。

## 3.5 D5 · 验证只对标"相对模式"，不对标"绝对高度"

绝对骨折高度（7–9 m）依赖 THUMS 特有的网格/材料/沙漏控制，**不可移植**。
可移植的是**相对结构**：聚类边界、OR 排序、阶跃部位集合、双路径。

## 3.6 D6 · OpenSim 版本与安装

- 项目 venv 是 **Python 3.12** → 支持 OpenSim **4.5.2（conda）/ 4.6（conda+pip）**。
- **首选**：`pip install opensim`（4.6，PyPI 有 cp312 wheel）。
- Windows 坑：wheel 的 DLL 注入在部分构建失效，脚本开头加
  `os.add_dll_directory(<site-packages>/opensim)`。
- **不要 pip 与 conda 混装 OpenSim**（Simbody/CasADi DLL 冲突）。
- Rajagopal2015 在 4.5+ 被改名 `Rajagopal2016.osim`；zip 里是 2015 版，4.6 仍可读。

---

# 五、分阶段计划（每阶段带 go/no-go 判据）

> 原则：**先端到端打通"一根骨 + 一个高度"的垂直切片**，再横向铺开。
> 切片不通过，加再多骨都没有意义。

## S0 · 环境（✅ 已完成）

- [x] S0.1 OpenSim 装到 `.venv` —— **`opensim 4.6-2026-06-22-85aaf64`**，
      `import opensim` + `osim.Model()` 正常。
- [x] S0.2 FEBio 定位与安全调用 —— 新增 `src/climbing/coupling/febio_run.py`
      （`find_febio()` / `probe_febio()` / `run_febio()` + `FebioNotFound`/`FebioRunError`）：
      探测到 `D:\Program\FEBioStudio\bin\febio4.exe` = **4.13.0**；
      **坏输入抛错而非静默 rc=1**（正是 pyfebio `run_model()` 的反例）。
      测试 `tests/test_febio_run.py` **8/8 通过**。
- [x] S0.3 模型已解压到 `model/opensim/FullBodyModel-4.0/` —— `Rajagopal2015.osim`
      加载成功：**22 bodies / 39 coords / 80 muscles / 75.34 kg**
      （模型名 `FullBodyModel_MuscleActuatedLowerLimb_TorqueActuatedUpperBody`）。
- **判据** ✅：`osim.Model(...)` 成功并打印维度；FEBio 冒烟 `rc=0`。

## S1 · 垂直切片：跟骨 @ 5 m 刚性（≈1–2 天）

最小闭环，验证**整条管线**能跑通：

- [x] S1.1 `pad.py` 刚性地面 → GRF(t)。**实施时修订**：`hard_surface()` 是
      **纯弹性**的，5 m 落地会把身体以 ~7–10 m/s 弹回（整窗 `∫F dt ≈ 2·m·v0`），
      而 [Y25] 是"停住"。故 GRF 截断在**全身动量首次归零**处（= 减速到静止），
      残余力线性降到 0。实现：`src/climbing/coupling/opensim_grf.py`。
- [x] S1.2 OpenSim：站立 + `pelvis_ty` 速度 `−√(2gh)`；GRF 作为时变外力。
      **实施时修订**：4.6 绑定的 `ExternalForce` **没有 `setForceFunction`**，
      改用 `PrescribedForce`（力/作用点写成时间的函数，全局系）。
      实现：`src/climbing/coupling/opensim_fall.py`。
- [x] S1.3 Forward Dynamics → 距下关节峰值 wrench。**实施时修订**：
      4.6 的 `JointReaction` 经 `Manager` 积分**不回填 `getStorageList()`**
      （实测 size=0），改用**牛顿自由体**（对 calcn+toes 远端子系统）直接算，
      天然可挂 T2。实现：`src/climbing/coupling/joint_loads.py`。
- [x] S1.4 跟骨几何 → tet → `.feb`。**实施时修订**：OpenSim 的 `r_foot.vtp`
      是 **1000 点、非流形、非封闭**的低模面（3 条开边、Euler=−330，trimesh
      判非水密），gmsh `classifySurfaces/createGeometry` 直接失败（几何自交）→
      回退**体素等值面 + scipy Delaunay 修复**。**已改用解剖级**
      `model/Model/Geometry/calcaneus_r.stl`（6718 面；trimesh 焊接后**水密**、
      体积 **57.7 cm³**，是真正的跟骨；而 `r_foot.vtp` 体积 121.7 cm³ 实为"整足"）。
      gmsh 对该 STL 仍卡在参数化（`Invalid exterior boundary mesh for
      parametrization`），但把 **classifySurfaces 角度设为 0.001 rad 可成体**
      （会保留全部面片 → 709k tet，过密；S2 应先抽面/decimate 再网格化）。S1 暂用
      体素修复网格。实现：`meshing.py`、`febio_model.py`。
- [x] S1.5 FEBio 求解 → von Mises。**实测**（`--height 5`，载荷 26.4 kN）：
      **FEBio 正常终止**。载荷改用**面压 p=F/A = 18.3 MPa**（非幽灵刚体）。实现：
      `fe_post.py`（新增 `von_mises_stats`）、入口 `meshing_driver.py` /
      `s1_vertical_slice.py`。

- **判据实测**：
  - **峰值 σ_vm 是边界奇异、不收敛**：同模型 4→3→2 mm 网格峰值
    **126.6 → 179.9 → 222.8 MPa**（单调增，h→0 无界）——面压只加到关节面边缘、
    跖面又全约束，棱边应力奇异。**峰值绝不可当强度判据**（这正是 [Y25] 只报
    相对模式的原因）。
  - **稳健判据（体分布）收敛**：p95 = **14.97**、p99 = 30.27、mean = 3.80 MPa
    （4/3/2 mm 下 p95 = 14.6 / 15.0 / 11.0，有界，O(名义压 18 MPa)）。
  - **结论：p95 = 15.0 MPa ≪ 150 MPa 阈值 → 5 m 刚性地面不骨折 ✅**，与预期一致；
    先前"576 MPa → 骨折"是伪影误判。
  - FEBio 收敛 ✅；T1 冲量 **+5.8% ✅**；T2 残差 **0.00 N ✅**。

> **S1 结论：全链路打通 ✅，稳健判据给出合理结论 ✅**（p95 = 15 MPa，不骨折）。
> 两个遗留（不阻塞 S1，交 S2）：(1) **网格质量**——体素修复网格 17.6% 单元质心
> Jacobian<0.3，gmsh 路线需先抽面再网格化；(2) **应力奇异**——后续统一以 p95/p99
> 或关节面下 ~5 mm 体积平均作强度判据，不用裸峰值。这正是方案 §六 D5
> "绝对量不可移植、只对相对模式"的体现。

> **⚠️ 实施中暴露的关键建模问题（比预期严重）**
> OpenSim 的关节是**理想无摩擦铰链、零被动刚度**。dead-drop 又关闭了肌肉，
> 于是一旦把 26 kN 的足底冲击施到脚上，小腿会被**直接甩飞** —— 实测踝关节
> 冲到 **−4.4 rad**（范围只有 [−0.7, 0.52]），膝关节 −1.55 rad，模型四肢乱飞。
> 这不是物理，是"零刚度关节 + 无肌肉"的**模型缺陷**。
>
> S1 的处理：把下肢关节（踝/膝/髋）**锁死**，近似 [Y25] 的"直立、腿绷直"
> （~20 ms 冲击内腿近似刚性柱体）。锁死后稳定：`pelvis_ty` 0.94→0.82 m，
> 距下关节纵向反力 **26.4 kN @ 20.1 ms**（右足施加峰值 26.5 kN）。
> **关节顺应性 / 肌肉（正是 [Y25] 双路径的成因）必须留到 S3** —— 届时要么给
> 下肢加真实被动刚度/Hill 肌肉，要么上真正的接触约束。

- **判据（全部满足才进 S2）**：
  - 峰值 σ_vm 非零且量级合理（跟骨皮质 ~100–150 MPa 量级）；
  - σ_vm < 跟骨强度阈值（5 m < 首次骨折高度，应"不骨折"）；
  - FEBio 收敛、网格 Jacobian > 0（不含负体积单元）；
  - 冲量不变量：`∫GRF dt ≈ m·√(2gh)`，相对误差 < 10%（§七 T1）。
    **实测 5 m = +5.8% ✅**（T2 自由落体残差 +0.00 N ✅）。

## S2 · 跟骨全高度扫描 + 过程区正则化（✅ 完成）

### S2.0 正则化（把发散的边界奇异峰值 → 有物理尺度的稳健应力）

问题：裸峰值是**边界奇异**（面压只加到关节面**边缘** + 跖面全约束），随网格
发散（4/3/2 mm → 127/180/223 MPa），不能作强度判据（LEARNINGS-023）。

**定义正则化长度 ℓ = 4 mm**（有物理依据，非调参）：松质骨 RVE ≈ 5 mm
（Harrigan 1988）；皮质骨断裂过程区 ≈ 1–5 mm。用它做两件事：

1. **网格单元尺寸取 ℓ** → 裸峰值 = "ℓ 尺度正则化峰值"（**操作判据**）；
2. **gauge 体积平均**（`fe_post.von_mises_gauge`，半径 R=ℓ、体积加权）→ 更强平滑的
   **下界**（交叉验证）。

> 实测：R=4–5 mm 时 cl=4 与 cl=3 两网格 gauge 一致（@5 m 约 25–26 MPa）；更细的
> cl=2 偏离（14 MPa）——因为 **voxel 修复网格不是有效加密族**，真正的收敛仍须先修
> 网格（decimate + gmsh，见 S3 前处理）。

### S2.1 / S2.2 扫描与结果

- [x] S2.1 扫 h = 1…50 m（刚性），每高度**真跑** OpenSim FD + 面压 FE。
      实现：`scripts/opensim_fe/s2_height_sweep.py`（网格算一次，逐高度仅换载荷）。
- [x] S2.2 σ_vm(h) 图 `scripts/opensim_fe/s2_plot.py`；产物
      `results/opensim_fe/s2_height_sweep_cl{3,4}.{json,csv}`、`s2_sigma_h.png`。

| 判据（ℓ=4 mm 网格，6 229 tet） | 首次 >150 MPa |
|---|---|
| **max（ℓ 正则化峰值，操作判据）** | **7.0 m ✅** |
| gauge(R=4 mm)（下界） | 未触发（≤50 m；@50 m 仅 122 MPa） |
| p99 / p95 | 未触发（≤50 m） |

- **结论**：
  - σ_vm(h) **单调递增 ✅**（T9 18/18 段非降）；**T1 全程 +5.6~+9.3% ✅**。
  - **ℓ 正则化峰值首次骨折 = 7.0 m**，落入 [Y25] 的 **7–9 m** ✅。
    ⚠️ **事后修正（S3.0）**：换成平滑 gmsh 网格后该数**不再成立**（gauge 26→49 MPa、
    raw 峰值高约 10 倍）——"7 m"是**粗 voxel 网格的正则化副作用**，不是对 [Y25] 的
    独立验证。可移植的是 σ_vm(h) 的**单调相对模式**，绝对高度不可移植（见 LRN-027）。
  - 裸峰值随网格漂移（cl=3 → 3 m），故**必须固定单元尺寸 = ℓ**（已写入
    `REG_LEN_MM`），绝对高度才可复现；gauge 给保守下界（≈1.5× 名义压）。
  - 与 E 无关；力控线弹性下 σ_vm 随**实测 F(h)** 走（**注意 F 快于 √h**，α≈0.66，
    非 0.5 —— 接触模型随冲击速度变硬）。

- **S2 判定：通过** —— 单调 ✅、ℓ 正则化峰值首次骨折 7 m ∈ [7,12] ✅。
  遗留（交 S3 前处理）：用 **decimate + gmsh** 得有效加密族，确认 gauge 在 ℓ 上
  mesh-independent；`REG_LEN_MM` 全局统一记录。

## S3.0 · 前处理：decimate + gmsh 平滑网格（✅ 完成，为 S3 必备）

**动机**：体素修复网格不是有效加密族（17.2% 单元质心 Jacobian<0.3），gauge 在
cl=2 上明显偏离（LEARNINGS-026），无法做网格收敛/正则化验证。

**做法**（`meshing.stl_to_tet_gmsh`，现为**主路径**）：
1. **抽面**：pyvista `decimate_pro(reduction=0.80)` → 6718 → **1342 面**（保持流形、
   体积 57161 mm³，较原始 −0.91%）；固定此面片 = **固定几何**。
2. **gmsh**：`classifySurfaces(0.001 rad)`（LEARNINGS-025）→ `generate(2)` 重网格化
   表面 → `generate(3)` Delaunay 四面体；`char_len` 只控制**内部**细化
   ⇒ 同一几何的有效 h-加密族。

**实测（@ h=5 m，载荷 26.4 kN）**：

| char_len | tet | 质心 J<0.3 | raw_max | **gauge(ℓ=4)** | gauge(R5) |
|---|---|---|---|---|---|
| 4.0 mm | 35 331 | 7.7 % | 1214.8 | **48.7** | 41.9 |
| 3.0 mm | 41 410 | 6.9 % | 1179.4 | **48.7** | 41.0 |
| 2.5 mm | 48 821 | 6.2 % | 1035.0 | **53.4** | 38.5 |

- **质量提升**：J<0.3 占比 **17.2 % → 5–8 %**；体积误差 **0.00 %**。
- **正则化成功** ✅：**gauge(ℓ) 在 2.5–4 mm 上收敛（48.7 / 48.7 / 53.4，~10 %）**，
  不再随网格漂移 —— 这正是 S3 需要的可信判据。
- **⚠️ 但绝对量变了**：平滑网格下应力场比体素网格高 2–3 倍（gauge 26→49 MPa），
  **raw 峰值高约 10 倍**。故 S2 里"首次骨折 7 m 恰好对上 [Y25] 7–9 m"应视为
  **粗网格的正则化副作用，而非独立验证**（见 LRN-027）。
>
> **gauge(ℓ) 判据下首次骨折 = 26.5 m**（**实测**，非外推）。力控线弹性使 σ ∝ F
> 严格成立（实测线性比 **0.99–1.00**）。5 m 时 gauge(ℓ)=48.7 MPa @ F=26.4 kN，
> 需 **F = 81.3 kN** 才到 150 MPa；而实测 F(h) **快于 √h**（5 m→26.4 kN、
> 50 m→121.4 kN，指数 α≈0.66 而非 0.5），故 F=81.3 kN 落在 **h≈26.5 m**
> （实测 25 m→144.8 MPa、30 m→162.5 MPa，内插）。
>
> 对照旧体素网格：gauge(ℓ)=26.3 MPa@5 m → 需 150 kN → **>50 m 永不触发**。
> 两套网格**都远高于** [Y25] 的 7–9 m（体素网格那个"7 m"来自裸峰值，非 gauge）。
> **结论：可移植的是 σ_vm(h) 的单调相对模式与部位间排序；绝对骨折高度不可移植**
> （与方案 §六 A5 一致）。

> **工程约束（重要）**：`pyfebio.xplt.to_hdf5` 把网格写成 HDF5 **attribute**，
> 超过约 **5 万 tet** 会 `OSError: object header message is too large`，读不回应力。
> 因此 FE 网格需 **≤ ~50 k tet**（S3 批量前处理必须遵守；必要时多站点拆批）。

## S3 · 载荷路径链全扫描（≈3–5 天，工作量最大）

- [ ] S3.1 批量网格化 9 部位骨 VTP → tet（复用 `meshing.stl_to_tet_gmsh`，
      固定 `reduction` 以保证各部位同 ℓ 尺度）。
- [ ] S3.2 实现 `load_transfer.py`：关节 wrench → 压力 + 力偶 + CSA 分配。
- [x] **S3.2b 跟腱力 —— ✅ 已实现并检验（结论：非主因）** —— 后足三点载荷
      （跖面↑ ＋ **跟腱↑** ＋ 距下关节↓）已补全。实现：`load_transfer.transfer` 给
      `F_ach = a·F_max(triceps)`（F_max = **10885 N**，`temp/opensim_fe/probe_achilles.py`
      实测）；`meshing` 选后上结节 `achilles` 面（x<−25, y>0）；`febio_model` 用
      `TractionLoad` 沿 ≈**(0.15, 0.99, 0)** 施加。
      **检验**（`scripts/opensim_fe/s3_achilles_sweep.py`，h=5/10/20，a=0 vs 1）：
      gauge_max **逐位相同**（26.31 / 42.58 / 68.39 MPa），仅体平均 +6~7%。
      诊断（`temp/opensim_fe/smoke_achilles_diag.py`）：后上区 σ_vm 随 a 显著升
      （a=0/1/5 → **0.2 / 20 / 101 MPa**），但**峰值恒在距下关节面棱边**
      （gauge argmax c≈(18,10,−2)），只由距下关节压力决定。
      **⇒ 跟腱在生理范围内（a≤1，≤10.9 kN）不改变跟骨峰值应力 → 它不是
      "26.5 m vs 7–9 m"差距的原因。** 真正主因指向下节。
- ⚠️ **S3.2c（新增，差距的真正主因）—— 距下关节载荷/BC 的边缘奇异** ——
      当前 BC = 关节面**均匀压力** + 跖面**全固定** → 峰值/ gauge 恒在关节面棱边，
      对网格发散（voxel gauge 26.3 vs gmsh gauge 48.7 MPa，差 1.9×）且对跟腱不敏感。
      下一步：用**接触/分布关节面载荷**替代均匀压力，或用**扩散部位（非棱边）**作判据。
- [ ] S3.3 关节链：距下/踝、膝、髋、腰骶、L4/L5、胸、颈、颅底。
- [ ] S3.4 9 部位 × 50 高度 → 峰值 σ_vm → **50×9 二值骨折矩阵**（对 `bone.py` 阈值）。
- **判据**：骨折矩阵中出现 **≥5 个阶跃部位**（完全分离）；
  全局自由体平衡 `ΣF_joint = m·a_com` 在每步成立（§七 T5）。
- **不通过**：无阶跃部位（阈值过低/过高）或全阶跃（载荷标定错）。

## S3.C · 方向 C：借 THUMS AM50 跟骨几何+材料（✅ 完成，未复现绝对高度）

**决策**：方向 C = **只借 THUMS AM50 跟骨的几何 + 材料**替换现有 IITD
`calcaneus_r.stl`；载荷仍走 OpenSim；修的是**几何保真 + 骨折判据**，
**不修动力学/BC**。

- **抽取**（`scripts/opensim_fe/kmesh_io.py` + `extract_thums_calcaneus.py`
  → `temp/opensim_fe/thums_calcaneus/`）：THUMS 甲板是 ASCII；右跟骨 = 两个独立 part：
  - `R_CALCANEUS_SPON` pid **81001200** —— 3323 单元（塌缩六面体记法的**四面体**），
    `*MAT_ELASTIC` RO=1e-9、E=**73.4** MPa、ν=**0.45**；
  - `R_CALCANEUS_CORT` pid **81001300** —— 634 **六面体**，
    `*MAT_PIECEWISE_LINEAR_PLASTICITY` E=**15000** MPa、ν=**0.3**、σy=**224** MPa。
  - 体积 **96.2 cm³**。**甲板内没有针对跟骨的 `*MAT_ADD_EROSION`**（仅有的侵蚀卡挂在
    `*MAT_SIMPLIFIED_RUBBER` 韧带且全为 0）；骨折 override 文件只改股/胫/腓。
- **配准（两版）**：
  - (a) `registration.py` ICP 到 OpenSim `r_foot.vtp` → **发现 `r_foot.vtp` 就是
    `calcn_r` 体挂的整只脚网格（193.6 mm），模型里没有独立跟骨几何**；ICP RMS 5.65 mm，
    但**原点落在骨后缘**。
  - (b) `registration_anatomical.py` 按解剖放置：**距下关节面质心→原点、跖面→−Y、
    前突→+X**；R 正交 det=1；于是载荷块 **222 mm²(14 面) → 2495 mm²(119 面)**。
    calcn 帧 local→ground 旋转 = **单位阵**（X前/Y上/Z内外、原点=距下中心）。
- **FE**（`scripts/opensim_fe/thums_feb.py`）：双材料域 cortical `"calcaneus"`（hex8
  E15000 / ν0.3）+ trabecular `"trabecular"`（tet4 E73.4 / ν0.45）；距下用
  `PressureLoad p=F/A`、跟腱 `TractionLoad`、跖面三向固定、`_fix_febio413`；
  FEBio 4.13 **rc=0**。判据 = 皮质 σ_vm ≥ **150** MPa（[Y25] 表）。
- **扫描**（`scripts/opensim_fe/s2_thums_sweep.py --mesh anat` →
  `results/opensim_fe/s2_thums_height_sweep_anat.{json,csv,png}`、`S2_THUMS_REPORT_ANAT.md`）：
  **38/38 rc=0**。**a=0 首次骨折 ≈2.0-3 m**（gauge：1m=114.9、2m=149.9、3m=156.6 MPa）；
  **a=1 时 gauge ≡ 裸 max**（跟腱 TractionLoad 边缘奇异主导）→ "1 m"非物理。
- **结论**：**方向 C 未复现 [Y25] 的足部骨折簇 7-9 m**（模型偏高 ~3-4×）；残余差距成因在
  **BC 层面**：跖面**全固定**固定端奇异、150 MPa 判据、跟腱 TractionLoad 边缘、多体
  F_subt 量级，均超出 C 的"不修动力学/BC"边界。

> **⚠️ caveat（重要）**：旧 S1-S3 的 `"calcaneus"` 域**实际是 `r_foot.vtp` 整足网格
> （121.7 cm³）**，不是单一跟骨；新模型 `"calcaneus"` 是真实 AM50 跟骨
> （CORT 18279 + SPON 77794 mm³）。**两者数值不可直接比较。**

## S4 · Logistic + 聚类对标 [Y25]（≈1 天）

- [ ] S4.1 对渐进部位做二元 Logistic 回归 → OR / 95% CI。
- [ ] S4.2 层次聚类（Jaccard 距离 + 平均联结）→ 树状图 / 群集边界。
- [ ] S4.3 双路径 PCA（各部位应力随高度协方差矩阵）。
- **验收判据（A1–A5，见 §六）**。

## S5 · 换软垫 → 抱石（≈1 天）

- [ ] S5.1 GRF 从刚性换成 `pad.py` 的 `CrashPad(0.20 m)` 软垫脉冲（**只换输入**）。
- [ ] S5.2 扫 0.5–4.5 m（[B25] 抱石量程），9 部位。
- [ ] S5.3 产出：**抱石软垫骨折矩阵 / 各部位受力分布 / 与刚性地面差异**。
- **判据**：≤4.5 m 落标准垫**基本无骨折**（这正是软垫的作用）；
  若 1–2 m 就骨折 → 软垫本构或载荷标定错。

---

# 六、验证判据（对表 [Y25]）

> **原则**：相对模式可移植，绝对高度不可移植。
> 与 v2 §五 V 表对齐：V1/V2/V4/V7 在本方案 S4 覆盖。

| 代号 | 判据 | 可移植？ | 目标 | Go |
|---|---|---|---|---|
| **A1** | 聚类边界 vs [Y25] 的 7 m / 16 m | ✅ | \|Δ\| < 2 m | 通过 |
| **A2** | 阶跃部位集合（足/胫骨两端/股骨颈/脊柱/颈椎） | ✅ | ≥4/5 命中 | 通过 |
| **A3** | 渐进部位 OR（腓骨 1.682 / 骨盆 1.236 / 颅骨 1.576） | ✅ | ≥2/3 在 ±30% 内 | 通过 |
| **A4** | 双路径结构（PCA 两个主成分：{胫/股/脊柱} vs {腓/骨盆/颅}） | ✅ | 可分 | 通过 |
| **A5** | 首次骨折高度（绝对） | ❌ | 放宽到 5–12 m | 参考（跖面 BC 放松后外推 **5.43 m**，落入区间内） |
| **A6** | 跖面 BC 敏感性诊断（fixed → 弹性）：gauge 降幅 | ✅ | ≥15% ⇒ BC 主导 | **通过（降 24%）** |

> **重要**：S4 通过**不代表** S5（软垫）自动通过 —— 换垫是**独立假设检验**。
> Logistic 遇完全分离会报 warning/拟出 ∞；必须**按 [Y25] 的规则剔除阶跃部位**再拟合。
> 现有 `bone.py` 命中 3/7，可作为阈值标定的锚点。

---

# 七、正交不变量（每个都对应一个静默失效）

> 沿用仓库铁律：物理 bug 不会自己暴露，必须靠不变量把关。

| # | 不变量 | 抓什么 |
|---|---|---|
| T1 | **冲量守恒**：`∫GRF dt = m·√(2gh)`，误差 <10% | pad.py 标定漂移 / 双边漏分 |
| T2 | **重力基线**：无冲击时 JointReaction ≈ 体重分配 | 符号约定漂移（压缩正） |
| T3 | **全局平衡**：`ΣF_joint = m·a_com` 每步成立 | 载荷映射漏项 |
| T4 | **载荷分配守恒**：共享关节 `Σ load_fraction = 1` | 胫/腓分配错 |
| T5 | **网格质量**：Jacobian > 0.3、长宽比 < 5 | gmsh 碎网格 → 伪应力 |
| T6 | **网格收敛**：加密一档 σ_vm 变化 < 20% | 网格未收敛 |
| T7 | **体积守恒**：FE 网格体积 vs 解剖参考（踝已 0.00%） | 单位/几何错 |
| T8 | **静-动一致性**：transient 峰值 ≤ 2× static 峰值 | 率效应失控 / 时间窗错 |
| T9 | **单调性**：σ_vm(h) 随高度单调 | 参数空间边界失效（类比 `ERR-...-001` t_max bug） |
| T10 | **单位自检**：.feb 显式单位、SI 一致 | mm/Pa 混用 |

---

# 八、风险与坑（已知）

| 风险 | 说明 | 缓解 |
|---|---|---|
| **OpenSim 不是 FEM** | 最根本的概念风险 | 架构显式拆分（§1.1） |
| Windows OpenSim DLL | PyPI wheel 缺 `opensim-cmd.exe`，DLL 注入失效 | 脚本加 `os.add_dll_directory` |
| VTP 仅支持 ASCII | 4.5+ 拒绝 binary/appended VTP | 用随包的 ASCII VTP；重导时选 ASCII |
| OpenSim 模型**未缩放** | Rajagopal 是 75 kg 通用模型 | 用 ScaleTool 或明确记录"未缩放" |
| JointReaction 在**关节原点**，非骨面 | 直接把 wrench 当骨面载荷会错位 | 在锚点处施加 follower force / 位移 |
| **BC 过约束** | 长骨两端刚性钳制虚高应力 | 用关节面压力，不用双端刚性 |
| **pyfebio 只能写不能读** | 0.3.0 往返 bug（`ERR-...-006`） | 单向：Python→.feb→febio4；读取走 xplt |
| **FEBio 不在 PATH** | `run_model()` 静默 rc=1 | 显式 `find_febio()` + 检查 rc |
| `.xplt` 布局反直觉 | 每 state 独立 dataset；位移是 (n,3) float32 | 按 `ERR-...-007` 的布局读 |
| **pad.py 输出是全身 GRF** | 施加到 OpenSim 需**双足拆分** | 按左右各半或解剖接触比拆分 |
| **冲击窗口** | 峰值只在"首触→首离"取 max | 复用 `pad.py` 的窗口逻辑，勿取全程 |
| **积分时长** | 类比 `t_max` 静默归零 bug | 高度扫描时加"未接触即抛错"守卫 |
| 阈值敏感性 | [Y25] 自认"阈值稍变结论就变" | 全部阈值写来源 + 灵敏度分析 |
| **跟骨=结构失效，非材料强度** | 用 `材料强度(135 MPa) × 截面` 高估 **17–34×**（`方案更新_v2` P0）：预测 **>40 m** vs [Y25] **7–9 m**（❌ 偏晚）。跟骨是足弓一环，由**非均匀小梁骨 + 足弓结构屈曲**主导 | 换 **Voo HIPC 接触力-骨折概率曲线**（Barnes, IRCOBI 2019, 全身 PMHS 验证），**不用** von Mises vs 135 MPa |
| **后足载荷缺跟腱力**（已检验，非主因） | S3.2b 已补跟腱牵引 + 参数化（`a·F_max`, F_max=10885 N）。**实测：生理范围 a≤1 内 gauge_max 逐位不变**；后上区局部 σ_vm 升高（a=1 达 20 MPa）但峰值恒在距下关节棱边 → **跟腱不是 26.5 m vs 7–9 m 差距的原因** | 保留跟腱项（a=0 基线 + a=1）；重心转向下条 |
| **跖面 BC 法向刚度 = 首折偏低的主因**（已量化 2026-10-04，原"判据被边界奇异支配"） | 跖面**三向全固定**的**法向（竖直）支承刚度**制造固定端奇异 → gauge 虚高。实验（`results/opensim_fe/PLANTAR_BC_REPORT.md`，方向 C anat @5 m，a=0）：fixed **186.8** → 三向 Winkler 弹簧 k≈1e3 **142.7 MPa（降 24%）**、峰值**离开跖面棱边**；首次骨折外推 **2.11 → 5.43 m**。**机制修正**：是**法向**固定，不是"三向 / 切向"（真放开切向反而 **+14%**）。关节面**均匀压力**角色见下条（③：±10% 内） | 换**弹性 / 接触**支承（`plantar_bc="spring"`，默认仍 `fixed`）。**即便最佳也只到 ~5.4 m，未达 7–9 m** → 属范式差（§3.5 D5） |
| **次级因素已排除**（③④，2026-10-04） | 距下**载荷表示**（均匀/梯度/节点）：gauge 变幅 **±10% 内**且**方向敏感**（+X 偏心 +2.7~9.9%，+Z −10.3%）；**STATIC vs DYNAMIC**：瞬态峰值仅 **+1.1%**（T8 通过）⇒ **STATIC 充分**；多体 **77 kg 校准**：F_subt +1.3% → 首折 2.02→1.99 m（**可忽略**）；**跖腱膜**：固定跖面下效应**精确为 0**（被反力吸收），弹性支承下 +0.28%@5 kN（**小**）。产物：`JOINT_LOAD_REPORT.md` / `MASS_CALIB_REPORT.md` / `FASCIA_REPORT.md` | 均**非**主因；主因仍是跖面**法向刚度**固定端奇异。峰值恒在跖面棱边（未因上述改变而迁移） |
| **绝对骨折高度不可移植** | 7–9 m 依赖 THUMS 特有网格/材料/沙漏；本项目实测：体素网格 → "7 m"（伪影）、平滑 gmsh 网格 gauge(ℓ) → **26.5 m** | 只对标**相对模式**（部位排序 / 聚类边界 / OR / 双路径）；绝对高度一律标"不可移植"（§3.5 D5、§六 A5） |
| **方向 C 未复现绝对高度** | 借 THUMS AM50 跟骨几何+材料后，首次骨折仍 ≈2-3 m；[Y25] 足部骨折簇 7-9 m **未复现**（模型偏高 ~3-4×）。残余差距在 **BC 层面**（跖面全固定奇异 / 150 MPa 判据 / 跟腱 TractionLoad 边缘 / 多体 F_subt 量级），超出方向 C 的"不修动力学/BC"边界 | 若必须追平绝对高度需**改 BC**（接触或分布关节面载荷、跖面非全固定、重标判据），属新方向，不在 C 范围 |
| **旧 S1-S3 `"calcaneus"` 是整足网格** | 旧域名为跟骨、实为 `r_foot.vtp` **整只脚**（121.7 cm³）；新模型才是真实 AM50 跟骨（CORT 18279 + SPON 77794 mm³）。**两者数值不可直接比较** | 引用旧结果必须标注"整足域"；跨模型比对先核几何身份（LRN-024 / LRN-027） |
| **gauge ≡ 裸 max = 边界奇异告警** | 方向 C 扫描 a=1 时 `gauge_max` 与裸 `max` 完全相等，峰值是**加载边界角点奇异**（4 mm 球内几乎只有单元本身），非骨内过程区峰值 | 报结果前识别并标注；见 LRN-029 |
| 停不下来 | 为好看调参 = 造假象 | 遵 `阶段总结.md §3.5`：有依据才标定，无依据先查结构 |

---

# 九、代码结构（建议）

```
src/climbing/coupling/            # 新增：OpenSim↔FE 耦合
├── __init__.py
├── opensim_fall.py     # 坠落 FD：装模型 / 初速度 / GRF 注入 / Forward Dynamics
├── joint_loads.py      # JointReaction + BodyKinematics 解析 → wrench 表
├── load_transfer.py    # wrench → FE BC（压力/力偶/CSA 分配/坐标符号）
├── meshing.py          # VTP → tet（包装 scripts/ankle_fe/geom.py）
├── febio_model.py      # pyfebio/XML 生成 .feb（材料/BC/步）
├── fe_post.py          # xplt → HDF5 → von Mises → 分部位峰值
└── y25_analysis.py     # 骨折矩阵 / Logistic / 聚类 / PCA

scripts/opensim_fe/                 # 分阶段入口（对应 S0–S5）
├── s0_env_check.py
├── s1_vertical_slice.py
├── s2_calcaneus_sweep.py
├── s3_chain_sweep.py
├── s4_validate_y25.py
└── s5_pad_bouldering.py

model/opensim/                      # 解压的 Rajagopal 模型 + Geometry/
temp/fe_meshes/                     # 各骨 tet 网格
results/opensim_fe/                 # 报告 / 矩阵 / 图
```

**测试**（`tests/test_coupling.py`）：把 §七 的 T1–T10 写成 pytest，
复用仓库"先写失败测试锁 bug"的纪律。

---

# 十、工作量与关键路径

| 阶段 | 工作量 | 依赖 |
|---|---|---|
| S0 环境 | 0.5 天 | — |
| S1 垂直切片 | 1–2 天 | S0 |
| S2 跟骨扫描 | 0.5 天 | S1 |
| S3 载荷链扫描 | 3–5 天 | S2（网格 + 映射是瓶颈） |
| S4 对标 [Y25] | 1 天 | S3 |
| S5 换软垫→抱石 | 1 天 | S4 |
| **合计** | **≈2 周（单人）** | 关键路径 = S3 网格与载荷映射 |

**关键路径**：S1 打通 → S3 的 **VTP→tet 批量网格化** 与 **wrench→BC 映射**
是风险最高的两步，建议 S1 就先跑通一根骨的全部环节。

---

# 十一、进度与下一步

**已完成**：

- **S0（环境）✅** —— OpenSim **4.6** 装入 `.venv`；FEBio **4.13.0** 已定位；
  Rajagopal2015.osim（22 bodies / 39 coords / 80 muscles / 75.34 kg）；
  `src/climbing/coupling/febio_run.py` 的 `find_febio()/probe_febio()/run_febio()`；
  `tests/test_febio_run.py` **8/8 通过**。
- **S1 多体半段（S1.1–S1.3）✅** —— `opensim_grf.py` / `opensim_fall.py` /
  `joint_loads.py`；`s1_vertical_slice.py --height 5`：T1 **+5.8%**、T2 残差
  **0.00 N**、距下关节 **26.4 kN @ 20.1 ms**；`tests/test_opensim_fe.py`
  **11 passed / 1 skipped**。
- **S1 全链路（S1.4–S1.5）✅** —— 解剖级跟骨 STL + 面压 `p=F/A` → FEBio 正常终止。
- **S2（跟骨全高度扫描 + 过程区正则化）✅** —— `s2_height_sweep.py` 扫 1…50 m：
  T9 单调 **18/18 ✅**、T1 全程 **+5.6~+9.3% ✅**；定义 **ℓ=4 mm**（RVE/过程区）
  作操作判据，`fe_post.von_mises_gauge` 提供体积加权下界。图
  `results/opensim_fe/s2_sigma_h.png`。
- **S3.0（前处理：decimate + gmsh 平滑网格）✅** —— `meshing.stl_to_tet_gmsh`
  （**主路径**）：抽面 `reduction=0.80`（1342 面，固定几何）→ gmsh(0.001 rad)
  重网格化。质心 J<0.3 占比 **17.2 % → 5–8 %**、体积误差 **0.00 %**；
  **gauge(ℓ) 在 cl=2.5–4 mm 上收敛（48.7 / 48.7 / 53.4 MPa）✅**。
  **⚠️ 结论修正**：平滑网格下 **首次骨折（gauge ℓ）= 26.5 m（实测）**，远晚于
  [Y25] 7–9 m；早先"7.0 m 合 [Y25]"是**粗体素网格的正则化伪影**（见 LRN-20261003-027）。
- **S3.C（方向 C：借 THUMS AM50 跟骨几何+材料）✅ 完成，未复现绝对高度** ——
  抽取 + 解剖配准 + 双材料 FE + 全高度扫描（**38/38 rc=0**）。新增脚本
  `scripts/opensim_fe/{kmesh_io.py, extract_thums_calcaneus.py, registration.py,
  register_calcaneus.py, registration_anatomical.py, thums_feb.py, thums_feb_run.py,
  s2_thums_sweep.py}`；产物 `temp/opensim_fe/thums_calcaneus/`、`results/opensim_fe/*anat*`。
  **a=0 首次骨折 ≈2.0-3 m**（gauge 1m=114.9 / 2m=149.9 / 3m=156.6 MPa）；**未复现
  [Y25] 足部骨折簇 7-9 m**（模型偏高 ~3-4×）。残余差距在 **BC 层面**（跖面全固定奇异 /
  150 MPa 判据 / 跟腱 TractionLoad 边缘 / 多体 F_subt 量级），超出 C 的"不修动力学/BC"边界。
  ⚠️ 旧 S1-S3 `"calcaneus"` 实为 `r_foot.vtp` 整足网格（121.7 cm³），与新 AM50 跟骨
  （96.2 cm³）**数值不可直接比较**。

- **S3.D（BC 诊断 + 次级因素排除）✅ 完成** —— 跖面 BC：fixed→三向 Winkler 弹簧(k≈1e3)
  gauge **186.8→142.7 MPa（降 24%）**、首折外推 **2.11→5.43 m**，机制=**法向刚度**（非切向）；
  ③ 距下表示 ±10%、**STATIC 充分**（DYNAMIC +1.1%）；④a 77 kg 校准 +1.3%（首折 2.02→1.99 m）；
  ④b 跖腱膜固定下 0。⇒ **残余差距主因唯一指向跖面法向刚度固定端奇异**。产物
  `results/opensim_fe/{PLANTAR_BC_REPORT,JOINT_LOAD_REPORT,MASS_CALIB_REPORT,FASCIA_REPORT}.md`。

- **坑①（gmsh 尺寸控制）✅ 解决** —— 根因 `classifySurfaces(默认 forReparametrization=False)`
  把输入面片钉成强制几何；修法 `classifySurfaces(ang,True,True)` + 关点/边界/曲率尺寸源。新 helper
  `scripts/opensim_fe/part_meshing.py`（mm-native，**未改 `meshing.py`**）。胫骨 cl=4 → **27670 tet /
  J<0.3=3.66%**（体素 17.4%）；跟骨 cl=3 → 17059/J2.90%。报告 `results/opensim_fe/GMSH_SIZE_REPORT.md`。

- **坑② 步①（壳侦察）✅ 完成** —— `scripts/opensim_fe/{shell_io.py,shell_scout.py}`：壳 **452920** 条；
  骨皮质厚度在 `*SECTION_SHELL`（非 `*ELEMENT_SHELL_THICKNESS`）；**壳↔实体节点 100% 共享**；
  骨盆 9 变体=**分区铺满**（须全取，只取 `_CORT_1.5` 仅盖 42%）；pyfebio 原生支持壳。报告 `SHELL_SCOUT_REPORT.md`。

- **坑② 步②（L3 shell+solid FE）✅ 完成** —— pyfebio `ShellDomain` 共节点**端到端跑通**
  （`l3_shell_feb.py`）；L3 SPON(tet4)+CORT(quad4,t=1.39) 100% 共节点、FEBio 4.13 rc=0；
  **去壳后松质 σ_vm 升 16–27×** ⇒ 纯松质严重失真、坑② 必须解决。报告 `L3_SHELL_FEB_REPORT.md`。

- **坑② 步③（泛化 FE builder）✅ 完成** —— `scripts/opensim_fe/bone_feb.py`：任意骨 1..n 实体域 + 0..n 壳域、
  按 PID 读材料、共节点；**4 骨 rc=0**（L3 / tibia[0壳] / parietal[2壳] / hip[9壳]）；发现 `*ELEMENT_SOLID`
  混 tet4+hex8（须按 distinct 节点数分类）。报告 `BONE_FEB_BUILDER_REPORT.md`。

- **各关节反力 ✅ 完成** —— `joint_reactions.py`：自由体法泛化到踝/膝/髋/脊柱；每关节 **T2≈1e-12**。@5 m：
  subtalar 25.4 / ankle 25.3 / knee 22.4 / hip 15.5 / **lumbar 22.2 kN**。修正 `_realized_states` **state 时间未写**
  的 bug（lumbar 旧 1.5 kN 作废；subtalar −4%）；刚性腿使近端 |F| 偏低（伪影）。报告 `JOINT_REACTIONS_REPORT.md`。

- **50×9 骨折矩阵骨架 ✅ 完成** —— `fracture_matrix.py`：9 代表骨 FE（`bone_feb`）+ 逐高度实测关节载荷 →
  9×50 矩阵（`fracture_matrix.{json,csv,png}` + `FRACTURE_MATRIX_REPORT.md`）。⚠️ 结果尚不可用：8/9 在 1 m 即骨折
  （**退化塌缩 tet4 点载荷使 σ 高 ~30×**，例股骨 15 kN→p95 1169 MPa）；3 骨高载发散须降载 5–25×。
  **下一步 = 修 FE 网格质量（`part_meshing`）后重跑矩阵。**

- **FE 网格质量修复 ✅ 完成** —— 4 骨 gmsh-reparam 干净网格；**证伪"退化 tet 是 σ 虚高主因"**：主因是 **BC**
  （斜向加载弯曲 + 端部夹持/点载荷奇异）。median 已物理；PCA 轴向后 p95 大降（跟骨 1.5×→32 m）。fibula 屈曲发散。
  报告 `MESH_FIX_REPORT.md`。**下一步 = 修 FE BC（轴向 + 端面分布力）。**

**下一步（🎯 清单⑤ = S3.1 + S4）**：

1. **S3.1 批量网格化 9 部位**（足/胫/腓/股/骨盆/腰/胸/颈/颅）—— 定几何来源（优先 THUMS）→
   `meshing.stl_to_tet_gmsh` 固定 `reduction`；单网格 **≤ ~5 万 tet**（HDF5 上限）。
2. **扩展载荷到 9 部位** —— `joint_loads` 现只出距下反力，需补踝/膝/髋/脊柱等各关节反力。
3. **泛化 FE builder** 到任意骨（现 `build_calcaneus_feb` / `thums_feb` 专用）。
4. **50×9 骨折矩阵**（1–50 m 扫描）→ **S4 Logistic + 层次聚类** → 对标 §六 **A1–A4**。
5. 补**下肢被动刚度/肌肉**（S1 锁腿是权宜，见 §五 S1 注）。
6. 对标**仅相对模式**（部位排序 / 聚类边界 / OR / 双路径）；绝对高度标"不可移植"。
7. **方向 C / S3.D 已封档** —— 绝对高度不追（D5）；跖面 BC 已诊断（法向刚度=主因）。

---

**参考**：
- [Y25] 袁红敏, 高树辉, 魏智彬. 中南大学学报(医学版). 2025;50(11):2051-2061.
- [T22] Li Z, et al. Forensic Sci Res. 2022;7(3):518-527.
- [B25] Beurienne E, et al. Front Sports Act Living. 2025;7:1609133.
- 本仓库：`docs/方案更新_v2.md`、`阶段总结.md`、`ANKLE_FE_TODO.md`、
  `PYFEBIO_TODO.md`、`.learnings/ERRORS.md`。

---

# 附 A · FEBio / pyfebio 落地要点（已核实，来源见下）

> 来源：FEBio User Manual 4.9（help.febio.org）、`febiosoftware/FEBio`、
> `febiosoftware/pyfebio` @ commit `e591eaf`（= 已装的 0.3.0）。
> 本节是 S1 直接可用的实现细节。

## A.1 pyfebio 0.3.0 能做什么

- **能**：Mesh / Material / Boundary / Loads / Rigid / Contact / Step / Control / Output 全部有类型封装。
- **能读 gmsh 网格**：`feb.mesh.translate_meshio(meshio.gmsh.read(".msh"))`——
  gmsh 的 Physical Group 自动变成具名 `Surface`/`ElementSet`（**施加载荷靠它**）。
- **MaterialParameter / Scale / Value** 用 `text="…"`（FEBio 是逗号分隔字符串，
  pydantic-xml 的已知限制）。

## A.2 载荷施加（对应 §3.2）

| 需求 | pyfebio 写法 |
|---|---|
| 关节 (F,M) | `add_simple_rigid_body` + `BCRigid` + `RigidForceLoad`/`RigidMomentLoad`（或 `RigidFollowerForceLoad` 带 insertion 点） |
| 面压力 | `loads.PressureLoad(surface=…, pressure=Scale(lc=…))` |
| 节点力 | `loads.NodalForce(node_set=…, value=Scale(lc=…, text="0,0,-150"))` |
| 指定位移 | `boundary.BCPrescribedDisplacement(node_set=…, dof="z", value=Value(lc=…))` |

## A.3 分析类型与载荷时程

- 冲击首选 `Control(analysis="DYNAMIC")` + `SolidSolver`（隐式）；不收敛再 `ExplicitSolver`。
- `TimeStepper(dtmax=…)` 要匹配载荷时程最快特征（~5 ms/100 ≈ 50 µs）。
- `OutputPlotfile(all_vars=[Var(type="displacement"), Var(type="stress")])` 才有应力输出。

## A.4 读结果 → von Mises

- `pyfebio.xplt.to_hdf5(xplt, hdf5)`；每个 state 是独立 dataset
  （沿用 `ERRORS.md [ERR-...-007]` 的布局）。
- 应力是 **per-element** 的 Voigt `[xx,yy,zz,xy,yz,xz]`：

```
σ_vM = sqrt( 0.5[(σxx−σyy)² + (σyy−σzz)² + (σzz−σxx)²] + 3(σxy²+σyz²+σxz²) )
```

- per-node 需自己把单元值平均到节点（meshio + 邻接表），再写 `.vtu` 供 ParaView。
- 按解剖区域筛选：gmsh Physical Surface → meshio `cell_sets` → pyfebio `ElementSet`。

## A.5 VTP 骨面 → 实体 tet（复用 `scripts/ankle_fe/geom.py`）

```
pyvista clean(tol) → fill_holes → trimesh 检查 is_watertight（Euler=2）
→ gmsh classifySurfaces(angle=30°, includeBoundary=True) → createGeometry()
→ addPhysicalGroup(3D=volume, 2D=具名关节面) → generate(3)
→ meshio.gmsh.read → translate_meshio
```

关键：`classifySurfaces` 重建拓扑（处理非流形边），`generate(3)` 限定体积做 Delaunay-tet。
关节面要在 gmsh 里加 `Physical Surface`，否则载荷没地方加。

## A.6 pyfebio 0.3.0 的硬限制（必须绕）

| 限制 | 后果 | 绕法 |
|---|---|---|
| `extra="forbid"` 只写不读 | 读回外部编辑的 .feb 会拒 | 只做单向；改动用 lxml/字符串 |
| **无 hyperfoam** | 软垫本构缺 | 存盘后字符串替换 `<material type="uncoupled foam">` |
| **无 rigid_wall contact** | 刚性地面接触缺 | 裸 XML 注入，或用 `SlidingElastic` + 刚体面 |
| `run_model()` 裸调 `febio4` | 不在 PATH → 静默 rc=1 | 显式 `subprocess.run([r"D:\Program\FEBioStudio\bin\febio4.exe","-i",…])` |
| `add_simple_rigid_body` 是 tiny tet | 只是参考点，非真刚体网格 | 正常，别当实体用 |

## A.7 单位

FEBio 与单位制无关。**MSK 骨 FE 惯例 = mm–N–MPa–s**：
应力 = MPa、密度 = tonne/mm³（1.7e-9）、力 = N。（也可 m–N–Pa，但**全程一致**。）
OpenSim 侧默认也是 mm，衔接时注意。**

> ⚠️ 单位不一致是本类项目最典型的静默失效（应力差 10⁶ 而不报错）——
> 见 §七 T10。
