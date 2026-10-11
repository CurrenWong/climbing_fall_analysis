# G4 损伤演化（最小路径）报告 —— 双线性弹塑 σ(h) 重算 + A4 诚实判定

> 预注册设计 `docs/methods/G4_损伤演化设计.md`（§3/§5/§6）。冻结输入只读：`results/opensim_fe/fracture_matrix_bc2.json`。**零 FE 重跑**（σ(h) 后处理层）。
> 本件**无墙钟字段**：固定平行分析种子 ⇒ 逐位可复现。

## 0. 一句话结论（如实）

- **A4 rank-2：FAIL** —— PC2 解释方差比 0.000242 → **0.008791**（×36.3），**数值上已脱离 rank-1 天花板**；但平行分析 95% null = 0.1381、断棍 null = 0.2417 ⇒ PC2 **未**超 null ⇒ rank-2 **不成立**。
- **A4 分离（LP）：FAIL** —— `(PC1,PC2)` 载荷上 G1/G2 不可线性分。
- **A2（股骨/胸椎）：未获区分度** —— 损伤律把脊柱 σ(h) 封顶在 σ_y=80 MPa < σ_c=150 ⇒ 腰椎/胸椎/颈椎由「阶跃/常值1」退化为「常值0」；股骨维持「常值0」（见 §4）。

> 结论：**双线性后处理律本身不足以**造出可分双路径结构（与 §4 预注册风险一致）。§7 时间域度量是文档化备选（本件未实现；其触发条件「PC2 超 null 但 LP 失败」也**未**被满足——PC2 连 null 都未超）。

## 1. 本构与参数（非调参）

`σ(ε)=E·ε`（ε≤σ_y/E）否则 `σ_y+E_t·(ε−σ_y/E)`；后处理等价映射 `σ=σ_lin`（σ_lin≤σ_y）否则 `σ_y+(E_t/E)(σ_lin−σ_y)`。

| 部位 | 代表骨 | 主承载域 | PID | E (MPa) | σ_y (MPa) | E_t (MPa) | E_t/E | 首屈服 h (m) | 分类 |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| 足部 | `calcaneus_r` | solid_calc_cort | 81001300 | 15000 | 224.00 | 1e-20 | 0.000 | 24 | saturation |
| 胫骨 | `tibia_r` | solid_tibia_cort | 81000700 | 18000 | 36.10 | 1.48e+04 | 0.822 | 14 | gradual |
| 腓骨 | `fibula_r` | solid_fib_cort | 81000900 | 18500 | 37.00 | 1.67e+04 | 0.903 | 1 | gradual |
| 股骨 | `femur_r` | solid_fem_cort | 81000100 | 17300 | 34.50 | 1.56e+04 | 0.902 | 22 | gradual |
| 骨盆 | `R_HIPBONE` | shell_hip_cort_83500205 | 83500205 | 17300 | 34.50 | 0 | 0.000 | 2 | saturation |
| 腰椎 | `L3` | shell_l3_cort | 89001501 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 胸椎 | `T6` | shell_t6_cort | 89000601 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 颈椎 | `C5` | shell_c5_cort | 87000501 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 颅骨 | `parietal_r` | shell_par_int | 88000006 | 14900 | 95.88 | 7.59e+03 | 0.510 | 1 | gradual |

> σ_y/E_t 逐 PID 取自 THUMS 卡（`read_material_damage`，卡 `nums[4]/nums[5]`），**非调参**。E_t=0 ⇒ 完美塑性平台（饱和端）；皮质骨 E_t/E≈0.82–0.90 ⇒ 膝拐（单调，非平台）。

## 2. A4 —— rank-2（数据驱动 null）

- A4 输入轨迹 = 与 `a4_from_matrix` **逐行同构**的 per-part σ(h)（长骨 ends；其余 median）。
- 有效变量 6/9（常数/饱和部位 ['腰椎', '胸椎', '颈椎'] 方差=0，从检验剔除 ⇒ null 门随之降低 0.2703→0.1381（对通过更有利、非保守）；即便如此 PC2 仍 < null ⇒ FAIL）：['足部', '胫骨', '腓骨', '股骨', '骨盆', '颅骨']。

| 量 | 基线（严格线弹性） | 损伤（双线性） |
|---|---:|---:|
| PC1 解释方差比 | 0.999753 | 0.991180 |
| **PC2 解释方差比** | 0.000242 | **0.008791** |
| 平行分析 95% null（PC2） | 0.2703 | 0.1381 |
| 断棍 null（PC2） | 0.2032 | 0.2417 |
| PC2 超 PA null？ | False | **False** |
| PC2 超断棍？ | False | **False** |

> 平行分析 null = Horn 法（同维/同方差/iid 正态，5000 trials，seed=20261011）；断棍 = 断棍 null（见上）。⇒ **rank-2 = FAIL**。

## 3. A4 —— 分离（LP，复用冻结规则）

- 基线（**直接调用**冻结 `s4_alignment_bc2.a4_from_matrix`）：PC2=0.000242，separable=False。
- 线性轨迹（转录规则，应逐值等于基线）：PC2=0.000242，separable=False。
- 损伤轨迹（转录规则）：PC1=0.991180、PC2=0.008791；LP **FAIL**（The problem is infeasible. (HiGHS Status 8: model_status is Infeasible; primal_status is None)）。

| 部位 | PC1 | PC2 | 组 |
|---|---:|---:|---|
| 足部 | -0.2234 | -0.9686 | 未标注 |
| 胫骨 | -0.0791 | -0.0134 | G1 |
| 腓骨 | -0.9212 | +0.1788 | G2 |
| 股骨 | -0.0552 | +0.0105 | G1 |
| 骨盆 | -0.0010 | -0.0085 | G2 |
| 腰椎 | +0.0000 | +0.0000 | G1 |
| 胸椎 | +0.0000 | +0.0000 | G1 |
| 颈椎 | +0.0000 | +0.0000 | G1 |
| 颅骨 | -0.3035 | +0.1717 | G2 |

## 4. A2（股骨/胸椎）—— 只报

| 部位 | 弱区 σ_c | 基线 状态 | 损伤 状态 | 基线首折 h | 损伤首折 h |
|---|---:|---|---|---:|---:|
| 足部 | 150 | complete_separation（阶跃） | complete_separation（阶跃） | 13.0 | 13.0 |
| 胫骨 | 70 | complete_separation（阶跃） | complete_separation（阶跃） | 38.0 | 45.0 |
| 腓骨 | 70 | all_one（常值1） | all_one（常值1） | 1.0 | 1.0 |
| 股骨 | 80 | all_zero（常值0） | all_zero（常值0） | — | — |
| 骨盆 | 180 | all_zero（常值0） | all_zero（常值0） | — | — |
| 腰椎 | 150 | complete_separation（阶跃） | all_zero（常值0） | 7.0 | — |
| 胸椎 | 150 | all_one（常值1） | all_zero（常值0） | 1.0 | — |
| 颈椎 | 150 | complete_separation（阶跃） | all_zero（常值0） | 7.0 | — |
| 颅骨 | 160 | complete_separation（阶跃） | complete_separation（阶跃） | 6.0 | 11.0 |

> 股骨/胸椎**未获区分度**；脊柱封顶 σ_y=80 < σ_c=150 ⇒ 全高度不判骨折。与 FE 载荷分配缺陷同根因（设计 §4/§5 预注册）。

## 5. 复现命令

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_damage_evolution.py
& .venv\Scripts\python.exe -m pytest tests\test_g4_damage.py -q
```

## 6. 产物

- `results/opensim_fe/g4_damage_matrix.json`（新；σ(h) 线性/损伤逐格 + A4/A2/rank2 判定）
- `results/opensim_fe/g4_damage_run.json`（新；判定摘要）
- `results/opensim_fe/G4_DAMAGE_REPORT.md`（本文件）
- `src/climbing/coupling/bone_damage.py`、`scripts/opensim_fe/g4_damage_evolution.py`（新）

> **载荷链顺应性杠杆（G4 第二杠杆，§3 行4 / §6 D6）已完成 → 见** `results/opensim_fe/G4_LOADCHAIN_REPORT.md`（驱动 `scripts/opensim_fe/g4_load_chain_compliant.py`；唯一冻结预设 `PassiveStiffness.nominal_landing()`）。三路对照：PC2 0.024%→0.88%(损伤)→**0.93%**(损伤+顺应链)，仍 < 平行分析 null 0.141 ⇒ A4 rank-2/LP **如实 FAIL**；A2 **股骨获阶跃、胸椎 all_zero**；§7 触发条件**未成立**。本件（最小路径）数值**未变**。

