# G4 载荷链顺应性杠杆报告 —— 顺应载荷链下 σ(h) 重算 + A4/A2 诚实判定

> 预注册设计 `docs/methods/G4_损伤演化设计.md`（§3 行4 / §5 / §6 D6 / §7）。冻结输入只读：`results/opensim_fe/fracture_matrix_bc2.json`。**零 FE 重跑**（σ(h) 后处理层）；OpenSim 载荷链在**唯一冻结预设** `PassiveStiffness.nominal_landing()`下逐高度重跑。本件**无墙钟字段**：固定平行分析种子 ⇒ 逐位可复现。

## 0. 一句话结论（如实）

- **A4 rank-2：FAIL** —— PC2 解释方差比（有效部位集）0.000242（rank-1）→ 0.008791（损伤-only）→ **0.009288**（损伤+顺应链）；平行分析 95% null = 0.1405、断棍 null = 0.2417 ⇒ PC2 **未超** null ⇒ rank-2 **FAIL**。
- **A4 分离（LP）：FAIL** —— `(PC1,PC2)` 载荷上 G1/G2不可线性分。
- **A2（股骨/胸椎）：部分改善（诚实报限）** —— 股骨 complete_separation（阶跃，首折 h=33.0）；胸椎 all_zero（常值0）。顺应链**直接抬高髋关节载荷**使**股骨获得阶跃区分度**；但**胸椎** σ(h) 取自冻结 median 轨迹且被损伤律 σ_y=80 < σ_c=150 **封顶** ⇒ 载荷链杠杆**单独不足以**修复胸椎（设计 §4 预注册的「与 FE 载荷分配缺陷同根因」对胸椎成立）。

> **§7 触发条件（PC2 超 null 但 LP 失败）：未成立**。本件**不实现**时间域度量（任务范围外）；如实记录触发状态。

## 1. 三路对照（rank-1 严格线弹性 vs 损伤-only vs 损伤+顺应载荷链）

| 量 | rank-1 基线（严格线性） | 损伤-only（最小路径） | **损伤+顺应链（本件）** |
|---|---:|---:|---:|
| PC2 解释方差比（有效集） | 0.000242 | 0.008791 | **0.009288** |
| PC1 解释方差比（有效集） | 0.999753 | 0.991180 | **0.990629** |
| 平行分析 95% null（PC2） | 0.2703 | 0.1381 | **0.1405** |
| 断棍 null（PC2） | 0.2032 | 0.2417 | **0.2417** |
| PC2 超 PA null？ | False | False | **False** |
| rank-2 成立？ | — | FAIL | **FAIL** |
| A4 分离（LP） | FAIL | FAIL | **FAIL** |
| A2 股骨 状态 | all_zero | all_zero | **complete_separation** |
| A2 胸椎 状态 | all_one | all_zero | **all_zero** |
| 有效部位数 | 9 | 6 | **6** |

> 平行分析 null = Horn 法（同维/同方差/iid 正态，5000 trials，seed=20261011）；断棍 = 断棍 null。**三路同规则、同种子**（转录/复用自最小路径）。

## 2. 顺应载荷链（唯一冻结预设，非调参）

- 预设 = `PassiveStiffness.nominal_landing()`；质量 = 75.337 kg；高度 = 1..50 m。
- 链路 = `ground_reaction(h, 75.337) -> run_dead_drop(grf, passive=预设) -> joint_reaction(fall, grf, joint=j)`，取 `peak_vertical_n`（纵向峰值，与 `fracture_matrix` 的 `load_ref` 同口径）。
- **只**用 `ankle_r`（tibia/fibula）/ `hip_r`（femur）；非长骨部位 σ(h) **逐位沿用**冻结 median 轨迹。

| 关节 | @1 m (kN) | @5 m (kN) | @50 m (kN) | 50/5 比 |
|---|---:|---:|---:|---:|
| `ankle_r` | 17.9 | 28.5 | 112.9 | 3.967 |
| `hip_r` | 14.7 | 25.7 | 115.0 | 4.466 |

> 与冻结（刚性腿）载荷链对照（`fracture_matrix_bc2.json` 的 `loads_n`）：顺应链改变各关节**峰值量级与沿高度形状**——这正是本杠杆的物理内容（`coupling/` **未改**，仅调用既有可选路径）。

## 3. 本构与参数（非调参）

| 部位 | 代表骨 | PID | E (MPa) | σ_y (MPa) | E_t (MPa) | E_t/E | 首屈服 h (m) | 分类 |
|---|---|---|---:|---:|---:|---:|---:|---|
| 足部 | `calcaneus_r` | 81001300 | 15000 | 224.00 | 1e-20 | 0.000 | 24 | saturation |
| 胫骨 | `tibia_r` | 81000700 | 18000 | 36.10 | 1.48e+04 | 0.822 | 14 | gradual |
| 腓骨 | `fibula_r` | 81000900 | 18500 | 37.00 | 1.67e+04 | 0.903 | 1 | gradual |
| 股骨 | `femur_r` | 81000100 | 17300 | 34.50 | 1.56e+04 | 0.902 | 9 | gradual |
| 骨盆 | `R_HIPBONE` | 83500205 | 17300 | 34.50 | 0 | 0.000 | 2 | saturation |
| 腰椎 | `L3` | 89001501 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 胸椎 | `T6` | 89000601 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 颈椎 | `C5` | 87000501 | 13020 | 80.00 | 0 | 0.000 | 1 | saturation |
| 颅骨 | `parietal_r` | 88000006 | 14900 | 95.88 | 7.59e+03 | 0.510 | 1 | gradual |

> σ_y/E_t 逐 PID 取自 THUMS 卡（`read_material_damage`，与非调参口径一致）。**唯一**旋钮 = 冻结被动预设，**不扫、不调**（§6 D6）。

## 4. A2（股骨/胸椎）—— 只报 + 三路状态

| 部位 | 弱区 σ_c | rank-1 状态 | 损伤-only 状态 | 损伤+顺应链 状态 | 顺应链首折 h |
|---|---:|---|---|---|---:|
| 足部 | 150 | complete_separation | complete_separation | **complete_separation** | 13.0 |
| 胫骨 | 70 | complete_separation | complete_separation | **complete_separation** | 44.0 |
| 腓骨 | 70 | all_one | all_one | **all_one** | 1.0 |
| 股骨 | 80 | all_zero | all_zero | **complete_separation** | 33.0 |
| 骨盆 | 180 | all_zero | all_zero | **all_zero** | — |
| 腰椎 | 150 | complete_separation | all_zero | **all_zero** | — |
| 胸椎 | 150 | all_one | all_zero | **all_zero** | — |
| 颈椎 | 150 | complete_separation | all_zero | **all_zero** | — |
| 颅骨 | 160 | complete_separation | complete_separation | **complete_separation** | 11.0 |

> 顺势链把**股骨**由「常值0」抬为「阶跃」（髋载荷随高度上移使股骨弱区 σ(h) 越过 σ_c=80）；**胸椎**仍为「常值0」——因其 σ(h) 取自冻结 median 轨迹且损伤律 σ_y=80 < σ_c=150 封顶。⇒ 载荷链杠杆**部分改善 A2（股骨✔ / 胸椎✘）**，但**未**修复 A4（rank-2 与 LP 双双 FAIL）。胸椎的修复**超出**本杠杆范围（设计 §4/§8 预注册）。

## 5. LP 载荷坐标（损伤+顺应链）

| 部位 | PC1 | PC2 | 组 |
|---|---:|---:|---|
| 足部 | -0.2246 | -0.9713 | 未标注 |
| 胫骨 | -0.0782 | -0.0104 | G1 |
| 腓骨 | -0.9146 | +0.1872 | G2 |
| 股骨 | -0.1148 | +0.0445 | G1 |
| 骨盆 | -0.0010 | -0.0086 | G2 |
| 腰椎 | +0.0000 | +0.0000 | G1 |
| 胸椎 | +0.0000 | +0.0000 | G1 |
| 颈椎 | +0.0000 | +0.0000 | G1 |
| 颅骨 | -0.3063 | +0.1393 | G2 |

## 6. 复现命令

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_load_chain_compliant.py
& .venv\Scripts\python.exe -m pytest tests\test_g4_load_chain.py -q
```

## 7. 产物

- `results/opensim_fe/load_chain_compliant.json`（新；顺应链逐关节逐高度 + meta）
- `results/opensim_fe/load_chain_compliant.npz`（新；设计 §3 命名的二进制缓存）
- `results/opensim_fe/g4_loadchain_matrix.json`（新；三路 σ(h) + A4/A2/rank2 判定）
- `results/opensim_fe/g4_loadchain_run.json`（新；判定摘要）
- `results/opensim_fe/G4_LOADCHAIN_REPORT.md`（本文件）
- `scripts/opensim_fe/g4_load_chain_compliant.py`、`tests/test_g4_load_chain.py`（新）

