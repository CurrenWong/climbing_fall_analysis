# S1 场景：垫子/地面刚度 → 首次骨折高度（抱石）

> 生成时间 2026-10-04 22:36:45；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s。
> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio，不重跑 OpenSim。
> 只读复用：冻结纯轴向基线 + 已有 FE 跖面 BC 实测（`plantar_bc_summary.json`）。

---

## 0. 一句话结论

- **垫子越软，首折越高**：FE 实测折减比最低点 R_min = **0.7642**（k = 1000 N/mm，即应力降 23.6%）。
- 刚性基线（R=1）：腓骨 7 m、足部 8 m、胫骨 8 m、股骨 19 m。
- 最好垫子（R_min，k≈1000 N/mm）：足部 8 → **11 m**、腓骨 7 → **10 m**、胫骨 8 → **11 m**、股骨 19 → **31 m**。
- **对比已有 FE 结果**：同一 k≈1e3，FE 单骨跟骨外推把首折从 2.11 m 抬到 5.43 m（×2.58）；冻结 1D（本脚本 M1）足部从 8 m 抬到 11 m。两者方向一致，但**绝对高度不同源**（见 §4）。

---

## 1. 地面/垫子当前如何建模（rigid vs compliant）

| 层 | 模块/文件 | 现状 | 刚性/柔性 |
|---|---|---|---|
| 1D 接触 | `src/climbing/pad.py` | `simulate_boulder_fall` 的渐进接触；`on_pad=False` 用 `hard_surface(k=1e7 N/m)` 线性弹簧（人体等效刚度，非真刚性墙） | **刚性等价**（基线） |
| GRF 生成 | `src/climbing/coupling/opensim_grf.py::ground_reaction` | `on_pad` 默认 False；把 1D 接触力翻译成左右脚的 `(t, F)` | **刚性地面** |
| 多体注射 | `src/climbing/coupling/opensim_fall.py::run_dead_drop` | `PrescribedForce` 把 GRF 施加到 `calcn`，腿关节锁死、肌肉关闭；**模型内没有地面接触单元** | 载荷边界（刚性） |
| FE 跖面 BC | `src/climbing/coupling/plantar_bc.py::apply_plantar_bc` | `fixed`/`roller`/`spring`(Winkler k)；`spring` 是可参数化的柔性支承 | **可柔性**（仅 FE 探针用） |
| 载荷链数据 | `temp/opensim_fe/fracture_matrix/load_chain.npz` | 50 高度的 `grf_peak_n` + 逐关节 `vert/abs` 峰值；全部由**刚性地面**链产生 | **刚性** |

**结论**：冻结基线（`risk_1d_loadshare` / `AXIAL_SUBREGION_REPORT` / `load_chain.npz`）是**刚性地面**口径，没有垫子刚度这一维。唯一已有的柔性地面证据是**独立的 FE 跖面 BC 探针**（`plantar_bc.py` + `plantar_bc_summary.json`），它给出 `k -> gauge 折减比`，正好可用作把刚性基线迁移到柔性地面的**实测转移函数**。

---

## 2. 垫子模型与公式

### 2.1 路线 M1（主口径）：FE 实测折减比转移

```
R_FE(k) = gauge_max(spring_k) / gauge_max(fixed)      # 来自 plantar_bc_summary.json，实测
sigma_pad(h; k) = R_FE(k) · sigma_rigid(h)            # sigma_rigid(h)=risk(h)·stored σ_c
首折 h*(k) = min h 使 sigma_pad(h;k) >= σ_c(弱子区域)
等价：risk(h) >= (σ_c_weak / σ_c_stored) / R_FE(k)
```

| k (N/mm) | gauge_max (MPa) | R_FE = gauge/fixed |
|---:|---:|---:|
| fixed | 186.76 | 1.0000 |
| 10 | 156.51 | 0.8380 |
| 100 | 151.90 | 0.8133 |
| 1000 | 142.73 | 0.7642 |
| 10000 | 173.44 | 0.9287 |
| 100000 | 184.40 | 0.9874 |
| 1e+06 | 186.44 | 0.9983 |

> R_FE(k) **非单调**：软端（k≤1e3）峰值移到内部/关节面边缘（跖面固定端奇异被消掉），最低 0.7642 @ k=1000；硬端峰值回到跖面边缘，k=1e6 残差相对 fixed 仅 0.17%（收敛回刚性）。

### 2.2 路线 M2（物理交叉校验）：串联线性弹簧

```
F_peak = v0 · sqrt(m · k_series)        # 能量 ½mv² = ½k_s x²，F = k_s x
k_series = k_pad · k_body / (k_pad + k_body)   # 垫子与人体组织串联
R_phys(k_pad) = F_pad / F_rigid = sqrt(k_pad / (k_pad + k_body))
k_body = 1.0e7 N/m  (pad.hard_surface 默认)
```

- 反推：FE 最佳折减比 R_min=0.7642 对应物理 **k_pad ≈ 1.4e+07 N/m**（≈ k_body 量级）；k=1000 N/mm 的 R=0.7642 对应 k_pad ≈ 1.4e+07 N/m。
- **选路说明**：主口径用 **M1（实测）**，因为它直接复用已有 FE 证据且与任务锚点（FE 跟骨 5.43 m @ k≈1e3）同源；M2 只做量级/方向校验（真实泡沫非线性，R_phys 会高估软垫收益）。

---

## 3. k → 首折高度曲线（M1 主口径，逐骨）

| 部位 | 弱子区域(σ_c) | 刚性 h* | 10 | 30 | 100 | 300 | 1000 | 3000 | 10000 | 30000 | 100000 | 1e+06 | ∞ (刚性) |
|---|---|---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---: |
| 足部 | calcaneus (150) | **8** | 10 | 10 | 10 | 11 | 11 | 10 | 9 | 8 | 8 | 8 | 8 |
| 腓骨 | fibula_ends (70) | **7** | 9 | 9 | 10 | 10 | 10 | 9 | 8 | 8 | 7 | 7 | 7 |
| 胫骨 | tibia_ends (70) | **8** | 10 | 10 | 11 | 11 | 11 | 10 | 9 | 9 | 8 | 8 | 8 |
| 股骨 | femoral_neck (80) | **19** | 27 | 27 | 28 | 29 | 31 | 26 | 22 | 21 | 20 | 19 | 19 |

> 表头 k 列（N/mm）：10、30、100、300、1000、3000、10000、30000、100000、1e+06、∞ (刚性)。
> 首折高度由离散的 1..50 m 扫描给出，故为整数米（分辨率 1 m）。

### 3.1 5 m 处的应力（MPa，随 k 变化）

| 部位 | 刚性 σ@5m | R_min σ@5m | 降幅 |
|---|---:|---:|---:|
| 足部 | 117.38 | 89.71 | 23.6% |
| 腓骨 | 57.59 | 44.01 | 23.6% |
| 胫骨 | 53.48 | 40.87 | 23.6% |
| 股骨 | 37.61 | 28.74 | 23.6% |

### 3.2 M2 物理曲线（k_pad N/m → 首折，逐骨）

| 部位 | 1e+05 | 3e+05 | 1e+06 | 3e+06 | 1e+07 | 3e+07 |
|---|:---:|:---:|:---:|:---:|:---:|:---:|
| 足部 | >50 | >50 | 46 | 22 | 12 | 9 |
| 腓骨 | >50 | >50 | 42 | 20 | 11 | 9 |
| 胫骨 | >50 | >50 | 48 | 23 | 13 | 10 |
| 股骨 | >50 | >50 | >50 | >50 | 35 | 25 |

---

## 4. 对比：刚性基线 / 软垫极限 / 已有 FE 结果

| 部位 | 刚性基线 h* (m) | 最好垫子 h* (m) @ R_min | Δ (m) |
|---|---:|---:|---:|
| 足部 | 8 | **11** | +3 |
| 腓骨 | 7 | **10** | +3 |
| 胫骨 | 8 | **11** | +3 |
| 股骨 | 19 | **31** | +12 |

**已有 FE 单骨跟骨结果（anat 尺度，实测外推）**：

| 量 | 值 |
|---|---:|
| FE fixed 首折高度 | 2.11 m |
| FE spring k=1000 N/mm 首折高度 | 5.43 m |
| 高度比 (k=1000 / fixed) | 2.58× |

> **尺度警告**：FE anat 尺度（fixed 2.11 m）与冻结 1D 尺度（足部 8 m）**绝对高度不同源**（FE 含弯曲/BC 伪影，1D 为纯轴向），二者**不可直接比较绝对高度**。可比较的是**变化方向与量级**：
> - FE：应力降 23.6% → 高度 ×2.58（2.11 → 5.43 m）；
> - 冻结 1D（本脚本 M1）：同一折减比使足部 8 → 11 m、腓骨 7 → 10 m。
> 两者方向一致（垫子抬高首折），1D 因 σ(h) 增长较缓，绝对抬高幅度更大。

---

## 5. 敏感性与不确定范围（避免单个误导数字）

1. **FE 转移非单调** ⇒ 应力折减比在 k∈[10, 1e6] 内落在 **[0.7642（k=1e3）, 1.0000（k=1e6 刚性极限）]**；对应足部首折的**有界范围** **[8 m（R=1）, 11 m（R_min）]**。不把“k=1e3 → 5.43 m”当作唯一答案。
2. **M1 折减比是跟骨局部实测**：外推到腓/胫/股骨为假设；若垫子主要消掉的是**跟骨固定端奇异**而非全局载荷，则其它骨的收益应更小（保守）。
3. **M2 物理路线**给的是方向/量级：R_phys(k_pad) 单调递减，k_pad∈[1e+05, 3e+07] N/m → R∈[0.100, 0.866]；真实泡沫非线性（平台+压实）会显著削弱软端收益，故 M2 的软端是**乐观上界**。
4. **名义应力 1D 局限**：σ 是纯轴向名义压缩，忽略弯曲/剪切/偏心/应力集中（对细长骨低估、对紧凑骨高估）；本报告只在**相对变化**上使用它。

---

## 6. 假设与诚实边界

1. 冻结基线（纯轴向 σ(h)=risk_by_height×stored σ_c）与子区域强度（足150/腓两端70/胫两端70/股骨颈80）逐位沿用 AXIAL_SUBREGION_REPORT.md；本脚本只加垫子这一维。
2. M1 假设垫子按同一折减比 R_FE(k) 折减整条载荷链峰值应力（1D 链线弹性 ⇒ 载荷与应力成比例）。R_FE 是**跟骨局部 gauge** 实测（单骨 THUMS AM50、h=5m、a=0、线弹性后处理、无单元删除），外推到腓/胫/股骨是**假设**，不是实测。
3. M1 的 R_FE(k) 非单调（软端峰值在内部/关节面边缘，硬端在跖面边缘；k=1e3 最低 0.764）；k>=1e6 视为刚性极限 R=1；k<10 平推软端 0.838。
4. M2 用线性弹簧；真实泡沫是平台+压实非线性，R_phys 只应视为量级/方向校验，若要定量应取相关压缩量下的**割线刚度**。
5. 绝对首折高度不可移植（[Y25] 范式）；本报告只给相对变化与排序。
6. 本脚本只读复用既有产物，不调用 FEBio，不重跑 OpenSim 正动力学。

- **测量 vs 建模**：M1 的 R_FE(k) 是**实测**（FE 探针）；σ(h) 映射与整链外推是**建模**；M2 全程**建模**。所有首折高度都是**线弹性比例外推**，不是重跑高度扫描。
- 本脚本不改变任何既有产物；`risk_1d_loadshare.json`、`plantar_bc_summary.json`、`fracture_matrix.json` 均为只读输入。

---

## 7. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\pad_stiffness.py
# 只读复用：
#   results/opensim_fe/risk_1d_loadshare.json    (冻结纯轴向 risk_by_height / stored σ_c)
#   results/opensim_fe/fracture_matrix.json      (heights_m 1..50)
#   results/opensim_fe/plantar_bc_summary.json   (FE 实测 gauge(k) / fixed / 外推高度)
#   src/climbing/bone.py MATERIAL_STRENGTH_MPA  (弱子区域 σ_c)
# 产物：
#   results/opensim_fe/pad_stiffness.json
#   results/opensim_fe/PAD_STIFFNESS_REPORT.md
#   results/opensim_fe/pad_stiffness_curve.png
```

---

## 8. 产物

- 本报告：`results/opensim_fe/PAD_STIFFNESS_REPORT.md`
- 新脚本：`scripts/opensim_fe/pad_stiffness.py`（additive）
- 机器可读：`results/opensim_fe/pad_stiffness.json`
- 曲线：`results/opensim_fe/pad_stiffness_curve.png`

