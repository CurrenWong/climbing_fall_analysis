# 全身各关节反力（自由体法）@ 5 m dead-drop

**目标**：把 `joint_loads.subtalar_reaction` 的自由体数学泛化到
**全身各关节**，为 50×9 骨折矩阵提供**载荷端**（每根骨/节段的受力与时刻）。

**方法**：对关节 `J` 的**远端刚体子树** `D`（child body + 全部子代）列牛顿第二定律：

```
R = Σᵢ mᵢ(aᵢ − g) − Σ F_ext
M = Σᵢ (pᵢ − p_J)×(mᵢ aᵢ) − Σᵢ (pᵢ − p_J)×(mᵢ g) − Σ (p_apply − p_J)×F_ext
```

关节中心 `p_J` = 该关节 **child body 的体坐标原点**；所有量在**全局系**。
与 `subtalar_reaction` **同一套数学**（`_mass_center` / `_as_vec3` 直接复用），
并修正了 `_realized_states` 的**状态时间**（见 §4）。只在 GRF 窗口内分析，
避开 `PrescribedForce` 的分段线性**窗口外推**（实测 +3 ms 处约 −88 kN）。

- 模型：`model\opensim\FullBodyModel-4.0\Rajagopal2015.osim`（总质量 75.337 kg，未缩放）
- 坠落：h=5.0 m；GRF peak=53.06 kN，窗口 0.00→21.10 ms，T1 err=+5.8%
- FD：1206 帧，末时 24.10 ms；腿关节锁定 + 肌肉禁用（刚性腿，S1 约定）

## 1. 关节 → 远端子刚体表（Rajagopal2015 joint 树）+ 子树质量

拓扑自动校验（声明集 vs 模型 joint 树）：**✅ 全部一致**；腿部子树质量近端→远端单调递增：**✅**。

| 关节 | 侧 | 远端子刚体（含 child body） | 子刚体数 | 子树质量 (kg) | 自重 (N) |
|---|---|---|---|---|---|
| `subtalar_r` | r | calcn_r;toes_r | 2 | 1.4666 | 14.38 |
| `ankle_r` | r | talus_r;calcn_r;toes_r | 3 | 1.5666 | 15.36 |
| `knee_r` | r | tibia_r;talus_r;calcn_r;toes_r | 4 | 5.2741 | 51.72 |
| `hip_r` | r | femur_r;tibia_r;patella_r;talus_r;calcn_r;toes_r | 6 | 14.6617 | 143.78 |
| `subtalar_l` | l | calcn_l;toes_l | 2 | 1.4666 | 14.38 |
| `ankle_l` | l | talus_l;calcn_l;toes_l | 3 | 1.5666 | 15.36 |
| `knee_l` | l | tibia_l;talus_l;calcn_l;toes_l | 4 | 5.2741 | 51.72 |
| `hip_l` | l | femur_l;tibia_l;patella_l;talus_l;calcn_l;toes_l | 6 | 14.6617 | 143.78 |
| `lumbar` | — | torso;humerus_r;ulna_r;radius_r;hand_r;humerus_l;ulna_l;radius_l;hand_l | 9 | 34.2366 | 335.75 |

> ⚠️ 与任务清单的两点差异（**均为模型 joint 树的必然结果，T2 已自检**）：
> - `hip_{r,l}` 子树**包含 `patella`**（`patellofemoral` 关节的 parent 是 `femur`，
>   patella 在髋切面远端）；任务列的 femur/tibia/talus/calcn/toes 是其真子集。
> - `knee_{r,l}` 子树**不含 patella**（patella 挂在 femur 侧 = 膝切面近端）。
> - 本模型**无 head/neck 刚体**：`torso` 是躯干链顶，`lumbar` = torso + 双上肢。

## 2. T2 自检（静止、无外力：R_y = 子树远端自重）

构造**全坐标锁死**的静止模型（T2 语义：身体被支撑住）。残差 `R_y − m_d·g`
应 << 冲击峰值。

| 关节 | 子树质量 (kg) | 自重 (N) | R_y (N) | 残差 (N) | 相对残差 |
|---|---|---|---|---|---|
| `subtalar_r` | 1.4666 | 14.38 | 14.38243 | +2.634e-11 | 1.83e-12 |
| `ankle_r` | 1.5666 | 15.36 | 15.36310 | +2.417e-11 | 1.57e-12 |
| `knee_r` | 5.2741 | 51.72 | 51.72125 | -3.526e-11 | 6.82e-13 |
| `hip_r` | 14.6617 | 143.78 | 143.78216 | -1.810e-10 | 1.26e-12 |
| `subtalar_l` | 1.4666 | 14.38 | 14.38243 | -2.311e-11 | 1.61e-12 |
| `ankle_l` | 1.5666 | 15.36 | 15.36310 | -2.834e-11 | 1.84e-12 |
| `knee_l` | 5.2741 | 51.72 | 51.72125 | -2.019e-10 | 3.90e-12 |
| `hip_l` | 14.6617 | 143.78 | 143.78216 | -6.403e-10 | 4.45e-12 |
| `lumbar` | 34.2366 | 335.75 | 335.74635 | -2.329e-09 | 6.94e-12 |

最大 T2 残差 = 2.329e-09 N，比 @5 m 最大关节峰值 25.38 kN **小 ~1e+13 倍** → ✅ 通过。

## 3. @5 m 各关节峰值（窗口内：0.00–21.10 ms）

主口径 = **修正状态时间**；`legacy` = 既有 `joint_loads._realized_states` 口径。

| 关节 | 侧 | 子树质量(kg) | peak\|F\| (kN) | t (ms) | 纵向 (kN) | t (ms) | peak\|M\| (N·m) | t (ms) | legacy peak\|F\| (kN) |
|---|---|---|---|---|---|---|---|---|---|
| `subtalar_r` | r | 1.467 | 25.380 | 20.10 | 25.354 | 20.10 | 156.72 | 20.20 | 26.383 |
| `ankle_r` | r | 1.567 | 25.295 | 20.10 | 25.266 | 20.10 | 1440.85 | 20.20 | 26.373 |
| `knee_r` | r | 5.274 | 22.402 | 20.10 | 22.200 | 20.10 | 1592.27 | 20.20 | 26.091 |
| `hip_r` | r | 14.662 | 15.489 | 20.10 | 15.228 | 20.10 | 989.40 | 20.20 | 25.806 |
| `subtalar_l` | l | 1.467 | 25.380 | 20.10 | 25.354 | 20.10 | 156.72 | 20.20 | 26.383 |
| `ankle_l` | l | 1.567 | 25.295 | 20.10 | 25.266 | 20.10 | 1440.85 | 20.20 | 26.373 |
| `knee_l` | l | 5.274 | 22.402 | 20.10 | 22.200 | 20.10 | 1592.27 | 20.20 | 26.091 |
| `hip_l` | l | 14.662 | 15.489 | 20.10 | 15.228 | 20.10 | 989.40 | 20.20 | 25.806 |
| `lumbar` | — | 34.237 | 22.168 | 20.10 | 22.088 | 20.10 | 1633.70 | 20.20 | 1.495 |

> **为何近端 |F| 反而小？** 本场景是**刚性锐峰**：峰值瞬间足底 GRF 向上 ~26.5 kN
> （单足），已**超过**远端腿本身加速所需；`R = m_d(a−g) − F_ext` 中 `−F_ext` 占主导，
> 而 `m_d(a−g)`（子树的惯性项，向上为正）随纳入的远端质量增多而增大，**抵消**掉一部分
> `F_ext`，于是净 `|R|` 从远端的 ~25 kN 递减到髋的 ~15 kN。这是**刚性腿 + 无肌肉**
> 近似的直接后果，不代表活体趋势。
>
> **lumbar** 修正口径 ~22 kN（旧口径 1.5 kN，被时间 bug 抹平）；代表腰-骨盆交界处
> 躯干+上肢子系统的合力，量级与全身减速一致。

## 4. 状态时间修正（关键发现）与距下回归对照

既有 `_realized_states` **从不写 state 的时间**（实测 `ws.getTime()==0`），
而 OpenSim `PrescribedForce` 按 `state.getTime()` 求值 → 所有帧都用 `GRF(0)`。
独立证据（骨盆 COM 的 `a_y` vs `pelvis_ty/speed` 的数值微分 `d(speed)/dt`）：

| 口径 | max\|a_y\| (m/s²) | max\|d(speed)/dt\| | max\|误差\| | 相关系数 |
|---|---|---|---|---|
| legacy | 7.7 | 2009.9 | 2002.2 | 0.1295 |
| timed | 2382.4 | 2009.9 | 372.5 | 0.9997 |

旧口径的加速度几乎是自由落体（≈ −g），与真实轨迹**不相关**；修正口径相关系数≈1。
对**含足底外力**的关节（subtalar/ankle/knee/hip），`R = m_d(a−g) − F_ext` 中 `−F_ext`
（用插值 GRF，本就正确）占主导，故旧口径看起来对；但对**子树内无外力**的
`lumbar`，`R` 完全由被污染的加速度决定 → 旧口径**不可信**。

- `subtalar_r` 修正口径峰值 = **25.380 kN**；legacy = 26.383 kN（既有 `JOINT_LOAD_REPORT.md` 记 F_subt ≈ 26.383 kN，落在两者之间/同量级）。差异来自修正后的真实足部惯性项。
- 本模块**未改动** `joint_loads.py`；legacy 值可由 `set_state_time=False` 复现。

## 5. 诚实边界

- **刚性腿**：S1 把下肢关节坐标**锁定**、肌肉禁用，故膝/踝/髋的反力是**刚体柱体**近似下的结果；真实关节有顺应性与主动肌力，会改变峰值与分布。
- **无肌肉**：关节反力完全来自**惯性 + 重力 + 足底 GRF**；活体肌肉会分担/重分配载荷（尤其髋、腰）。
- **矩的口径**：`M` 用现有自由体的动量矩率口径（对关节中心，未做移动参考点的`−v_J×P` 修正），与 `subtalar_reaction` 完全一致；T2 静止验证成立，冲击期的绝对矩应作**同口径相对比较**，非临床预测。
- **lumbar 语义**：模型只有单个 `back`（pelvis→torso）关节，无分节腰椎；`lumbar` 反力代表**腰-骨盆交界**处躯干+上肢子系统的合力/矩。
- **外力近似**：GRF 来自 `pad.py` 的 1D 等效双质点模型、施加于 calcn 原点，是**载荷边界近似**（方案 §3.1）。

## 6. 产物

- `results/opensim_fe/JOINT_REACTIONS_REPORT.md`（本文件）
- `results/opensim_fe/joint_reactions_summary.json` / `.csv`
- `temp/opensim_fe/joint_reactions/joint_reactions_timeseries.npz`（逐帧 F/M）
- `temp/opensim_fe/joint_reactions/subtree_mass_table.csv` / `t2_residuals.csv` / `run.log`

