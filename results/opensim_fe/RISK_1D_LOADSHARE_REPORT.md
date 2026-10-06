# Route-1 1D 骨折风险 —— 载荷分配修订版（S4 清单⑤ 输入）

> 生成时间 2026-10-04 22:13:51；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s。
> 本报告**新增**，不改任何既有模块/默认值，不覆盖既有产物，不调用 FEBio。
> 修订对象：`BC_ROBUST_METRIC_REPORT.md` §2 的 Route-1 1D 名义利用率 `risk(h)=F_joint(h)/(A_section·σ_c)`。

---

## 0. 一句话结论

Route-1 的两个载荷映射伪影（胫/腓共用整条踝反力；上段中轴骨共用整条腰椎反力）已用**载荷分配**修订：

1. **胫/腓并联分流**：按轴向刚度 `k=E·A/L` 并联分配，得 **tibia 76.8% / fibula 23.2%**。纯并联刚度给腓骨 ~23%，**高于**文献的 5–10%——原因是真实腓骨经胫腓联合（骨间膜/下胫腓）串联合加载，并非等位移并联柱；报告同时给出文献锚定 10% / 5% 两种分流。
2. **中轴骨按“其上方质量”分配**：`F_level(h) = (m_above(level)/m_above(L5/S1))·F_lumbar(h)`，参考质量 = 49.7% BM = 37.44 kg。L3/T6/C5/parietal 的力分别降到腰椎反力的 89.5%、52.9%、15.1%、9.5%。
3. **排序显著改变**：无分配时腓骨 #1、C5 #2；修订后（纯并联）`calcaneus_r` 升到 #1，腓骨退到 #2，C5 退到 #7，T6 升到 #3/4。

**重要范围声明**：本修订修的是**载荷账目（load bookkeeping）**，不改变 1D 名义应力的根本局限（忽略弯曲/剪切/偏心）。它把 Route-1 从“载荷明显错误”变成“载荷口径正确但仍是 1D”，是当前 S4 可用的**较不坏**输入。

---

## 1. 方法

### 1.1 基础定义（沿用 Route 1，未改）

```
risk(h) = F_model(h) / (A_section × σ_c)
first_fracture_height = 最小的 h 使 risk(h) ≥ 1（否则 > 50 m）
```

| 量 | 来源 | 性质 |
|---|---|---|
| `A_section` | `results/opensim_fe/bc_robust_metric_route1.json` | 4 长骨 measured CORT；5 中轴骨 estimated body |
| `σ_c` | `src/climbing/bone.py` `MATERIAL_STRENGTH_MPA`（压缩） | assumed（[Y25] 表1） |
| `F_joint(h)` | `results/opensim_fe/fracture_matrix.json` → `loads_n` | measured（逐高度） |

载荷来源映射（@5 m 实测参考值：`ankle_r`=25,266.2 N，`subtalar_r`=25,354.4 N，`hip_r`=15,227.7 N，`lumbar`=22,087.6 N）。

### 1.2 修正 a —— 胫/腓并联分流（k = E·A/L）

把胫、腓视为**并联轴向弹簧**（同一位移，分担总踝反力），按刚度分流：

```
k_i = E · A_i / L_i ,   f_i = k_i / (k_tibia + k_fibula) ,   F_i = f_i · F_ankle
E = 18000 MPa（皮质骨，assumed）
```

| 骨 | A_section (mm²) | L (mm) | k = E·A/L (N/mm) | **分流 f** |
|---|---:|---:|---:|---:|
| `tibia_r` | 362.65 | 327.6 | 19,924.4 | **0.7677** |
| `fibula_r` | 101.94 | 304.3 | 6,030.5 | **0.2323** |

> **结果：腓骨分担 23.2%**（胫骨 76.8%）。任务预期腓骨 ~5–10%；纯并联公式给出的 23% 偏高。物理原因：
> - 真实小腿的腓骨经**胫腓联合 + 骨间膜**间接加载，存在**串联柔度**（位移不等），不是理想并联柱；
> - 腓骨还承受弯曲/张力带作用，轴向柱刚度并非全部用于分流；
> - `bone.py` 明确记载腓骨单侧仅分担 ~5%（双侧 ~10%），即 **5–10%** 才是可信区间。

因此本脚本**同时**给出三档腓骨分流用于敏感性（§5）：纯并联 23.2%、文献锚定 10%、文献锚定 5%。
**缺少胫腓联合刚度的实测值，故无法从第一性原理推出 5–10%；此处明确标注为 assumed。**

### 1.3 修正 b —— 中轴骨/颅骨按“上方质量”分配

**第一步（rigorous 尝试）：用 `joint_reactions.py` 算各上位节段的自由体关节反力。**

实测：**不可行**。Rajagopal2015 模型只有一个刚性 `torso` 体（无 head/neck/椎体刚体），`joint_reactions.JOINT_SPECS` 仅定义固定关节：

- 支持关节：`['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']`；
- `lumbar` 远端子刚体 = `['torso', 'humerus_r', 'ulna_r', 'radius_r', 'hand_r', 'humerus_l', 'ulna_l', 'radius_l', 'hand_l']`（整体躯干+双上肢，无法再按椎体切分）；
- 逐级探测（C5/T6/L3/neck/head/skull）全部无法解析。

<details><summary>探测明细（点击展开）</summary>

- `C5` → **ValueError: 未知关节 'C5'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**
- `T6` → **ValueError: 未知关节 'T6'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**
- `L3` → **ValueError: 未知关节 'L3'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**
- `neck` → **ValueError: 未知关节 'neck'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**
- `head` → **ValueError: 未知关节 'head'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**
- `skull` → **ValueError: 未知关节 'skull'（side=None）；可选：['ankle_l', 'ankle_r', 'hip_l', 'hip_r', 'knee_l', 'knee_r', 'lumbar', 'subtalar_l', 'subtalar_r']**

</details>

**第二步（fallback）：mass-above 模型。**

对上位节段，隔离体为其**上方全部质量**；刚性下降期内各段加速度近似一致，由自由体恒等式 `F_level = m_above·(a−g)` 与 `F_lumbar = m_ref·(a−g)` 相除得：

```
F_level(h) = ( m_above(level) / m_above(L5/S1) ) × F_lumbar(h)
```

`m_above` 取自公开人体节段质量分数（Pearsall 1996 胸腰椎 + Ivancic 2006 头/颈 + Dempster/Winter 上肢），按“经过该节段下位椎间盘”的累计口径（含该节段自身）。

| 节段 | 来源 | 质量 (% BM) |
|---|---|---:|
| head | Ivancic 2006 | 4.7 |
| C1–C7 | Ivancic 2006 | 0.6/0.7/0.5/0.5/0.5/0.6/0.7 |
| T1–T12 | Pearsall 1996 | 1.1/1.1/1.4/1.3/1.3/1.3/1.4/1.5/1.6/2.0/2.1/2.5 |
| L1–L5 | Pearsall 1996 | 2.4/2.4/2.3/2.6/2.6 |
| 双上肢 | Dempster/Winter | 10.0（挂于肩/T1） |

**参考（分母）** = 经过 L5/S1 椎间盘以上质量 = 头+全颈+全胸+全腰+双上肢 = **49.7% BM = 37.44 kg**（模型 `lumbar` 子树实测质量 34.24 kg，同量级，差异 <10%）。

**上肢归属**：双上肢载荷在肩部（约 T1）传入脊柱 → 计入 T6/L3/L5 的上方质量，**不计入 C5 及颅骨**（其位于肩以上）。

| 骨 | 对应节段 | 上方质量 (% BM) | 上方质量 (kg) | **力比 = m/m_ref** |
|---|---|---:|---:|---:|
| `L3` | L3 | 44.5 | 33.52 | **0.8954** |
| `T6` | T6 | 26.3 | 19.81 | **0.5292** |
| `C5` | C5 | 7.5 | 5.65 | **0.1509** |
| `parietal_r` | head | 4.7 | 3.54 | **0.0946** |

**颅骨口径说明**：`parietal_r` 是终端骨，其“上方质量”按**头部自身质量**（4.7% BM）解释（颅骨承受头部惯性载荷），而非“颅顶之上的身体质量”（后者≈0，不合理）。

### 1.4 逐骨力映射（修订后）

| 骨 | 关节 | 力倍率（模型） |
|---|---|---:|
| `calcaneus_r` | `subtalar_r` | 1.0000 ×（单骨，无分配） |
| `fibula_r` | `ankle_r` | 0.2323 ×（并联分流） |
| `tibia_r` | `ankle_r` | 0.7677 ×（并联分流） |
| `T6` | `lumbar` | 0.5292 ×（上方质量比） |
| `L3` | `lumbar` | 0.8954 ×（上方质量比） |
| `femur_r` | `hip_r` | 1.0000 ×（单骨，无分配） |
| `C5` | `lumbar` | 0.1509 ×（上方质量比） |
| `R_HIPBONE` | `hip_r` | 1.0000 ×（单骨，无分配） |
| `parietal_r` | `lumbar` | 0.0946 ×（上方质量比） |

---

## 2. 修订后逐骨表（主要：纯并联分流，按修订排名排序）

| # | 骨 | 部位 | 关节 | σ_c (MPa) | A_section (mm²) | F_model@5m (N) | 力倍率 | **risk@5m** | **首次骨折 h (m)** |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
| 1 | `calcaneus_r` | 足部 | `subtalar_r` | 150 | 216.01 | 25,354.4 | 1.0000 | **0.7825** | **8** |
| 2 | `fibula_r` | 腓骨 | `ankle_r` | 160 | 101.94 | 5,870.5 | 0.2323 | **0.3599** | **23** |
| 3 | `tibia_r` | 胫骨 | `ankle_r` | 200 | 362.65 | 19,395.8 | 0.7677 | **0.2674** | **37** |
| 4 | `T6` | 胸椎 | `lumbar` | 150 | 291.20 | 11,688.2 | 0.5292 | **0.2676** | **47** |
| 5 | `L3` | 腰椎 | `lumbar` | 150 | 692.60 | 19,776.6 | 0.8954 | **0.1904** | **>50** |
| 6 | `femur_r` | 股骨 | `hip_r` | 220 | 404.87 | 15,227.7 | 1.0000 | **0.1710** | **>50** |
| 7 | `C5` | 颈椎 | `lumbar` | 150 | 178.83 | 3,333.1 | 0.1509 | **0.1243** | **>50** |
| 8 | `R_HIPBONE` | 骨盆 | `hip_r` | 180 | 1947.68 | 15,227.7 | 1.0000 | **0.0434** | **>50** |
| 9 | `parietal_r` | 颅骨 | `lumbar` | 160 | 877.53 | 2,088.8 | 0.0946 | **0.0149** | **>50** |

**修订排序（纯并联，早骨折优先，其次 risk@5m）**：

`calcaneus_r(#1, 8m) > fibula_r(#2, 23m) > tibia_r(#3, 37m) > T6(#4, 47m) > L3(#5, >50m) > femur_r(#6, >50m) > C5(#7, >50m) > R_HIPBONE(#8, >50m) > parietal_r(#9, >50m)`

---

## 3. with-share vs without-share 对比

| 骨 | risk@5m 无分配 | risk@5m 并联 | Δ倍数 | 首次骨折 无分配 | 首次骨折 并联 | 排名 无→并联 |
|---|---:|---:|---:|---|---|---|
| `fibula_r` | 1.5491 | 0.3599 | 0.232× | 2 | 23 | #1 → **#2** |
| `C5` | 0.8234 | 0.1243 | 0.151× | 8 | >50 | #2 → **#7** |
| `calcaneus_r` | 0.7825 | 0.7825 | 1.000× | 8 | 8 | #3 → **#1** |
| `T6` | 0.5057 | 0.2676 | 0.529× | 17 | 47 | #4 → **#4** |
| `tibia_r` | 0.3484 | 0.2674 | 0.768× | 25 | 37 | #5 → **#3** |
| `L3` | 0.2126 | 0.1904 | 0.895× | >50 | >50 | #6 → **#5** |
| `femur_r` | 0.1710 | 0.1710 | 1.000× | >50 | >50 | #7 → **#6** |
| `parietal_r` | 0.1573 | 0.0149 | 0.095× | >50 | >50 | #8 → **#9** |
| `R_HIPBONE` | 0.0434 | 0.0434 | 1.000× | >50 | >50 | #9 → **#8** |

### 3.1 腓骨如何移动（关键）

| 口径 | 腓骨分流 | 腓骨 risk@5m | 腓骨首次骨折 | 腓骨排名 |
|---|---:|---:|---|---:|
| 无分配 | 1.000 | 1.5491 | 2 | #1 |
| **纯并联 k=E·A/L** | 0.232 | 0.3599 | 23 | **#2** |
| 文献锚定 10% | 0.100 | 0.1549 | >50 | #6 |
| 文献锚定 5% | 0.050 | 0.0775 | >50 | #7 |

> 腓骨从**无分配的 #1（risk@5m=1.55，2 m 即“骨折”）**降到**纯并联 #2（0.36）**；若用文献 10%，进一步降到 **#6（0.15，>50 m 不骨折）**。这正说明旧 Route-1 的“腓骨最危险”是**载荷未分配的人为放大**（`BC_ROBUST_METRIC_REPORT.md` §7 边界 1）。

---

## 4. 与其他排序的对比

| 方法 | 排序（1→9） |
|---|---|
| **本报告（无分配）** | `fibula_r > C5 > calcaneus_r > T6 > tibia_r > L3 > femur_r > parietal_r > R_HIPBONE` |
| **本报告（并联分流）** | `calcaneus_r > fibula_r > tibia_r > T6 > L3 > femur_r > C5 > R_HIPBONE > parietal_r` |
| **本报告（文献 10%）** | `calcaneus_r > tibia_r > T6 > L3 > femur_r > fibula_r > C5 > R_HIPBONE > parietal_r` |
| **FE 矩阵（修正后 p95/σ_c, util@5m）** | `fibula_r > tibia_r > T6 > femur_r > C5 > parietal_r > L3 > calcaneus_r > R_HIPBONE` |

Spearman 秩相关（9 骨，1 = 同序）：

| 对比 | ρ | 解读 |
|---|---|---|
| 并联分流 vs 无分配 | **0.700** | 载荷分配把腓骨/中轴骨位置显著改写 |
| 并联分流 vs FE 矩阵 | **0.367** | 中等 —— FE p95 受弯曲 BC 伪影污染 |
| 无分配 vs FE 矩阵 | 0.517 | （与 `BC_ROBUST_METRIC_REPORT.md` §4 的 0.217 同量级） |
| 文献10% vs FE 矩阵 | 0.183 | |

> **注意**：FE 矩阵排序（`fibula_r > tibia_r > T6 > femur_r > C5 > parietal_r > L3 > calcaneus_r > R_HIPBONE`）这里用的是**修正后的**`fracture_matrix.json` 的 `utilization_at_5m`（p95/σ_c）。但 `BC_ROBUST_METRIC_REPORT.md` §3–4 已证明 FE 的 p95 被边界条件（斜置骨沿 bbox 轴加载 → 中段弯曲）主导，故其 9 骨相对排序**不可信**；此处仅作**交叉参照**，不作为真值。

---

## 5. 敏感性

### 5.1 腓骨分流档位

| 腓骨分流 | risk@5m | 首次骨折 | 备注 |
|---|---:|---|---|
| 纯并联 k=E·A/L (0.232) | 0.3599 | 23 | 本报告主口径 |
| 文献 10% | 0.1549 | >50 | `bone.py` 双侧分担 |
| 文献 5% | 0.0775 | >50 | `bone.py` 单侧分担 |

### 5.2 中轴骨上方质量比 ±20%（S4 排序稳定性）

| 骨 | 力比 | risk@5m ×0.8 | risk@5m ×1.0 | risk@5m ×1.2 |
|---|---:|---:|---:|---:|
| `L3` | 0.8954 | 0.1523 | 0.1904 | 0.2284 |
| `T6` | 0.5292 | 0.2141 | 0.2676 | 0.3211 |
| `C5` | 0.1509 | 0.0994 | 0.1243 | 0.1491 |
| `parietal_r` | 0.0946 | 0.0119 | 0.0149 | 0.0179 |

> ±20% 的质量分数扰动不改变中轴骨“全部远低于 1”的结论（最高的 T6 在 ×1.2 下仍 0.321 < 1）→ **中轴骨排序对质量分数不敏感**，主要受 `F_lumbar(h)` 时间历程控制。

---

## 6. 假设（ASSUMPTIONS）与诚实边界

**A. 胫/腓分流规则**

- 并联轴向弹簧 `k=E·A/L`，E=18000 MPa（皮质骨，**assumed**）；A=`A_section`（measured CORT）；L=bbox 最长轴（measured）。
- 纯并联给腓骨 ~23%，**高于**文献 5–10%；差异归因于胫腓联合串联合柔度 + 腓骨弯曲，脚本另给 10% / 5% 两档。**无实测胫腓联合刚度，故 5–10% 为 assumed 而非推导。**

**B. 中轴骨/颅骨力规则 + 质量分数**

- `F_level = (m_above(level)/m_above(L5/S1))·F_lumbar`（刚性上体，同加速度）。
- 质量分数：Pearsall 1996（T/L）+ Ivancic 2006（head/C）+ Dempster/Winter（arms）；参考 = 49.7% BM。上肢挂于肩（计入 T6/L3/L5，不计 C5/颅骨）。
- 颅骨按头部自身质量（4.7% BM）为终端惯性载荷。
- 模型无 head/neck/椎体 → **无法用 `joint_reactions.py` 做真正的逐级自由体**；mass-above 是给定模型下的严格刚体近似（非“任意节段自由体”）。
- 质量分数为人群均值，个体差异 ±20% 已由 §5.2 覆盖。

**C. 1D 极限（未修）**

- `risk = 压缩力/(截面×σ_c)` **忽略弯曲、剪切、偏心、应力集中**：对细长骨（胫/腓）低估风险，对紧凑骨（跟骨）高估 2–3×（`bone.py` 明确记载）。
- 中轴骨 `A_section` 是 THUMS 实体的 **estimate**（THUMS CORT 为壳，无实体），与真实皮质截面可能差 ~2×。
- `σ_c` 是**材料**强度而非整体骨失效载荷，绝对值不可信；可用的是**相对排序**。

**D. 未做什么**

- 不改任何既有模块默认/签名；不覆盖既有产物；不调用 FEBio；不重跑 OpenSim 正动力学。
- 未把本修订回灌进 `bc_robust_metric_route1.json`（保持只读）。

---

## 7. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\risk_1d_loadshare.py
# 仅只读复用：
#   results/opensim_fe/fracture_matrix.json           (loads_n / heights_m)
#   results/opensim_fe/bc_robust_metric_route1.json   (A_section / sigma_c)
#   src/climbing/coupling/joint_reactions.py          (逐级探测，未跑正动力学)
# 产物：
#   results/opensim_fe/risk_1d_loadshare.json
#   results/opensim_fe/RISK_1D_LOADSHARE_REPORT.md
```

---

## 8. 产物

- 本报告：`results/opensim_fe/RISK_1D_LOADSHARE_REPORT.md`
- 新脚本：`scripts/opensim_fe/risk_1d_loadshare.py`（additive）
- 机器可读：`results/opensim_fe/risk_1d_loadshare.json`

