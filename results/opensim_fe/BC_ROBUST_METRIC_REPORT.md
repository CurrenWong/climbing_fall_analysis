# BC-健壮性交叉验证报告 —— 9 骨相对骨折排序是否可信？

> 任务：在不重跑全量 FEBio 的前提下，用**两条独立路线**检验 `results/opensim_fe/fracture_matrix.json`（9 骨 × 50 m）的**相对**骨折排序是否可信，并回答「S4（Logistic + 层次聚类）该用什么指标」。
> 前置诊断：`results/opensim_fe/MESH_FIX_REPORT.md`（σ_vm 不是网格受限，而是**边界条件受限**）。
> 本报告数字除「estimate/assumed」标注外均为本机实测；生成时间 2026-10-04；仓库 `D:\Project\climbing_fall_analysis`；FEBio 4.13.0。
> 复现命令见 §8；回归结果见 §9。

---

## 0. 一句话结论（verdict）

**当前 50×9 矩阵的相对排序不可信，不能直接喂给 S4；必须重做。** 三条独立证据：

1. **无 FE 的 1D 名义应力排序（Route 1）与矩阵 p95 排序只有弱一致**（Spearman ρ=0.22）。矩阵把 `femur_r/T6/tibia_r` 排在最前，而 BC-free 路线把 `fibula_r` 排第一。
2. **端带剔除（Route 2，剔除顶/底 30% 夹持带）不但没有降低 p95，反而使每根骨的 p95 都升高**（tibia 2224→4308 MPa）。→ 热点在**中段弯曲**，不是端部夹持奇异；"end-clamp singularity 是 p95 主控" 这一假设被证伪。FE 的 p95 无法救回。
3. **额外发现（可复现缺陷）**：矩阵对 3 根降载骨（tibia/fibula/T6）的 `sigma_vm_ref` **没有乘回降载因子**（`fracture_matrix.py` 第 570 行），导致同一张矩阵里混入了**不同载荷下的原始应力**（tibia 少乘 5.05×、fibula 少乘 25.27×、T6 少乘 2.21×）。报告正文却声称已折算（§6）——**报告与代码不一致**。

**该喂 S4 的指标**：目前唯一不带 BC 伪影的是 **Route 1 的 1D 名义利用率 `risk(h)=F_joint(h)/(A_section·σ_c)`**（BC-free，9 骨区分度良好）。但 Route 1 也有明确假设（未建模胫/腓载荷分配、未建模弯曲），故它是「**当前最不坏**」的 S4 输入，而不是已验证的真值。**FE 的 p95/volavg/median 都不应作为 S4 的相对排序输入**（p95 被弯曲主导；volavg/median 与 Route 1 仅中等一致，见 §4）。

---

## 1. 背景与判据

- 清单⑤ 交付物是 **9 骨 × 50 m 的相对模式**（Logistic + 层次聚类 S4），不是绝对临界高度。
- 矩阵现有主口径 = 主承载域 σ_vm 的 **p95**；判据 `σ_vm(h) ≥ σ_c`，`σ_vm(h)=σ_vm_ref·load(h)/load_ref`（线弹性）。
- `MESH_FIX_REPORT` 已实测：p95 由 **BC**（斜置骨沿 bbox 轴加载 → 弯曲 + 端部夹持 + 点载荷）决定，而该伪影**骨形状相关**（倾斜角/细长比），因此会污染相对排序 → 必须用 BC-无关路线复核。

**Route 1（必需，FE-free）**：1D 截面强度估算。
**Route 2（尽力，复用现有干净网格）**：端带剔除 + 体积加权 FE 应力；不重跑 FEBio。
**交叉对比（必需）**：Route-1 排序 vs Route-2 排序 vs 当前矩阵排序。

---

## 2. Route 1 —— 1D 强度截面估计（FE-free）

### 2.1 定义与数据来源

```
risk(h) = F_joint(h) / (A_section × σ_c)
first_fracture_height = 最小的 h 使 risk(h) ≥ 1（否则 > 50 m）
```

| 量 | 来源 | 性质 |
|---|---|---|
| `F_joint(h)` | `results/opensim_fe/fracture_matrix.json` → `loads_n[<joint>]`（逐高度实测） | **measured** |
| `σ_c` | `src/climbing/bone.py` `MATERIAL_STRENGTH_MPA` 的**压缩**强度（元组第 2 项） | **assumed（[Y25] 表1）** |
| `A_section` | `temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5` 的实体域体积 ÷ bbox 最长轴长度 | 4 长骨 **measured**；5 中轴骨 **estimated** |

- **σ_c 映射**：`calcaneus_r→calcaneus(150)`、`tibia_r→tibia_shaft(200)`、`fibula_r→fibula_shaft(160)`、`femur_r→femur_shaft(220)`、`R_HIPBONE→pelvis(180)`、`L3/T6/C5→spine(150)`、`parietal_r→skull(160)`。与矩阵 JSON 的 `sigma_c_mpa` 逐骨一致。
- **A_section 来源（关键）**：4 长骨的 CORT 是 THUMS 实体 hex，V/L **精确复现** `MESH_FIX_REPORT` §2.3 的实测值（calc 216.0 / tibia 362.6 / fibula 101.9 / femur 404.9 mm²，误差 <0.05%）。5 中轴骨的 CORT 在 THUMS 中是**壳单元（0 实体）**，无法取 V/L；故取**实体体（SPON 椎体/髋骨体、颅骨 diploë）的 V/L** 作为承载截面 —— 明确标注 **estimate**。
- **F_joint 逐骨映射（必须）**：`calcaneus_r←subtalar_r`、`tibia_r/fibula_r←ankle_r`、`femur_r/R_HIPBONE←hip_r`、`L3/T6/C5/parietal_r←lumbar`。⚠️ **同一关节的多根骨共享同一个整体关节反力**（Route 1 未建模胫/腓载荷分配，见 §7 边界 1）。

### 2.2 Route 1 结果表（实测 V/L；σ_c assumed）

| 骨 | 部位 | 关节 | σ_key | σ_c (MPa) | V_solid (mm³) | L (mm) | **A_section (mm²)** | 来源 | F@5m (N) | 阈值 (kN) | **risk@5m** | **首次骨折 h (m)** |
|---|---|---|---|---|---|---|---|---|---:|---:|---:|---:|
| `fibula_r` | 腓骨 | ankle_r | fibula_shaft | 160 | 31,015.3 | 304.3 | **101.94** | measured CORT | 25,266 | 16.31 | **1.549** | **2** |
| `C5` | 颈椎 | lumbar | spine | 150 | 9,719.7 | 54.4 | **178.83** | est. body | 22,088 | 26.82 | **0.823** | **8** |
| `calcaneus_r` | 足部 | subtalar_r | calcaneus | 150 | 18,279.1 | 84.6 | **216.01** | measured CORT | 25,354 | 32.40 | **0.783** | **8** |
| `T6` | 胸椎 | lumbar | spine | 150 | 16,999.3 | 58.4 | **291.20** | est. body | 22,088 | 43.68 | **0.506** | **17** |
| `tibia_r` | 胫骨 | ankle_r | tibia_shaft | 200 | 118,811.5 | 327.6 | **362.65** | measured CORT | 25,266 | 72.53 | **0.348** | **25** |
| `L3` | 腰椎 | lumbar | spine | 150 | 50,901.6 | 73.5 | **692.60** | est. body | 22,088 | 103.89 | **0.213** | >50 |
| `femur_r` | 股骨 | hip_r | femur_shaft | 220 | 182,990.3 | 452.0 | **404.87** | measured CORT | 15,228 | 89.07 | **0.171** | >50 |
| `parietal_r` | 颅骨 | lumbar | skull | 160 | 119,376.8 | 136.0 | **877.53** | est. diploë | 22,088 | 140.40 | **0.157** | >50 |
| `R_HIPBONE` | 骨盆 | hip_r | pelvis | 180 | 400,380.3 | 205.6 | **1,947.68** | est. body | 15,228 | 350.58 | **0.043** | >50 |

**Route 1 排序（早骨折优先，其次 risk@5m）**：
`fibula_r(#1, 2 m) > C5(#2, 8 m) > calcaneus_r(#3, 8 m) > T6(#4, 17 m) > tibia_r(#5, 25 m) > L3(#6) > femur_r(#7) > parietal_r(#8) > R_HIPBONE(#9)`。

> 注：Route 1 只在 50 m 内给出 1 根骨（腓骨）"骨折"；其余为渐进的利用率梯度。这与矩阵"8/9 在 1 m 骨折"的退化图景形成强烈反差（§4）。

---

## 3. Route 2 —— 现有 FE 输出的 BC-健壮应力（端带剔除 + 体积平均）

### 3.1 方法

对每个存储的 hdf5：读末 state 的 per-element Voigt 应力 → von Mises → **剔除沿加载轴（bbox 最长轴）顶/底各 30% 的夹持带**（`bone_feb` 的 `top_frac=bot_frac=0.30`），并计算**体积加权平均**（实体按 tet/hex 体积，壳按面积）。**不重跑 FEBio**。

### 3.2 干净 remesh（4 长骨，主 BC = 斜向加载；`temp/opensim_fe/fracture_matrix_remesh/`）

| 骨 | 状态 | σ_nom(A_cort) (MPa) | 全域 p95 | **端带剔除 p95 (mid)** | **体积平均 volavg** | 全域 median | mid median |
|---|---|---:|---:|---:|---:|---:|---:|
| `calcaneus_r` | ok | 117.37 | 43.65 | **51.96** | **15.15** | 13.09 | 18.01 |
| `tibia_r` | ok（1万N×2.5266） | 69.67 | 2224.15 | **4308.38** | **303.82** | 38.17 | 995.17 |
| `femur_r` | ok | 37.61 | 707.60 | **1012.41** | **112.31** | 3.86 | 406.41 |
| `fibula_r` | **FAILED** | 247.9 | — | — | — | — | — |

> **关键实测：端带剔除后 p95 反而升高（每根骨都是）。** `fibula_r` 所有候选载荷（25266/10000/5000/2500/1000 N）均发散（`negative jacobians`，细长斜置杆屈曲），无 BC-健壮 FE 应力可用。

### 3.3 PCA 轴向加载诊断（隔离"斜向加载"这一 BC 分量）

| 骨（轴向） | σ_nom(A_fe) (MPa) | 全域 p95 | 端带剔除 p95 | 体积平均 volavg | 全域 median | mid median |
|---|---:|---:|---:|---:|---:|---:|
| `calcaneus_r_axial` | 23.79 | 34.63 | 39.62 | 13.36 | 12.27 | 19.32 |
| `tibia_r_axial` | 28.98 | 758.04 | 1785.06 | 122.38 | 21.27 | 358.42 |
| `femur_r_axial` | 12.67 | 213.77 | 305.16 | 35.30 | 3.99 | 127.07 |

> 即便把载荷摆正到骨轴，胫骨 p95 仍 26× 名义，且**端带剔除后更高**（758→1785）——说明残余高应力在**中段**（弯曲/截面形状），端部夹持不是 p95 主控。`calcaneus` 是唯一紧凑、无斜置、p95 已物理者（0.37× 名义）。

### 3.4 原始 9 骨重分析（`fracture_matrix/`，与矩阵同网格，去掉端带）

按矩阵规则（有壳取跨壳最大）对每骨候选主域重算；`all_p95 = 原始解应力 × (ref/solve)`（**正确折算**）。

| 骨 | 折回因子 scale | all_p95 | 端带剔除 mid_p95 | volavg | all_median | **矩阵存的 stored_p95** | volavg/σ_c | median/σ_c |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `calcaneus_r` | 1.000 | 257.15 | 311.61 | 87.83 | 65.46 | 257.15 | 0.59 | 0.44 |
| `tibia_r` | **5.053** | **4402.78** | 5651.88 | 1016.16 | 192.48 | 871.28 | 5.08 | 0.96 |
| `fibula_r` | **25.266** | **16596.45** | 20471.07 | 4786.66 | 2032.52 | 656.87 | **29.92** | **12.70** |
| `femur_r` | 1.000 | 1168.84 | 1295.84 | 352.58 | 37.91 | 1168.84 | 1.60 | 0.17 |
| `R_HIPBONE` | 1.000 | 145.68 | 176.36 | 51.51 | 44.16 | 145.68 | 0.29 | 0.25 |
| `L3` | 1.000 | 303.14 | 340.88 | 102.41 | 77.91 | 303.14 | 0.68 | 0.52 |
| `T6` | **2.209** | **1716.77** | 2377.36 | 575.31 | 541.81 | 777.24 | 3.84 | 3.61 |
| `C5` | 1.000 | 641.14 | 677.06 | 198.75 | 131.58 | 641.14 | 1.32 | 0.88 |
| `parietal_r` | 1.000 | 423.41 | 515.44 | 198.16 | 183.27 | 423.41 | 1.24 | 1.15 |

> 同样：**每根骨 `mid_p95 > all_p95`**（端带剔除恶化）。`stored_p95` 与 `all_p95` 在 3 根降载骨上不同，正是 §5 的折算缺陷。

---

## 4. 交叉对比（REQUIRED）

按「应力/σ_c 利用率」排序（1 = 最危险）：

| 方法 | 排序（1→9） |
|---|---|
| **当前矩阵 p95**（stored） | `femur_r > T6 > tibia_r > C5 > fibula_r > parietal_r > L3 > calcaneus_r > R_HIPBONE` |
| **Route 1（1D 名义，BC-free）** | `fibula_r > C5 > calcaneus_r > T6 > tibia_r > L3 > femur_r > parietal_r > R_HIPBONE` |
| **Route 2（FE volavg/σ_c）** | `fibula_r > tibia_r > T6 > femur_r > C5 > parietal_r > L3 > calcaneus_r > R_HIPBONE` |
| **Route 2（FE median/σ_c）** | `fibula_r > T6 > parietal_r > tibia_r > C5 > L3 > calcaneus_r > R_HIPBONE > femur_r` |
| **Route 2（FE scaled-p95/σ_c）** | `fibula_r > tibia_r > T6 > femur_r > C5 > parietal_r > L3 > calcaneus_r > R_HIPBONE` |

Spearman 秩相关（ρ，9 骨）：

| 对比 | ρ | 解读 |
|---|---|---|
| Route 1 vs **矩阵 p95** | **0.217** | 弱 —— 矩阵排序与 BC-free 路线基本不一致 |
| Route 1 vs Route 2 volavg | **0.517** | 中等 |
| Route 1 vs Route 2 median | **0.500** | 中等 |
| Route 1 vs Route 2 scaled-p95 | **0.517** | 中等 |
| 矩阵 p95 vs Route 2 scaled-p95 | **0.767** | 高 —— 矩阵 ≈ 折算错误版 FE-p95 排序 |

**唯一稳定的跨路线信号**：`fibula_r` 在两条 BC-健壮路线里都最危险；`R_HIPBONE` 都最安全。但 `calcaneus_r`（Route1 #3 vs Route2 #8）、`C5`（#2 vs #5）、`femur_r`（#7 vs #3/#4）等**相互矛盾**。

**Answer (a) 是否一致？** **不一致。** Route 1（BC-free）与矩阵 p95 只有 ρ=0.22；与 Route 2（BC-健壮 FE）仅 ρ≈0.5（中等）。排序随指标显著翻转 → 相对模式**不稳健**。

**Answer (b) 哪个指标可信？** 目前唯一不带 BC 伪影的是 **Route 1 的 1D 名义利用率 `risk(h)=F_joint(h)/(A_section·σ_c)`**，它可作为 S4 的输入。**FE 的 p95 明确不可信**（弯曲主导、端带剔除反而恶化）；FE 的 volavg/median 量级物理，但其 9 骨排序与 Route 1 仅中等一致，且两者都被斜向加载弯曲塑形，故**不足以**作为 S4 的排序依据。

**Answer (c) 当前矩阵可否直接用？** **不可以，必须重做。** 三条独立理由：(i) p95 相对排序由骨倾斜/细长（BC 伪影）塑形（§3.3）；(ii) 3 根降载骨应力未折算（§5），矩阵混入不同载荷应力；(iii) 二值矩阵退化为「8/9 在 1 m」——对 S4 无任何可聚类梯度。

---

## 5. 附加发现：降载应力未折算（可复现缺陷）

`FRACTURE_MATRIX_REPORT.md` §3/§6 声称 `σ_vm_ref = 收敛载荷下的统计量 × 线性折回因子`。实测**代码未乘折回因子**：

- `scripts/opensim_fe/fracture_matrix.py:570`：`row["sigma_vm_ref_mpa"] = {s: main["stats"][s] ...}`（直接用 solved 统计量）；
- 第 576 行 `sig[s][pi] = main["stats"][s] * ratio`，其中 `ratio = load_h / ref_load`——把 solved 应力当成了 ref 应力。

**实测证据**：`tibia_r` 的 `fe_solve_load_n=5000`、`load_reduction_factor=5.0532`，但其 `tibia_r.hdf5` 末 state 的 CORT p95=**871.285 MPa**，与矩阵存的 `stored_p95=871.28` **完全相等**；`tibia_r.feb` 的节点力 `<value>-0,-0,-1.311303436</value>`（≈5 kN 量级，非 25 kN）。→ 矩阵存的确实是**降载原始应力**。正确折算后应为 871.285×5.0532=**4402.8 MPa**。

| 骨 | fe_solve_load_n | 未折算（矩阵 stored p95） | 正确折算 (×lrf) | 少乘 |
|---|---:|---:|---:|---:|
| `tibia_r` | 5,000 | 871.28 | 4402.8 | **5.05×** |
| `fibula_r` | 1,000 | 656.87 | 16596.5 | **25.27×** |
| `T6` | 10,000 | 777.24 | 1716.8 | **2.21×** |

> 二值判据（≥σ_c）不因此改变（它们本就远超阈值），但**利用率/相对排序被扭曲**；这是矩阵不可信的独立原因。

---

## 6. 结论

1. **相对排序可信度：不可信。** 三条路线（BC-free 1D / BC-健壮 FE / 现有矩阵）两两 Spearman 仅 0.22–0.77，且 `calcaneus/C5/femur` 位置互换；矩阵 p95 排序与 BC-free 路线弱相关（0.22）。
2. **端部夹持假设被证伪**：剔除端带后 p95 **升高**（中段弯曲主导），说明 `MESH_FIX_REPORT` 的残余伪影主要是**斜向加载 → 中段弯曲**，而非端部奇异。
3. **S4 输入建议**：改用 **Route 1 的 1D 名义 `risk(h)`**（BC-free、9 骨区分度好）；或先修 FE（沿骨轴加载 + 修正折算）再谈用 FE 指标。**当前矩阵的 p95/二值结果不应直接进入 S4。**
4. **矩阵必须重做**：修正降载折算（§5）、修正载荷方向（沿骨轴而非 bbox 轴）、并用 BC-健壮统计量（避免 p95）。

---

## 7. 诚实边界

1. **Route 1 未建模载荷分配**：同一关节的多根骨共享整体关节反力，未按刚度分流。胫/腓共用 `ankle_r` 的 25.3 kN，但腓骨实际仅分担约 10%（`bone.py` 注：单侧 5%）——故 `fibula_r` 的 Route-1 风险是**上界**，其"最危险"含此人为放大。
2. **Route 1 是 1D 压缩**：忽略弯曲/剪切/偏心；对细长骨（胫/腓）会低估风险，对紧凑骨（跟骨）会高估（`bone.py` 明确记载材料强度×截面会高估跟骨 2–3×）。
3. **中轴骨 A_section 是 estimate**：THUMS 的中轴骨 CORT 是壳（0 实体），用 SPON/diploë 实体 V/L 近似承载截面，与真实皮质截面可能差 2× 量级。
4. **Route 2 端带剔除的副作用**：剔除 60% 长度后中段元素变少（尤其薄壳/小网格），mid 统计量样本变小；对小域（部分髋壳）会退化到全域（脚本已记录 `mid_n`）。
5. **`fibula_r` 无 Route-2 值**：干净 remesh 全载荷发散；原始 mesh 的 FE 也不收敛（矩阵靠 1000 N 折回）。
6. **原始 mesh Route-2 是"重分析"**：用的是产生矩阵的同一（较粗、可能退化的）网格，只做了端带剔除与正确折算；它验证的是"BC/折算"而非"网格质量"。
7. **本报告未改任何既有模块/脚本/默认值**，未覆盖任何既有产物；未重跑 FEBio（4 长骨+3 轴向诊断全部复用现有 hdf5）。

---

## 8. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\bc_robust_metric.py
# 产物：
#   results/opensim_fe/bc_robust_metric_route1.json
#   results/opensim_fe/bc_robust_metric_route2.json
#   results/opensim_fe/bc_robust_metric_result.json
```

数据来源（只读复用，未修改）：
- `results/opensim_fe/fracture_matrix.json`（`loads_n`、`sigma_vm_ref_mpa`、`first_fracture_*`）
- `temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5` + `<b>_fe.json`（原始 9 骨 per-element 应力 / 主域 / 折回因子）
- `temp/opensim_fe/fracture_matrix_remesh/<b>/*.hdf5`（干净 remesh：calc/tibia/femur + 3 轴向诊断；fibula 发散）
- `src/climbing/bone.py::MATERIAL_STRENGTH_MPA`

---

## 9. 回归（保持绿）

```powershell
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
# 实测：19 passed, 1 skipped, 1 warning in 96.19s（warning = .pytest_cache 目录权限，与本任务无关）
```

---

## 10. 产物

- 本报告：`results/opensim_fe/BC_ROBUST_METRIC_REPORT.md`
- 新脚本：`scripts/opensim_fe/bc_robust_metric.py`（additive，不改既有模块）
- 机器可读：`results/opensim_fe/bc_robust_metric_route1.json`、`bc_robust_metric_route2.json`、`bc_robust_metric_result.json`
