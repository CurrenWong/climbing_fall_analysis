# 1D 骨折风险 —— 弯曲感知修订（梁理论，Euler–Bernoulli）

> 生成时间 2026-10-04 22:19:52；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s。
> 本报告**新增**，不改任何既有模块/默认值，不覆盖既有产物，不调用 FEBio。
> 修订对象：`RISK_1D_LOADSHARE_REPORT.md` 的纯轴向 1D `risk(h)=F_axial(h)/(A_section·σ_c)`，补上其 §6 明确缺失的**弯曲项**。

---

## 0. 一句话结论

把关节反力沿**骨轴斜置角 θ** 分解后补上梁弯曲项 `σ = F_axial/A + M·c/I`（`M = F_perp·L/β`，主口径 β=4）：

- **纯轴向下不骨折**（≤4.5 m 有 0 根骨折；最早 8.0 m）。
- **加上弯曲后**，≤4.5 m 有 **4 根**骨折，最早 **1.0 m** —— 与 FE 矩阵的 「8/9 在 1 m 骨折」落在**同一区间**。
- **缺口被闭合**：弯曲项把 1D 从「严重低估」推到「与 FE 同量级甚至更高」，证明此前 8× 缺口的主因**就是被忽略的弯曲**（与 `MESH_FIX_REPORT.md` §5 的诊断一致）。
- **代价是可能过冲**：弯曲是**上界式**估计（等效圆截面 + 全骨长 + 简单支承），对细长骨（腓骨/胫骨）给出的风险**高于** FE；需按 §6 诚实边界解读。

---

## 1. 方法

### 1.1 梁公式

```
σ(h) = F_axial/A_section + M_max·c/I
F_axial = F(h)·cosθ      F_perp = F(h)·sinθ
M_max   = F_perp·L/β     β = 4（简支，主口径）；敏感性 β ∈ {2,4,8}
risk(h) = σ(h)/σ_c       first_fracture = 最小的 h 使 risk ≥ 1（否则 >50 m）
```

| 量 | 来源 | 性质 |
|---|---|---|
| `F(h)` | `risk_1d_loadshare.json` 的 `force_multiplier` × `fracture_matrix.json.loads_n` | measured（逐高度）+ 载荷分配 |
| `A_section` | `bc_robust_metric_route1.json` | 4 长骨 measured CORT；5 中轴骨 estimated body |
| `L` | `bc_robust_metric_route1.json` `L_bbox_mm` | measured（bbox 最长轴） |
| `σ_c` | `src/climbing/bone.py` `MATERIAL_STRENGTH_MPA`（压缩） | assumed（[Y25] 表1） |
| `θ` | 4 长骨：`MESH_FIX_REPORT.md` §5 measured；其余 5 骨：FE 网格 PCA 计算 | 见表 |

### 1.2 截面模型

**DEFAULT —— 等效圆截面**（由 `A_section` 反推）：

```
r = sqrt(A_section/π) ;  I = π·r⁴/4 ;  c = r
```

**替代 —— 真实中段皮质截面**（4 长骨，若可行）：从干净 CORT remesh `temp/opensim_fe/fracture_matrix_remesh/meshes/<b>_cort/<b>_cort_mesh.npz` 取 PCA 长轴中点 ±slab 内的四面体，以**体积加权质心矩**求截面 A、I₁/I₂、极值纤维距离 c，取弱轴 I 与其对应 c（偏保守）。方法同为 Riemann 切片，对皮质壳网格采样的是**皮质环**。

### 1.3 θ 的取法

`θ = 骨 PCA 长轴 与 载荷轴（= FE bbox 最长轴，即关节反力施加方向）的夹角`。

| 骨 | θ (°) | 来源 | PCA 长轴 | bbox 跨度 (mm) | 载荷轴 |
|---|---:|---|---|---|---|
| `fibula_r` | **42.8** | MESH_FIX_REPORT.md sec.5 (measured) | [-0.657, -0.169, -0.734] | [275.4, 81.1, 304.3] | Z |
| `tibia_r` | **42.4** | MESH_FIX_REPORT.md sec.5 (measured) | [-0.664, -0.136, -0.735] | [299.4, 115.5, 327.6] | Z |
| `calcaneus_r` | **25.0** | MESH_FIX_REPORT.md sec.5 (measured) | [-0.362, 0.106, 0.926] | [71.5, 51.0, 84.6] | Z |
| `femur_r` | **23.5** | MESH_FIX_REPORT.md sec.5 (measured) | [-0.919, 0.080, 0.385] | [452.0, 118.4, 222.0] | X |
| `L3` | **90.0** | computed from FE mesh node PCA | [-0.718, -0.000, 0.696] | [66.2, 73.5, 72.3] | Y |
| `T6` | **8.0** | computed from FE mesh node PCA | [-0.990, -0.000, 0.139] | [58.4, 43.6, 51.5] | X |
| `R_HIPBONE` | **24.7** | computed from FE mesh node PCA | [-0.909, -0.206, -0.363] | [205.6, 139.5, 162.2] | X |
| `C5` | **0.0** | computed from FE mesh node PCA | [0.000, 1.000, 0.000] | [43.7, 54.4, 19.9] | Y |
| `parietal_r` | **10.6** | computed from FE mesh node PCA | [-0.983, 0.178, 0.045] | [136.0, 70.1, 105.7] | X |

**方法验证（4 长骨：FE 网格 PCA 计算值 vs §5 measured）**：

| 骨 | §5 measured (°) | 本脚本网格计算 (°) | 差 (°) |
|---|---:|---:|---:|
| `calcaneus_r` | 25.0 | 22.2 | 2.8 |
| `tibia_r` | 42.4 | 42.7 | 0.3 |
| `fibula_r` | 42.8 | 42.7 | 0.1 |
| `femur_r` | 23.5 | 23.2 | 0.3 |

> 计算值与 §5 一致到 1–3°，说明「PCA 长轴 vs bbox 载荷轴」的 θ 定义可复现；4 长骨按任务要求直接采用 §5 的 measured 值。
> ⚠️ 对**近各向同性**的中轴骨（C5、L3），PCA 长轴本身病态，其 θ（0° / 90°）是数值伪影，不代表真实骨-载荷夹角 —— 见 §6 边界。

---

## 2. 主口径逐骨表（β=4，等效圆截面，按 risk@5m 排序）

| # | 骨 | θ (°) | A_section (mm²) | r (mm) | I (mm⁴) | M_max@5m (N·mm) | σ_axial@5m (MPa) | σ_bend@5m (MPa) | σ_total (MPa) | σ_c (MPa) | **risk@5m** | **首次骨折 h (m)** |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `fibula_r` | 42.8 | 101.94 | 5.70 | 827 | 303,399 | 42.3 | 2090.0 | 2132.3 | 160 | **13.33** | **1** |
| 2 | `tibia_r` | 42.4 | 362.65 | 10.74 | 10466 | 1,071,208 | 39.5 | 1099.7 | 1139.2 | 200 | **5.70** | **1** |
| 3 | `calcaneus_r` | 25.0 | 216.01 | 8.29 | 3713 | 226,687 | 106.4 | 506.2 | 612.6 | 150 | **4.08** | **1** |
| 4 | `femur_r` | 23.5 | 404.87 | 11.35 | 13044 | 686,104 | 34.5 | 597.1 | 631.6 | 220 | **2.87** | **1** |
| 5 | `L3` | 90.0 | 692.60 | 14.85 | 38173 | 363,365 | 0.0 | 141.3 | 141.3 | 150 | **0.94** | **6** |
| 6 | `T6` | 8.0 | 291.20 | 9.63 | 6748 | 23,702 | 39.7 | 33.8 | 73.6 | 150 | **0.49** | **18** |
| 7 | `R_HIPBONE` | 24.7 | 1947.68 | 24.90 | 301873 | 326,617 | 7.1 | 26.9 | 34.0 | 180 | **0.19** | **>50** |
| 8 | `C5` | 0.0 | 178.83 | 7.54 | 2545 | 0 | 18.6 | 0.0 | 18.6 | 150 | **0.12** | **>50** |
| 9 | `parietal_r` | 10.6 | 877.53 | 16.71 | 61280 | 13,070 | 2.3 | 3.6 | 5.9 | 160 | **0.04** | **>50** |

**修订排序（β=4，早骨折优先，其次 risk@5m）**：

`fibula_r(#1, 1m) > tibia_r(#2, 1m) > calcaneus_r(#3, 1m) > femur_r(#4, 1m) > L3(#5, 6m) > T6(#6, 18m) > R_HIPBONE(#7, >50m) > C5(#8, >50m) > parietal_r(#9, >50m)`

> `σ_bend` 全面大于 `σ_axial`（细长骨尤甚），说明**弯曲项主导**了本次修订 —— 这正是纯轴向 1D 低估的根源。

---

## 3. 对比表：纯轴向 vs 轴向+弯曲 vs FE 矩阵

| 骨 | 纯轴向 risk@5m | 纯轴向 首次骨折 | **+弯曲 risk@5m** | **+弯曲 首次骨折** | FE util@5m | FE 首次骨折 p95 |
|---|---:|---:|---:|---:|---:|---:|
| `fibula_r` | 0.360 | 23 | **13.33** | **1** | 103.73 | 1 |
| `tibia_r` | 0.267 | 37 | **5.70** | **1** | 22.01 | 1 |
| `calcaneus_r` | 0.783 | 8 | **4.08** | **1** | 1.71 | 1 |
| `femur_r` | 0.171 | >50 | **2.87** | **1** | 5.31 | 1 |
| `L3` | 0.190 | >50 | **0.94** | **6** | 2.02 | 1 |
| `T6` | 0.268 | 47 | **0.49** | **18** | 11.45 | 1 |
| `R_HIPBONE` | 0.043 | >50 | **0.19** | **>50** | 0.81 | 8 |
| `C5` | 0.124 | >50 | **0.12** | **>50** | 4.27 | 1 |
| `parietal_r` | 0.015 | >50 | **0.04** | **>50** | 2.65 | 1 |

排序对照：

| 方法 | 排序（1→9） |
|---|---|
| **纯轴向（载荷分配）** | `calcaneus_r > fibula_r > tibia_r > T6 > L3 > femur_r > C5 > R_HIPBONE > parietal_r` |
| **轴向+弯曲（β=4）** | `fibula_r > tibia_r > calcaneus_r > femur_r > L3 > T6 > R_HIPBONE > C5 > parietal_r` |
| **FE 矩阵（util@5m）** | `fibula_r > tibia_r > T6 > femur_r > C5 > parietal_r > L3 > calcaneus_r > R_HIPBONE` |

Spearman 秩相关（9 骨）：

| 对比 | ρ | 解读 |
|---|---|---|
| 轴向+弯曲 vs 纯轴向 | 0.800 | 弯曲项改变了相对排序 |
| 轴向+弯曲 vs FE util@5m | **0.500** | （FE util 本身受 BC 伪影污染，仅作交叉参照） |
| 纯轴向 vs FE util@5m | 0.367 | 原有缺口 |

### 3.1 缺口是否闭合

| 口径 | ≤4.5 m 骨折骨数 | 最早首次骨折 | 与 FE 的关系 |
|---|---:|---|---|
| 纯轴向 1D | 0 | 8 m | **低估**（FE 8/9@1m） |
| **轴向+弯曲 1D** | **4** | **1 m** | **同量级 / 略过冲** |
| FE 矩阵 | — | 1 m（8/9） | — |

**结论**：

1. 纯轴向 1D 的「≤4.5 m 无骨折」是**结构性的**：它完全忽略偏心弯矩，而真实骨轴与竖直关节反力夹角达 24–43°（4 长骨），弯矩不可忽略。
2. 补上梁弯曲后，4 根长骨的**首骨折高度从 ≥8 m 骤降到 1 m**，与 FE 矩阵的 1 m 落在同一量级 → **~8× 的缺口被弯曲项解释并闭合**。
3. 由于 FE 的 p95 本身被**边界条件伪影**（斜置骨沿 bbox 轴加载 → 中段弯曲）主导（`BC_ROBUST_METRIC_REPORT.md` §3、`MESH_FIX_REPORT.md` §5），「闭合」不能被解读为「FE 是金标准」；更准确的说法是：**两条路线在'弯曲主导'这一点上收敛**。
4. 弯曲 1D 与 FE 的**相对排序**仍只有中等一致（ρ≈0.50）—— 因为 FE 的 p95 排序被骨形状相关的 BC 伪影塑形。

---

## 4. 截面模型对比（等效圆 vs 真实中段皮质截面，4 长骨）

| 骨 | 等效圆 A (mm²) | 等效圆 I (mm⁴) | 真实中段 A (mm²) | 真实 I_weak (mm⁴) | 真实 c_weak (mm) | 真实 S_weak (mm³) | 等效圆 S (mm³) | 真实/等效 S |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `calcaneus_r` | 216.0 | 3713 | 201.5 | 29758 | 23.56 | 1263 | 448 | 2.82 |
| `tibia_r` | 362.6 | 10466 | 277.9 | 8951 | 11.50 | 778 | 974 | 0.80 |
| `fibula_r` | 101.9 | 827 | 82.6 | 682 | 6.32 | 108 | 145 | 0.74 |
| `femur_r` | 404.9 | 13044 | 362.6 | 21825 | 14.94 | 1460 | 1149 | 1.27 |

**真实截面口径的逐骨结果（β=4，弱轴）：**

| 骨 | 真实截面 risk@5m | 真实截面 首次骨折 | 等效圆 risk@5m | 等效圆 首次骨折 |
|---|---:|---:|---:|---:|
| `calcaneus_r` | 1.91 | 1 | 4.08 | 1 |
| `tibia_r` | 7.08 | 1 | 5.70 | 1 |
| `fibula_r` | 17.83 | 1 | 13.33 | 1 |
| `femur_r` | 2.29 | 1 | 2.87 | 1 |

> 真实皮质中段的截面模量 S 与等效圆**同量级**：tibia 0.80×、fibula 0.74×、femur 1.27×（弱轴偏保守），说明**等效圆假设并非主要误差源**；主因是 β（弯矩分布）与载荷方向，而非截面。跟骨为紧凑块状（非线性长骨），其中段截面口径不适用（S 2.82×），该行仅作参考。

---

## 5. 敏感性（β ∈ {2,4,8}）

| 骨 | β=2 risk@5m | 首次 | β=4 risk@5m | 首次 | β=8 risk@5m | 首次 |
|---|---:|---:|---:|---:|---:|---:|
| `fibula_r` | 26.39 | 1 | 13.33 | 1 | 6.80 | 1 |
| `tibia_r` | 11.19 | 1 | 5.70 | 1 | 2.95 | 1 |
| `calcaneus_r` | 7.46 | 1 | 4.08 | 1 | 2.40 | 1 |
| `femur_r` | 5.59 | 1 | 2.87 | 1 | 1.51 | 2 |
| `L3` | 1.88 | 1 | 0.94 | 6 | 0.47 | 19 |
| `T6` | 0.72 | 10 | 0.49 | 18 | 0.38 | 27 |
| `R_HIPBONE` | 0.34 | 34 | 0.19 | >50 | 0.11 | >50 |
| `C5` | 0.12 | >50 | 0.12 | >50 | 0.12 | >50 |
| `parietal_r` | 0.06 | >50 | 0.04 | >50 | 0.03 | >50 |

> β 越小（弯矩越集中）风险越高。三档下 **4 根长骨首骨折高度都 ≤2 m**，结论对 β 不敏感（弯曲项主导的定性不变）；中轴骨的风险仍远低于 1。

---

## 6. 假设（ASSUMPTIONS）与诚实边界

**A. 梁模型**

- Euler–Bernoulli，`σ = F_axial/A + M·c/I`；**不含剪切、不含应力集中、不含轴-弯交互失效准则**（简单线性叠加，可能高估）。
- `M_max = F_perp·L/β`，主口径 β=4（简支中载）；这是**端部支承假设**，真实关节端更接近固定端（弯矩反而更小）→ 4 是中间偏保守。
- `L` = bbox 最长轴（measured），非骨的真实力学跨距；对斜置骨这是外接跨度上界。

**B. θ**

- 4 长骨用 `MESH_FIX_REPORT.md` §5 的 measured 值（calc 25.0、tibia 42.4、fibula 42.8、femur 23.5°）；与「从 FE 网格 PCA 重算」一致到 1–3°。
- 其余 5 骨由 FE 网格 node PCA 计算（`temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5`）。**C5（0°）与 L3（90°）是近各向同性骨的数值伪影**，其弯曲项不可信，中轴骨结论应以纯轴向为准。

**C. 截面**

- DEFAULT：等效圆（`A_section` → r, I, c）。4 长骨另给真实中段皮质截面（§4），两者同量级，故等效圆可用。
- `A_section`：4 长骨 measured CORT；5 中轴骨 estimated body（THUMS CORT 为壳），与真实皮质截面可能差 ~2×。
- 真实中段截面为**四面体体积加权的 Riemann 切片**，非严格几何截面；取弱轴 I 与其 c 偏保守。

**D. 诚实边界 / 未做什么**

- **弯曲是上界式估计**：全骨长 + 等效圆 + 简单支承共同把细长骨（腓/胫）推到 FE 之上；本文**如实报告过冲**，不粉饰。
- FE 矩阵的 p95 本身被 BC 伪影污染（`BC_ROBUST_METRIC_REPORT.md`），「与 FE 同量级」≠「FE 已证实」。
- `σ_c` 是材料强度而非整体骨失效载荷；绝对值不可信，可用的仍是**相对排序**。
- 不改任何既有模块默认/签名；不覆盖既有产物；不调用 FEBio；不重跑 OpenSim。

---

## 7. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\risk_1d_bending.py
# 只读复用：
#   results/opensim_fe/fracture_matrix.json          (loads_n / heights / FE util)
#   results/opensim_fe/risk_1d_loadshare.json        (force_multiplier, 纯轴向)
#   results/opensim_fe/bc_robust_metric_route1.json  (A_section / L / sigma_c)
#   src/climbing/bone.py MATERIAL_STRENGTH_MPA
#   temp/opensim_fe/fracture_matrix_remesh/meshes/<b>_cort/<b>_cort_mesh.npz
#   temp/opensim_fe/fracture_matrix/<b>/<b>.hdf5    (PCA)
# 产物：
#   results/opensim_fe/risk_1d_bending.json
#   results/opensim_fe/RISK_1D_BENDING_REPORT.md
```

---

## 8. 产物

- 本报告：`results/opensim_fe/RISK_1D_BENDING_REPORT.md`
- 新脚本：`scripts/opensim_fe/risk_1d_bending.py`（additive）
- 机器可读：`results/opensim_fe/risk_1d_bending.json`

