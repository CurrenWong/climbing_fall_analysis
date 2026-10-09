# B13_DIFF — 新 BC（align_pca, median） vs 旧矩阵

## 1. 首次骨折高度对照（m；None=≤50m 未折）

| 部位 | 旧 p95 | 旧 median | **新 median(BC2)** | Δ(新med−旧med) |
|---|---|---|---|---|
| 足部 (`calcaneus_r`) | 1 | 17 | **13** | -4 |
| 胫骨 (`tibia_r`) | 1 | 6 | **42** | 36 |
| 腓骨 (`fibula_r`) | 1 | 1 | **1** | 0 |
| 股骨 (`femur_r`) | 1 | None | **47** | None |
| 骨盆 (`R_HIPBONE`) | 8 | None | **None** | None |
| 腰椎 (`L3`) | 1 | 16 | **7** | -9 |
| 胸椎 (`T6`) | 1 | 1 | **1** | 0 |
| 颈椎 (`C5`) | 1 | 7 | **7** | 0 |
| 颅骨 (`parietal_r`) | 1 | 4 | **6** | 2 |

## 2. σ_vm_ref 排序 Spearman ρ（≥3 对）

- 新median vs **旧median**（同统计量，纯 BC 效应）：ρ = **0.700** （n=9）
- 新median vs **旧p95**（旧生产口径）：ρ = **0.433** （n=9）

## 3. 首折一致性

- 新median vs 旧p95 完全同值：2/9
- 首折高度（双方均非 None，n=7）Spearman ρ = **0.771**

_生成：rerun_matrix_bc2.py --diff；数据源为两份 JSON 实测。_
