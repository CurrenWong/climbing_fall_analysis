# 载荷折算修复报告 —— fracture_matrix 降载骨 σ_vm_ref 未乘 load_reduction_factor

> 任务：修复 `scripts/opensim_fe/fracture_matrix.py` 中 **降载骨（FEBio 发散后降载求解）的 FE 应力未折算回参考载荷** 的已证实缺陷，并在**不重跑 FEBio / 不重跑 OpenSim** 的前提下，从缓存记录重放 9×50 矩阵。
> 仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13.0；生成时间 2026-10-04。
> 前置：`results/opensim_fe/BC_ROBUST_METRIC_REPORT.md`（§5 首次报告该缺陷）、`results/opensim_fe/MESH_FIX_REPORT.md`（σ_vm p95 实为边界条件受限）。
> 本报告除标注外全部为本机实测；命令原文见 §3。

---

## 0. 一句话结论

`compute_matrix` 把**降载求解得到的 σ_vm 直接当作参考载荷（ref-load）应力**，漏乘 `load_reduction_factor`（lrf），导致同一张矩阵混入**不同载荷下**的应力。修复后（`sigma_vm_ref[s] = solved_stats[s] × lrf`，`sig[s](h) = solved_stats[s] × lrf × load_h/load_ref`），3 根降载骨的 `sigma_vm_ref`/利用率按各自 lrf 放大，**数值与代码/报告口径一致**。

- 本修复**只修正载荷账目（load bookkeeping）**：`tibia_r ×5.0532`、`fibula_r ×25.266`、`T6 ×2.2088`；其余 6 骨 lrf=1.0，**逐字段不变**。
- 本修复**不改变相对骨折排序的结论**：p95 由边界条件伪影（斜置骨沿 bbox 轴加载 → 中段弯曲 + 端部夹持）主导，**排序仍不可信**（见 `BC_ROBUST_METRIC_REPORT.md`）；二值矩阵仍退化为「8/9 在 1 m 骨折」。
- 核心数值验证：`tibia_r sigma_vm_ref.p95` = 871.285 × 5.0532 = **4402.777 MPa ≈ 4402.8 MPa** ✅（与任务要求一致）。

---

## 1. 缺陷与数学

`extract_main_domain`（第 492–530 行）返回的是 **FE 实际求解载荷 `fe_solve_load_n` 下的原始 σ_vm 统计量，且不做任何缩放**。当该骨在 @5 m 参考载荷下发散、被迫降载求解时，`fe_solve_load_n < load_ref_n`，`load_reduction_factor = load_ref_n / fe_solve_load_n > 1`。

原 `compute_matrix`（修复前第 570、576 行）直接用 `main["stats"]` 当 ref-load 应力：

```
row["sigma_vm_ref_mpa"] = {s: main["stats"][s] ...}      # 缺 × lrf
sig[s][pi]             = main["stats"][s] * ratio        # ratio = load_h/ref_load，仍缺 × lrf
row["utilization_at_5m"] = main["stats"][primary] / σ_c  # 缺 × lrf
```

**正确数学**（严格线弹性）：

```
sigma_vm_ref[s] = solved_stats[s] × lrf
sig[s](h)       = solved_stats[s] × lrf × (load_h / ref_load)
                = solved_stats[s] × (load_h / fe_solve_load_n)
```

其中 `ratio = load_h/ref_load` 已存在，只需再乘 `lrf` 即可修正（因为 `lrf × ratio = load_h/fe_solve_load_n`）。

受影响骨（3 根，均在本机缓存中状态 `ok`）：

| 骨 | `fe_solve_load_n` | `load_ref_n` | **lrf** |
|---|---:|---:|---:|
| `tibia_r` | 5,000 | 25,266.249 | **5.0532** |
| `fibula_r` | 1,000 | 25,266.249 | **25.266** |
| `T6` | 10,000 | 22,087.580 | **2.2088** |
| 其余 6 骨 | = ref | = ref | 1.0000 |

---

## 2. 代码改动（最小、局部）

仅改 `scripts/opensim_fe/fracture_matrix.py::compute_matrix`（第 570–592、596–597 行），**不改默认值、CLI、主统计量（仍为 p95）、其他脚本**：

```python
# FE-solved stats live at ``fe_solve_load_n`` (<= ref load whenever FEBio
# had to reduce the load).  Scale them back to the reference load with the
# load-reduction factor before using them as "ref-load" stress; otherwise
# the matrix mixes stresses computed at different loads.
lrf = float(rec.get("load_reduction_factor") or 1.0)
row["sigma_vm_ref_mpa"] = {s: main["stats"][s] * lrf for s in STATS}
row["sigma_vm_ref"] = main["stats"][primary_stat] * lrf
# raw per-domain stats are left as solved (at fe_solve_load_n); the
# ref-load-scaled values are in ``sigma_vm_ref_mpa`` above.
row["main_domain_detail"] = main.get("domains", {})

ratio = load_h / ref_load if ref_load > 0 else np.full_like(load_h, np.nan)
for s in STATS:
    sig[s][pi] = main["stats"][s] * lrf * ratio        # ← 加 lrf
    frac[s][pi] = (sig[s][pi] >= part.sigma_c_mpa).astype(int)
    ...
row["utilization_at_5m"] = (main["stats"][primary_stat] * lrf
                            / part.sigma_c_mpa)         # ← 加 lrf
```

- `extract_main_domain` **未改**（按任务要求，缩放只放 `compute_matrix`）。
- `main_domain_detail` 保留为**求解载荷下的原始逐域统计量**（未缩放），供审计；ref-load 口径见 `sigma_vm_ref_mpa`。
- 循环内日志也改为打印缩放后值（同一缩放操作），保持日志与 JSON 一致。
- `frac` 由 `sig` 派生，自动随之更新；`load_scale_1_to_50` 等字段不变。

---

## 3. 重放方式与命令（**不重跑 FEBio / OpenSim**）

正常入口 `fracture_matrix.main` 虽复用缓存 FE，但仍会调用 `find_febio(verify=True)` + `probe_febio` 去**启动 FEBio 可执行文件**取版本号。为严格遵守「不调用 FEBio」，新增**纯缓存重放驱动** `scripts/opensim_fe/regenerate_fracture_matrix.py`（additive，不改既有 CLI 语义）：

- 载荷：直接复用缓存 `temp/opensim_fe/fracture_matrix/load_chain.npz`（逐高度实测关节纵向峰值，**不重跑 OpenSim**）。
- FE：直接读缓存 `temp/opensim_fe/fracture_matrix/<bone>/<bone>_fe.json`（**不重跑 FEBio**）。
- 求解器元数据（`febio_exe` / `febio_version` / `time_steps`）从既有 `fracture_matrix.json` 的 `meta` 复用，**不启动 FEBio**。
- 随后调用 `compute_matrix(...)` → `write_outputs(...)` → `write_report(...)`，写出与主入口相同的产物。

**执行命令（仓库根）：**

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\regenerate_fracture_matrix.py
# 实测耗时 ~3.1 s；日志：temp/opensim_fe/fracture_matrix/regen_run.log
```

**载荷来源**：`load_chain.npz`（逐高度 `ground_reaction → run_dead_drop → joint_reaction` 的实测关节纵向峰值；`ref_idx`=h=5 m）。@5 m 参考载荷：`ankle_r`=25,266.25 N、`subtalar_r`=25,354.40 N、`hip_r`=15,227.65 N、`lumbar`=22,087.58 N。

---

## 4. 3 根受影响骨 before / after

前值 = `*_prebugfix_backup`；后值 = 修复并重放后。**因子 = 后/前 = lrf**。

| 骨 | lrf | `sigma_vm_ref.p95` before → after | `sigma_vm_ref.median` before → after | `utilization_at_5m` before → after | 首次骨折 p95 before → after |
|---|---:|---|---|---|---|
| `tibia_r` | **5.0532** | 871.29 → **4402.78** (×5.0532) | 38.09 → 192.48 (×5.0532) | 4.356 → **22.014** | 1 m → 1 m |
| `fibula_r` | **25.266** | 656.87 → **16596.45** (×25.266) | 80.44 → 2032.52 (×25.266) | 4.105 → **103.728** | 1 m → 1 m |
| `T6` | **2.2088** | 777.24 → **1716.77** (×2.2088) | 245.30 → 541.81 (×2.2088) | 5.182 → **11.445** | 1 m → 1 m |

补充：四统计量全部按同一 lrf 缩放（`max`：tibia 1352.90→6836.47、fibula 1061.45→26818.65、T6 1419.72→3135.87）。其余 6 骨（lrf=1.0）逐字段与备份**完全一致**（脚本比对：`only 3 bones changed`）。

**首次骨折高度（四口径）before → after：**

| 骨 | 口径 | max | p95 | mean | median |
|---|---|---|---|---|---|
| `tibia_r` | before | 1 | 1 | 7 | >50 |
| | after | 1 | 1 | **1** | **6** |
| `fibula_r` | before | 1 | 1 | 3 | 14 |
| | after | 1 | 1 | **1** | **1** |
| `T6` | before | 1 | 1 | 1 | 1 |
| | after | 1 | 1 | 1 | 1 |

> p95 主口径上 3 骨的首次骨折高度**不变**（本就远超阈值）——这与 `BC_ROBUST_METRIC_REPORT.md` §5 的预见一致；受影响的是 `mean`/`median` 与非二元利用率。

---

## 5. 全部 9 骨修正后首次骨折高度（主口径 p95 加粗）

| # | 部位 | 代表骨 | σ_c (MPa) | 修正后 σ_vm_ref **p95** (MPa) | 首次骨折 **p95** | max | mean | median |
|---|---|---|---|---:|---|---|---|---|
| 1 | 足部 | `calcaneus_r` | 150 | 257.154 | **1** | 1 | 11 | 17 |
| 2 | 胫骨 | `tibia_r` | 200 | **4402.777** | **1** | 1 | 1 | 6 |
| 3 | 腓骨 | `fibula_r` | 160 | **16596.454** | **1** | 1 | 1 | 1 |
| 4 | 股骨 | `femur_r` | 220 | 1168.840 | **1** | 1 | 4 | >50 |
| 5 | 骨盆 | `R_HIPBONE` | 180 | 145.678 | **8** | 2 | 44 | >50 |
| 6 | 腰椎 | `L3` | 150 | 303.144 | **1** | 1 | 10 | 16 |
| 7 | 胸椎 | `T6` | 150 | **1716.765** | **1** | 1 | 1 | 1 |
| 8 | 颈椎 | `C5` | 150 | 641.141 | **1** | 1 | 2 | 7 |
| 9 | 颅骨 | `parietal_r` | 160 | 423.414 | **1** | 1 | 3 | 4 |

> 主口径 p95：**9 骨中 8 骨在 1 m 即判骨折**，仅 `R_HIPBONE` 为 8 m。修正后仍高度退化 —— 说明「p95 + 材料强度」判据在单骨点载荷 FE 上严重过判（见 §6 与 `BC_ROBUST_METRIC_REPORT.md` §4）。

---

## 6. 本修复做了什么 / 没做什么

**做了什么（仅此一项）：**

- 修正**载荷账目**：降载骨的 `sigma_vm_ref`、`sig[s](h)`、`utilization_at_5m` 现在都乘上 `load_reduction_factor`，使矩阵内所有应力统一到**参考载荷 @5 m**口径，与代码注释/既有报告声称的「σ_vm_ref = 收敛载荷统计量 × 线性折回因子」一致。

**没做什么（明确不变）：**

- **不挽救相对骨折排序。** p95 的主控项仍是**边界条件伪影**（斜置骨沿 bbox 长轴加载 → 中段弯曲 + 端部 30% 夹持 + 点节点力奇异），该伪影骨形状相关：`BC_ROBUST_METRIC_REPORT.md` 实测端带剔除后 p95 **不降反升**（tibia 2224→4308），且矩阵 p95 排序与 BC-free 的 1D 名义排序 Spearman **ρ≈0.22**（弱）。本次修复**不改变**该结论。
- 不改变二值矩阵的退化图景（仍 8/9 在 1 m），不改变主统计量（仍 p95），不改变映射假设、σ_c、FE 边界、载荷方向。
- 不重跑 FEBio、不重跑 OpenSim、不换网格、不修 FE 的 BC。
- 不改 `extract_main_domain` 的行为（其返回值仍为求解载荷下的原始统计量）。
- 若要让 FE 指标可用作相对排序，仍需按 `BC_ROBUST_METRIC_REPORT.md` §6 的建议：**先修载荷方向（沿骨轴而非 bbox 轴）+ 采用 BC-健壮统计量**，或改用 Route-1 的 1D 名义 `risk(h)`。本修复**不替代**这些工作。

---

## 7. 产物：覆盖 vs 新增 vs 备份

| 文件 | 状态 | 说明 |
|---|---|---|
| `results/opensim_fe/fracture_matrix.json` | **被覆盖（重放）** | 修正后 9×50 矩阵 + 逐高度载荷/σ_vm/四统计量 |
| `results/opensim_fe/fracture_matrix.csv` | **被覆盖（重放）** | 宽表（行=部位，列=h1..h50，主口径 p95） |
| `results/opensim_fe/fracture_matrix.png` | **被覆盖（重放）** | 二值矩阵 + 利用率热图 |
| `results/opensim_fe/FRACTURE_MATRIX_REPORT.md` | **被覆盖（重放）** | §3/§4 表格现显示折算后数值（如 tibia 4402.78） |
| `results/opensim_fe/fracture_matrix_prebugfix_backup.json` | **新增（备份）** | 修复前 JSON 原件 |
| `results/opensim_fe/fracture_matrix_prebugfix_backup.csv` | **新增（备份）** | 修复前 CSV 原件 |
| `results/opensim_fe/fracture_matrix_prebugfix_backup.png` | **新增（备份）** | 修复前 PNG 原件 |
| `results/opensim_fe/FRACTURE_MATRIX_REPORT_prebugfix_backup.md` | **新增（备份）** | 修复前报告原件 |
| `results/opensim_fe/MATRIX_SCALING_FIX_REPORT.md` | **新增** | 本报告 |
| `scripts/opensim_fe/fracture_matrix.py` | **被编辑** | 仅 `compute_matrix` 的缩放（§2） |
| `scripts/opensim_fe/regenerate_fracture_matrix.py` | **新增** | 纯缓存重放驱动（无 FEBio） |
| `temp/opensim_fe/fracture_matrix/regen_run.log` | **新增** | 本次重放日志（未覆盖既有 `run.log`） |

> 备份在任何写入之前完成（`Copy-Item` 到 `*_prebugfix_backup.*`，同目录）。未触碰 `scripts/ankle_fe/`、`temp/pyfebio_demo/`、`temp/ankle_*`、`FE_PIPELINE_TODO.md`。

---

## 8. 回归（保持绿）

```powershell
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
# 实测：19 passed, 1 skipped, 1 warning in 43.66s
# （warning = .pytest_cache 目录权限拒绝访问，与本任务无关）
```

`tests/` 中**无**任何用例引用 `fracture_matrix` / `compute_matrix` / `sigma_vm_ref`，故本次改动不触及既有测试面；回归结果与任务要求一致（19 passed / 1 skipped）。

---

## 9. 诚实边界

1. **线性折算本身有误差**：3 根降载骨的双载荷自检实测 σ(2L)/σ(L)=**1.84–2.67**（严格线性=2.0，见 `FRACTURE_MATRIX_REPORT.md` §7 边界 10），故 lrf 折算带约 ±8–33% 量级误差。本次修复只保证「用了 lrf」，不消除 FE 的非线性。
2. **本修复不提高 FE 结果可信度**：p95 的 BC 伪影与排序不可信问题原样保留，必须配合 `BC_ROBUST_METRIC_REPORT.md` §6 的建议使用。
3. **`main_domain_detail` 未缩放**：JSON 的 `parts[].main_domain_detail` 仍为求解载荷下的原始逐域量；引用时须知其口径为 `fe_solve_load_n`，ref-load 口径请用 `sigma_vm_ref_mpa`。
4. **重放不代表重算**：矩阵完全来自既有缓存，未重跑 FEBio/OpenSim；若缓存被清理，须先恢复缓存或重跑。
5. **未改任何模块默认/签名**，未覆盖未备份的既有产物。
