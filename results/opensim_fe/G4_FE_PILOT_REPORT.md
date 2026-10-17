# G4 (2) P1 PILOT — fibula_r 逐高度 reactive-damage FE σ(h)（V7 BC）

> **单骨 pilot（非 9 部件矩阵、非冻结预注册运行）**。用 **`V7` BC**（顶/底 50% 带 `NodalForce` 沿骨 PCA 长轴）对 `fibula_r` 在**每个高度**`load(h)=ankle_r__vert(h)×load_fraction` 上**真正重解损伤 FE**，替代冻结矩阵的 `σ_vm(h)=σ_ref·load(h)/load_ref` 线性外推。THUMS σ_y/E_t/c 与 Dmax 逐字保持。

_生成时间：2026-10-11 23:36:31；σ_c(fibula_shaft)=160 MPa；主统计=p95；参考 h=5 m，参考解 = V7 弹性 @1000 N。_

## 0. 结论（先看这里）

- **V7 弹性参考** @1000 N：σ_vm max=107.5 MPa，p95=79.0 MPa（交付 |R|=1000 N）。（对比冻结矩阵的点载荷基线 p95≈657 MPa @1000 N。）
- **物理分配 load_fraction=0.07**（文献腓骨份额 ~0.05–0.10）：@5 m load=1769 N；@50 m load=7905 N。
- **冻结分配 load_fraction=1.0**：@5 m load=25266 N；@50 m load=112934 N（全部远超模型承载，见下）。
- **可达高度（物理分配，逐高度损伤重解）**：D0.9:h∈[3, 4, 11, 14, 35]；D1:h∈[6, 7]。最大可达 h=35。**收敛是零散的** —— 多数高度以负 Jacobian 发散（损伤软化在点载荷量级下失去正定切线），即 **损伤-FE 在操作载荷量级根本无法给出 σ(h)**；这本身相对线性外推就是关键差异。
- **首次断裂高度**：损伤-FE（σ_vm max/p95 ≥ σ_c） = D0.9/p95:h=14, D0.9/max:h=14；线性基线（V7-弹性外推，σ_lin = σ_ref·load/load_ref）= max:h=4, p95:h=7。
- 冻结矩阵（点载荷 BC，σ_ref×25.27 线性外推）首次断裂 = h=1 m（σ_vm_ref p95=16596 MPa ≫ σ_c）。

## 1. σ_vm(h) 表：损伤-FE vs 线性基线

`σ_lin(h) = σ_ref · load(h)/load_ref`（σ_ref = V7 弹性 @1000 N，主统计=p95；等价于 σ_ref/1000 · load(h)）。损伤-FE = 在 load(h) 上重解 `elastic damage`+`CDF Simo`+`DC Drucker` 的 σ_vm。

| h (m) | ankle_r (N) | load_frozen=1.0 (N) | load_phys=0.07 (N) | σ_vm max (D0.9) | σ_vm p95 (D0.9) | σ_lin p95 | ratio p95 | 断裂 | 状态 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 15454 | 15454 | 1082 | — | — | 85.5 | — |  | diverged |
| 2 | 20239 | 20239 | 1417 | — | — | 112.0 | — |  | diverged |
| 3 | 21150 | 21150 | 1481 | 67.2 | 49.4 | 117.0 | 0.422 | ok | converged |
| 4 | 22790 | 22790 | 1595 | 136.9 | 87.8 | 126.1 | 0.696 | ok | converged |
| 5 | 25266 | 25266 | 1769 | — | — | 139.8 | — |  | diverged |
| 6 | 28218 | 28218 | 1975 | — | — | 156.1 | — |  | diverged |
| 7 | 31355 | 31355 | 2195 | — | — | 173.4 | — |  | diverged |
| 8 | 34524 | 34524 | 2417 | — | — | 191.0 | — |  | diverged |
| 9 | 37655 | 37655 | 2636 | — | — | 208.3 | — |  | diverged |
| 10 | 40644 | 40644 | 2845 | — | — | 224.8 | — |  | diverged |
| 11 | 43520 | 43520 | 3046 | 93.5 | 67.8 | 240.7 | 0.282 | ok | converged |
| 12 | 46305 | 46305 | 3241 | — | — | 256.1 | — |  | diverged |
| 13 | 48932 | 48932 | 3425 | — | — | 270.7 | — |  | diverged |
| 14 | 51448 | 51448 | 3601 | 261.7 | 171.4 | 284.6 | 0.602 | **FRACT** | converged |
| 15 | 53819 | 53819 | 3767 | — | — | 297.7 | — |  | diverged |
| 16 | 56091 | 56091 | 3926 | — | — | 310.3 | — |  | diverged |
| 17 | 58293 | 58293 | 4080 | — | — | 322.5 | — |  | diverged |
| 18 | 60443 | 60443 | 4231 | — | — | 334.3 | — |  | diverged |
| 19 | 62547 | 62547 | 4378 | — | — | 346.0 | — |  | diverged |
| 20 | 64606 | 64606 | 4522 | — | — | 357.4 | — |  | diverged |
| 21 | 66623 | 66623 | 4664 | — | — | 368.5 | — |  | diverged |
| 22 | 68603 | 68603 | 4802 | — | — | 379.5 | — |  | diverged |
| 23 | 70547 | 70547 | 4938 | — | — | 390.2 | — |  | diverged |
| 24 | 72458 | 72458 | 5072 | — | — | 400.8 | — |  | diverged |
| 25 | 74335 | 74335 | 5203 | — | — | 411.2 | — |  | diverged |
| 26 | 76182 | 76182 | 5333 | — | — | 421.4 | — |  | diverged |
| 27 | 78001 | 78001 | 5460 | — | — | 431.5 | — |  | diverged |
| 28 | 79791 | 79791 | 5585 | — | — | 441.4 | — |  | diverged |
| 29 | 81556 | 81556 | 5709 | — | — | 451.1 | — |  | diverged |
| 30 | 83297 | 83297 | 5831 | — | — | 460.8 | — |  | diverged |
| 31 | 85015 | 85015 | 5951 | — | — | 470.3 | — |  | diverged |
| 32 | 86709 | 86709 | 6070 | — | — | 479.6 | — |  | diverged |
| 33 | 88379 | 88379 | 6187 | — | — | 488.9 | — |  | diverged |
| 34 | 90030 | 90030 | 6302 | — | — | 498.0 | — |  | diverged |
| 35 | 91630 | 91630 | 6414 | 170.7 | 112.2 | 506.9 | 0.221 | ok | converged |
| 36 | 93147 | 93147 | 6520 | — | — | 515.3 | — |  | diverged |
| 37 | 94650 | 94650 | 6625 | — | — | 523.6 | — |  | diverged |
| 38 | 96137 | 96137 | 6730 | — | — | 531.8 | — |  | diverged |
| 39 | 97609 | 97609 | 6833 | — | — | 539.9 | — |  | diverged |
| 40 | 99067 | 99067 | 6935 | — | — | 548.0 | — |  | diverged |
| 41 | 100510 | 100510 | 7036 | — | — | 556.0 | — |  | diverged |
| 42 | 101942 | 101942 | 7136 | — | — | 563.9 | — |  | diverged |
| 43 | 103358 | 103358 | 7235 | — | — | 571.7 | — |  | diverged |
| 44 | 104761 | 104761 | 7333 | — | — | 579.5 | — |  | diverged |
| 45 | 106155 | 106155 | 7431 | — | — | 587.2 | — |  | diverged |
| 46 | 107532 | 107532 | 7527 | — | — | 594.8 | — |  | diverged |
| 47 | 108902 | 108902 | 7623 | — | — | 602.4 | — |  | diverged |
| 48 | 110258 | 110258 | 7718 | — | — | 609.9 | — |  | diverged |
| 49 | 111600 | 111600 | 7812 | — | — | 617.3 | — |  | diverged |
| 50 | 112934 | 112934 | 7905 | — | — | 624.7 | — |  | diverged |

> σ_vm 与 σ_lin 均为 fib_cort 域统计；`—` = 该格发散（无解）。

## 2. 损伤-FE 曲线是否非线性？（P1 payoff）

- 在**可达**的 5 个高度上，比值 σ_vm(damage)/σ_lin ∈ [0.221, 0.696]。
- **是非线性的**：若损伤-FE 与线性基线一致，比值应恒为 1.000；实测比值跨 0.475（≈48%），说明损伤软化改变了 σ 对载荷的比例关系 —— **这就是 P1 相对线性外推的差异**。

## 3. 可达性 / 非确定性 / 发散格

- **物理分配最大可达高度 = h=35**；可达格零散分布（见 §0），其余为负 Jacobian 发散（原文见下）。
- 发散格数（全部 frac/Dmax 组合）：**143/150**。
- **顶高度重复（非确定性）**：

| h | Dmax | load (N) | 收敛 | σ_p95 | 错误 |
|---|---|---|---|---|---|
| 35 | 0.9 | 6414 | yes | 112.2 |  |
| 35 | 0.9 | 6414 | yes | 112.2 |  |
| 35 | 0.9 | 6414 | yes | 112.2 |  |
| 35 | 1 | 6414 | no | — | Negative jacobian detected. | failed to converge a |
| 35 | 1 | 6414 | no | — | Negative jacobian detected. | failed to converge a |
| 35 | 1 | 6414 | no | — | Negative jacobian detected. | failed to converge a |

**发散原文（示例）**：

```
h=1 phys D0.9 load=1082N : 8 negative jacobians detected. | failed to converge at time : 1
h=2 phys D0.9 load=1417N : 96 negative jacobians detected. | 192 negative jacobians detected. | 30 negative jacobians detected.
h=5 phys D0.9 load=1769N : 49 negative jacobians detected. | failed to converge at time : 1
h=6 phys D0.9 load=1975N : 93 negative jacobians detected. | failed to converge at time : 1
h=7 phys D0.9 load=2195N : 16 negative jacobians detected. | failed to converge at time : 1
h=8 phys D0.9 load=2417N : 71 negative jacobians detected. | 77 negative jacobians detected. | 68 negative jacobians detected.
```

## 4. 复现与产物

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_fe_pilot.py
```

- `results/opensim_fe/G4_FE_PILOT_REPORT.md`（本文件）
- `temp/fe_out/g4_fe_pilot/`（各 `.feb`/`.xplt`/`.hdf5`/`.log` + `g4_fe_pilot.json`）
- `scripts/opensim_fe/g4_fe_pilot.py`（本脚本；只 import 既有模块）

**未改任何既有模块**：`g4_fe_p1b` / `g4_fe_spike` / `g4_fe_stab` / `bone_feb` 均只 import 复用。

