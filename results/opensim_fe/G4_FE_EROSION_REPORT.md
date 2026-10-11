# G4 P1-b — FEBio FEAMR 侵蚀自适应（erosion adaptor）能否稳定损伤骨模型（fibula_r）

> 决定性实验：把 FEBio 4.13 的 **FEAMR 侵蚀 mesh adaptor**（`<MeshAdaptor><mesh_adaptor type="erosion">`，删除损伤达到阈值的单元）注入 `fibula_r` 反应损伤模型（**固定步长**），看它能否让模型在**力控制下达到 1000 N**、**位移控制下突破 1088 N**，乃至操作载荷。**不调参**：σ_y/E_t/c 与 Dmax（0.9 / 1.0）逐字保持；只有 **erosion adaptor + 步进** 变化。

_生成时间：2026-10-11 22:31:13；FEBio `D:\Program\FEBioStudio\bin\febio4.exe` (probe=4.13.0)。fib_cort 单元数 = 842；σ_y=37.0 MPa，Simo a=21.362，D1(弹性标定)=198.9289 mm。_

## 0. FEBio 接受 MeshAdaptor？侵蚀是否触发？（先看这里）

**接受（ACCEPT）。** FEBio 4.13 完整解析了注入的 `<MeshAdaptor>` 段（`feamr.dll` 存在）；日志回显 adaptor 参数块（原文，来自一个实际运行）：

```
 MESH ADAPTOR DATA
===========================================================================

mesh adaptor 1 - Type: erosion
	max_iters ......................................... : 1
	max_elems ......................................... : 20
	sort .............................................. : 2
	remove_islands .................................... : yes (1)
	erode_surfaces .................................... : no (0)


 LOAD CONTROLLER DATA
===========================================================================


 OUTPUT DATA
===========================================================================
 LINEAR SOLVER DATA
===========================================================================
	Default linear solver ............................. : pardiso
```

**侵蚀已触发。** 在 2 个 cell 中删除了单元。单步 adaptor 应用原文：

```
=== Applying mesh adaptors: iteration 1
*mesh adaptor 1 (erosion):
	Deactivated elements: 8


 *************************************************************************
 * Selecting linear solver pardiso                                       *
 *************************************************************************
```

> **关键机制发现**：CDF Simo 损伤在 `Dmax` 处**从下方渐近**——即 D(t)→Dmax 但永不等于。实测一个 1×1×2 的 toy：`Dmax=0.9` 时 D 收敛到 **0.8988 < 0.9**。因此方案给定的阈值 `<min>0.9</min>` 与 `Dmax=0.9` **不相容**（阈值永远够不到）；只有 `Dmax=1.0`（D 可越过 0.9）或把阈值降到 0.9 以下（本文另测 `min=0.85`）时侵蚀才会真正删除单元。

## 1. 对照表：{FC 250/1000 N, DC δ*} × {erosion off vs on}

全部 cell = **固定步长**（已删除 pyfebio 自动 `<time_stepper>`）；on/off 之间**只差** 一个 `<MeshAdaptor>` 段。DC 的 δ 为 `-198.929 / -397.858 / -795.716 mm`（=1×/2×/4× 弹性标定 D1）。

| cell | ctrl | load/disp | Dmax | steps | erosion | thr | rc | wall(s) | steps done | **eroded** | **R_peak (N)** | σ_vm p95 | Dmax reached | FEBio error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `fc250_D0.9_off_s20` | FC | 250 N | 0.9 | 20 | off | — | 1 | 1.4 | 0 | **0** | — | — | — | Max nr of reformations reached. | failed to converge at time : 0.05 |
| `fc250_D0.9_on0.90_s20` | FC | 250 N | 0.9 | 20 | on | 0.90 | 1 | 3.6 | 5 | **0** | — | — | — | Max nr of reformations reached. | failed to converge at time : 0.3 |
| `fc250_D1_off_s20` | FC | 250 N | 1 | 20 | off | — | 1 | 0.6 | 0 | **0** | — | — | — | 265 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc250_D1_on0.90_s20` | FC | 250 N | 1 | 20 | on | 0.90 | 1 | 0.6 | 0 | **0** | — | — | — | 265 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D0.9_off_s20` | FC | 1000 N | 0.9 | 20 | off | — | 1 | 0.6 | 0 | **0** | — | — | — | 304 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D0.9_on0.90_s20` | FC | 1000 N | 0.9 | 20 | on | 0.90 | 1 | 0.5 | 0 | **0** | — | — | — | 304 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D1_off_s20` | FC | 1000 N | 1 | 20 | off | — | 1 | 0.3 | 0 | **0** | — | — | — | 98 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D1_on0.90_s20` | FC | 1000 N | 1 | 20 | on | 0.90 | 1 | 0.3 | 0 | **0** | — | — | — | 98 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D1.0_on0.90_s1` | FC | 1000 N | 1 | 1 | on | 0.90 | 1 | 0.3 | 0 | **0** | — | — | — | 98 negative jacobians detected. | failed to converge at time : 1 |
| `dc1x_D0.9_off` | DC | -198.93 mm | 0.9 | 20 | off | — | 0 | 3.9 | 20 | **0** | 1088.3 | 256.2 | 0.8925 |  |
| `dc1x_D0.9_on0.90` | DC | -198.93 mm | 0.9 | 20 | on | 0.90 | 0 | 3.8 | 20 | **0** | 1088.3 | 256.2 | 0.8925 |  |
| `dc2x_D0.9_off` | DC | -397.86 mm | 0.9 | 20 | off | — | 1 | 0.9 | 1 | **0** | — | — | — | 31 negative jacobians detected. | failed to converge at time : 0.1 |
| `dc2x_D0.9_on0.90` | DC | -397.86 mm | 0.9 | 20 | on | 0.90 | 1 | 0.8 | 1 | **0** | — | — | — | 31 negative jacobians detected. | failed to converge at time : 0.1 |
| `dc4x_D0.9_on0.90` | DC | -795.72 mm | 0.9 | 20 | on | 0.90 | 1 | 0.3 | 0 | **0** | — | — | — | 21 negative jacobians detected. | failed to converge at time : 0.05 |
| `dc1x_D1.0_off` | DC | -198.93 mm | 1 | 20 | off | — | 1 | 4.9 | 16 | **0** | — | — | — | Negative jacobian detected. | failed to converge at time : 0.85 |
| `dc1x_D1.0_on0.90` | DC | -198.93 mm | 1 | 20 | on | 0.90 | 1 | 1.1 | 2 | **27** | — | — | — | 5 negative jacobians detected. | failed to converge at time : 0.15 |
| `dc2x_D1.0_off` | DC | -397.86 mm | 1 | 20 | off | — | 1 | 0.5 | 0 | **0** | — | — | — | 32 negative jacobians detected. | failed to converge at time : 0.05 |
| `dc2x_D1.0_on0.90` | DC | -397.86 mm | 1 | 20 | on | 0.90 | 1 | 0.5 | 0 | **0** | — | — | — | 32 negative jacobians detected. | failed to converge at time : 0.05 |
| `fc1000_D0.9_on0.85_s20` | FC | 1000 N | 0.9 | 20 | on | 0.85 | 1 | 0.6 | 0 | **0** | — | — | — | 304 negative jacobians detected. | failed to converge at time : 0.05 |
| `dc1x_D0.9_on0.85` | DC | -198.93 mm | 0.9 | 20 | on | 0.85 | 1 | 1.8 | 5 | **30** | — | — | — | 8 negative jacobians detected. | failed to converge at time : 0.3 |
| `dc2x_D0.9_on0.85` | DC | -397.86 mm | 0.9 | 20 | on | 0.85 | 1 | 0.9 | 1 | **0** | — | — | — | 31 negative jacobians detected. | failed to converge at time : 0.1 |

> FC 行的 “R_peak” 列显示收敛载荷（N）；未收敛显示 `—`。

### 1.1 已知的 erosion-off 参照（P1-a，自动步长）

| 配置 | 结果（G4_FE_STAB_REPORT.md） |
|---|---|
| FC 250 N, Dmax=0.99 | 收敛（2/3）|
| FC 250 N, Dmax=1.0 | 失败 `265 negative jacobians … failed to converge at time : 1` |
| FC 500/1000/25266 N（任意 Dmax） | 全部失败 `negative jacobians … failed to converge` |
| DC δ=-198.929, Dmax=0.9 | 收敛，**R_peak = 1088.3 N**（σ_vm p95=256 MPa，Dmax reached=0.8925）|
| DC δ=-397.858, Dmax=0.9 | 失败 `32 negative jacobians …` |

## 2. 位移控制 R(δ) 轨迹（侵蚀开/关对照）

### `dc1x_D0.9_off`  disp=-198.929 mm, Dmax=0.9, erosion=off  →  R_peak=1088.3 N, eroded=0

| step | R_axis (N) | u_z top (mm) |
|---|---|---|
| 0 | +0.00 | +0.00000 |
| 1 | +210.65 | -9.94645 |
| 2 | +265.28 | -19.89289 |
| 3 | +304.31 | -29.83933 |
| 4 | +336.33 | -39.78578 |
| 5 | +364.75 | -49.73222 |
| 6 | +391.60 | -59.67867 |
| 7 | +417.86 | -69.62511 |
| 8 | +444.24 | -79.57156 |
| 9 | +471.52 | -89.51801 |
| 10 | +500.25 | -99.46445 |
| 11 | +531.00 | -109.41090 |
| 12 | +564.33 | -119.35734 |
| 13 | +601.17 | -129.30379 |
| 14 | +642.37 | -139.25023 |
| 15 | +689.12 | -149.19667 |
| 16 | +743.09 | -159.14313 |
| 17 | +806.31 | -169.08957 |
| 18 | +881.79 | -179.03601 |
| 19 | +973.71 | -188.98245 |
| 20 | +1088.34 | -198.92889 |

### `dc1x_D0.9_on0.90`  disp=-198.929 mm, Dmax=0.9, erosion=on, thr=0.90  →  R_peak=1088.3 N, eroded=0

| step | R_axis (N) | u_z top (mm) |
|---|---|---|
| 0 | +0.00 | +0.00000 |
| 1 | +210.65 | -9.94645 |
| 2 | +265.28 | -19.89289 |
| 3 | +304.31 | -29.83933 |
| 4 | +336.33 | -39.78578 |
| 5 | +364.75 | -49.73222 |
| 6 | +391.60 | -59.67867 |
| 7 | +417.86 | -69.62511 |
| 8 | +444.24 | -79.57156 |
| 9 | +471.52 | -89.51801 |
| 10 | +500.25 | -99.46445 |
| 11 | +531.00 | -109.41090 |
| 12 | +564.33 | -119.35734 |
| 13 | +601.17 | -129.30379 |
| 14 | +642.37 | -139.25023 |
| 15 | +689.12 | -149.19667 |
| 16 | +743.09 | -159.14313 |
| 17 | +806.31 | -169.08957 |
| 18 | +881.79 | -179.03601 |
| 19 | +973.71 | -188.98245 |
| 20 | +1088.34 | -198.92889 |

## 3. 显式 PASS/FAIL

| gate | 阈值 | 本次实测最佳 | 判定 | 依据 |
|---|---|---|---|---|
| (a) FC 1000 N 收敛 | 1000 N | 未收敛 | **FAIL** | 所有 FC 1000 N cell 都在第 1 个载荷增量的 full-load 点载荷奇异处失败；侵蚀无法作用（无收敛步） |
| (b) DC 突破 1088.3 N（需侵蚀真正删除单元） | > 1088.3 N | **1088.3 N**（无侵蚀触发且收敛的 cell） | **FAIL** | 最佳 R_peak=1088.34 N == 已知 1088.3 N，且来自**未删任何单元**的 no-op 格；侵蚀一旦真正触发，R_peak 反而更低、更早发散 |
| (c) 操作载荷 25266.2 N | 25266.2 N | **1088.3 N** | **FAIL** | 距操作载荷 23.2× |

**结论**：侵蚀让模型达到的最高载荷 = **1088.3 N**；操作载荷 = 25266.2 N；**差 23.2×**。

### 3.1 机制发现（为什么侵蚀救不了）

1. **侵蚀只在“某个时间步收敛之后”才作用**（手册 3.9）。力控制下第 1 个增量就把**满载荷**点载荷施加下去（实测：`steps=20` 与 `steps=1` 的第 1 步 residual 完全相同 `1.526718e+03 / 1.986168e+08`，即 FEBio 不默认 ramp 常数节点力）→ 立即负 Jacobian，**没有任何步收敛 → erosion 从未触发**（所有 FC cell `eroded=0`）。
2. **方案阈值 `min=0.9` 与 `Dmax=0.9` 不相容**：CDF Simo 损伤在 `Dmax` 处**从下方渐近**（toy 实测 Dmax=0.9 → D 收敛到 0.8988）。恰好是 1088 N 所在的 `Dmax=0.9` 配置里，侵蚀是**空操作**（`dc1x_D0.9_on0.90` 与 off 完全同解：R_peak=1088.34 N，eroded=0）。
3. **一旦侵蚀真正删除单元，结果更差**：把阈值降到 0.85 或改用 `Dmax=1.0` 让侵蚀触发后，它删除的是**承力单元**，模型**更早**发散：
   - `dc1x_D1.0`：off 走 **16** 步（t=0.85）才失败；on（删 27 单元）**2** 步（t=0.15）就失败。
   - `dc1x_D0.9`：off **20** 步收敛 R=1088 N；on（thr=0.85，删 30 单元）**5** 步（t=0.3）就失败。
4. **FEAMR 强制固定步长，本身削弱鲁棒性**：手册 3.9.3 明确 FEAMR 不能配 auto-time-stepper，而正是 auto-stepper 的 dt 回退（cutback）在近奇异的点载荷步上救回了已知的那些运行。去掉 auto-stepper 后，`dc1x_D1.0_off`（固定步长）在 t=0.85 失败，而 P1-a 的 auto-stepper 版本能在 Dmax=1.0 下收敛（虽仅 R_peak=142 N）——**“用 FEAMR”这件事本身就把模型的求稳能力拿掉了**。

**根因未变**：`fib_cort` 弹性 σ_vm 在 1000 N 已 ~1054 MPa ≈ 29×σ_y（655 节点点载荷 + 薄皮质环），是**点载荷奇异**。侵蚀既不能在奇异处收敛第一步（FC），也不能在收敛后把载荷推高（DC）——删除单元只会挪动/削弱传力路径。

## 4. 阻塞/结论与复现

**阻塞/限制（verbatim 依据见上）**：

- FC 1000 N：第 1 步满载荷即 `98 negative jacobians detected. | failed to converge at time : 0.05`（`steps=1` 时 `failed to converge at time : 1`）。侵蚀**从未触发**。
- DC 突破 1088 N：最高 R_peak=1088.34 N（未删单元的 no-op 格）；侵蚀真正触发时更早失败。
- 操作载荷 25266 N：**不可达**（差 23.2×）。

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_fe_erosion.py
```

- `results/opensim_fe/G4_FE_EROSION_REPORT.md`（本文件）
- `temp/opensim_fe/g4_fe_erosion/`（各 `.feb` / `.xplt` / `.hdf5` / `.log` + `g4_fe_erosion.json`）
- `scripts/opensim_fe/g4_fe_erosion.py`（本脚本；只 import 既有模块）

**未改任何既有模块**：`bone_feb` / `g4_fe_spike` / `g4_fe_stab` / `febio_run` 均只 import 复用。

