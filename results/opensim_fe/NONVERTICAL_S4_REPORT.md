# Phase-S4 非垂直落地 · 3D wrench 上 FE + 跖面弹性支撑 报告

> 生成时间 2026-10-04 23:09:51；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。
> 方案 `docs/非垂直落地扩展方案.md` §3 L4 / §4 Phase-S4 / §7。
> 本报告**新增**，不改任何模块默认、不覆盖既有产物。

## 0. 一句话结论

1. **轴向回归（硬门槛）: PASS** —— 默认 `build_thums_feb`（fixed 跖面、`use_rigid=False`、F=26383.1 N）复现改前缓存，最大相对误差 **0.00e+00**（gauge_max 186.756369 vs 186.756369 MPa，逐位一致）。
2. **S4 新增 API 全部 opt-in**：`subtalar_force` / `subtalar_moment`（默认 `None`）把完整 3D wrench 送到 FE；跖面弹性支撑 `plantar_bc="spring"`。
3. **力矩两条路**：`use_rigid=True` 的 `RigidMomentLoad` 精确（须放开承载力矩的刚体转动 DOF）；`use_rigid=False` 的**关节面线性压力梯度**在本解剖曲面关节面上**病态**（需 ±193 MPa 带压、含拉伸）→ 负结果。
4. **A7 边缘奇异**：轴向**均匀面压**下弹性支撑把 gauge 降 18.7%（复现 PLANAR_BC 结论）；但在 **3D wrench**（剪切/力矩把热点移到关节面）下**不缓解**（Δ-7.4%）。
5. **3D 载荷（方向效应，同口径同 BC、均无力矩）**：20° 倾斜相对轴向 `gauge` ×0.87；叠加力矩后 `gauge` ×0.68（含转动 BC 影响）。

---

## 1. 本阶段 API（全部 opt-in、默认逐位复现旧行为）

| 文件 / 符号 | 新增参数（默认） | 作用 |
|---|---|---|
| `febio_model.build_calcaneus_feb` | `subtalar_force=None` | 3 分量力矢量；None→`load_n×load_dir` |
| `febio_model.build_calcaneus_feb` | `subtalar_moment=None` | 3 分量力矩 (N·m)；None→零力偶 |
| `febio_model.build_calcaneus_feb` | `moment_bands=6` | 压力梯度分带数 |
| `thums_feb.build_thums_feb` | 同上三参数 | 同构路径（本演示所用） |
| `febio_model.joint_face_frame` | 新 | 关节面面积/形心/外法向/二阶矩 |
| `febio_model.moment_to_pressure_gradient` | 新 | 力矩→线性压力梯度（逐面片精确）|
| `febio_model.linear_pressure_bands` | 新 | 线性压力场→分带常压（合力守恒） |

**默认不变**：`subtalar_force=None, subtalar_moment=None` 时走的是改前的旧代码分支
（不新增任何 wrench 元素——结构测试锁定；FE 结果逐位复现缓存——§2 数值门槛）；
`plantar_bc` 默认仍 `"fixed"`（三向全固定，A7 原状）。

## 2. 轴向回归（硬门槛，方案 §7）

| 量 | 本次（默认） | 改前缓存 | 相对误差 |
|---|---:|---:|---:|
| gauge_max (MPa) | 186.756369 | 186.756369 | 0.00e+00 |
| raw max (MPa) | 213.502779 | 213.502779 | 0.00e+00 |
| p95 (MPa) | 136.296503 | 136.296503 | 0.00e+00 |

- 最大相对误差 **0.00e+00**（门槛 1e-06）→ **PASS**。
- 缓存来源：`results/opensim_fe/plantar_bc_summary.json` 的 `fixed` 行（PLANAR_BC，h=5 m a=0，F_subt=26383.13 N，`use_rigid=False`）。

## 3. S1 20° 倾斜 wrench 的 FE 表达

- 关节载荷面：238 面片，面积 2495.3 mm²，形心 (-1.43,2.22,0.57) mm，外法向 (-0.138,0.891,-0.432)。
- 力 `F=(-5538.7,-20764.4,-0.0) N`：相对面法向分解为 法向 -17741 N + 切向 |F_t|=12128 N。
- 力矩 `M=(3.69,-6.10,76.46) N·m`（|M|=76.79 N·m）：
  - **精确路径（本演示主口径）**：`RigidMomentLoad`（值按 N·mm 写），给力矩时放开承载该力矩的刚体转动 DOF（否则被固定 DOF 吸收、力矩无效）。
  - **压力梯度路径（负结果）**：线性压力场对面片求矩 `M(g)=J g` 精确，但解出
    坡度 12.332 MPa/mm（远大于平面估计 ~0.25），带压范围 [-93.5, 192.9] MPa，**含大幅拉伸**（不物理——真实关节面不能受拉）。
    - 目标面内分量 `(-1.69,28.64,59.62) N·m`；均匀法向压力 p0=7.1098 MPa
      自带曲率力矩 `(-33.94,-6.59,-10.62) N·m`（曲面效应）。
    - → 该几何/力矩下压力梯度表达**病态**（`pressure_path_ill_conditioned=True`）；结论：**本曲面关节面上力矩必须走刚体路径**。

## 4. σ_vm 结果（域 `calcaneus`/CORT，gauge R=4 mm）

| 工况 | 载荷路径 | 跖面 BC | k (N/mm) | max (MPa) | p95 (MPa) | gauge_max (MPa) | 峰值位置 |
|---|---|---|---:|---:|---:|---:|---|
| 轴向基线（均匀面压 F/A） | pressure | fixed | — | 213.50 | 136.30 | 186.76 | 跖面边缘 (plantar rim) |
| 轴向基线（均匀面压 F/A） | pressure | spring | 100 | 151.90 | 121.96 | 151.90 | 关节面边缘 (joint rim) |
| wrench 轴向（F 仅） | rigid | fixed | — | 319.37 | 87.69 | 319.37 | 关节面边缘 (joint rim) |
| wrench 轴向（F 仅） | rigid | spring | 100 | 256.98 | 129.16 | 256.98 | 关节面边缘 (joint rim) |
| S1 20° 倾斜 wrench（F+M） | rigid | fixed | — | 190.05 | 115.98 | 189.05 | 跖面边缘 (plantar rim) |
| S1 20° 倾斜 wrench（F+M） | rigid | spring | 100 | 202.95 | 148.36 | 202.95 | achilles 边缘 |
| S1 20° 倾斜（F 仅，M=0） | rigid | fixed | — | 278.75 | 81.97 | 278.75 | 关节面边缘 (joint rim) |
| S1 20° 倾斜（压力梯度，负结果） | pressure | fixed | — | 1272.80 | 806.14 | 1063.82 | 跖面边缘 (plantar rim) |

### 4.1 A7 边缘奇异：弹性支撑缓解多少？

| 工况 | gauge_max fixed | gauge_max spring | Δgauge | raw max fixed | raw max spring | Δmax |
|---|---:|---:|---:|---:|---:|---:|
| 轴向基线（均匀面压） | 186.76 | 151.90 | +18.7% | 213.50 | 151.90 | +28.9% |
| wrench 轴向（幽灵刚体） | 319.37 | 256.98 | +19.5% | 319.37 | 256.98 | +19.5% |
| S1 20° 倾斜 wrench（幽灵刚体） | 189.05 | 202.95 | -7.4% | 190.05 | 202.95 | -6.8% |

**判定**：弹性支撑对 **A7 跖面固定端奇异**有效（轴向均匀面压下 peak 在跖面棱边、弹性地基把 gauge 降 18.7%）；但 3D wrench 把峰值移到**关节面棱边**（剪切/力矩驱动），跖面弹性支撑**不缓解**该新奇异。

### 4.2 3D 载荷（20° 倾斜）改变多少？

- **方向效应（同口径+同 BC，仅力矢量不同、都无力矩）**：轴向 gauge=319.37 vs 20° 倾斜 gauge=278.75（×0.87）。倾斜把力沿关节面分解，法向分量下降 + 新增横向剪切。
- **力矩效应**（刚体、fixed）：带力矩 gauge=189.05 vs 无力矩 gauge=278.75（×0.68）。⚠ 带力矩必须放开转动 DOF，与无力矩（锁转动）的 BC 不同，故此差值含 BC 影响。
- 压力梯度路径 gauge=1063.82 vs 刚体 gauge=189.05（×5.63）：病态梯度把应力炒高，**不可用**。

## 5. 诚实边界

- **载荷参考帧**：S1 的 3D wrench 在 OpenSim 全局帧；本演示按既有轴向口径**直接把矢量施加在 FE 帧**，未做额外的帧旋转变换（与 `SUBTALAR_DIR=(0,−1,0)` 的既有假设一致）。
- **压力梯度路径是负结果**：在本解剖配准的**曲面**关节面上，线性压力梯度表达该力矩需要病态的大坡度 → 带压出现大幅**拉伸**（不物理）。因此 §4 的主口径用幽灵刚体`RigidMomentLoad`。这是本阶段最有价值的**诚实结论**（方案 §3 L4 的 (a) 路在本几何失败）。
- **刚体力矩的正确性依赖放开转动 DOF**：给力矩时必须放开承载它的刚体转动自由度，否则 `RigidMomentLoad` 被 `RigidFixed` 吸收而静默无效（本实现已处理）。
- **弹性支撑是参数**：k=100 N/mm 取自 PLANAR_BC 最优点，非标定值；三向 Winkler 地基不是真实接触/摩擦。
- **幽灵刚体载荷口径 ≠ 均匀面压口径**：两者把同一关节力以不同方式作用（刚体 → 关节面整体刚性平移；面压 → 分布压力），绝对应力不可直接互比；§4.2 用**同口径**（刚体）做轴向 vs 倾斜的方向比较。
- 本演示仍是**筛选性**研究；σ_vm≥150 MPa 仅后处理判据，无单元删除。

## 6. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s4.py
# 回归：
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
& .venv\Scripts\python.exe -m pytest tests\test_nonvertical_s1.py tests\test_nonvertical_s4.py -q
```

## 7. 产物

- `results\opensim_fe\NONVERTICAL_S4_REPORT.md`（本报告）
- `results\opensim_fe\nonvertical_s4.json`（机器可读）
- `temp\opensim_fe\s4_nonvertical/s4_*.feb/.xplt`（每次求解输入/结果）

