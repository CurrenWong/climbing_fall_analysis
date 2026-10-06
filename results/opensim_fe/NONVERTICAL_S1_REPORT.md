# Phase-S1 非垂直落地（矢量通路）报告

> 生成时间 2026-10-04 22:43:41；仓库 `D:\Project\climbing_fall_analysis`；单位 mm-N-MPa-s。
> 方案 `docs/非垂直落地扩展方案.md` §3 L1–L3 / §4 Phase-S1 / §7。
> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。

---

## 0. 一句话结论

1. **轴向回归（硬门槛）: PASS** —— 默认（竖直）路径与改前缓存 `s1_h5_opensim.json` 的最大相对误差 **0.00e+00**；显式 `tilt_deg=0` 与默认的 GRF 数组**逐位一致** (True)。
2. **倾角 20° 演示**：距下峰值力 25.38 kN（轴向）→ 21.49 kN（倾斜），其中新增横向分量 5.54 kN；`load_transfer` 现保留完整 3D 力矢量 + 力矩。
3. 逐骨纯轴向 1D 风险见 §3；倾斜是否改变排序见 §3.3。

---

## 1. 本阶段实现的 API（全部 opt-in、默认复现旧行为）

| 层 | 文件 / 符号 | 新增（默认） |
|---|---|---|
| L1 GRF | `opensim_grf.ground_reaction` | `ground_normal=(0,1,0)`、`tilt_deg=None`（度）|
| L1 GRF | `GroundReaction.ground_normal` / `.is_vertical` / `.force_vec_n(side)` | 默认竖直，`F⃗=f(t)·n̂` |
| L2 FD | `opensim_fall._add_foot_force` | `f_x=None, f_z=None, point=None`（默认 0 / 原点）|
| L2 FD | `opensim_fall.run_dead_drop` | `vx0_ms=0.0, vz0_ms=0.0` |
| L3 传递 | `load_transfer.FeLoadSpec` | `subtalar_force=None, subtalar_moment=None`；`force_vector_n` / `moment_vector_nm` |
| L3 传递 | `load_transfer.transfer` | `subtalar_force=None, subtalar_moment=None` |
| 消费 | `joint_reactions.joint_reaction` | 自由体按 `grf.force_vec_n` 计入完整 3D 外力（默认竖直 = 旧行为）|

> 中段 `joint_loads.subtalar_reaction` **未改**（本就 3D）。

---

## 2. 轴向回归（硬门槛，方案 §7）

默认链路重跑 vs 缓存：

| 量 | 本次默认 | 缓存 | 相对误差 |
|---|---:|---:|---:|
| `impulse_ns` | 789.328187 | 789.328187 | 0.00e+00 |
| `grf_peak_n` | 53064.8332 | 53064.8332 | 0.00e+00 |
| `subtalar_peak_vertical_n` | 26383.1338 | 26383.1338 | 0.00e+00 |
| `subtalar_peak_force_n` | 26383.134 | 26383.134 | 0.00e+00 |
| `subtalar_peak_moment_nm` | 20.0250922 | 20.0250922 | 0.00e+00 |

- 最大相对误差 **0.00e+00**（阈值 1e-6）。
- 显式 `tilt_deg=0` vs 默认：GRF 数组逐位一致 = **True**；距下纵向峰值相对差 = 0.00e+00。
- **判定：PASS**

---

## 3. 倾角演示（载荷链层面）

地面法向 `n̂ = (0.3420, 0.9397, 0.0000)`（竖直占比 0.9397）。

### 3.1 距下关节（subtalar_r）wrench

| 量 | 轴向 | 倾斜 |
|---|---:|---:|
| 峰值 |F| (N) | 25,379.57 | 21,490.38 |
| 峰值 纵向 y (N) | 25,354.40 | 20,764.38 |
| 峰值 横向 x (N) | 1,139.10 | 5,542.64 |
| 峰值 |M| (N·m) | 156.72 | 273.64 |

### 3.2 `load_transfer` 保留的 FE 载荷规格

- `subtalar_n`（纵向标量）= 20,764.38 N
- `force_vector_n` = (-5,538.68, -20,764.38, -0.00) N，|F| = 21,490.38 N
- `moment_vector_nm` = (3.6912, -6.0970, 76.4603) N·m
- 旧口径等价物 = `subtalar_n × SUBTALAR_DIR`（仅纵向）；现完整矢量与力矩不再被丢弃。

### 3.3 逐骨纯轴向 1D 风险（Route-1 + 载荷分配）

`risk = (mult · F_joint) / (A_section · σ_c)`；`F_joint` 取关节反力的**纵向分量**（Route-1 口径）为主，另列**合成 |F|** 作为 3D 载荷代理。mult 来自 `risk_1d_loadshare`（胫/腓并联 k=E·A/L；中轴骨 mass-above）。

| 骨 | 关节 | F_y 轴向 (N) | F_y 倾斜 (N) | risk 轴向 | risk 倾斜(纵向) | risk 倾斜(合成) | 偏移(合成) |
|---|---|---:|---:|---:|---:|---:|---:|
| `calcaneus_r` | `subtalar_r` | 25,354.4 | 20,764.4 | 0.7825 | 0.6408 | 0.6633 | ×0.848 |
| `fibula_r` | `ankle_r` | 25,266.2 | 20,694.6 | 0.3599 | 0.2948 | 0.3046 | ×0.846 |
| `T6` | `lumbar` | 22,087.6 | 17,670.2 | 0.2676 | 0.2141 | 0.2182 | ×0.816 |
| `tibia_r` | `ankle_r` | 25,266.2 | 20,694.6 | 0.2674 | 0.2190 | 0.2263 | ×0.846 |
| `L3` | `lumbar` | 22,087.6 | 17,670.2 | 0.1904 | 0.1523 | 0.1553 | ×0.816 |
| `femur_r` | `hip_r` | 15,227.7 | 12,215.9 | 0.1710 | 0.1371 | 0.1401 | ×0.819 |
| `C5` | `lumbar` | 22,087.6 | 17,670.2 | 0.1243 | 0.0994 | 0.1013 | ×0.816 |
| `R_HIPBONE` | `hip_r` | 15,227.7 | 12,215.9 | 0.0434 | 0.0348 | 0.0356 | ×0.819 |
| `parietal_r` | `lumbar` | 22,087.6 | 17,670.2 | 0.0149 | 0.0119 | 0.0121 | ×0.816 |

- 轴向 @5m 已越阈（risk≥1）：`无`。
- 倾斜（合成 |F|）@5m 越阈：`无`。
- 说明：倾角把法向大小按方向拆分，**纵向分量下降、横向分量新增**；纯轴向 1D 指标只吃纵向，因此倾斜下纵向风险普遍**下降**，而 3D 合成载荷（含横向剪切）才是倾角真实加载强度 —— 两者之差即"纯轴向指标的盲区"。

---

## 4. 诚实边界

- **矢量通路是方向化的最小版**：GRF 标量大小不变，只旋转方向；接触本构、摩擦锥、CoP 偏移、单脚均未实现（方案 §3 L1 完全版 / S2 / S4）。
- **未上 FE**：`FeLoadSpec` 已保留完整 3D 力 + 力矩，但当前 FE 适配器仍只消费大小×方向；力矩/偏心压力映射属 S4，本阶段不宣称已实现。
- 逐骨风险是 **1D 名义应力 + 载荷分配**，忽略弯曲/剪切/偏心与应力集中；`σ_c`、`A_section`、质量分数均为既有 assumed/estimated 口径，绝对值不可信，只用**相对**。
- 倾角动力学为**两足对称、刚性腿**（方案 S1 范围），非单脚。

## 5. 复现命令

```powershell
cd D:\Project\climbing_fall_analysis
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s1.py --tilt-deg 20
# 回归：
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
& .venv\Scripts\python.exe -m pytest tests\test_nonvertical_s1.py -q
```

## 6. 产物

- `results\opensim_fe\NONVERTICAL_S1_REPORT.md`（本报告）
- `results\opensim_fe\nonvertical_s1.json`（机器可读）

