import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑬ ①独立能量审计 ⇒ **A 判定成立**；②非单调机制定案（2026-10-16，用户「按你的建议」，提交 `db6c29e`）

### ① 第 1 步：立有效性判据并重判 A —— **A 判定成立**（撤回 → 未定 → **成立**）
新增 **`scripts/pad_energy_check.py::audit_energy`**：沿轨迹
`E_mech(t) = KE + PE_grav + E_elastic` 必须非增（`E_mech(0)=0` 取 t=0 静止参考）。
**结果：A 的 12 个格全部 `violation_j = 0.00 J`、`E_mech` 非增**（容差 1%·mgh = 3.30 J）
⇒ **`pad.py` 动力学能量守恒、无数值缺陷 ⇒ A 的数值物理有效 ⇒ 判定恢复成立**
（A-J0 PASS / A-J1 FAIL / A-J2 PASS / A-J3 如实）。

**★ 审计器口径坑（我踩了两次，已钉死入库 `PAD_GOTCHAS.md` §6c）**：
1. `result.t_s` 从 **t=0（释放）** 起，**不是**接触时刻（`meta["impact_window_s"]` 才是冲击窗）；
2. 窗首速度是 **0**，**不是** `impact_speed_ms` —— 拿它当偏移会得到**系统性 ≈2·mgh 的假 FAIL**
   （实测 violation 恰好 = 2×329.5 = **659.6 J**，一眼可辨）；
3. 接触压缩量 = `z_foot − result.height_m`，**不是** `z_foot − z_foot[0]`
   （后者把整个下落量当压缩 ⇒ 报 **435.7 mm** 而真值 **15.7 mm**）。
   钉死手段：`temp/v1_arrdump.py` 逐项打印数组首尾 + `meta` 窗口。
4. `f_flex_peak_n` **只在冲击窗内**取 max ⇒ 与全轨迹 max 不可比（全轨迹峰可大 4–8×）。

### ② 第 2 步：软接触 `F_flex(h)` 非单调 —— **机制定案：进 `k_stop` 止段**
`anchored-self`、`k_n=2.1e6`（`k_flex=8.2 kN/m`、`x_flex_max=0.40 m`、`k_stop=3.0e5`）：

| h | 冲击窗 | `F_flex`峰 | 弹簧峰 | 阻尼峰 | **`xf@峰`** | `>x_flex_max`? |
|---|---:|---:|---:|---:|---:|---|
| 1 m | 143 ms | 3.548 | 2.710 | 2.613 | 0.1298 | 否 |
| **2 m** | **252 ms** | **9.309** | **9.165** | 3.711 | **0.4191** | **是** |
| 3 m | 43 ms | 6.051 | 2.463 | 4.568 | 0.2052 | 否 |
| 5 m | 40 ms | 7.781 | 2.966 | 5.949 | 0.2532 | 否 |

**机制（二重证据）**：软接触 ⇒ 腿压得更深 ⇒ 恰在 h=2 m **跨过 `x_flex_max=0.40` 进止段**
⇒ 峰值被顶起：`3.3 kN + 3.0e5×(0.4191−0.40) = 9.03 kN` ≈ 实测 **9.165 kN** ✓。
同档 `k_n=None` 窗口仅 **17–18 ms**、`xf ≤ 0.140 m` ⇒ **永不进止段 ⇒ 单调** ✓。
⇒ **是未标定字段（`x_flex_max` / `k_stop`，继承自 `controlled-drop`）的产物**，
非积分误差（三段解算器设置逐位一致）、非能量缺陷。
**⇒ 软接触区 `h*` 仍不可引用**，须先标定/冻结这两个字段。

### 回归 / 交付
**771 passed / 8 skipped / 0 failed**（+3 条 `TestEnergyAudit`：守恒 / **压缩量口径锁死** / 软接触深行程）。
新增 `scripts/pad_energy_check.py`；`PAD_GOTCHAS.md` §6b 重写 + 新增 §6c；
`FP_CONTACT_REPORT.md` 与 `总体进度.md` §九之二 状态同步为**成立**。

### ⇒ 下一步 B 的**前提**（新发现，必须一起处理）
`x_flex_max=0.40` / `k_stop=3.0e5` 是**未标定**字段，且已被证实会**制造首折台阶**
⇒ **B（膝扭簧）的预注册必须同时冻结/标定这两个字段**，否则软接触区依旧不可用。

**锚点**：`scripts/pad_energy_check.py`、`tests/test_climbing.py::TestEnergyAudit`、
`results/FP_CONTACT_REPORT.md`、`docs/notes/PAD_GOTCHAS.md` §6b/§6c、
`temp/v1_arrdump.py`、`temp/v1_window_decomp.py`、提交 `db6c29e`
"""

PTR = """

## ①独立能量审计 ⇒ A 判定成立；②非单调机制定案（提交 `db6c29e`，用户「按你的建议」）

**① 第 1 步**：新增 `scripts/pad_energy_check.py::audit_energy`（沿轨迹 `E_mech = KE+PE+E_elastic` 非增）
⇒ **A 的 12 个格全部 `violation_j = 0.00 J`** ⇒ **`pad.py` 能量守恒无数值缺陷 ⇒ A 判定恢复成立**
（【撤回】→【未定】→**【成立】**）。
**★ 审计器口径坑（踩两次，已入库 §6c）**：`t_s` 从 **t=0（释放）** 起；窗首速度是 **0** 不是
`impact_speed_ms`（用错 ⇒ **系统性 ≈2·mgh 假 FAIL**，实测恰 659.6 J = 2×329.5）；
压缩量 = `z_foot − height_m` 不是 `z_foot − z_foot[0]`（后者报 435.7 mm vs 真值 15.7 mm）。
**② 第 2 步**：软接触 `F_flex(h)` 非单调**机制定案** —— **进 `k_stop` 止段**：
self@2.1e6 的 h=2 m **`xf@峰`=0.4191 > `x_flex_max`=0.40** ⇒ 弹簧峰 9.165 kN ≈ `3.3+3e5×0.0191` ✓；
h=1/3/5 m 未进止段 ⇒ 3.548/6.051/7.781；同档 `None` 窗口仅 17–18 ms、`xf≤0.14` ⇒ 单调。
⇒ **是未标定字段 `x_flex_max`/`k_stop`（继承 `controlled-drop`）的产物**，非数值误差/非能量缺陷；
**软接触区 `h*` 仍不可引用**。
回归 **771 passed / 8 skipped**；**⇒ B 的前提：预注册须同时冻结/标定 `x_flex_max` 与 `k_stop`。**
详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑬、`PAD_GOTCHAS.md` §6b/§6c。
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
