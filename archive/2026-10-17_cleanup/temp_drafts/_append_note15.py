import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑮ B2a 两段腿等效单关节（2026-10-16，用户「继续」，提交 `fe64e4c`）

### ★ 设计定理（零新增自由度）
两段腿 `L(θ)=2l·cos(θ/2)` ⇒ `θ(s)=2·arccos(1−s/(2l))`；`M=k_θ·θ`
⇒ **等效轴向力 `F(s) = k_θ·θ/(l·sin(θ/2))`** ⇒ **准静态关节下与单自由度框架完全等价**。
自检：`dF/ds` 随 s 上升 **649→1015 N/m**（硬化）✓。
**冻结值**：`l=0.425 m`（2l=0.85=`leg_ext0`）｜`k_θ=350 N·m/rad`（0.08 N·m/kg/deg×76 kg；Shamaei ~400）
｜`θ_max` 由实测 COM 深度反解 ⇒ **120.8 / 111.2 / 86.3°**。门 `TWO_SEG_LEG` 默认 `None`。

### 判定
**B2a-J0 逐位 PASS｜B2a-J1（软接触单调）三档全 PASS｜B2a-J2 FAIL（stiff 不可定）｜
B2a-J3 只报 0.804/0.762/0.636｜B2a-J5 遵守**

| 档 | **B1（线性）** | **B2a（两段腿）** |
|---|---|---|
| soft | 3.657 m | **3.707 m** |
| **self** | **不可定** | **3.302 m** ✅ 修好 |
| stiff | 不可定 | 不可定（`θ_max`=86.3° **@1 m 用尽**）|

### 结论
- **B2a 修好了单调性**（线性模型在软接触下非单调，因 `k_stop` 假台阶；`θ≤θ_max` 是**几何约束**），
  **并让 `self` 档可定**；
- **`stiff` 仍不可定** ⇒ **与 B1 的 `k_stop` 瓶颈是同一缺口的两种语言表达**
  ⇒ **缺的是「行程耗尽之后」的物理机制（肌肉离心吸能 / 骨性-软组织接触），不是再加几何自由度。**

### ⚠️ 两个必须记住的读数
1. `F(θ→0) → 2k_θ/l ≈ 1.65 kN` **有限**（线性力矩律的产物）⇒ 小角段偏「恒力弹簧」，
   **与真实直腿的无穷抗屈曲刚度极限不符**（量级对、极限不对）；
2. 两段腿的**弹性**力量级偏低：全程仅 **1.65–1.99 kN**（stiff 档线性弹簧 0.42 m 是 7.77 kN）
   ⇒ **软接触口径下 `F_flex` 由阻尼项主导**（本件未预期但必须记录的读数）。
⚠️ 算力：`h=0.20 m` + 实测跟垫 **~186 s/格**（软化非线性律 ⇒ 长时振铃、`settled` 不触发）
⇒ 探针网格**跑前**缩为 `{0.50…9.00}`（跑前选择，非事后裁剪）。

### 回归 / 交付
**777 passed / 8 skipped / 0 failed**（+3 条 `TestTwoSegmentLeg`）。
新增 `scripts/leg2seg_probe.py`、`results/leg2seg_probe.json`、`results/LEG2SEG_REPORT.md`、
`docs/methods/P3_两段腿等效关节_预注册.md`；`PAD_GOTCHAS.md` §6e；`总体进度.md` §九之四。

**锚点**：`results/LEG2SEG_REPORT.md`、`docs/methods/P3_两段腿等效关节_预注册.md`、
`scripts/leg2seg_probe.py`、`pad.py`（`TWO_SEG_LEG`）、`docs/notes/PAD_GOTCHAS.md` §6e、提交 `fe64e4c`
"""

PTR = """

## B2a 两段腿等效单关节完成（提交 `fe64e4c`，用户「继续」）

**★ 设计定理（零新增自由度）**：两段腿 ⇒ `F(s)=k_θ·θ/(l·sin(θ/2))`、`θ(s)=2·arccos(1−s/(2l))`
⇒ 准静态关节下与单自由度框架**完全等价**；`dF/ds` 649→1015 N/m（硬化 ✓）。
**冻结**：`l=0.425`（2l=0.85=`leg_ext0`）｜`k_θ=350 N·m/rad`｜`θ_max` 由实测 COM 深度反解 ⇒ **120.8/111.2/86.3°**。
**判定**：**J0 逐位 PASS｜J1 软接触单调三档全 PASS｜J2 FAIL（stiff 不可定）｜J3 只报 0.804/0.762/0.636**。
**★ 对比 B1**：soft 3.657→**3.707**；**self 不可定→3.302 ✅（修好）**；stiff 仍不可定（`θ_max`=86.3° @1 m 用尽）。
**⇒ 结论：B2a 修好单调性与 self 档；stiff 与 B1 的 `k_stop` 瓶颈是同一缺口的两种语言表达
⇒ 缺的是「行程耗尽之后」的物理机制（肌肉离心吸能/骨性接触），不是再加几何自由度。**
⚠️ `F(θ→0)→2k_θ/l≈1.65 kN` **有限**（小角段偏恒力，极限与直腿不符）；两段腿**弹性力偏低**
（全程 1.65–1.99 kN vs stiff 线性 7.77 kN）⇒ 软接触下 `F_flex` **阻尼主导**。
⚠️ 算力：h=0.20 m + 实测跟垫 **~186 s/格** ⇒ 网格跑前缩为 {0.50…9.00}。回归 **777 passed / 8 skipped**。
详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑮、`PAD_GOTCHAS.md` §6e。
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
