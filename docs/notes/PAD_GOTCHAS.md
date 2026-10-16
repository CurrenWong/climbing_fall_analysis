# `pad.py` / 1D-2D 链模型 · GOTCHAS

> 与 `FE_BC_GOTCHAS.md`（FE/BC 专用）同类，但范围是 **`src/climbing/pad.py`、`pad2d.py`、
> `chain2d*.py`、`bodies.py`** 这条落地动力学线。
> 新增条目请附**实测数字**与**复现脚本**。

---

## 1. 🔴 `peak_force_n`（`F_pad`）对 `k_flex` **不敏感** —— 不得当 GRF / 腿刚度观测量

**现象（硬地面 `on_pad=False`，@0.42 m、80 kg，`anchored-*` 三档）**：

| 档 | `k_flex` (kN/m) | **`F_pad` (kN)** | `F_flex` (kN) |
|---|---:|---:|---:|
| soft | 6.3 | **8.61** | 2.27 |
| self | 8.2 | **8.61** | 2.34 |
| stiff | 18.5 | **8.63** | 2.78 |

**⇒ `k_flex` 变 2.9 倍，`F_pad` 只动 0.2% ⇒ 它由 `hard_surface(k_n=1e7)` × `m_low`
的**刚性接触尖峰**决定，不是腿刚度的响应量。**

- ✅ **腿刚度的观测量 = `result.meta["f_flex_peak_n"]`**（对 k 敏感：2.27→2.78，+22%）。
- ⛔ 用 `peak_force_n` 做的「模型 vs 力台 GRF」对比**会得到无意义的结论**
  （2026-10-16 曾据此误判「模型 @0.42 m 高 3.3×、缺屈服吸能」——已作废）。
- 复现：`scripts/fp_anchor_probe.py`、`tests/test_climbing.py::TestAnchoredPostures::
  test_f_pad_insensitive_to_k_flex`。

---

## 2. 🔴 **L2 刚性门按姿势的「键」判定 ⇒ 传 `Posture` 对象会静默绕门**

`pad.py` 的 `RIGID_STANDING_CHAIN = {"kf_mult": 2.0, "x_max": 0.015, "k_stop": None}`
（2026-10-08 L2 标定采纳）**只在 `on_pad=False` 且姿势的*键* == `"feet-first-stiff"`** 时生效
⇒ 把 `k_flex` **×2.0**（1.1e6 → 2.2e6）、行程限到 **15 mm**。

实现（`simulate_boulder_fall`）：

```python
post = POSTURES[posture] if isinstance(posture, str) else posture
# L2 门用 POSTURES 的**键**判姿势（post.name 是中文显示名）
posture_key = posture if isinstance(posture, str) else None
```

**⇒ 传 `replace(POSTURES["feet-first-stiff"], ...)` 这类**对象**时 `posture_key=None`
⇒ 门不生效 ⇒ 同名姿势两种调用结果不同。**

**实测（`feet-first-stiff` @1 m，80 kg，硬地面）**：

| 调用 | `F_pad` (kN) | `F_flex` (kN) |
|---|---:|---:|
| `posture="feet-first-stiff"`（**键**，含门） | **37.46** | **31.10** |
| `posture=<Posture 对象>`（**无门**） | **33.65** | **27.43** |

**⇒ 差 13.4%。** 后果与用法：

1. **引用既有数字必须标「含门/无门」**：`V1_WHY_20KN_AND_GAP_LOCUS.md` §3 表里标着
   「1.1e6（模型 `feet-first-stiff`）」的 27.43 kN 是**无门**值；真实 `feet-first-stiff`
   在刚性口径下是 **2.2e6 / 15 mm**。
2. **新姿势的键若不是 `"feet-first-stiff"` ⇒ 门不生效**（设计如此：「pad 路径与其他姿势
   逐位不变」）。**照抄姿势字段造新档时，务必确认你要的是有门还是无门口径**
   （`anchored-*` 就是要**无门**，因为比的是 `k_flex` 本身 —— 预注册 §4-O5）。
3. 守卫测试：`tests/test_climbing.py::TestAnchoredPostures::test_l2_rigid_gate_key_vs_object_differs`
   （**特征化**测试：锁定这不是 bug、锁死 ratio 区间）。

---

## 3. 🟡 `POSTURES` 的键 ≠ `Posture.name`

`Posture.name` 是**中文显示名**（如「脚先落但腿绷直」），**不是**键。
任何按名字分支的逻辑都应在**键**上做（同 §2）；日志/报告里出现中文名时不要反查成键。

---

## 4. 🟡 `x_flex_max` 不是主控参数，行程由 `x = F/k` 决定

`x_flex_max` 只是**行程上限**。在触到上限之前，实际行程 `x = F/k` 由刚度决定
（例：`k=1.1e6` ⇒ `x≈25 mm` ≪ 30 mm 上限 ⇒ 上限根本没被触及）。
**⇒ 对 `x_flex_max` 单变量做灵敏度往往会得到"无影响"的假阴性**
（2026-10-16 实测：`feet-first-stiff` 基座 0.03→0.20 m，`F_flex@4 m` 45.98→60.30 kN **不降反升**）。

---

## 5. 🔴 `hard_surface` 的 `k_n` **不是**总刚度 —— 有效刚度 = `k_n × A(κ)`

`CrashPad.force_n(compression, area) = stress_pa(xi) × area`，而
`hard_surface(k_n_per_m=...)` 把该形参喂给 **`linear_modulus_pa`（Pa）**、`thickness=1`
⇒ **有效总刚度**

```
k_eff = k_n × A(κ)          # A(κ) = Posture.contact_schedule 的接触面积
```

**实测（`anchored-stiff` @0.42 m、80 kg）**：

| `k_n` | `F_pad` (kN) | `x` (mm) | **实测 `k_eff` = F/x** | `k_eff/k_n` |
|---|---:|---:|---:|---:|
| 1.0e7 | 8.625 | 15.728 | **548.4 kN/m** | 0.055 |
| 3.0e6 | 5.194 | 29.333 | 177.1 kN/m | 0.059 |
| 1.0e6 | 3.870 | 57.236 | 67.6 kN/m | 0.068 |

**⇒ `hard_surface` docstring 声称的「默认 1.0e7 N/m（动态人体等效刚度）」比**实际**硬 ~18×。**

**用法纪律**：
1. **引用硬地面接触硬度时必须报实测 `k_eff`（= `F_pad/max_compression_m`），不得引用 `k_n` 本身。**
2. 想调接触硬度就调 `RIGID_CONTACT["k_n_per_m"]`（**L3 门**，默认 `None` = 1.0e7）——
   ⚠️ **`k_eff` 是 `k_n` 的非线性函数**（`A(κ)` 随压缩变化，`k_eff/k_n` 在 0.055–0.16 间漂移），
   所以**不能**按 `k_eff` 直接线性反解；须实测。
3. 特征化/守卫测试：`tests/test_climbing.py::TestRigidContactGate`。
4. 锚点与实测见 `results/FP_CONTACT_REPORT.md`、`docs/methods/P3_足部接触柔度_预注册.md`。

---

## 6. 🔴 `absorbed_fraction` / `energy_residual_j` 是**口径混乱量** —— 不得当「能量守恒判据」

**现象**：`absorbed_fraction = (energy_into_pad_j + energy_into_flex_j) / impact_energy_j`
在软接触下可达 **1.583**（>1），表面上"违反能量守恒"。

**根因（读码 `pad.py:733-741` 定案：是账的口径，不是物理缺陷）**：

- **分子** = `cumtrapz(f_flex · d(xf)/dt, t)`，而 `f_flex = k·xf（+`k_stop` 段） + c_flex·rate`
  ⇒ **含不可恢复的阻尼耗散**（`c_flex = 400 N·s/m`），且符号变号的弹性项被积掉
  ⇒ 分子 ≈ **阻尼耗散累积 + 净弹性功**，随身体**继续下沉**而增长，**不受 `mgh` 约束**；
- **分母** `impact_energy_j = m·g·h` 只对应**首次接触那一刻**的能量。

⇒ **不同物理量 + 不同时间基准 ⇒ `>1` 完全正常。**（`pad.py` 注释自承历史上有过
「`E_pad + E_flex > m·g·h`（吸能超过输入能量）这种自相矛盾的结果」的读法 —— 他们改过实现，
但**根本口径问题仍在**：阻尼混入 + 时间基准不一致。）

**实测（@0.42 m、`anchored-*`）**：`None` 档 0.276 / 0.289 / 0.329；`k_n=2.1e6` 档 **1.568 / 1.513 / 1.376**。

**⇒ 用法纪律**：
1. ⛔ **不得**用 `absorbed_fraction > 1` 判「数值不物理」（2026-10-16 曾据此误撤 `FP_CONTACT` 判定）；
2. ⛔ **不得**把 `energy_residual_j` 当「必须由组织损伤承担的能量」（**同源污染**；越界区被夹到 0）；
3. ✅ 要判能量守恒，须**沿轨迹**算 `E_mech(t) = KE + PE_grav + E_elastic` 并检查**非增**：
   参照 `E_mech(0)=0`，`PE_grav = −[m_low·g·z_foot + m_up·g·(z_torso + leg_ext0)]`，
   `E_elastic = ∫F_pad dx + 分段弹簧能`。
   ⚠️ **`result.t_s` 只覆盖冲击窗，窗起点速度是 `v_impact = √(2gh)`（不是 0）**
   ⇒ 速度重建**必须带这个偏移**（`result.impact_speed_ms`），否则检验无效（已踩过一次）。

---

## 6b. 🟡 软接触下 `F_flex(h)` **收敛地非单调** —— 机制已定案：**进 `k_stop` 止段**

`anchored-self`、`k_n = 2.1e6`（`k_flex = 8.2 kN/m`、`x_flex_max = 0.40 m`、`k_stop = 3.0e5 N/m`）：

| h | **冲击窗长度** | `F_flex` 峰 | 弹簧峰 | 阻尼峰 | **`xf@峰`** | **`> x_flex_max`?** |
|---|---:|---:|---:|---:|---:|---|
| 1 m | 143 ms | 3.548 kN | 2.710 | 2.613 | 0.1298 m | 否 |
| **2 m** | **252 ms** | **9.309** | **9.165** | 3.711 | **0.4191 m** | **是** |
| 3 m | 43 ms | 6.051 | 2.463 | 4.568 | 0.2052 m | 否 |
| 5 m | 40 ms | 7.781 | 2.966 | 5.949 | 0.2532 m | 否 |

**⇒ 机制（不可翻转的二重证据）**：软接触 ⇒ 腿压得更深 ⇒ **恰在 h=2 m 跨过 `x_flex_max = 0.40 m`
进入 `k_stop` 止段** ⇒ 峰值被顶起：`3.3 kN + 3.0e5×(0.4191−0.40) = 9.03 kN` ≈ 实测弹簧峰 **9.165 kN** ✓。
同档 `k_n = None`（刚性）窗口只有 17–18 ms、`xf ≤ 0.140 m` ⇒ **永不进止段 ⇒ 单调** ✓。

**⇒ 定性**：**不是积分误差**（三段解算器设置逐位一致）、**不是能量缺陷**
（沿轨迹 `E_mech` 审计 `violation = 0.00 J`，见 §6），而是
**未标定字段 `x_flex_max` / `k_stop`（从 `controlled-drop` 继承）**的产物。

**⇒ 后果与纪律**：
1. ⛔ **软接触区的 `h*`（首折内插）不可直接引用** —— 止段会制造假台阶；
2. ⚠️ 冲击窗长度随接触刚度剧变（40 / 143 / 252 ms）⇒ 峰值对窗口定义敏感；
3. ✅ 要修的是**先标定/冻结 `x_flex_max` 与 `k_stop`**（或把窗口限定在 `xf ≤ x_flex_max` 内），
   而不是改接触层。

---

## 6c. 🔴 `pad.py` 结果数组的**口径**（审计必需，2026-10-16 实测钉死）

| 量 | 口径 |
|---|---|
| `result.t_s` | 从 **t = 0（释放瞬间）** 起，**不是**接触时刻；`meta["impact_window_s"]` 才是冲击窗 |
| 窗首速度 | **0**（`t=0` 静止）；`impact_speed_ms` 是**接触瞬间**的速度 —— 拿它当窗首偏移会给系统性 ≈2·mgh 的假 FAIL |
| 接触压缩量 | `max(0, z_foot − result.height_m)`；**不是** `z_foot − z_foot[0]`（那样会把整个下落量当压缩，报 435.7 mm 而真值 15.7 mm） |
| `z_torso` 初值 | `−leg_ext0 = −0.85` ⇒ 重力势参考点取 `t=0` 时 `PE=0` |
| 加速度符号 | `accel_*_g` 以**向上减速为正** ⇒ `az_down = −g·accel_*_g` |
| `f_flex_peak_n` | 只在**冲击窗内**取 `max` ⇒ 与全轨迹 `max` 不可比（全轨迹峰值可大 4–8×） |

**⇒ 现成工具**：`scripts/pad_energy_check.audit_energy(result, posture, mass_kg)`
（回归守卫：`tests/test_climbing.py::TestEnergyAudit`）。

---

## 6d. 🔧 L4/B1 门 `STROKE_ANCHOR` + `h*` 内插的**跨无效格**陷阱

**门**（`pad.py`，默认 `None` ⇒ 既有路径逐位不变；§6b 的「要修的是先标定 `x_flex_max`/`k_stop`」即此）：

```python
STROKE_ANCHOR = {"x_flex_max_m": None, "k_stop_n_per_m": None}
```
用于把 `x_flex_max` 锚定到**实测 COM 深度**（`PMC9465758`「Drop jump depth」列，与三档 `k_flex` 同源一一对应）：
**soft 0.43 / self 0.37 / stiff 0.23 m**；`k_stop` 是**未标定**的硬编码值，可覆盖以便标定性探针。
⚠️ **L2 门排在 L4 门之后** ⇒ 对 `feet-first-stiff` 键由 L2 覆盖；`anchored-*` 不受 L2 影响。

**实测结论**（`results/FP_STROKE_REPORT.md`）：
- **刚性口径：止段从不触发 ⇒ B1 对 `LEG_K_ANCHOR` 的 `h*`（3.418 / 3.216 / 2.244 m）逐位零影响**；
- **实测跟垫口径：止段在 `stiff`@1 m、`self`@2 m 就触发** ⇒ 峰值被顶起（`stiff`@1 m 3.64 → **6.51 kN**，假台阶）
  ⇒ `self`/`stiff` 在**有效范围内跨不过 6.2 kN ⇒ `h*` 不可定**；仅 `soft` 可定 = **3.657 m**。
  ⇒ **`k_stop`（未标定）已是软接触区的主控瓶颈。**

**⇒ ⇒ ★ 内插陷阱（我踩了一次）**：算 `h*` 时**不能**「先滤掉止段污染格、再在剩下的格之间内插」——
那会**跨过**无效格：实测因此把 `anchored-self` 报成 **3.506 m 的假首折**
（真正的跨越段是 1 m → 5 m，中间夹着被污染的 2 m 格）。
**正确做法**：只在**原网格上相邻、且两端都未触发止段**的格之间内插
（`scripts/fp_stroke_probe.py::_hstar`）。

---

## 6e. 🔧 B2a 门 `TWO_SEG_LEG`（两段腿 = 非线性硬化弹簧；零新增自由度）

```python
TWO_SEG_LEG = {"k_theta_nm_per_rad": None, "segment_l_m": None, "theta_max_rad": None}
```

**设计定理**：两段腿 `L(θ)=2l·cos(θ/2)` ⇒ `θ(s)=2·arccos(1−s/(2l))`；`M=k_θ·θ`
⇒ **等效轴向力 `F(s) = k_θ·θ/(l·sin(θ/2))`** ⇒ **准静态关节下与单自由度框架完全等价**（不用加自由度）。
自检：`dF/ds` 随 s 上升 **649 → 1015 N/m**（硬化）✓。

**冻结值**：`l = 0.425 m`（2l = 0.85 = `leg_ext0`）｜`k_θ = 350 N·m/rad`
（0.08 N·m/kg/deg × 76 kg；Shamaei 2013 膝 ~400）｜`θ_max` 由实测 COM 深度反解 ⇒ **120.8 / 111.2 / 86.3°**。

**实测收益/局限**（`results/LEG2SEG_REPORT.md`）：
- ✅ 软接触下 **`F_flex(h)` 三档全部单调**（线性模型在此**非单调**，因 `k_stop` 假台阶）；
  `θ ≤ θ_max` 是**几何约束**，不再需要 `k_stop` 造力；
- ✅ **`self` 档从「不可定」→ 可定（3.302 m）**（B1 下不可定）；
- ❌ **`stiff` 仍不可定**（`θ_max` = 86.3° 在 h=1 m 就被用尽）⇒ **「行程耗尽之后」仍缺物理机制**。

**⚠️ 两个必须记住的读数**：
1. **`F(θ→0) → 2k_θ/l ≈ 1.65 kN` 是有限的** —— 线性力矩律的产物；小角段偏「恒力弹簧」，
   **与真实直腿的无穷抗屈曲刚度极限不一致**（量级对、极限不对）；
2. **两段腿的**弹性**力量级偏低**：全程仅 **1.65–1.99 kN**（stiff 档线性弹簧在 0.42 m 是 7.77 kN）
   ⇒ 软接触口径下输出的 `F_flex` **由阻尼项主导**。

**⚠️ 算力**：`h=0.20 m` + 实测跟垫口径 **~186 s/格**（软化非线性律 ⇒ 长时振铃、`settled` 不触发）
⇒ 探针网格跑前缩为 `{0.50…9.00}`（跑前选择，非事后裁剪）。

---

## 6f. ⚠️ 「门」是**临时脚手架** —— 已实际咬过人，须定期收编

`pad.py` 现有四个「门」（模块级 dict / 按姿势键判定 ⇒ **运行时开关**，默认态 = 原逻辑）：

| 门 | 档位 | 作用 | 默认识别 |
|---|---|---|---|
| `RIGID_STANDING_CHAIN` | L2 | 刚性+直立的标定链（`k_flex`×2、行程限 15 mm） | **按姿势键** == `feet-first-stiff` |
| `RIGID_CONTACT` | L3 | 刚性地面接触刚度 `k_n` | `None` → `1.0e7` |
| `STROKE_ANCHOR` | L4 | 覆盖 `x_flex_max` / `k_stop` | `None` → 沿用姿势值 |
| `TWO_SEG_LEG` | B2a/C1 | 两段腿非线性弹簧 / 离心饱和律 | `None` → 线性 `k_flex` + 粘性阻尼 |

**风险（已实际咬过）**：**L2 门按姿势的「键」判定** ⇒ 传 `Posture` **对象**会**静默绕门**，
同一姿势两种调用差 **13.4%**（§2）；且 `pad.py` 自己的 docstring 曾把 `hard_surface` 的
`k_n` 写成总刚度（实为「模量 × 面积」，差 18×，§5）。

**纪律**：
1. 门**只用于标定性探针**；探针必须 `finally` 恢复，且**「门关态逐位正对照」是提交的前提**；
2. **新门必须在 gotcha 里登记**（档位/作用/默认/生效条件）；
3. **标定定稿后应收编** —— 把值写进 `Posture` 字段并删门；否则门越积越多，
   而「按键判定」这类**隐式**行为会随重构**静默改变数值**。

---

## 6g. 🔴 纪律：标为「外部锚」的数字**必须带可点击来源**；凭印象的必须显式标 `【未核实】`

**为什么**：2026-10-16 同一会话里我犯了**至少三次**同类错，每次都差点把结论带偏：

| # | 我用过的说法 | 实情 | 后果 |
|---|---|---|---|
| 1 | `PMC9465758` 表头「Stiffness, N/m」当成 kN/m | 量纲实为 **N·kg⁻¹·m⁻¹**（须 × 体重） | 模型「过刚」从 **12×** 误判成 0.90× |
| 2 | `f_ecc_max = M_ecc/(l·sin(θ/2))`（**肢体**力臂） | 关节力矩的臂是**肌腱**的（`r(θ)` ≈ 3–4 cm） | θ→0 **发散** ⇒ 积分器崩；小角段上限表全错 |
| 3 | `s_hard` =「深蹲 COM 下降 ≈ 0.28–0.32·H」 | **搜不到来源**（未核实且未标注） | 改用**解剖 ROM**（深蹲膝屈曲 110–135°）才站住 |

**⇒ 纪律**：
1. **凡标「外部锚 / 实测 / 文献值」的数字，必须写清**：**出处**（可检索 ID）、**量纲**、
   以及**被试/工况**（谁、多少、什么条件）；
2. **凭印象的系数必须显式标 `【未核实】`**，且**不得**用于判据或结论；
3. **量纲单独校验**：把文献值代回公式，检查**量纲是否闭合**、**数量级是否合理**（#1 即由此露馅）；
4. **换算前先问「分子分母是不是同一个物理量、同一时刻」**（#1 的根因；`absorbed_fraction` 同源）；
5. **检查器本身要有已知正对照**（本会话的「逐位正对照」两次拦住「拿错口径当基准」，
   见 `results/PASSIVE_PATH_REPORT.md` §3）。

---

**最近更新**：2026-10-16（`LEG_K_ANCHOR` → `FP_CONTACT` → 账的口径 → 审计 → B1 → B2a → C1 → **C2 全 PASS** → **纪律 §6g**）。
**锚点**：`src/climbing/pad.py`（四个门 + `flex_force` + 能力账 `:733-741`）、
`scripts/pad_energy_check.py`、`scripts/fp_stroke_probe.py`、`scripts/leg2seg_probe.py`、
`scripts/ecc_muscle_probe.py`、`results/{FP_STROKE,LEG2SEG,ECC_MUSCLE}_REPORT.md`
