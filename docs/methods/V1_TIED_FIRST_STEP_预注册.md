# V1 跖面 `tied`/`sticky` 首步 —— 预注册（2026-10-17）

> **性质**：预注册。**以下判据在跑任何新算例前写死**，跑后只做判定，不改判据、不改参数。
> **上游事实**：`V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` §2.3/§2.4、`V1_OPTIMIZATION_OPTIONS.md` §1/§2.3。

---

## §0 · 问题与已排除项

**要解决的事**：in-vivo 跟垫总刚度 **100–175 N/mm**（per-node ≈ **0.203 N/mm**，粗网格 479 节点）
在本模型内**不收敛**（`V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` §2.3：k = 0.05 / 0.2028 / 0.5 / 2.0 全灭）。

**真因（已定）**：**欠约束刚体模态** —— 弹簧太软 ⇒ 跖面几乎无支撑 ⇒ 切向/法向零刚度方向
⇒ 负 Jacobian。手册 §8.1.3 / §8.3.1.5 明示。**不是**「物理上跟垫不该这么软」。

**已排除（不要重试）**：

| 路 | 结果 | 出处 |
|---|---|---|
| 切向释放 `roller_free`（k_t=1 N/mm） | gauge **+14%（更差）** | `PLANTAR_BC_REPORT.md` |
| 跟垫模量 `pad_e=0.185` / 厚 15.4 / 板可变 `bottom_only` | **逐位相同（板模量不参与）** | `V1_OPTIMIZATION_OPTIONS.md` §1 |
| `penalty=5000` / `auto_penalty=1` | **破坏收敛（负 J）** | 同上 |
| 手册五条处方（初始接触 / 两步法 / False Transient / `symmetric_stiffness` / `auto_penalty`） | 均未让 `hold=None` 收敛 | `V1_UNDERCONSTRAINT_METHODS_SURVEY.md` |
| 跖面 BC 族全成员 | ＝Winkler 弹簧床及其 k→∞ 极限，无接触承载 | `V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md` |

---

## §1 · 假设（跑前写死）

- **H1**：`tied`（或 `sticky`）**首步**在 t=0 强制跖面-支撑面 tie，可消除刚体模态
  ⇒ **k = 0.203 N/mm 收敛**（rc=0，0 个 negative Jacobian）。
- **H2**：若 H1 成立，tie 释放后由弹簧承载 ⇒ `gauge_max(k)` 在可解区间**单调**
  （即不再是 §3 记录的 U 形洼地）。
- **H0（待否证）**：tie 只是把刚体模态推迟到释放步 ⇒ 释放后仍不收敛。

---

## §2 · 判据（跑前写死，不可事后调整）

| 编号 | 判据 | 通过条件 |
|---|---|---|
| **J0（正对照）** | k=1000 复现既有值 | `gauge_max = 142.729 MPa ± 1%`（`V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` §2.3）|
| **J1（可解性）** | k ∈ {0.2028, 1.0, 10.0} 三档 | 全部 `rc=0` 且 **negative jacobian 计数 = 0** |
| **J2（单调性）** | 若 J1 过，k ∈ {0.2028, 1, 10, 100, 1000} 的 `gauge_max` | **单调不增**（软 ⇒ 高应力），无内部极小 |
| **J3（可复现）** | k=0.2028 连跑 2 次 | `gauge_max` 逐位一致（±1e-6 相对）|

**判定规则**：J0 不过 ⇒ **工具链作废，全部结果不采用**。
J0 过 ∧ J1 过 ⇒ **H1 成立**，继续 J2。
J0 过 ∧ J1 不过 ⇒ **H0 成立**，本路线记为「不可解」，**不再重试**（避免无限试错）。

---

## §3 · 红线

1. ⛔ **不得**为让首折进 6–9 m 而调 `tied` 刚度 / 脱开阈值 / `k` / `gauge` 阈值。
2. ⛔ **不得**只报成功档、隐藏失败档（三档全报）。
3. ⛔ **不得**把本路线结果与 `hold_k` 族混链比较（口径必须同链，闸门 10）。
4. `tied` 的 tie 刚度若需取值，**必须取自手册 §3.14.9 的默认/明示范围**并在跑前记录，
   不得「试到收敛为止」。
5. 单次运行前须核对：deck 内实际写入的参数（防「传的参数 ≠ deck 里写的 ≠ FEBio 读到的」，闸门 3）。

---

## §4 · 实现要点

- 现状：`grep tied/sticky` 在 `src/climbing/*.py` 与 `scripts/opensim_fe/*.py` **均为空**
  ⇒ 属**新功能**，需在 deck 构造层新增。
- 参考位置：`docs/methods/踝关节FE建模交接.md:859`（「换回 force 驱动：需先解决刚体模态（弱弹簧 / 初始 tied 接触）」）。
- 需改：`src/climbing/plantar_bc.py`（弹性地基 → 可选 `tied` 接触面）或 `thums_feb.py` deck builder。
- 口径固定（与既有扫描同）：`h=5 m`、`a=0`、`anat` 网格、`F_sub=26383.1 N`、`PressureLoad`、ℓ=4 mm。

---

## §5 · 产物

- `results/opensim_fe/V1_TIED_FIRST_STEP.md` + `.json`
- 表列：k / rc / negative-Jac 计数 / `gauge_max` / 是否触发释放
- **不改任何既有判定**（V1 仍为 ❌【口径不可复现】＋【绝对量不可引用】），本件只回答「H1 是否成立」。
