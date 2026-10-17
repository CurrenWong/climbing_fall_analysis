# CANONICAL —— 唯一真相源清单（Single Source of Truth）

> **规则**
> 1. 引用任何数字前**先查本表**；表中「作废」的数字**禁止引用**。
> 2. 每次跑完探针 ⇒ **更新本表** + 把被替代的产物移进 `archive/<日期>_<事由>/`。
> 3. 本表与 `pad.py` 注释 / `docs/notes/PAD_GOTCHAS.md` 三处必须一致。

---

## 0. ★ 接口契约（违反 ⇒ 静默算错，务必先读）

| 接口 | **契约** | 违反后果 | 守卫 |
|---|---|---|---|
| `TWO_SEG_LEG["m_theta_table"]` | `[[θ_deg, M_Nm], ...]`：**2 维、两列、θ 用度数、严格升序** | 传 `(x_rad_array, y_array)` 元组 ⇒ 变 2×N ⇒ 角度二次 `radians` ⇒ 插值撞末端 ⇒ **力矩被钉成常数** | `pad.py` 硬校验（2026-10-17）+ `tests/test_pad_gates_guard.py` |
| 腿力观测量 | `result.meta["f_flex_peak_n"]`，且必须与 `rhs` 用**同一个** `flex_force` | 旧版用线性 `k_flex·xf` ⇒ 与动力学不符 | `pad.py`（已修 `e19e758`） |
| `fp_anchor_probe.json` | **7 个脚本的输入锚**（J0 正对照的期望值） | 移动/删除 ⇒ 断链 | `scripts/tools/archive_cleanup.py` 引用守卫**拒移** |
| `Posture` 传参 | 传**字符串键**才触发 `RIGID_STANDING_CHAIN` 等按「键」判定的门 | 传 `Posture` 对象 ⇒ **静默绕门** | `pad.py` 注释（已知风险） |

---

## 1. pad 线（1D 落地）—— 现行数字

| 数字 | 现行来源 | 备注 |
|---|---|---|
| 力台锚定三档 `k_flex` = 6.3 / 8.2 / 18.5 kN/m | `scripts/1d/fp_anchor_probe.py` → `results/fp_anchor_probe.json` + `results/FP_ANCHOR_REPORT.md` | `PMC9465758` |
| 跟垫柔度 / 硬地面接触 | `scripts/1d/fp_contact_probe.py` → `results/FP_CONTACT_REPORT.md` | 🟡 仅跟垫几何 |
| 行程锚 `θ_max`（COM 深度 0.43/0.37/0.23 m） | `scripts/1d/fp_stroke_probe.py` → `results/FP_STROKE_REPORT.md` | 外部锚实测 |
| 卸载过程 `D1`（力-长 + 率）| `scripts/1d/unload_process_probe.py` → `results/unload_process_probe.json` | **表契约已修正** |
| **首折 `h*`** | **`0.402 m`** —— `scripts/1d/recompute_hstar.py` → `results/hstar_canonical.json`（自检①正面 `interp(1 rad)=254.7`、②反面拒错格式） | ⚠️ **soft/self/stiff 三档同值** ⇒ 门路径下「力台锚定」已失效（`k_flex` 被 `M(θ)` 取代）⇒ 三档只在 `θ > θ_max` 后分化 |
| 二次接触可达性 / 能量-行程账 | `results/SECOND_CONTACT_REPORT.md` §1/§2 | J1 PASS、Q2 缺口 0.46–0.71 m |

---

## 2. 已作废（禁引）

| 作废陈述 | 作废原因 | 记录处 |
|---|---|---|
| `h*` = 1.34–3.66 m | 后处理用线性观测量（与 ODE 力律不符） | `SECOND_CONTACT_REPORT.md` §0b |
| `h*` = >9 / 7.17 / 2.09 m | `m_theta_table` 传错格式（力矩被钉成常数） | §0c/§0d |
| 「与流行病学 1.1–4 m 一致」 | 随上两条一起作废 | §0b/§0d |
| 旧「波杆佐证接触压溃主导」 | 已撤回 | `MEMORY` |
| V1「实测 >40 m」/ 绝对量 | 口径不可复现 | `docs/reports/最终结论_2026-10-10.md` |

---

## 3. 归档约定

- 目录：`archive/YYYY-MM-DD_<事由>/`，内部按 `results_invalid/`、`scripts_invalid/`、`docs_stale/`、`temp_drafts/` 分。
- **只 move，不删除**（可恢复）。
- 每次归档必须生成 `MANIFEST.md`（原路径 → 新路径 + 事由 + 被替代者）。
- 归档前自动跑**引用守卫**：若某文件名在 `archive/` 之外仍被引用 ⇒ **拒移**并打印引用方。
