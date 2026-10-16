import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑩ A 件 `FP_CONTACT` 执行完成（2026-10-16，用户「继续A」）

### 交付（提交 `95c2bc7`，7 文件 / +782 −3）
`docs/methods/P3_足部接触柔度_预注册.md`｜`scripts/fp_contact_probe.py`｜
`results/fp_contact_probe.json`＋`results/FP_CONTACT_REPORT.md`｜
`tests/test_climbing.py::TestRigidContactGate`（3 条）｜`docs/notes/PAD_GOTCHAS.md` §5｜`docs/reports/总体进度.md` §九之二。
**★ `src/climbing/pad.py` 零改动**（`git status` 证）—— 只用既有门 `RIGID_CONTACT["k_n_per_m"]`（默认 `None`）。

### ★★★★ 又一条语义 gotcha：`hard_surface` 的 `k_n` **不是**总刚度
`CrashPad.force_n = stress_pa(xi) × area`，而 `hard_surface` 把 `k_n_per_m` 喂给
`linear_modulus_pa`（**Pa**，`thickness=1`）⇒ **有效刚度 `k_eff = k_n × A(κ)`**。
实测 `k_n=1.0e7` ⇒ `F_pad` 8.625 kN、`x` 15.728 mm ⇒ **`k_eff` = 548.4 kN/m**（`k_eff/k_n` = 0.055–0.16）
⇒ **docstring 声称的「1.0e7 N/m」（18×）比实际硬。引用必须报实测 `k_eff`，不得引 `k_n`。**

### 判定（`results/FP_CONTACT_REPORT.md`）
**A-J0 PASS（逐位复现上游 `fp_anchor_probe.json`）｜A-J1 FAIL｜A-J2 PASS｜A-J4 PASS**
| 判据 | 结果 |
|---|---|
| **A-J1**（锚带 `k_eff`∈[100,175]，判在 `k_n`=2.1e6 ⇒ `k_eff`=127.7 kN/m） | soft **1.691❌** / self **1.547❌** / **stiff 1.077✅** |
| **A-J2** `F_pad/F_flex` 单调 + 上端 ≤2.50 | **PASS**（3.79/3.67/3.10 → 2.30/2.19/1.49，三档单调） |
| **A-J3** `F_flex` 如实报（无阈值） | soft 2.89→2.83、self 2.99→3.11、**stiff 3.54→4.31（J1a 0.661→0.805✅）** |

### 结论（本件的主要收获）
1. **★ 诊断成立（`stiff` 档）**：换实测跟垫后 stiff 档 `J1b` **2.051→1.077✅**、`J1a` **0.661→0.805✅**、
   `F_pad/F_flex` **3.10→1.49** ⇒ **「`F_pad` 超 105–225% ＝ 刚性接触律尖峰伪影」被证实**。
2. **⇒ 但 `soft`/`self` 仍偏高 1.55–1.69× ⇒ 剩余缺口不在接触层** —— 我们三档**只差 `k_flex`**，
   而真人 soft/self 的差别是**行程/运动学**（实测 COM 深度 **0.43/0.37/0.23 m**）
   ⇒ **模型类别不足 ⇒ 下一件 B（运动学自由度）的必要性坐实**。
3. ⛔ **反面读数（不可采纳）**：`k_n` 降到 1.0e6（`k_eff`=66，**锚带外**）三档全 PASS ⇒ 调参凑靶
   （A-J5）。**锚带内三档全 PASS 不可达。**
4. 不改变 V1 判定；不改变 `LEG_K_ANCHOR` 判定（后者是**门关**口径）。

### 回归
**768 passed / 8 skipped / 0 failed**（-n auto --dist load；+3 条）。

**锚点**：`results/FP_CONTACT_REPORT.md`、`results/fp_contact_probe.json`、`scripts/fp_contact_probe.py`、
`docs/methods/P3_足部接触柔度_预注册.md`、`docs/notes/PAD_GOTCHAS.md` §5、
`pad.py`（`hard_surface` / `RIGID_CONTACT`）、`tests/test_climbing.py::TestRigidContactGate`、提交 `95c2bc7`
"""

PTR = """

## A 件 `FP_CONTACT` 完成（2026-10-16，用户「继续A」，提交 `95c2bc7`）

把刚性地面**有效**接触刚度（实测 **548 kN/m**，注意 `hard_surface` 的 `k_n` 被**乘面积**：
`k_eff = k_n × A(κ)`，docstring 声称的 1e7 N/m 比实际硬 ~18×）换成**实测跟垫 in-vivo 100–175 N/mm**
⇒ 判定 **A-J0 PASS（逐位）｜A-J1 FAIL｜A-J2 PASS**。
★ **`stiff` 档两个判据同时修好**（J1b 2.051→**1.077✅**、J1a 0.661→**0.805✅**，比值 3.10→1.49）
⇒ **「`F_pad` 超 105–225% ＝ 刚性接触律伪影」被证实**；
**但 `soft`/`self` 仍偏高 1.55–1.69×** ⇒ 因为三档**只差 `k_flex`**，真人 soft/self 的差别是**运动学**
（COM 深度 0.43/0.37/0.23 m）⇒ **下一件 B（运动学自由度/膝扭簧）的必要性坐实**。
**`pad.py` 零改动**（只用既有门 `RIGID_CONTACT`）；回归 **768 passed / 8 skipped**。
⛔ 反例（不可采纳）：`k_n`=1.0e6（`k_eff`=66，**锚带外**）三档全 PASS ＝ 凑靶。
**详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑩。**
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
