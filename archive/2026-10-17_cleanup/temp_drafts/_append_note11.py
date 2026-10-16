import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑪ 🔴🔴 A 件判定**撤回** —— 自检发现软接触下 `pad.py` 能量不守恒（2026-10-16，提交 `0fd4b8d`）

用户问「现在高度到多少了」时，为算带接触层的 `h*` 而扫高度网格，**当场发现**：

### 缺陷
`absorbed_fraction = (energy_into_pad_j + energy_into_flex_j) / impact_energy_j`
（`pad.py:479-482`）**物理上不可能 > 1**，但实测（@0.42 m、80 kg、`anchored-*`）：

| 接触设定 | `k_eff` | soft | self | stiff |
|---|---:|---:|---:|---:|
| **None**（=1e7，现行刚性） | 548 kN/m | 0.276 ✅ | 0.289 ✅ | 0.329 ✅ |
| `k_n=3.0e6` | 177 | 0.584 ✅ | 0.624 ✅ | **1.358 ❌** |
| **`k_n=2.1e6`（= A-J1 判定档）** | 128 | **1.568 ❌** | **1.513 ❌** | **1.376 ❌** |
| `k_n=1.6e6` | 100 | **1.583 ❌** | **1.526 ❌** | **1.391 ❌** |

**⇒ 接触层越软越易越界；同设定下随 `k_flex` 变。**（h=0.5 m 亦见：`into_flex` 541.6 J vs 输入 392 J）

**症状二（同源）**：越界区 `F_flex(h)` **非单调** —— `anchored-self`@`k_n=2.1e6`：
**2 m 给 9.31 kN，5 m 给 7.78 kN**（同档 `k_n=None` 下 5.03→7.91 单调 ✓）⇒ **`h*` 不可引用**。

### 撤回范围（提交 `0fd4b8d`，替换式 banner）
`results/FP_CONTACT_REPORT.md` 的 **A-J1/A-J2/A-J3 与 §4 全部结论**、`总体进度.md` §九之二 对应判定。
⛔ **缺陷修复前，「接触层软化」的一切后果（含 `stiff` 档「两判据修好」）不得引用。**

### 未受影响（重要）
- **`LEG_K_ANCHOR`（`FP_ANCHOR_REPORT.md`，门关 `k_n=None`）全部有效**：
  `absorbed` = 0.276/0.289/0.329 ≤1 ✓、`energy_residual_j` = 238.5/234.3/221.3 J >0 ✓、
  `F_flex(h)` 单调 ✓ ⇒ **首折 2.244/3.216/3.418 m 仍可引用**。
- **`hard_surface` 语义发现**（`k_eff = k_n × A(κ)`，实测 548 kN/m，docstring 声称 1e7 比实际硬 ~18×）
  —— 与能量账无关，**仍然有效**。

### 嫌疑机制（**均未定案，待查**）
① `A(κ)` 渐进接触面积随压缩增长；② `k_stop=3.0e5` 止挡；③ 初始速度与落高是否**双重计入**；
④ `energy_into_flex_j` 的积分口径。
**⇒ 下一步优先级已改变：先定位并修复能量不守恒，之后才重跑 A、再谈 B。**

**锚点**：`temp/v1_contactA_selfcheck.py`、`temp/v1_contactA_energy.py`、`temp/v1_contactA_hstar.py`、
`results/FP_CONTACT_REPORT.md`（撤回）、`docs/notes/PAD_GOTCHAS.md` §6、提交 `0fd4b8d`
"""

PTR = """

## 🔴 A 件判定撤回（2026-10-16，提交 `0fd4b8d`）

回答「现在高度到多少了」时扫高度网格，**当场发现 `pad.py` 软接触下能量不守恒**：
`absorbed_fraction`（物理上必 ≤1）实测 **1.358–1.583**（`k_n=2.1e6/1.6e6` 三档全越界；
`None` 档 0.276–0.329 正常）。**⇒ A-J1 判定档正落在越界区 ⇒ A 的判定全部撤回**
（`results/FP_CONTACT_REPORT.md` + `总体进度.md` §九之二 已加替换式 banner）。
同源症状：越界区 `F_flex(h)` **非单调**（self@2.1e6：2 m 9.31 kN 而 5 m 7.78 kN）⇒ `h*` 不可引用。
**未受影响**：① **`LEG_K_ANCHOR`（门关）全部有效 ⇒ 首折 2.244/3.216/3.418 m 仍可引用**
（`absorbed` 0.276–0.329、`resid` >0、`F_flex(h)` 单调）；② `hard_surface` 语义发现
（`k_eff = k_n × A(κ)` = 548 kN/m）仍有效。
嫌疑机制未定案：`A(κ)` 渐进面积 / `k_stop` 止挡 / 初始速度与落高双重计入 / 能量积分口径。
**⇒ 下一步改为：先修能量不守恒，再重跑 A，之后才谈 B。**
详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑪、`docs/notes/PAD_GOTCHAS.md` §6。
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
