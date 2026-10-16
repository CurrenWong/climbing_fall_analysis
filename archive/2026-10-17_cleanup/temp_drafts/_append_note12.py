import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑫ 🔴 §⑪ **作废**（替换式）—— 不是能量不守恒，是 `pad.py` 能量账的**口径混乱**（2026-10-16，提交 `067e8b9`）

**§⑪ 的结论「软接触下 `pad.py` 能量不守恒」错误，已作废。** 真机制（**读码 `pad.py:733-741` 定案**）：

```python
p_flex = f_flex[sl] * np.gradient(xf, t_out)[sl]     # f_flex = k·xf(+k_stop 段) + c_flex·rate
e_flex_t = cumtrapz(p_flex, t)                       # 符号变号 ⇒ 弹性项被积掉，剩≈阻尼耗散
k = argmax(e_pad_t + e_flex_t); e_flex = e_flex_t[k]
e_impact = m * g * height_m                          # 分母只对应"首次接触那一刻"
```
**⇒ 分子含不可恢复的阻尼耗散**（`c_flex=400 N·s/m`），随身体**继续下沉**累积、**不受 `mgh` 约束**；
**分母 `mgh` 只对应首次接触** ⇒ **不同物理量 + 不同时间基准** ⇒
`absorbed_fraction > 1`（实测 **1.376–1.583**）**是预期行为，不是缺陷**。
（`pad.py` 注释自承历史上有过「吸能超过输入能量」的自相矛盾读法 —— 实现改了，**口径问题仍在**。）

### 三项独立检验（`temp/` 只读脚本）
1. `v1_energy_verdict.py`：用**同时刻**能量平衡重算 ⇒ 峰值 `xf` 时刻可用 PE **595.8 J ≥ 弹性 519 J** ⇒ **不越界**（推翻 §⑪ 依据）。
2. `v1_emechan.py`：沿轨迹 `E_mech(t)` —— **该版检验因速度偏移无效**
   （`result.t_s` **只覆盖冲击窗**，窗起点速度是 `v_impact=√(2gh)` 而非 0 ⇒ 必须先偏移，否则漏检）⇒ **已写进 gotcha**。
3. `v1_conv_check.py`：非单调格的收敛检验 —— 三段解算器设置
   （`1e-8/1e-11/1e-4`、`1e-9/1e-12/3e-5`、`1e-10/1e-13/1e-5`）**逐位一致**
   （`anchored-self`@`2.1e6`：2 m **9.309** / 5 m **7.781** / 9 m **25.979 kN**）
   ⇒ **`F_flex(h)` 的非单调是收敛的模型行为，不是积分精度问题**（**未定案**：可疑 `A(κ)` 增长 / `k_stop` / `c_flex·rate` 占比）。

### 状态更新
- **`FP_CONTACT` 判定：【撤回】→【未定 · 判据口径待立】**（数值未被证明不物理，但也未被正面验证；
  重判前须先立**不依赖 `absorbed_fraction`** 的有效性判据：沿轨迹 `E_mech` 非增 + **速度须带 `v_impact` 偏移**）。
- **未受影响**：**`LEG_K_ANCHOR`（门关）全部有效 ⇒ 首折 2.244 / 3.216 / 3.418 m 仍可引用**；
  `hard_surface` 语义发现（`k_eff = k_n × A(κ)` = 548 kN/m）仍有效。

### ⚠️ 本会话自错模式（值得记）
**同日三次翻转**：A 首判 → 撤回 → 更正。**两次错都出在「用比值/单点判物理」，而不是用同时刻守恒量。**
⇒ 教训：**判「数值是否物理」时，必须比较同一个时间基准的量，并明确写出参照系**（与项目既有的
「方向/符号结论须有不可被约定翻转的证据」同源）。

**锚点**：`docs/notes/PAD_GOTCHAS.md` §6/§6b、`results/FP_CONTACT_REPORT.md`（未定）、
`temp/v1_energy_verdict.py`、`temp/v1_emechan.py`、`temp/v1_conv_check.py`、`pad.py:733-741`、提交 `067e8b9`
"""

PTR = """

## 🔴 §⑪ 作废（替换式）：不是能量不守恒，是 `pad.py` 能量账的**口径混乱**（提交 `067e8b9`）

**读码 `pad.py:733-741` 定案**：`absorbed_fraction` 的分子 `cumtrapz(f_flex·d(xf)/dt, t)` 里
`f_flex` **含阻尼项 `c_flex·rate`**（不可恢复、随继续下沉累积、**不受 mgh 约束**），
而分母 `impact_energy_j = m·g·h` 只对应**首次接触那一刻** ⇒ **不同量 + 不同时间基准**
⇒ `>1`（实测 1.376–1.583）**正常**。⛔ 不得据此判「数值不物理」；⛔ `energy_residual_j` 同源污染。
**三项检验**：① 同时刻能量平衡 ⇒ 不越界（推翻撤回依据）；② 沿轨迹 `E_mech` 检验**无效**（`result.t_s`
只覆盖冲击窗，窗起点速度是 `v_impact` 而非 0 —— 陷阱已入库）；③ 收敛检验：三段解算器设置**逐位一致**
⇒ `F_flex(h)` 非单调（self@2.1e6：2 m **9.309** / 5 m **7.781 kN**）是**收敛的模型行为**，未定案。
**状态**：`FP_CONTACT` 判定 【撤回】→**【未定 · 判据口径待立】**；**未受影响**：`LEG_K_ANCHOR`（门关）
全部有效 ⇒ **首折 2.244/3.216/3.418 m 仍可引用**；`hard_surface` 语义（`k_eff = k_n × A(κ)`=548 kN/m）仍有效。
⚠️ **本会话自错模式**：**同日三次翻转**（A 首判 → 撤回 → 更正），两次错都出在**用比值/单点判物理**，
而非同时刻守恒量 ⇒ **判「数值是否物理」必须同时间基准 + 明写参照系**。
详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑫、`docs/notes/PAD_GOTCHAS.md` §6/§6b。
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
