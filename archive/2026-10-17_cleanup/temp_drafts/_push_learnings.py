import io
import os

P = r'D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md'
ENTRY = """
## 2026-10-17 · pad 线同一处返工三次（接口契约未核实）

### [LRN-20261017-130] 改「接口 / 诊断量」前必须读**契约方**的真实代码；A/B 对拍要测**代码行为**而不是**我的假设**

* **类别**：correction（自我更正，返工三次）
* **情境**：`pad.py` 的 `TWO_SEG_LEG["m_theta_table"]` 契约是 `[[θ_deg, M_Nm], ...]`
  （2 维两列、**度数**、升序；`pad.py` 自己 `np.radians(_mt[:, 0])`）。我**没读这一段**
  就断定「探针传的 `_TAB` 是废表」，并自己写脚本传 `(x_rad, y)` 元组 ⇒ 元组被 `np.asarray`
  变成 2×N ⇒ `_mt[:, 0]` 取到的是**整行 x** ⇒ 再 `radians` 一次 ⇒ 角度被压到 0~0.04 rad
  ⇒ `np.interp` 对真实 θ(0.5–1.8 rad) **全部撞到末端** ⇒ 力矩被钉成常数 **222.4 N·m**
  ⇒ 我据此报的 `h* = >9 / 7.17 / 2.09 m` 是**废值**。
* **更糟的一步**：我做了「A/B 对拍」来验证，但**测的是我假设的语义**
  （`np.interp(tab[0], tab[1])`），而不是**代码的真实行为**（`np.radians(arr[:, 0])`）
  ⇒ 对拍「通过」，结论却是反的。
* **后果**：`h*` 在 1.34–3.66 m → <0.5 m → >9/7.17/2.09 m → **0.402 m** 之间来回；
  同一处返工三次，用户三次质疑。
* **处置（写死）**：
  1. **改任何接口/观测量前，先 `read_file` 契约方的实现** —— 不是文档、不是注释、不是我的记忆。
  2. **A/B 对拍必须复刻契约方的代码路径**；不得用「我以为的语义」当对照。
  3. **给易错接口加硬校验**：`pad.py` 现对 `m_theta_table` 的形状/升序/量纲**当场 raise**；
     守卫测试 `tests/test_pad_gates_guard.py`（4 passed）。
  4. **唯一真相源** `CANONICAL.md`（接口契约表 + 现行数字 + 已作废清单 + 归档约定）。
* **一般化**：**「静默算错」比「报错」危险得多。** 凡是被 `np.asarray` **隐式转形**的输入
  （元组 vs 列表对、行 vs 列、度 vs 弧度），**都必须有形状/量纲校验**；
  否则错值会一路流进论文。
* **锚点**：`src/climbing/pad.py` 的 `m_theta_table` 段；`CANONICAL.md` §0；
  `tests/test_pad_gates_guard.py`；`results/SECOND_CONTACT_REPORT.md` §0b/§0c/§0d。

### [LRN-20261017-131] 「诊断副本未与实现同步」—— 同一类错误的**第二次**

* **类别**：correction（**复发**）
* **情境**：① 本次会话早前：改了 `pad.py` 弹性项（加 `r(θ)` 归零因子），**忘了同步探针里的
  重建公式** ⇒ O1 报出 **673.42 kN** 的荒谬值 ⇒ 判定作废一次。
  ② 本次：`meta["f_flex_peak_n"]` 的后处理用**线性** `k_flex·xf`，而 ODE 用 `M(θ)`/`k_θ·θ`
  ⇒ **观测量与动力学不是同一个力律** ⇒ `h*` 全错（已修 `e19e758`：改为调用 `rhs` 同一个
  `flex_force`，单一真相源）。
* **根因**：**同一个物理量在两处各写一份公式。**
* **处置**：诊断量 / 观测量**必须与实现共用同一个函数**；无法共用时**加断言比对**。
* **一般化**：**「复制公式」是技术债的种子** —— 复制的当天就对，改实现的那天就错，
  而错误是**静默**的。
* **锚点**：第一次的锚点 + 本次 `e19e758`；`results/SECOND_CONTACT_REPORT.md` §4。
"""

cur = io.open(P, encoding='utf-8').read() if os.path.isfile(P) else ''
if 'LRN-20261017-130' in cur or 'LRN-20261017-131' in cur:
    print('已存在，跳过')
else:
    io.open(P, 'a', encoding='utf-8').write(ENTRY)
    print('LEARNINGS.md += LRN-20261017-130 / -131')
    print('新大小', os.path.getsize(P))
