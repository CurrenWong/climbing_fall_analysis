"""落盘：ERR-20261005-008（lc 挂错）+ LRN-071 更新。"""
import io
import re

R = r"D:\Project\climbing_fall_analysis"
ERR = R + r"\.learnings\ERRORS.md"
LEARN = R + r"\.learnings\LEARNINGS.md"

SNIP_E = """
---

## [ERR-20261005-008] FEBio `<traction lc="...">` 被静默忽略 —— 载荷曲线必须挂在 `<scale>`

- **症状**：多步分析 step2（力控制）**无论总载荷多大**（100 / 600 / 5000 N 实测），
  都在第一个增量报 42~48 个负 Jacobian，且**看似**"卡在 11~12 N"。
  更早的表现：`DRIVER=force` 历史上一律失败，看起来像"力控制病态"。
- **根因**：写成 `<traction lc="2">-0,-0,-9.93</traction>`。
  **FEBio 不认 `<traction>` 上的 `lc`，静默忽略 ⇒ 满载荷在第一增量瞬加。**
- **正解**（手册 §3.13.2.2 / p121 原文示例）：
  ```xml
  <surface_load type="traction" surface="tibia_top">
    <scale lc="2">1.0</scale>      <!-- ★ 载荷曲线挂在 scale 上 -->
    <traction>0,0,-9.93</traction>
  </surface_load>
  ```
  原文："An optional load curve can be defined for the **scale** element using the `lc` attribute."
- **诊断信号**（有用）：**失败点对总载荷量级不敏感** ⇒ 一定是"满值瞬加"而不是"斜坡"。
  用 `载荷 = FORCE_N × (t-1)` 反推会得到假象的"绝对值天花板"。
- **影响面**：此前所有 `DRIVER=force` 实验都未真正测过"斜坡力控制"。
- **修法**：`ankle_lig_feb.py` 的 traction 输出改为 `<scale lc=..>1.0</scale>` + `<traction>`；
  单步路径同步加 `<scale lc="1">1.0</scale>`。
"""

SNIP_L = """
---

## [LRN-20261005-071] 多步落位第二轮：`lc` 挂错点已修；确认"撤约束回弹"只是部分原因

- **承接** `[LRN-20261005-070]`（step1 落位成功；step2 失败）

### 已验证正确的东西（不是问题）

- **多步机制完全正确**：`FORCE_N≈0` 时模型 **NORMAL 跑完 43 个时间点**，
  且 step2 的顶面支反力 ≈ 0 ⇒ **step1 的指定位效约束在 step2 被正确解除** ✓
- **落位量与实际反力成对可测**（SEAT_MM → step1 末胫骨顶 Rz）：
  `0.02→31.2`、`0.05→43.7`、`0.10→65.2`、`0.20→109.0`、`0.30→155.0`、`0.50→245.1` N
  （腓骨顶另有约 37% 于胫骨的值）⇒ **0.02 mm 就有 31 N ⇒ 关节在 t=0 几乎已接触**，
  不需要改几何（推翻 `[LRN-070]` 里"间隙 1 mm 压不闭"的推测）。
- **诊断口径修正**：FEBio 的 `reac_*.txt` 是**每节点一行**，
  文件 = 节点数 × 时间点；须 `reshape(nstep, nnode, 3).sum(axis=1)` 才是合力。
  第一列是**节点 ID 不是时间**。

### 本轮结论

- **`lc` 挂错是主因之一**（`[ERR-20261005-008]`）：已改为 `<scale lc="2">1.0</scale>`，
  .feb 已确认写出正确形式。
- **`SEAT_HOLD=1`（step2 保持落位位移、不撤约束）确实有效果**：
  step2 从 `t=1.01818`（10.9 N）推进到 **`t=1.055`（33.0 N）** ⇒
  **"撤约束后弹性回弹、关节张开"是一个真实因素**，但它**不是全部**。
- **仍未通过**：step2 在 ~33 N 处仍失败（`Negative jacobian detected`）。
  下一步候选（按优先级）：
  1. **按 Anderson 原设计：把落位/约束加在「距骨」而不是胫骨顶** ——
     距骨垂直位移作为**常驻运动学把手**，胫骨顶**从不被位移约束**、只受载荷 ⇒
     不存在"撤约束回弹"，且与 Anderson §2.2 一致。
  2. 查 step1 收敛态是否本身含翻负单元（打开 `ls_check_jacobians`）。
  3. 软骨网格加密（`cart_patch` m≥2 目前建网格即失败）。
"""

for path, snip, key in ((ERR, SNIP_E, "ERR-20261005-008"),
                        (LEARN, SNIP_L, "LRN-20261005-071")):
    txt = io.open(path, encoding="utf-8").read()
    n = len(re.findall(re.escape(key), txt))
    if n:
        print(f"★ {key} 已存在 {n} 次 ⇒ 跳过")
        continue
    io.open(path, "a", encoding="utf-8").write(snip)
    after = io.open(path, encoding="utf-8").read()
    print(f"{path.split(chr(92))[-1]}: {len(txt)} -> {len(after)} 字符；"
          f"{key} × {after.count(key)}")
