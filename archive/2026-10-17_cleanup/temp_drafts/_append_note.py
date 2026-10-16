import pathlib

p = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")
NOTE = """

## 刚度文献标定（2026-10-16，用户问「有没有文献可以标定真实刚度」）

**按惯例：结论先落主题件** → `memory/2026-10-16/v1-stiffness-literature-calibration.md`；
脚本 `temp/v1_kcalib.py`（只读 `climbing.pad`，零 FEBio，**无仓库改动**）。

三句：① **有文献** —— `pad.py:194` 的 `k_flex` 定义（踝+膝+髋串联等效刚度）与文献 leg stiffness
**同类量**；实测 **82.4 / 107.7 / 243.3 kN/m**（`PMC9465758`，24 名运动员 0.42 m 落箱力台，soft/self/stiff）
⇒ **`controlled-drop` 220 kN/m = 0.90×**（✅ 已在带内），**`feet-first-stiff` 1100 kN/m = 4.52×**（❌ 超带）。
② **但标定把绝对量推更远**：用文献 k 反算首折 **0.14–0.57 m**（比 6–9 m 低 **11–64×**）；
同高度（0.42 m）GRF 高 **3.3×**（`controlled-drop` 17.68 BW vs 实测 stiff 5.36 BW）。
③ **★★★★ 真机制**：文献的「刚度」是**早期瞬态**刚度（峰值力在 ~50 ms、COM 才走 2–3 cm），
真腿随后**屈服**（肌肉离心收缩）在 **0.23–0.43 m** 行程上继续吸能 —— 线性 220 kN/m 弹簧走 230 mm
需要 **50 kN**，所以给不出 5.36 BW。**缺口不是 k，是「未含主动肌肉」**（论文关键局限 (4)）。
量级指示（待验）：3.3× 力差 ⇒ `h* ∝ F²` ⇒ 首折 **×11～42**（~1 m → 11–42 m）
—— 唯一有外部锚定、且**不是**凑靶的路线（参数由力台数据外部定）。
外部锚定：坠落伤高度 **1.1–4 m 占 46.1%**（`PMC5717375`）、**≤5 ft（1.52 m）亦可致跟骨骨折**
⇒ 流行病学站「首折 ~1 m」侧，**与 [Y25]「群集 1（1~6 m）＝无骨折」冲突**
（后者要求足部耐受 30–50 kN，比尸体 6.2 kN 高 5–8×，项目已记）。

**不动任何判定**：V1 仍 ❌【口径不可复现】＋【绝对量不可引用】；
⛔ 不得用「把 `k_flex` 调到文献值」使 V1 进 6–9 m（**方向相反**）。

**锚点**：`memory/2026-10-16/v1-stiffness-literature-calibration.md`、`temp/v1_kcalib.py`、
`PMC9465758`、`PMC7375`、`pad.py:194`、`pad.py:290-298`。
"""

with p.open("a", encoding="utf-8") as f:
    f.write(NOTE)
print("appended", len(NOTE.encode("utf-8")), "bytes;", p.stat().st_size, "total")
