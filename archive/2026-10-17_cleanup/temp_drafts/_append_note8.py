import pathlib

p = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                 r"\v1-stiffness-literature-calibration.md")
NOTE = """

## ⑧ 「加」= 加膝段 → 实际落成「力台锚定预注册」（2026-10-16，用户答「加」）

**用户「加」的是我在上一条末尾提的「膝扭簧 1D 两段腿预注册提案」。执行中发现两处自错 + 一个结构判决，
提案的**内容**因此被替换（**加的是「外部锚定」，不是自由度**）。**

### 落盘产物
- `docs/methods/P3_腿刚度力台锚定_预注册.md`（新，**v2 替换式**；结构照抄 `G1_多段链设计.md`/`G4_压缩相本构设计.md`）
- `LEG_K_ANCHOR_TODO.md`（session 笔记，`*_TODO.md` 已 gitignore）
- `temp/v1_kcalib.py`、`temp/v1_fpanchor_predry.py`（只读）
- **仓库零修改**（`git status --short` 只有 1 个未跟踪新件）

### ★★ 自错一：单位读错（替换式）
`PMC9465758` 原文公式逐字 = **`Stiffness = GRFmax / Depth`**，而 `GRFmax` 以 **N·kg⁻¹** 报 ⇒
表值 82.36/107.71/243.28 的量纲是 **N·kg⁻¹·m⁻¹**（复算 77.05/98.00/228.48，偏差 7–10% 吻合）
⇒ **换算物理刚度须 ×76 kg ⇒ 6.3 / 8.2 / 18.5 kN/m**（能量法 `k=F_peak²/2E` 独立反推 17.4–25.5 ✓）。
⇒ **`controlled-drop` 220 kN/m = 12×、`feet-first-stiff` 1100 = 59×**（vs 18.5）。
⛔ 同日早前「`controlled-drop` 0.90× 已合格 / `feet-first-stiff` 4.52×」**作废**。

### ★★ 自错二：观测量选错（替换式）—— **新 gotcha**
预演实测（@0.42 m、80 kg、硬地面）：`F_pad` = **8.61 / 8.61 / 8.63 kN**（k 变 2.9× 而力只动 0.2%）
⇒ **硬地面下 `peak_force_n` 由 `hard_surface(k_n=1e7)` × `m_low` 的刚性接触尖峰决定，对 `k_flex` 不敏感**
⇒ **不能当 GRF / 腿刚度观测量**；观测量须用 `meta["f_flex_peak_n"]`。
⛔ 早前「模型 @0.42 m 17.68 BW vs 实测 5.36 BW（3.3×）⇒ 缺屈服吸能」**不得再引**（该数出自伪影）。
判据方向也写反过：**`h* ∝ 1/k` ⇒ k↑ ⇒ h*↓**（实测 3.418/3.216/2.244 m ✓）。

### ★★★ 结构判决（本件最有价值的部分）
| 档 | k(kN/m) | `F_pad`/BW | **`F_flex`/BW** | 实测靶 | ratio | h*(m) |
|---|---:|---:|---:|---:|---:|---:|
| soft | 6.3 | 10.97 | 2.90 | 3.38 | 0.86 | 3.418 |
| self | 8.2 | 10.98 | 2.98 | 3.70 | 0.81 | 3.216 |
| stiff | 18.5 | 10.99 | **3.54** | 5.36 | **0.66** | 2.244 |

- 两个候选观测量**骑跨**实测（`F_flex` 低 14–34%、`F_pad` 高 105%）⇒
  **模型「力的绝对值」不可靠，根因在「地面-足部接触层」（非 k）** ⇒ 这是「加」的正当理由（非凑靶）。
- **★ 正面**：锚定 k 把力判据首折从旧 k 的 **~0.25 m** 移到 **2.24–3.42 m** ⇒
  **落进流行病学 1.1–4 m 带**（首次）；对靶 6–9 m 仍差 **~2–3×**。
- 预注册 §5 预判：**J1 FAIL**（J1a stiff 0.66 出界、J1b 2.05×）、**J2 FAIL**（2.24–3.42 > 1.52 m，超 1.5–2.3×）、J3 PASS。
- §9 下一件：**A 足部接触柔度（跟垫 E=185 kPa/t=15.4 mm，复用 G5/`chain2d_tissue.py`）优先于 B 两段腿+膝扭簧**
  （膝 quasi-stiffness 锚 0.08 Nm·kg⁻¹·deg⁻¹ ≈ 350 Nm/rad，Shamaei 2013 ~400 Nm/rad）。

### 验证脚本现状（跑过）
`verify_paper_links` **PASS**（78/60/39 checked, 0 dead）｜`verify_paper_chain` **PASS**｜
`verify_paper_crossloc` **PASS**｜`verify_index_linehits` **FAIL —— 既有、与本次无关**：
`论文v2_出处索引.md` 记录的**代码行号**已漂移（**13 mismatch / 7 matched**，涉及 `pad.py`/`pad2d.py`），
且该脚本**不属于**项目已验证的四条（links/chain/numbers/crossloc）⇒ 报告为既有技术债，**未擅自修改**。

### 待用户
签 §5（判据带 ±30%、J2 ≤1.52 m）+ §7（k 锚点/观测量/网格）⇒ 才实现 `pad.py` 门控增量与正式判定件。

**锚点**：`docs/methods/P3_腿刚度力台锚定_预注册.md`、`LEG_K_ANCHOR_TODO.md`、
`temp/v1_fpanchor_predry.py`、`verify_index_linehits.py`、`pad.py:194`/`:566-579`
"""
with p.open("a", encoding="utf-8") as f:
    f.write(NOTE)
print("ok", p.stat().st_size)
