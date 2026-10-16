import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑨ `LEG_K_ANCHOR` 执行完成（2026-10-16，用户「同意，按你说的做」）

### 交付
| 文件 | 内容 |
|---|---|
| `docs/methods/P3_腿刚度力台锚定_预注册.md` | 预注册 v2（判据/锚点/观测量/网格**跑前冻结**）＋ O5 口径守卫 ＋ J0 正对照 |
| `src/climbing/pad.py` | **纯增量**：`_anchored_posture` + `POSTURES["anchored-soft/self/stiff"]`（其余字段 `replace` 自 `controlled-drop`） |
| `tests/test_climbing.py` | `TestAnchoredPostures` **5 条**（冻结 k / 字段不搭车 / 门不生效 / 门语义特征化 / `F_pad` 不敏感） |
| `scripts/fp_anchor_probe.py` | 探针（含**双正对照**，不过则作废） |
| `results/fp_anchor_probe.json`、`results/FP_ANCHOR_REPORT.md` | 产物 + 判定 |
| `docs/notes/PAD_GOTCHAS.md` | **新建**：`pad.py`/链模型 gotcha（4 条） |
| `docs/reports/总体进度.md` | §九 新增 |

### ★★★★ 执行中挖出的最重要一条：**L2 刚性门按「键」判定 ⇒ 传对象静默绕门**
`pad.py` 的 `RIGID_STANDING_CHAIN={"kf_mult":2.0,"x_max":0.015}`（2026-10-08 L2 标定采纳）
**只在 `on_pad=False` 且姿势的键 == `"feet-first-stiff"`** 时生效（k ×2.0、行程限 15 mm）。
门用**键**判定：`posture_key = posture if isinstance(posture, str) else None`
⇒ **传 `Posture` 对象（`replace(...)` 造的）会绕门**。实测同姿势两种调用差 **13.4%**：
`feet-first-stiff` @1 m **传键** `F_pad`/**`F_flex` = 37.46/31.10 kN**；**传对象 = 33.65/27.43 kN**。
⇒ ① `V1_WHY_20KN_AND_GAP_LOCUS.md` §3 里标「1.1e6（模型 `feet-first-stiff`）」的 27.43 是**无门**值
（真实 `feet-first-stiff` 刚性口径是 2.2e6/15 mm）；② 这正是我先前 dry-run 与项目记录不一致的根因；
③ 锚定档键 `anchored-*` ⇒ 门不生效（= 想要的口径）。

### 判定（`results/FP_ANCHOR_REPORT.md`，与 §5 预判**逐位一致**）
**J0 PASS（逐位）｜J1a soft 0.856✅ self 0.808✅ stiff 0.661❌｜J1b 三档全❌(2.05–3.25×)｜
J1 FAIL｜J2 FAIL(2.244/3.216/3.418 m > 1.52)｜J3 PASS(3.418>3.216>2.244)｜J4 PASS**
- `F_pad` 三档 **8.61/8.61/8.63 kN**（对 k 不敏感，再证 §1.5）
- `F_flex` 2.27/2.34/2.78 kN = **2.89/2.99/3.54 BW**（实测靶 3.38/3.70/5.36）
- **结论：k 锚定成功（方向+量级），但力的绝对值**骑跨**实测（`F_flex` −14~−34%、`F_pad` +105~+225%）
  ⇒ 缺口在「地面-足部接触层」，不在腿刚度**
- **★ 正面**：力判据首折 **~0.25 m → 2.24–3.42 m**，**首次落进流行病学 1.1–4 m 带**；靶 6–9 m 仍差 ~2–3×
- 全套回归 **765 passed / 8 skipped / 0 failed**（本件 +5；据此变动前基线 760，记录里 759 略陈旧）

### 未做 / 待办
- `MEMORY.md` **原始字节 10238 / 上限 10240（仅剩 2 B）** ⇒ **未硬塞**（按项目规则需先走
  「备份 → 搬不删 → 验尺寸+无死链+信息未丢」的压缩工序）。**待办：MEMORY.md 压缩后再加本件指针。**
- `verify_index_linehits` **既有 FAIL**（索引记录的代码行号漂移 **13 mismatch / 7 matched**，
  涉及 `pad.py`/`pad2d.py`；该脚本不属项目认证的四条）⇒ **未擅自修**，`fix_index_linehits.py` 可用。
- 下一步（预注册 §9 触发已成立）：**A 足部接触柔度（跟垫，复用 G5/`chain2d_tissue.py`）≻ B 膝扭簧**；
  两者均须另行预注册。

**锚点**：`results/FP_ANCHOR_REPORT.md`、`results/fp_anchor_probe.json`、`scripts/fp_anchor_probe.py`、
`docs/methods/P3_腿刚度力台锚定_预注册.md`、`docs/notes/PAD_GOTCHAS.md`、`pad.py`（`RIGID_STANDING_CHAIN`）、
`tests/test_climbing.py::TestAnchoredPostures`
"""

PTR = """

## `LEG_K_ANCHOR` 完成（2026-10-16，用户「同意，按你说的做」）

力台锚定预注册**已实现并跑完**：`k_flex` 220/1100 kN/m → **6.3/8.2/18.5 kN/m**（`PMC9465758` 表值
量纲是 N·kg⁻¹·m⁻¹，×76 kg 换算）。判定 **J0 PASS / J1 FAIL / J2 FAIL / J3 PASS / J4 PASS**
（与预注册预判**逐位一致**）；`F_pad` 对 k 不敏感（8.61/8.61/8.63 kN）而 `F_flex` 2.89/2.99/3.54 BW
vs 实测 3.38/3.70/5.36 ⇒ **力的绝对值骑跨实测 ⇒ 缺口在「地面-足部接触层」**。
★力判据首折 **~0.25 m → 2.24–3.42 m**，**首次落进流行病学 1.1–4 m 带**。
★★★★ 执行中挖出：**`pad.py` L2 刚性门按「键」判定 ⇒ 传 `Posture` 对象静默绕门**
（同姿势差 13.4%：传键 37.46/31.10 kN vs 传对象 33.65/27.43 kN）。
回归 **765 passed / 8 skipped**；新 gotcha 件 `docs/notes/PAD_GOTCHAS.md`。
⚠️ `MEMORY.md` 仅剩 2 B ⇒ 未塞，待压缩；`verify_index_linehits` 既有 FAIL（行号漂移 13/7）未修。
**详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑨。**
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
