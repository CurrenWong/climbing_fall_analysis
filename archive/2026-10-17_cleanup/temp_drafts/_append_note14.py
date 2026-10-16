import pathlib

TOPIC = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16"
                     r"\v1-stiffness-literature-calibration.md")
DAILY = pathlib.Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-16.md")

SEC = """

## ⑭ B1 屈曲行程锚定（2026-10-16，用户「告诉我加完这个条件高度是多少，然后按你的建议推进」，提交 `2fc4851`）

### ★ 用户问的「高度」（B1 加完后）
| 口径 | soft | self | stiff |
|---|---|---|---|
| **刚性（= `LEG_K_ANCHOR`）** | **3.418 m** | **3.216 m** | **2.244 m**（**逐位一致 ⇒ B1 零影响**） |
| **实测跟垫（`k_n=2.1e6`）** | **3.657 m** | **不可定** | **不可定** |

**⇒ 现在可引用的首折 = 2.24–3.66 m**；对 `[Y25]` 靶 7~9 m 仍差 ~2–3×；落流行病学 1.1–4 m 内 ✓。

### 实现
门 `STROKE_ANCHOR = {"x_flex_max_m": None, "k_stop_n_per_m": None}`（`pad.py`；**默认 None ⇒ 既有路径逐位不变**）。
**锚** = `PMC9465758`「Drop jump depth」列（实测 COM 深度），与三档 `k_flex` **同源一一对应**：
soft **0.43** / self **0.37** / stiff **0.23 m**。
⚠️ L2 门排在 L4 门之后 ⇒ `feet-first-stiff` 由 L2 覆盖、`anchored-*` 不受影响。

### 判定
**B1-J0 正对照逐位 PASS｜B1-J1（刚性 `h*` 与 `LEG_K_ANCHOR` 逐位一致）PASS｜
B1-J2（实测跟垫口径，只采信相邻且未触发止段的格）soft 可定 / self·stiff 不可定｜B1-J3/J4 PASS。**

### ★ 结论：`k_stop`（未标定）是软接触区**主控瓶颈**
`stiff`@1 m 3.64 → **6.51 kN** 是假台阶（= 3.3 + 3e5×0.0191 ✓）、`self`@2 m 同理
⇒ **`self`/`stiff` 在有效范围内跨不过 6.2 kN**。
**B1 对现行口径零影响** —— 其价值是**把止段污染格识别出来**。

### ★★ 新 gotcha（我踩了一次）
算 `h*` 内插**不能**「先滤掉止段污染格、再在剩下的格之间内插」—— 会**跨过**无效格：
实测把 `anchored-self` 报成 **3.506 m 的假首折**（真正的跨越段是 1 m→5 m，中间夹着被污染的 2 m 格）。
⇒ 只在**原网格上相邻、且两端都未触发止段**的格之间内插（`fp_stroke_probe.py::_hstar`）。

### 回归 / 交付
**774 passed / 8 skipped / 0 failed**（+3 条 `TestStrokeAnchor`）。
新增 `scripts/fp_stroke_probe.py`、`results/fp_stroke_probe.json`、`results/FP_STROKE_REPORT.md`；
`PAD_GOTCHAS.md` §6d；`总体进度.md` §九之三。

### ⇒ B2（膝扭簧）的**设计分叉** —— 需用户定，不可替他选
把实测 COM 深度（0.43/0.37/0.23 m）当**两段腿的几何缩短量**反解等效屈曲角，会得到
**θ ≈ 120° / 104° / 86°**（取 l=0.43 m：`ΔL = 2l(1−cos(θ/2))`）—— **「stiff landing」也要 86° 膝屈曲，
说明实测 COM 深度并非纯膝**（含踝背屈 + 髋屈曲）⇒ 单膝模型**在几何上装不下**。
**⇒ 两条路**：**(B2a) 等效单关节**（θ 为复合屈曲角，θ_max 由 COM 深度反解 = 几何恒等式；膝 quasi-stiffness
**350 Nm/rad** 仅作量级锚）vs **(B2b) 膝 + 踝两关节**（解剖更真，但参数与不变量都翻倍）。
**我建议 B2a**（先看「加一个几何屈曲自由度」是否就够），B2b 留作备选。

**锚点**：`results/FP_STROKE_REPORT.md`、`scripts/fp_stroke_probe.py`、`results/fp_stroke_probe.json`、
`pad.py`（`STROKE_ANCHOR`）、`docs/notes/PAD_GOTCHAS.md` §6d、提交 `2fc4851`
"""

PTR = """

## B1 屈曲行程锚定完成（提交 `2fc4851`，用户「告诉我加完这个条件高度是多少，然后按你的建议推进」）

**★ 高度（B1 加完后）**：**刚性口径 3.418 / 3.216 / 2.244 m（与 `LEG_K_ANCHOR` 逐位一致 ⇒ B1 零影响）；
实测跟垫口径 soft 3.657 m、self/stiff 不可定（止段在 2.0/1.0 m 触发、有效范围内跨不过 6.2 kN）
⇒ 现在可引用的首折 = 2.24–3.66 m**（对靶 7~9 m 差 ~2–3×）。
**实现**：门 `STROKE_ANCHOR`（默认 None，既有路径逐位不变）；**锚** = `PMC9465758`「Drop jump depth」
实测 COM 深度 soft 0.43 / self 0.37 / stiff 0.23 m（与三档 `k_flex` 同源一一对应）。
**判定**：B1-J0 逐位正对照 PASS｜B1-J1 PASS｜B1-J2 soft 可定、self/stiff 不可定｜B1-J3/J4 PASS。
**★ 结论：`k_stop`（未标定）是软接触区主控瓶颈**（`stiff`@1 m 3.64→**6.51 kN** 假台阶 = 3.3+3e5×0.0191）。
**★★ 新 gotcha**：`h*` 内插**不能先滤污染格再内插** —— 会跨过无效格（把 self 报成 **3.506 m 假首折**）；
只在原网格相邻且两端都有效的格间内插。
回归 **774 passed / 8 skipped**。
**⇒ B2 设计分叉待用户定**：**(B2a) 等效单关节**（θ 由 COM 深度反解 = 几何恒等式；膝 350 Nm/rad 仅作量级锚）
vs **(B2b) 膝+踝两关节**。反解出 stiff 也要 **86°**、soft **120°** ⇒ 实测 COM 深度**非纯膝** ⇒ 单膝几何装不下；
**我建议 B2a**。详见 `memory/2026-10-16/v1-stiffness-literature-calibration.md` §⑭、`PAD_GOTCHAS.md` §6d。
"""

with TOPIC.open("a", encoding="utf-8") as f:
    f.write(SEC)
with DAILY.open("a", encoding="utf-8") as f:
    f.write(PTR)
print("topic", TOPIC.stat().st_size, "| daily", DAILY.stat().st_size)
