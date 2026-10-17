"""清理：删 3 文件里的旧错横幅 + 修使用说明 8 处（含 §1"必开"）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"

# ── 1) 三个文件：删掉残留的旧横幅 ──
for rel in ("ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md", "TIE_PORT_TODO.md"):
    p = os.path.join(R, rel)
    t = io.open(p, encoding="utf-8").read()
    m = re.search(r"^> 🚨 \*\*2026-10-05 重大更正.*?^---\n", t, re.S | re.M)
    if m:
        t = t[:m.start()] + t[m.end():]
        io.open(p, "w", encoding="utf-8").write(t)
        print(f"✓ {rel} 旧横幅已删除 -> {len(t.encode('utf-8'))} B")
    else:
        print(f"⚠ {rel} 未找到旧横幅")

# ── 2) 使用说明 ──
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
d = io.open(D, encoding="utf-8").read()
n0 = len(d.encode("utf-8"))
REPS = [
    # §1 环境变量表：LOAD_SIGN 由"必开"改为"无效/别设"
    ("| **`LOAD_SIGN=+1`** | **必开** | **压实关节**。默认 `-1` 是**拉伸**（把关节拉开）⇒ "
     "反力由韧带承担、量级只有 ~10~200 N。见 `[ERR-20261005-009]` |",
     "| ~~`LOAD_SIGN`~~ | **对 `dispx` 无效，不用设** | 实测写入的 z 只由 `DISP_MM` 的符号决定。"
     "**压实 = `DISP_MM > 0`**（见 §2 符号说明与 `[ERR-20261006-003]`）|"),
    # §5 符号表
    ("| `LOAD_SIGN=+1` / `DISP_MM<0` | **压实关节** ✅ | `−0.25 → −132.56 N`（压缩）；`−1.0 → −608.82 N` |",
     "| **`DISP_MM > 0`** | **压实关节** ✅ | `+1.0 → 418.10 N`；**`+1.6 → 592.60 N`** |"),
    ("**⇒ 默认值是\"拉\"，不设 `LOAD_SIGN=+1` 会静默跑出虚低的力**（构建日志会打 ⚠️ 告警）。\n"
     "**⇒ 历史上\"验证过的\" `187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N` 全是拉伸值。**",
     "**⇒ 取 `DISP_MM < 0` 会静默跑成\"拉开\"**（反力由韧带承担、量级虚低且随韧带几何变化）。\n"
     "**⇒ 而 `187.77 / 187.90 / 202.29 / 209.5589 / 592.61 N`（`DISP_MM>0`）是**压实**值、有效。**\n"
     "**⇒ `DISP_MM=−0.25…−1.6 → −132…−1017 N` 才是**拉伸**值，已作废。**"),
    # §6 坑表
    ("| **不设 `LOAD_SIGN=+1`（或 `DISP_MM` 取正）** | **静默跑成\"拉伸\"** ⇒ 力量级虚低（~10~200 N）、"
     "看着像\"收敛不了\"。构建日志有 ⚠️ 告警 |",
     "| **`DISP_MM` 取负** | **静默跑成\"拉开（拉伸）\"** ⇒ 反力由韧带承担、量级虚低，"
     "且**会随韧带几何变化**（实测改韧带 ⇒ 力变 23%）⇒ 极易误判。**压实要 `DISP_MM > 0`** |"),
    # 旧结果段的标签
    ("**实测（`DISP_MM=-1.0`，压缩侧，**全部 `NORMAL`**）**：",
     "**实测（`DISP_MM=-1.0`＝**拉伸**方向，全部 `NORMAL`）** —— ⚠️ 旧口径，见下方 §4 已修正的三口径表："),
    ("**载荷随姿势变化**（`DISP_MM=-1.0`，变形距骨）：",
     "**载荷随姿势变化**（`DISP_MM=-1.0`＝**拉伸**方向，变形距骨）—— ⚠️ 旧口径；"
     "**压缩方向的重测见 §4 ③ 行**："),
]
n = 0
for old, new in REPS:
    if old in d:
        d = d.replace(old, new)
        n += 1
    else:
        print(f"  ⚠ 未命中：{old[:56]}")
io.open(D, "w", encoding="utf-8").write(d)
print(f"✓ 使用说明改写 {n} 处  {n0} -> {len(d.encode('utf-8'))} B")

# ── 3) 复查残留 ──
print("\n" + "=" * 86)
print("复查：仍声称「DISP<0 是压缩」或「历史值是拉伸」的位置")
print("=" * 86)
PAT = re.compile(r"(DISP_MM<0.{0,20}压缩|全是拉伸|全是\*\*拉伸|拉伸值。|LOAD_SIGN=\+1\*\* \| \*\*必开)")
for rel in ("docs/踝FE可用模型使用说明.md", "ANKLE_FIX_TODO.md",
            "docs/总体进度.md", "FE_PIPELINE_TODO.md", "ANKLE_GOAL_TODO.md",
            "ANKLE_600N_TODO.md", "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md"):
    p = os.path.join(R, rel)
    if not os.path.exists(p):
        continue
    t = io.open(p, encoding="utf-8").read()
    hits = [(i + 1, l.strip()) for i, l in enumerate(t.splitlines()) if PAT.search(l)]
    if hits:
        print(f"\n--- {rel}（{len(hits)}）")
        for i, l in hits[:8]:
            print(f"  {i:>4}: {l[:140]}")
