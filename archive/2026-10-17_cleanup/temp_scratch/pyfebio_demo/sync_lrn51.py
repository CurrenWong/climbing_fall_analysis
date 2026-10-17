"""把 LRN-051 的更正同步到 ANKLE_FIX_TODO / FE_PIPELINE_TODO / CONV_TODO。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")

# ---- ① ANKLE_FIX_TODO.md ----
P = ROOT / "ANKLE_FIX_TODO.md"
s = P.read_text(encoding="utf-8")
OLD = ("      `cart_refine.py` 的加密本身正确（wedge6 加密 0 翻负、L=0 体积偏差 +0.0000%），\n"
       "      但 **hex8 加密恒有 ~30% 子单元翻负**")
NEW = ("      ⚠️ **已更正（LRN-051）**：原写「wedge6 加密 0 翻负」是误判 —— 三角形细分\n"
       "      我写成 1→3（**漏中央三角形**），**静默丢 25% 体积**；L=1 的 −14~−18% 偏差\n"
       "      正是它（补上中央后 → **−0.77~−2.05%**，单元数 93→**372** ✓）。\n"
       "      但 **hex8 加密恒有 ~28% 子单元翻负** 这条仍然成立")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    P.write_text(s, encoding="utf-8")
print(f"ANKLE_FIX_TODO.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-051' in P.read_text(encoding='utf-8') else 'XX'}]")

# ---- ② FE_PIPELINE_TODO.md：入口「新通路」行补一句 ----
P = ROOT / "FE_PIPELINE_TODO.md"
s = P.read_text(encoding="utf-8")
OLD = "⚠️ 但**软骨面内红加密本身仍被母单元扭曲挡住**（hex8 加密恒 30% 子单元翻负，`[LRN-050]`）\n⇒ **机制通、收益尚未兑现**"
NEW = ("⚠️ 但**软骨面内红加密的 hex8 侧仍被母单元扭曲挡住**（**28%** 子单元翻负，`[LRN-050]`；\n"
       "注意 `[LRN-051]` 更正了旁支：三角形细分漏中央三角形曾静默丢 25% 体积，已修）\n"
       "⇒ **机制通、收益尚未兑现**")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    P.write_text(s, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-051' in P.read_text(encoding='utf-8') else 'XX'}]")

# ---- ③ CONV_TODO.md ----
P = ROOT / "CONV_TODO.md"
s = P.read_text(encoding="utf-8")
OLD = "但**软骨面内加密本身**目前被母单元扭曲挡住（P1d2b），故尚未能兑现。"
NEW = ("但**软骨面内加密的 hex8 侧**目前被母单元扭曲挡住（P1d2b），故尚未能兑现。\n"
       "（`[LRN-051]`：三角形细分曾因漏中央三角形静默丢 25% 体积，已修；\n"
       "hex8 侧 28% 子单元翻负的原因不变。）")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    P.write_text(s, encoding="utf-8")
print(f"CONV_TODO.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if 'LRN-051' in P.read_text(encoding='utf-8') else 'XX'}]")
