from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 子单元翻负，`[LRN-046~050]`）",
     "  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 30% 翻负，`[LRN-046~050]`）"),
    ("L=0 罚法**单调收敛**\n  到共形解 436.817 N",
     "L=0 罚法**单调收敛**到共形解\n  436.817 N"),
    ("  坑：tie 面须是**真实 facet**（tet4 只三角形）",
     "  坑：tie 面须是**真实 facet**（tet4 只三角）"),
    ("（`cart_refine.py`+`TIE=1`）：", "（`cart_refine.py`+`TIE=1`）:"),
]
for a, b in CUTS:
    s = s.replace(a, b, 1)
M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} 余量 {L - sz:+,}")
