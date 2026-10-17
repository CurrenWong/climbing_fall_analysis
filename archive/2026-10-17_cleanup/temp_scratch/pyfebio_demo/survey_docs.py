"""盘点：哪些文档已同步 tied-elastic 结果，哪些还没。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

DOCS = [
    ROOT / "FE_PIPELINE_TODO.md",
    ROOT / "CONV_TODO.md",
    ROOT / "ANKLE_FIX_TODO.md",
    ROOT / "ANKLE_600N_TODO.md",
    ROOT / "ANKLE_LIG_TODO.md",
    ROOT / "docs/Anderson双论文_关键信息.md",
    ROOT / "docs/踝关节FE建模交接.md",
    ROOT / "docs/吴恺论文_重解析关键信息.md",
    ROOT / ".learnings/LEARNINGS.md",
    ROOT / ".learnings/ERRORS.md",
    WS / "MEMORY.md",
]
MEM = ROOT / "memory"
for d in sorted(MEM.glob("2026-10-0*")):
    if d.is_dir():
        DOCS += sorted(d.glob("*.md"))

KEYS = {"tied-elastic": "tied-elastic", "TIE=1": "TIE=1",
        "432.4": "432.4", "penalty": "penalty", "cart_refine": "cart_refine",
        "LRN-046": "LRN-046", "2.2×": "2.2×"}

print(f"{'文件':52s} {'大小':>9s}  命中")
for p in DOCS:
    if not p.exists():
        print(f"{str(p.relative_to(p.parents[1])):52s} {'MISS':>9s}")
        continue
    t = p.read_text(encoding="utf-8", errors="replace")
    hit = [k for k, v in KEYS.items() if v in t]
    try:
        rel = p.relative_to(ROOT)
    except ValueError:
        rel = p.name
    print(f"{str(rel):52s} {p.stat().st_size:>9,}  {hit if hit else '--'}")
