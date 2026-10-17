"""盘点记忆文件现状（为压缩上下文前收口做准备）。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
print("=== memory/ 下的文件 ===")
for p in sorted(WS.glob("memory/**/*.md")):
    if "2026-10" in str(p):
        print(f"  {p.stat().st_size:>9,}  {p.relative_to(WS)}")

p = WS / "memory/2026-10-03/ankle-fe-modeling.md"
if p.exists():
    s = p.read_text(encoding="utf-8")
    print(f"\nankle-fe-modeling.md: {len(s.encode('utf-8')):,} bytes")
    print("章节:")
    for l in s.split("\n"):
        if l.startswith("## "):
            print("   ", l.strip()[:60])

m = WS / "MEMORY.md"
if m.exists():
    print(f"\nMEMORY.md: {len(m.read_text(encoding='utf-8').encode('utf-8')):,} bytes "
          f"(上限 10,240)")
