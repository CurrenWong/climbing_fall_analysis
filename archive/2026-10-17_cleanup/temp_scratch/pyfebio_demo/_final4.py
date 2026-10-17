import io
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-03\ankle-fe-modeling.md")
t = io.open(M, encoding="utf-8").read()
print(f"记忆文件 {len(t.encode('utf-8')):,} bytes")
print("frontmatter:", t.startswith("---"), "| 有结尾 ---:",
      t.split("---", 2)[2][:1] if t.count("---") >= 2 else "N/A")
for k in ("## 九", "§九", "A✅"):
    print(f"  含 {k!r}: {k in t}")
print("章节:", [l.strip() for l in t.split("\n") if l.startswith("## ")])

R = Path(r"C:\Users\curre\.copaw\workspaces\default\references\fe-modeling-gotchas.md")
print(f"\nreferences/fe-modeling-gotchas.md {R.stat().st_size:,} bytes")
m = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md").read_text(encoding="utf-8")
print(f"MEMORY.md {len(m.encode('utf-8')):,} bytes（自设上限 10240，⚠️ 待整理）")
print("  gotcha #21 指向 references:", "fe-modeling-gotchas" in m)
