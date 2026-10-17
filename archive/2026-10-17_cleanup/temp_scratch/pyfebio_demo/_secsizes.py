"""按段统计 MEMORY.md 体积，找出还值得压的地方。"""
from pathlib import Path

P = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
lines = P.read_text(encoding="utf-8").split("\n")
cur, acc, size = "(frontmatter)", [], 0
rows = []
for l in lines:
    if l.startswith("## ") or l.startswith("### "):
        rows.append((size, cur))
        cur, size = l.strip("# ").strip(), 0
    size += len(l.encode("utf-8")) + 1
rows.append((size, cur))
rows.sort(reverse=True)
tot = sum(r[0] for r in rows)
print(f"总计 {tot:,} bytes（上限 10240，需再砍 {tot-10240:,}）\n")
for s, n in rows[:12]:
    print(f"  {s:>6,}  {n[:60]}")
print("\n--- 最长的单行 top8 ---")
for s, n in sorted(((len(l.encode('utf-8')) + 1, l) for l in lines),
                   reverse=True)[:8]:
    print(f"  {s:>5,}  {n[:88]}")
