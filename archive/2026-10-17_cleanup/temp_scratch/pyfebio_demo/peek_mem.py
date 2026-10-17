from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L2 = M.read_text(encoding="utf-8").split("\n")
print(f"disk={M.stat().st_size:,} 行数={len(L2)}")
print("--- 最后 6 行 ---")
for ln in L2[-6:]:
    print(f"  [{len(ln.encode())}] {ln[:150]}")
print("--- 最长的 8 行 ---")
for ln in sorted(L2, key=lambda x: -len(x.encode()))[:8]:
    print(f"  [{len(ln.encode()):4d}] {ln[:120]}")
