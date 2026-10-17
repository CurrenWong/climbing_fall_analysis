"""显示 MEMORY.md 踝关节段 + 计算可裁剪处。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
L = M.read_text(encoding="utf-8").split("\n")
i0 = next(i for i, l in enumerate(L) if l.startswith("- **接手入口**"))
i1 = next(i for i in range(i0, len(L)) if L[i].startswith("- **详见**"))
print("=" * 80)
for i in range(i0 - 1, i1 + 1):
    print(f"{i:4d}| {L[i]}")
    print(f"      ({len(L[i].encode('utf-8'))} bytes)")
print("=" * 80)
print(f"段字节 {sum(len(x.encode('utf-8'))+2 for x in L[i0:i1+1]):,}")
print(f"全文磁盘 {M.stat().st_size:,} / 10240")
