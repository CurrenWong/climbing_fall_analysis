from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
s = (WS / "MEMORY.md").read_text(encoding="utf-8")
lines = s.split("\n")
i0 = next(i for i, l in enumerate(lines) if "踝关节 FE 子线" in l)
i1 = next(i for i in range(i0 + 1, len(lines))
          if lines[i].startswith("### ") or lines[i].startswith("## "))
print("=" * 74)
for i in range(i0 - 2, i1 + 1):
    print(f"{i:4d}| {lines[i]}")
print("=" * 74)
print(f"该段字节数: {len(chr(10).join(lines[i0-2:i1+1]).encode('utf-8')):,}")
print(f"MEMORY.md 总计: {len(s.encode('utf-8')):,} / 10240")
