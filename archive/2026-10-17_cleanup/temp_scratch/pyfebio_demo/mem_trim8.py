from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
L = 10240
s = M.read_text(encoding="utf-8")
CUTS = [
    ("| 21 | **FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法** | `references/fe-modeling-gotchas.md`。**一句话**：",
     "| 21 | **FE 接触四铁律 + tet 翻转奇偶性 + `.feb` 语法** | `references/fe-modeling-gotchas.md`。**一句**："),
    ("④ **优化前先确认它是不是瓶颈**（我曾把\"网格质量\"当瓶颈，结果质量提升 42 倍反而更差）|",
     "④ **优化前先确认它是不是瓶颈**（曾把\"网格质量\"当瓶颈，提升 42 倍反而更差）|"),
    ("| 3 | **copaw MCP 是死代码** —— 直接用 skill/npx，别等内置工具 | `references/mcp-bypass.md` + `[LRN-20260604-001]` |",
     "| 3 | **copaw MCP 是死代码** —— 直接用 skill/npx，别等内置工具 | `references/mcp-bypass.md` |"),
]
for a, b in CUTS:
    if M.stat().st_size <= L:
        break
    if a in s and a != b:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        print(f"  -> {M.stat().st_size:,} 余量 {L - M.stat().st_size:+,}")
print(f"disk = {M.stat().st_size:,} / {L}  余量 {L - M.stat().st_size:+,}")
