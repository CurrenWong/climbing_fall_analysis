"""MEMORY 再削 ~37 B。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")
REPS = [
    ("""| `references/accounts.md` | **Neon 项目 + 账户/分类 ID 速查 + 记账规则** |
| `references/neon-mcp.md` | Neon MCP stdio 接管方案 |""",
     """| `references/accounts.md` / `neon-mcp.md` | **Neon 项目 + 账户/分类 ID + 记账规则** / stdio 接管 |"""),
    ("- **中文优先**（专有名词留英文）；**先行动后报告**；**宁缺毋滥**（搜索硬 3-5 条）；**保留原文**",
     "- **中文优先**（专有名词留英文）；**先行动后报告**；**宁缺毋滥**（硬 3-5 条）；**保留原文**"),
    ("：k=F/δ 趋近物理值 **49/151/451/880** N/mm", "：k=F/δ 趋近 **49/151/451/880** N/mm"),
]
for a, b in REPS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配: {a[:36]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

import re
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ROOT = Path(r"D:\Project\climbing_fall_analysis")
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
print(f"字数 {len(ms):,} chars / {ms.count(chr(10))+1} 行")
