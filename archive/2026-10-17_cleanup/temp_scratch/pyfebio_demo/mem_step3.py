"""MEMORY 再削 ~115 B。"""
from pathlib import Path

M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")

REPS = [
    ("| 5 | **cmd `set X=Y && Z` 尾追 space**（40013）| `references/wechat-publish.md` + "
     "`[ERR-20260629-001]` + `[ERR-20260717-001]`（wenyan.cmd set 链会丢 env → 必须多行 set 或 PowerShell）|",
     "| 5 | **cmd `set X=Y && Z` 尾追 space**（40013）| `references/wechat-publish.md`；"
     "wenyan.cmd **set 链会丢 env** ⇒ 多行 set 或 PowerShell |"),
    ("④ **优化前先确认它是不是瓶颈**（曾把\"网格质量\"当瓶颈，提升 42 倍反而更差）|",
     "④ **优化前先确认它是不是瓶颈**（曾误判\"网格质量\"，提升 42 倍反而更差）|"),
    ("""### 数据库 / 财务
- **记账/查账**：Neon MCP bypass → `references/accounts.md` + `references/neon-mcp.md`
- 65 个 skill 清单（按加载方式分 3 类）→ `references/tools-cheatsheet.md`""",
     """### 数据库 / 财务
- **记账/查账**：Neon MCP bypass → `references/accounts.md` + `neon-mcp.md`；skill 清单 → `references/tools-cheatsheet.md`"""),
]
for a, b in REPS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配: {a[:40]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

# 死链复检
import re
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ROOT = Path(r"D:\Project\climbing_fall_analysis")
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
for k in ("ALIGN_AXIS", "LRN-060", "3.2×", "12.18"):
    print(f"  [{'OK' if k in ms else 'XX'}] {k}")
