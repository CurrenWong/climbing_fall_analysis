"""MEMORY 削 ~110 B + 修 TODO 更新。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
M = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
s = M.read_text(encoding="utf-8")

REPS = [
    ("""- ⛔ **600 N 到不了**（`[LRN-061]`）：稳定窗仅轴向 **0.1~0.3mm**（≤235 N），硬天花板 **~340 N**
  （`DET_MIN=0.95` 抬到 0.6mm；叠加厚软骨 ~380 N）；`CART_E=1000` 可达 561 N ⇒ **限流=软骨单元质量**""",
     """- ⛔ **600 N 到不了**（`[LRN-061]`）：稳定窗仅轴向 **0.1~0.3mm**（≤235 N），天花板 **~340 N**
  （`DET_MIN=0.95`→0.6mm）；`CART_E=1000` 达 561 N ⇒ **限流=软骨单元质量**"""),
    ("② **接触面要与网格真实面重合且是合法 2-流形**（否则接触不检测，却报 `rc=0` 无警告、**只是应力全 0**）",
     "② **接触面要与网格真实面重合且是合法 2-流形**（否则静默不检测，只报 `rc=0`、应力全 0）"),
    ("- **emoji 用在合适处**（Discord/Slack 表情回应）；正文中不用",
     "- **emoji 用在合适处**（Discord/Slack 回应）；正文中不用"),
    ("""### 其他
- MCP 死代码 → 自接管：`references/mcp-bypass.md`（Tavily / Flint 同款）；定时任务 → `copaw cron`；
  Server 酱 → `references/operational-gotchas.md`""",
     """### 其他
- MCP 死代码 → 自接管：`references/mcp-bypass.md`（Tavily / Flint 同款）
- 定时任务 → `copaw cron`；Server 酱 → `references/operational-gotchas.md`"""),
]
for a, b in REPS:
    if a in s:
        s = s.replace(a, b, 1)
    else:
        print(f"⚠ 未匹配: {a[:34]}...")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

# TODO（上次未匹配，改用更短锚点）
TD = ROOT / "ANKLE_GOAL_TODO.md"
t = TD.read_text(encoding="utf-8")
ANCH = "- [ ] Step 3 `ankle_lig_feb.py` 加 `TAL_SPRING`"
ADD = ("- [x] **Step 2-1** 扫出稳定窗：**轴向 0.1~0.3mm ⇒ ≤235 N**；天花板 **~340 N**（`DET_MIN=0.95`），\n"
       "  叠加厚软骨 ~380 N ⇒ **600 N 到不了**。`[LRN-20261005-061]`\n"
       "- [x] **Step 2-2** 否证两假设：**间隙说**（`T_TARGET=3.0` 单独更差、穿透过 0）与"
       "**距骨自由说**（`FIX_TALUS=all` 无效）；`NODE_RELOC=1` 更差；`LIG_MODE=none` 发散\n"
       "- [x] **Step 2-3** 最强杠杆 = **`DET_MIN`**（0.3→0.6mm）⇒ **限流 = 近退化单元**；"
       "`CART_E=1000` 达 561 N 佐证\n"
       "- [ ] **Step 3a ⭐ 要一个质量好的软骨网格**（唯一剩下的强候选）：`cart_patch.py` 现成但输出\n"
       "  prism6+`tie_faces` 新格式 ⇒ **需给 `ankle_lig_feb.py` 移植 TIE**（`build_feb.py` 可参照）；\n"
       "  或在共享节点格式下重做 2D 参数化网格\n")
if ANCH in t:
    t = t.replace(ANCH, ADD + ANCH, 1)
    print("✓ TODO 已插入 Step 2-1~3a")
else:
    # 兜底：追加到「开工顺序」节末尾
    t = t.rstrip() + "\n\n## Step 2 (实测补充，2026-10-05)\n" + ADD
    print("✓ TODO 已追加（兜底）")
TD.write_text(t, encoding="utf-8")
print(f"ANKLE_GOAL_TODO.md -> {TD.stat().st_size:,} B")

import re
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
ms = M.read_text(encoding="utf-8")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
for k in ("LRN-061", "600 N 到不了", "ALIGN_AXIS", "DET_MIN"):
    print(f"  [{'OK' if k in ms else 'XX'}] MEMORY 含 {k}")
