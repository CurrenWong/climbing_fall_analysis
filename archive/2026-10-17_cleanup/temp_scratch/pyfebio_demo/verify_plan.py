"""给 MEMORY 留余量 + 全量校验。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

s = M.read_text(encoding="utf-8")
TRIMS = [
    ("- **🗺️总览**：`docs/总体进度.md`（两条线+子线）", "- **🗺️总览**：`docs/总体进度.md`"),
    ("- ⚠️ **刚度 2 倍未归因**：tet4 259.3 vs penta6 116.9 N（恒 2.2~2.35，五假设已排除，`CONV_TODO.md`）",
     "- ⚠️ **刚度 2 倍未归因**：tet4 259.3 vs penta6 116.9 N（恒 2.2~2.35，五假设排除，`CONV_TODO.md`）"),
    ("- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜`docs/Anderson双论文_关键信息.md`",
     "- **入口**：`FE_PIPELINE_TODO.md`「🧭接手」｜`ANKLE_FIX_TODO.md`｜`docs/Anderson双论文_关键信息.md`"),
]
for a, b in TRIMS:
    if a in s and a != b:
        s = s.replace(a, b, 1)
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'✓' if sz<=10240 else '✗'}")

print("\n" + "=" * 76)
print("计划落盘校验")
print("=" * 76)
lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
AF = (ROOT / "ANKLE_FIX_TODO.md").read_text(encoding="utf-8")
OP = (ROOT / "docs/总体进度.md").read_text(encoding="utf-8")
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
DAILY = WS / "memory/2026-10-04.md"
sess = SESS.read_text(encoding="utf-8") if SESS.exists() else ""
daily = DAILY.read_text(encoding="utf-8") if DAILY.exists() else ""
ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", lt, re.M)
dups = sorted({x for x in ids if ids.count(x) > 1})

CHK = [
    ("LRN-057 入库", "LRN-20261005-057" in lt),
    ("LRN-057 含力控制卡口", "DRIVER=force`（600 N 面力）**会失败**" in lt),
    ("LRN-057 含循环论证纠正", "不是验证，是同义反复" in lt),
    ("LRN-057 含双端夹持", "双端夹持" in lt),
    ("LRN-057 含 Anderson §2.2 处方", "距骨平移弱弹簧** `k = 500 N/mm`" in lt),
    ("LRN-057 含 5 条验收", "姿势排序对罚因子鲁棒" in lt),
    ("LRN-057 含不做清单", "明确不做" in lt),
    ("ANKLE_FIX 有计划节", "下一步（已定，2026-10-05" in AF),
    ("ANKLE_FIX 含验收 5 条", "姿势排序对罚因子鲁棒" in AF),
    ("总体进度 有子线行", "LRN-20261005-057" in OP),
    ("session §十六", "十六、踝子线的下一个动作" in sess),
    ("daily 有索引", "踝子线的下一个动作" in daily),
    ("MEMORY 有指针", "LRN-057" in ms),
    ("MEMORY 已纠正韧带线", "循环论证" in ms),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")

print(f"\n  LEARNINGS {len(lt.encode()):,} B / {len(ids)} 条 / 重复 ID {len(dups)} {dups}")
for p in (ROOT / "ANKLE_FIX_TODO.md", ROOT / "docs/总体进度.md",
          ROOT / "ANKLE_LOCAL_TODO.md", ROOT / "TIE_DIVERGE_TODO.md",
          ROOT / "CART_HEX_TODO.md", ROOT / "FE_PIPELINE_TODO.md",
          SESS, DAILY, M):
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {p.name}")

# 死链（只查 MEMORY 引用）
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if not ((WS / p).exists() or (ROOT / p).exists())
        and "/" in p]
print(f"\n  MEMORY 引用 {len(refs)}，带路径的死链 {len(dead)} {dead}")
