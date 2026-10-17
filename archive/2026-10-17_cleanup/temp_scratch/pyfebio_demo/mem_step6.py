"""MEMORY 加"唯一运行入口"指针（复用现有 bullet，省字节）+ 终检。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
s = M.read_text(encoding="utf-8")

OLD = """- ⛔ **「不用软骨」≠解法**（`[LRN-062]`）：力变**罚因子伪影**（`CART_E=5000`+`FIX_TALUS=all`→**59,590 N**；
  罚因子 1e-1/2/3→**36,156/8,813/1,248 N**）。**通用判据：力随罚因子显著变化 ⇒ 不可信**"""
NEW = """- ⛔ **「不用软骨」≠解法**（`[LRN-062]`）：力变**罚因子伪影**（`CART_E=5000`+`FIX_TALUS=all`→**59,590 N**；
  罚因子 1e-1/2/3→**36,156/8,813/1,248 N**）。**通用判据：力随罚因子显著变化 ⇒ 不可信**
- 📘 **唯一运行入口**（`[LRN-063]`）= `docs/踝FE可用模型使用说明.md`；**姿势=`POSE_TILT`/`POSE_TILT_ML`**（倾斜式，
  ±10° 剪切 0→55~118 N 全 NORMAL；几何旋转式 `POSE_DF` 因关节面相切而受阻）"""
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    print("✓ MEMORY 已加运行入口")
else:
    print("⚠ 未匹配 MEMORY 锚点")
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md -> {sz:,}/10240  余量 {10240-sz:+,}  {'OK' if sz <= 10240 else 'XX 超'}")

print("\n" + "=" * 74)
LT = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
DOC = ROOT / "docs/踝FE可用模型使用说明.md"
doc = DOC.read_text(encoding="utf-8") if DOC.exists() else ""
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
DAILY = WS / "memory/2026-10-04.md"
ids = re.findall(r"^##\s*\[LRN-(\d{8}-\d{3})\]", LT, re.M)
dups = sorted({x for x in ids if ids.count(x) > 1})
CHK = [
    ("LRN-057~063 全在档", all(f"LRN-20261005-{k:03d}" in LT for k in (57, 58, 59, 60, 61, 62, 63))),
    ("LRN-063 含倾斜式姿势", "POSE_TILT" in LT and "倾斜式" in LT),
    ("使用说明存在", DOC.exists()),
    ("使用说明含一条命令", "acceptance.py" in doc and "ALIGN_AXIS=1" in doc),
    ("使用说明含已知限制", "载荷上限 ~235 N" in doc and "罚因子伪影" in doc),
    ("使用说明含避坑表", "别踩的坑" in doc),
    ("session §二十一", "二十一(10-05)" in SESS.read_text(encoding="utf-8")),
    ("daily 含 LRN-063", "LRN-20261005-063" in DAILY.read_text(encoding="utf-8")),
    ("MEMORY 含入口指针", "踝FE可用模型使用说明.md" in ms),
    ("MEMORY 含姿势开关", "POSE_TILT" in ms),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")
print(f"\n  LEARNINGS {LT and (ROOT/'.learnings/LEARNINGS.md').stat().st_size:,} B / {len(ids)} 条 / 重复 {dups}")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json))`", ms)))
dead = [p for p in refs if "/" in p and not ((WS / p).exists() or (ROOT / p).exists())]
print(f"  MEMORY 引用 {len(refs)} 带路径死链 {len(dead)} {dead}")
for p in [ROOT / "docs/踝FE可用模型使用说明.md", ROOT / "ANKLE_GOAL_TODO.md",
          ROOT / "scripts/ankle_fe/acceptance.py", ROOT / "scripts/ankle_fe/ankle_lig_feb.py",
          SESS, DAILY]:
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {p.name}")
