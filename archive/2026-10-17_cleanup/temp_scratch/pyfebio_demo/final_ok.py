"""补回踝关节段的详见指针，并做最终全量校验。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240
s = M.read_text(encoding="utf-8")

if "ankle-fe-modeling.md" not in s:
    # 插在"韧带线已降级"那一行之后
    L = s.split("\n")
    i = next(k for k, l in enumerate(L) if "韧带线已降级" in l)
    L.insert(i + 1, "- **详见**：`memory/2026-10-03/ankle-fe-modeling.md`")
    s = "\n".join(L)
    M.write_text(s, encoding="utf-8")
sz = len(s.encode())
print(f"disk = {sz:,} / {LIMIT}   余量 {LIMIT - sz:+,}")
assert sz <= LIMIT, "超限"
print("详见指针:", "OK" if "ankle-fe-modeling.md" in s else "缺")

# ---- 最终校验 ----
print("\n" + "=" * 70)
DOCS = [
    ("docs/Anderson双论文_关键信息.md", ROOT / "docs/Anderson双论文_关键信息.md"),
    ("docs/吴恺论文_重解析关键信息.md", ROOT / "docs/吴恺论文_重解析关键信息.md"),
    ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
    ("CONV_TODO.md", ROOT / "CONV_TODO.md"),
    ("FE_PIPELINE_TODO.md", ROOT / "FE_PIPELINE_TODO.md"),
    (".learnings/LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md"),
    ("scripts/ankle_fe/*.py (3)", ROOT / "scripts/ankle_fe/ankle_lig_feb.py"),
    ("temp/thums/r1/cartilage.json", ROOT / "temp/thums/r1/cartilage.json"),
]
for nm, p in DOCS:
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {nm}")

paths = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb|k))`", s)))
dead = [p for p in paths
        if not ((WS / p).exists() or (ROOT / p).exists()
                or p.startswith(("memory/", "references/", "digest/", "docs/")))
        and not list(ROOT.glob(f"**/{Path(p).name}"))]
print(f"\n  MEMORY 引用 {len(paths)} / 死链 {len(dead)} {dead}")

lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
print("  本轮 LRN:", " ".join(f"{n}{'✓' if f'LRN-20261004-0{n}' in lrn else '✗'}"
                              for n in (39, 40, 41, 42, 43, 44, 45)))
print(f"\n  ALL = {'PASS' if not dead else 'CHECK'}")
