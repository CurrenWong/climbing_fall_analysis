"""补全路径前缀 + 最终验证（输出纯 ASCII，避免 cmd 编码吞行）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
F = ROOT / "FE_PIPELINE_TODO.md"
t = F.read_text(encoding="utf-8")

B = chr(96)          # `
ARROW = chr(8594)    # →
n = 0
for name in ("bone_materials.json", "thums_lowerlimb.npz"):
    old = f"{ARROW} {B}{name}{B}"
    new = f"{ARROW} {B}temp/thums/{name}{B}"
    n += t.count(old)
    t = t.replace(old, new)
F.write_text(t, encoding="utf-8")
print(f"prefixed {n} bare output names")

head = t[t.index("# \U0001F9ED"):t.index("## C. ")]
paths = sorted(set(re.findall(r"`([\w\-./\\]+\.(?:py|md|json|npz|zip|k|fe?b))`", head)))
miss = [p for p in paths
        if not (ROOT / p.replace("\\", "/")).exists() and not (WS / p).exists()]
print(f"entry paths: {len(paths)}, unresolvable: {len(miss)} {miss}")

doc = (ROOT / "docs/踝关节FE建模交接.md").read_text(encoding="utf-8")
ank = (ROOT / "ANKLE_FE_TODO.md").read_text(encoding="utf-8")
mem = (WS / "MEMORY.md").read_text(encoding="utf-8")
lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")

checks = [
    ("HANDOVER entry section exists", "# \U0001F9ED" in t),
    ("entry commands marked verified", "has been run-verified" in t
     or "本节命令已实跑验证" in t),
    ("stale status cleared", "A ⚠️ 进行中" not in t and "卡在 Jacobian" not in t),
    ("no bare file cells in switch table",
     "| `build_feb.py` |" not in t and "| `cartilage.py` |" not in t),
    ("ANKLE_FE marked superseded", "IITD" in ank and "已停止" in ank),
    ("handover doc cross-refs entry", "接手入口" in doc),
    ("MEMORY.md <= 10KB", len(mem.encode("utf-8")) <= 10240),
    ("LRN-016..020 present",
     all(f"LRN-20261003-0{k}" in lrn for k in (16, 17, 18, 19, 20))),
    ("all entry paths resolve", len(miss) == 0),
]
print()
for k, v in checks:
    print(f"  [{'OK' if v else 'XX'}] {k}")
print(f"\nFE_PIPELINE_TODO.md = {len(t.encode('utf-8')):,} bytes")
print(f"ALL PASS = {all(v for _, v in checks)}")
