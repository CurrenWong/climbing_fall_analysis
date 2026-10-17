import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
t = (P / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
head = t[t.index("# 🧭 接手入口"):t.index("## C. 解出 THUMS")]
for name in ("bone_materials.json", "build_feb.py", "cartilage.py",
             "thums_lowerlimb.npz"):
    print(f"--- {name} 在入口节的引用 ---")
    for line in head.split("\n"):
        if name in line:
            for h in re.findall(r"`([^`]+)`", line):
                if name in h:
                    p = P / h.replace("\\", "/")
                    print(f"    {'OK  ' if p.exists() else 'MISS'} {h}")
