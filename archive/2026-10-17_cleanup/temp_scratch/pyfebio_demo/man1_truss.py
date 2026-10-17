"""查 FEBio 手册：truss 单元写法、discrete domain、tension-only 材料。"""
import re
from pathlib import Path

m = Path(r"D:\Project\climbing_fall_analysis\temp\thums/febio_manual.txt")
t = m.read_text(encoding="utf-8", errors="ignore")
print(f"手册 {len(t):,} 字符\n")

for kw in ("truss3", "truss", "DiscreteDomain", "discrete domain",
           "tension-only", "tension only", "rigid body", "shell4",
           "ShellDomain"):
    idx = [mm.start() for mm in re.finditer(re.escape(kw), t, re.I)]
    print(f"{kw!r}: {len(idx)} 处")
    for i in idx[:4]:
        seg = t[max(0, i - 180):i + 320].replace("\n", " ")
        print(f"    ...{seg}...")
    print()
