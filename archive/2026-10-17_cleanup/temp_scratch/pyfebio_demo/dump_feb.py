"""导出自前 .feb 的 Boundary / Contact / Rigid / Discrete 段，看真实约束。"""
import io
import re
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
p = W / "ankle_lig_WK.feb"
t = p.read_text(encoding="latin-1", errors="replace")
print(f"=== {p.name}  {p.stat().st_size/1024:.0f} KB  行数 {t.count(chr(10))+1}")

for sec in ("MeshDomains", "Boundary", "Contact", "Constraint", "Discrete",
            "Rigid", "Loads", "Material"):
    for m in re.finditer(rf"<{sec}>(.*?)</{sec}>", t, re.S):
        body = m.group(1)
        lines = [l for l in body.split("\n") if l.strip()]
        print(f"\n----- <{sec}>  {len(lines)} 行")
        for l in lines[:70]:
            print("   " + l.strip()[:150])
        if len(lines) > 70:
            print(f"   ... 余 {len(lines)-70} 行")
