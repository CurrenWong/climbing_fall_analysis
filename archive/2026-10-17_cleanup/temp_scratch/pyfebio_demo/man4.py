import re
from pathlib import Path

m = Path(r"D:\Project\climbing_fall_analysis\temp\thums/febio_manual.txt")
t = m.read_text(encoding="utf-8", errors="ignore")
out = []

# 找 beam/truss 的材料族
for kw in ("Beam Material", "beam material", "truss material",
           "BeamMaterial", "linear-truss", "elastic-truss",
           "axial stiffness", "cross-sectional area", "area"):
    idx = [mm.start() for mm in re.finditer(re.escape(kw), t, re.I)]
    out.append(f"\n{'='*25} {kw!r}: {len(idx)} 处 {'='*25}")
    for i in idx[:3]:
        out.append(f"\n-- @{i} --\n" + t[max(0, i - 350):i + 900])

Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_man4.txt"
     ).write_text("\n".join(out), encoding="utf-8")
print("wrote", len("\n".join(out)))
