import re
from pathlib import Path

m = Path(r"D:\Project\climbing_fall_analysis\temp\thums/febio_manual.txt")
t = m.read_text(encoding="utf-8", errors="ignore")

out = []
for kw in ("BeamDomain Section", "truss and beam", "type=\"truss\"",
           "linear elastic\""):
    for mm in re.finditer(re.escape(kw), t):
        i = mm.start()
        out.append(f"\n{'='*30} {kw!r} @ {i} {'='*30}")
        out.append(t[max(0, i - 300):i + 1800])
        break

# BeamDomain 附近
i = t.find("3.7.3 BeamDomain")
if i < 0:
    i = t.find("BeamDomain Section", 5000)
out.append(f"\n{'='*30} BeamDomain 区段 {'='*30}")
out.append(t[i:i + 2500])

# 找 truss 材料：搜 "area" 与 truss 同现
out.append(f"\n{'='*30} 含 area 的材料参数 {'='*30}")
for mm in list(re.finditer(r"area", t, re.I))[:400]:
    i = mm.start()
    seg = t[max(0, i - 400):i + 400]
    if "truss" in seg.lower() or "beam" in seg.lower():
        out.append(f"\n-- @{i} --\n{seg}")
        break

Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_man2.txt"
     ).write_text("\n".join(out), encoding="utf-8")
print("wrote", len("\n".join(out)))
