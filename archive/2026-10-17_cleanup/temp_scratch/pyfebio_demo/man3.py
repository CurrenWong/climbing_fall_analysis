import re
from pathlib import Path

m = Path(r"D:\Project\climbing_fall_analysis\temp\thums/febio_manual.txt")
t = m.read_text(encoding="utf-8", errors="ignore")

out = []
out.append(f"{'='*30} 全部 12 处 'truss' 上下文 {'='*30}")
for k, mm in enumerate(re.finditer(r"truss", t, re.I)):
    i = mm.start()
    out.append(f"\n-- #{k} @{i} --\n" + t[max(0, i - 500):i + 700])

out.append(f"\n\n{'='*30} 搜 'nonlinear elastic' 材料段 {'='*30}")
for mm in list(re.finditer(r"[Nn]onlinear [Ee]lastic", t))[:3]:
    i = mm.start()
    out.append(f"\n-- @{i} --\n" + t[i:i + 2000])

Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\_man3.txt"
     ).write_text("\n".join(out), encoding="utf-8")
print("wrote", len("\n".join(out)))
