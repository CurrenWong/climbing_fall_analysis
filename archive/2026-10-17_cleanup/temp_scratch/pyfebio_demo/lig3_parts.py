"""在 219MB 主 .k 里列出所有 *PART 标题，找足部骨（舟/骰/楔）与全部韧带。"""
import re
import zipfile
from pathlib import Path

ZIP = Path(r"D:\Project\climbing_fall_analysis\paper\AM50_V71_Occupant.zip")
INNER = "AM50_V71_Occupant/THUMS_model/THUMS_AM50_V71_Occupant_202411.k"

pat = re.compile(r"^\*PART\s*$", re.I)
parts = []
pend = None
lastkw = None

with zipfile.ZipFile(ZIP) as z, z.open(INNER) as f:
    for raw in f:
        line = raw.decode("latin-1").rstrip("\r\n")
        if line.startswith("*"):
            kw = line.split()[0].upper() if line.split() else ""
            lastkw = kw
            pend = kw == "*PART"
            continue
        if pend:
            if line.strip():
                parts.append(line.strip())
            pend = False

print(f"共 {len(parts)} 个 *PART 标题\n")
print("=== 足部/踝部相关（舟 navicular / 骰 cuboid / 楔 cuneiform）===")
kw = ("NAVIC", "CUBOID", "CUNEIF", "METATARS", "PHALAN", "TALUS", "CALCAN",
      "LIG", "TIBIA", "FIBULA", "PLANTAR", "FOOT")
for i, p in enumerate(parts):
    if any(k in p.upper() for k in kw):
        print(f"  {i:3d}  {p}")
print("\n=== 全部 PART（前 80）===")
for i, p in enumerate(parts[:80]):
    print(f"  {i:3d}  {p}")
