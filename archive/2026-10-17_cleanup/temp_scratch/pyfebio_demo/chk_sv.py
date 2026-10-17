import re
from pathlib import Path

s = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
         ).read_text(encoding="utf-8").split("\n")
print(f"{len(s)} 行")
for i, l in enumerate(s, 1):
    if re.search(r"def detJ_|def sv_|^FN =|sv_penta6|\"penta6\":" , l):
        print(f"{i:5d}| {l[:110]}")
