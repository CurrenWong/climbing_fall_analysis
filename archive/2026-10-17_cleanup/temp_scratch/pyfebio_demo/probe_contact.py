"""定位 Contact 段 / all_surfs / SurfacePair 的写法，便于插 tie。"""
import io
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = io.open(P, encoding="utf-8").read().split("\n")
for i, l in enumerate(L, 1):
    if re.search(r"all_surfs|SurfacePair|surface_pair|<Contact>|</Contact>|CONTACT_PAIRS|PAIRS|tibiotalar_pair|blocks = |def ns_of", l):
        print(f"{i:>5}|{l[:118]}")
