"""定位 </Boundary> 与 _RIGID_LINES。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if "</Boundary>" in l or "_RIGID_LINES" in l or "<Output>" in l.replace("w('", ""):
        print(f"{i:>5}|{l[:100]}")
