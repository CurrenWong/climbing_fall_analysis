"""找 </MeshDomains> 与 Material 段的位置，确定 <Discrete> 该插在哪。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if any(k in l for k in ("</MeshDomains>", "<MeshDomains>", "<Material>",
                            "<Rigid>", "<Constraints>", "<LoadData>",
                            "<Boundary>", "TAL_GROUND")):
        print(f"{i:>5}|{l[:105]}")
