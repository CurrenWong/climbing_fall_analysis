"""核对段序 + TAL_ROTLOCK 是否已定义。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
print("=== 关键段序 ===")
for i, l in enumerate(L, 1):
    if any(k in l for k in ("<MeshDomains>", "</MeshDomains>", "<Rigid>", "</Rigid>",
                            "<Discrete>", "</Mesh>", "<Boundary>", "</Boundary>",
                            "<LoadData>", "<Material>", "TAL_ROTLOCK",
                            "<DiscreteSet", "<NodeSet name=\"tal")):
        print(f"{i:>5}|{l[:100]}")
