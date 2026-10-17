"""核对 <Nodes> 写出位置 vs 地面节点追加位置。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if any(k in l for k in ("<Mesh>", "Nodes name=", "coords.append",
                            "GROUND_IDS.append", "tal_faces = extreme",
                            "tal_ground_ns", "TAL_GROUND:", "if TALUS_RIGID:",
                            "talus_rigid_ns")):
        print(f"{i:>5}|{l[:100]}")
