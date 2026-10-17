"""定位 NodeSet / </Mesh> / BC 的关键行号。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
KEYS = ("talus_low_ns", "navicular_all_ns", "</Mesh>", "tibiotalar_pair",
        "<Boundary>", "fix_talus", "tal_dofs =", "lig_only_ns", "<Discrete>")
for i, l in enumerate(L, 1):
    if any(k in l for k in KEYS):
        print(f"{i:>5}|{l[:105]}")
