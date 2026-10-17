"""健壮性守卫（按行定位，避开文本匹配）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")


def find(pred, start=0):
    return next(i for i in range(start, len(L)) if pred(L[i]))


def indent(i0, i1):
    for k in range(i0, i1 + 1):
        if L[k].strip():
            L[k] = "    " + L[k]


# ① calcaneus BC 加守卫
i0 = find(lambda l: "fix_calcaneus" in l)
i1 = find(lambda l: l.strip() == "w('    </bc>')", i0)
L.insert(i0, "if bot_calc:")
indent(i0 + 1, i1 + 1)

# ② NodeSet 守卫（三处）
for nm in ("tibia_top_ns", "fibula_top_ns", "calcaneus_bottom_ns"):
    i = find(lambda l, nm=nm: nm in l and "NodeSet" in l)
    var = {"tibia_top_ns": "top_tibia", "fibula_top_ns": "top_fibula",
           "calcaneus_bottom_ns": "bot_calc"}[nm]
    # 该语句可能跨 2 行
    j = i
    while L[j].count("'") % 2 == 1 or L[j].count('"') % 2 == 1:
        j += 1
    L.insert(i, f"if {var}:")
    indent(i + 1, j + 1)

p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
