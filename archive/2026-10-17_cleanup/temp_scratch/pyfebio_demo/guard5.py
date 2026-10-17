"""补：空 Surface 不写 + 不写零 traction 的载荷（按行定位）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")

# ① fibula_top / tibia_top Surface 块加守卫
for nm, var in (("fibula_top", "top_fibula"), ("tibia_top", "top_tibia")):
    i0 = next(k for k, l in enumerate(L) if l.strip() == f"""w('    <Surface name="{nm}">')""")
    i1 = next(k for k in range(i0, len(L)) if L[k].strip() == "w('    </Surface>')")
    for k in range(i0, i1 + 1):
        if L[k].strip():
            L[k] = "    " + L[k]
    L.insert(i0, f"if {var}:")

# ② 载荷：跳过空面
i = next(k for k, l in enumerate(L) if "for surf, F in ((" in l)
i_end = i
while i_end + 1 < len(L) and (L[i_end + 1].startswith("    ")
                              or not L[i_end + 1].strip()):
    i_end += 1
body = ["    " + x if x.strip() else x for x in L[i:i_end + 1]]
L[i:i_end + 1] = ['for surf, F, _faces in (("tibia_top", F_t, top_tibia),',
                  '                      ("fibula_top", F_f, top_fibula)):',
                  '    if not _faces or F <= 0:',
                  '        continue'] + body[1:]

p.write_text("\n".join(L), encoding="utf-8")
s = p.read_text(encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
