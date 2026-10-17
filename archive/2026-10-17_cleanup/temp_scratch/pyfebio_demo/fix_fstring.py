"""给含 {VAR} 占位符的 w('...') 行补 f 前缀。"""
import ast
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
L = p.read_text(encoding="utf-8").split("\n")
n = 0
for i, l in enumerate(L):
    if ("w('" in l and "{" in l and "}" in l and "w(f'" not in l
            and any(v in l for v in ("ETOL", "DTOL", "MAX_REFS", "MAX_UPS"))):
        L[i] = l.replace("w('", "w(f'", 1)
        n += 1
        print(f"  {i+1}: {L[i].strip()[:78]}")
p.write_text("\n".join(L), encoding="utf-8")
ast.parse("\n".join(L))
print(f"\n修了 {n} 行，syntax OK")
