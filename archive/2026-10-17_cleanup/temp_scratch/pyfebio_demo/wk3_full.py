import io
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\paper/_text/吴恺_足踝FE.txt")
t = io.open(p, encoding="utf-8", errors="ignore").read()
L = t.split("\n")
compact = [l.strip() for l in L if l.strip()]

print("=== 结果段（含韧带位移）完整拼接 ===")
print(" ".join(compact[160:178]))

print("\n\n=== 边界条件段 ===")
print(" ".join(compact[150:165]))

print("\n\n=== 材料表 ===")
print(" ".join(compact[112:132]))
