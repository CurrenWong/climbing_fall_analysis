"""打印 <Nodes> 写出代码。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
idx = [k for k, l in enumerate(L) if 'Nodes name=' in l]
print("Nodes 写出位置:", [k + 1 for k in idx])
for k in idx:
    for j in range(k, min(k + 14, len(L))):
        print(f"{j+1:>5}|{L[j]}")
