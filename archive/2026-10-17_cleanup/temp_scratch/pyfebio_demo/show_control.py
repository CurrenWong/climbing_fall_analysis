"""打印 ankle_lig_feb.py 的 <Control> 段。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
idx = [k for k, l in enumerate(L) if "<Control>" in l]
print("Control 行:", [k + 1 for k in idx])
for k in idx:
    for j in range(k, min(k + 30, len(L))):
        print(f"{j+1:>5}|{L[j]}")
    print("-" * 40)
