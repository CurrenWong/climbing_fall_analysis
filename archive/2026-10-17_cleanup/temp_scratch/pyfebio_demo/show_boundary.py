"""打印 Boundary / Loads 段的写出代码。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
idx = [k for k, l in enumerate(L) if "<Boundary>" in l or "<Loads>" in l
       or "</Loads>" in l or "DRIV" in l]
print("命中行:", [k + 1 for k in idx])
for k in idx:
    s = max(0, k - 3)
    for j in range(s, min(k + 12, len(L))):
        print(f"{j+1:>5}|{L[j][:105]}")
    print("-" * 50)
