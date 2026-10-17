"""打印 ankle_lig_feb.py 里所有 <bc ...> 的完整写法（照抄用）。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for k, l in enumerate(L):
    if "prescribed" in l or "zero displacement" in l:
        print(f"--- 行 {k+1} ---")
        for j in range(max(0, k - 2), min(k + 10, len(L))):
            print(f"{j+1:>5}|{L[j]}")
        print()
