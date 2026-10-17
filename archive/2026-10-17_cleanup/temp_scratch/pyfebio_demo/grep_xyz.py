"""列出 XYZ / coords 在 ankle_lig_feb.py 里的所有使用位置。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
print("=== XYZ 使用 ===")
for i, l in enumerate(L, 1):
    if "XYZ" in l:
        print(f"{i:>5}|{l[:110]}")
print("\n=== coords 使用 ===")
for i, l in enumerate(L, 1):
    if "coords" in l:
        print(f"{i:>5}|{l[:110]}")
