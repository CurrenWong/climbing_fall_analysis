"""查 ns_of 的定义位置与缩进（是否在条件块内）。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if "def ns_of" in l:
        print(f"--- 定义在 行 {i}，缩进 {len(l) - len(l.lstrip())} ---")
        for j in range(max(0, i - 12), min(i + 8, len(L))):
            mark = ">>" if j == i - 1 else "  "
            print(f"{mark}{j+1:>5}|{L[j]}")
