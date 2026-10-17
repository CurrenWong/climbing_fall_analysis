"""定位 TALUS_RIGID 定义/使用位置 + 顶部 env 块。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if "TALUS_RIGID" in l or "TIE_SEC_TOL" in l or "TALUS_SUPPORT" in l or \
       "TAL_ROTLOCK" in l:
        print(f"{i:>5}|{l[:100]}")
