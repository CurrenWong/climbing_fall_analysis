"""查韧带相关的 env 开关（LIG_* / NO_LIG）在脚本里怎么用。"""
import io
import re

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
L = io.open(P, encoding="utf-8").read().splitlines()
for i, l in enumerate(L, 1):
    if re.search(r"NO_LIG|LIG_MODE|LIG_SUB|LIG_ONLY|environ.*LIG", l):
        print(f"{i:>5}|{l[:105]}")
