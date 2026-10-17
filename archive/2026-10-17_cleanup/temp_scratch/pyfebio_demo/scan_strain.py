"""找逐韧带应变的计算代码。"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
ls = t.splitlines()
for i, l in enumerate(ls):
    if re.search(r"应变|strain|LIG_NAMES|lig_names", l):
        print(f"{i+1:>5}: {l.rstrip()[:130]}")
