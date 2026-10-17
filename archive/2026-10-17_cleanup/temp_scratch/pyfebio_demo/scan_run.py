"""定位 FEBio 子进程调用点，为加 NO_RUN 开关做准备。"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
ls = t.splitlines()
for i, l in enumerate(ls):
    if re.search(r"subprocess|\.feb.{0,20}write|write_text|FEBIO|rc\s*=|Popen|run\(", l):
        print(f"{i+1:>5}: {l.strip()[:126]}")
