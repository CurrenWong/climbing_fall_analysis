"""查韧带定义：名称 / 起止骨 / 内外侧标签是否与解剖一致。"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
ls = t.splitlines()
print("--- 命中韧带的定义行 ---")
for i, l in enumerate(ls):
    if re.search(r"ATFL|CFL|PTFL|MD\d|PTCL|LIG_DEF|LIGS|lig_def|deltoid|deltoid", l, re.I):
        print(f"{i+1:>5}: {l.strip()[:126]}")
