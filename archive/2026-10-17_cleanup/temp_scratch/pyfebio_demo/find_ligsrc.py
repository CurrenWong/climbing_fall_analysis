"""ankle_lig_feb.py 的韧带数据来源：读 ligaments.npz 还是自己从 THUMS 重建？"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
for i, l in enumerate(t.splitlines()):
    if re.search(r"ligaments|lig_npz|_NPZ|load\(|npz", l) and re.search(
            r"ligaments|_NPZ|_LIG|load\(|np\.load", l):
        print(f"{i+1:>5}: {l.rstrip()[:130]}")
print("\n--- 含 LIG_DEF / BONES_OF / truss 的行 ---")
for i, l in enumerate(t.splitlines()):
    if re.search(r"LIG_DEF|BONES_OF|truss|a_bone|b_bone", l):
        print(f"{i+1:>5}: {l.rstrip()[:130]}")
