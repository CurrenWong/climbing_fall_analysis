"""看 LOG_DATA 相关代码：写出哪些文件、记录的列是什么。"""
import io

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
ls = t.splitlines()
hits = [i for i, l in enumerate(ls) if "LOG_DATA" in l or "lig_str" in l
        or "LOG_FILE" in l or "<logfile" in l]
for i in hits:
    lo, hi = max(0, i - 2), min(len(ls), i + 12)
    print(f"--- 行 {i+1} ---")
    for j in range(lo, hi):
        print(f"{j+1:>5}: {ls[j].rstrip()[:126]}")
    print()
