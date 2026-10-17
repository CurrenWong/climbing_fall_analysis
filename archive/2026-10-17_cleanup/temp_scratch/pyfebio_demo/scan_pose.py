"""定位 POSE_IE / POSE_DF 的施加代码 + 有无 no-run 开关。"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
print("--- 出现 POSE 的行 ---")
for i, l in enumerate(t.splitlines()):
    if "POSE" in l:
        print(f"{i+1:>5}: {l.strip()[:120]}")
print("\n--- no-run / build only ---")
for i, l in enumerate(t.splitlines()):
    if re.search(r"NO_RUN|BUILD_ONLY|--no-run|DRY", l):
        print(f"{i+1:>5}: {l.strip()[:120]}")
