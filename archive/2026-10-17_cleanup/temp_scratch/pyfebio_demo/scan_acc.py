"""看 acceptance.py / ankle_lig_feb.py 里已有的 Kabsch / 节点集 / 姿势提取能力。"""
import io
import re

for f in (r"scripts\ankle_fe\acceptance.py", r"scripts\ankle_fe\ankle_lig_feb.py"):
    t = io.open(f, encoding="utf-8").read()
    print("=" * 26, f, len(t.encode("utf-8")), "B", len(t.splitlines()), "行")
    for i, l in enumerate(t.splitlines()):
        m = re.match(r"^(def |class |[A-Z_]{3,}\s*=|#\s*---)", l)
        if m:
            print(f"{i+1:>5}: {l.strip()[:118]}")
