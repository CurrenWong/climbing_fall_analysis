"""查双靶验证所需的现成工具：接触面积 / 软骨应力。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
print("=== scripts/ankle_fe 下的脚本 ===")
d = os.path.join(R, "scripts", "ankle_fe")
for f in sorted(os.listdir(d)):
    if f.endswith(".py"):
        print(f"  {f:<30} {os.path.getsize(os.path.join(d, f)):>8} B")

print("\n=== 全项目里提到 接触面积 / gap_area / contact_area / 应力 的脚本 ===")
for root, dirs, files in os.walk(R):
    dirs[:] = [x for x in dirs if x not in (".venv", ".git", "__pycache__")]
    for f in files:
        if not f.endswith(".py"):
            continue
        p = os.path.join(root, f)
        try:
            t = io.open(p, encoding="utf-8", errors="ignore").read()
        except OSError:
            continue
        if re.search(r"gap_area|contact_area|接触面[积]|von ?mises|stress|应力", t):
            rel = os.path.relpath(p, R)
            ks = []
            for k in ("gap_area", "contact_area", "接触面积", "von_mises", "stress",
                      "应力"):
                if re.search(k, t, re.I):
                    ks.append(k)
            print(f"  {rel:<52} {ks}")
