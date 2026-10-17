"""定位 Jacobian 自检代码块 + 姿势段的行号关系。"""
import io
import re

t = io.open(r"scripts\ankle_fe\ankle_lig_feb.py", encoding="utf-8").read()
ls = t.splitlines()
keys = ["Jacobian 自检", "DET_MIN", "朝向修正", "韧带壳退化剔除",
        "if POSE_DF or POSE_IE or POSE_TILT", "def _min_detJ_penta6",
        "def detJ_penta6", "已应用（输出坐标"]
for i, l in enumerate(ls):
    if any(k in l for k in keys):
        print(f"{i+1:>5}: {l.rstrip()[:126]}")
