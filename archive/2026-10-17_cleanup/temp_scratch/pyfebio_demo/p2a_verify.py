"""用正确 env 跑 builder，只看面积自检（验证修复）。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
env = dict(os.environ)
env.update({"DRIVER": "dispx", "DISP_MM": "+1.6", "DET_MIN": "0.0",
            "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
            "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
            "CART_DIR": str(R / "temp/thums/cart_m1"), "POSE_IE": "0",
            "PYTHONIOENCODING": "utf-8"})
r = subprocess.run([str(R / ".venv/Scripts/python.exe"), "-u",
                    str(R / "scripts/ankle_fe/ankle_lig_feb.py")],
                   cwd=str(R), env=env, capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
o = (r.stdout or "") + (r.stderr or "")
(Path(R) / "temp/pyfebio_demo/_bld2.txt").write_text(o, encoding="utf-8")
for l in o.splitlines():
    if re.search(r"area自检|顶面投影面积|外表面顶|FEBio rc|两顶面合力|剔除", l):
        print("  " + l.strip())
