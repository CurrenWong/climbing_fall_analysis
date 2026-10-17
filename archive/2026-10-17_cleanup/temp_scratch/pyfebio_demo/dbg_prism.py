"""跑单个配置并显示错误（诊断 pconv 崩在 prism 上的原因）。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
env = dict(os.environ)
env.update(dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
                TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20",
                DISP_MM="0.1472", CART_ELEM="prism"))
env["PYTHONIOENCODING"] = "utf-8"
r = subprocess.run([str(ROOT / ".venv/Scripts/python.exe"), "-u",
                    str(ROOT / "scripts/ankle_fe/ankle_lig_feb.py")],
                   cwd=str(ROOT), env=env, capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
o = (r.stdout or "") + (r.stderr or "")
print("=== 尾部 2000 字符 ===")
print(o[-2000:])
