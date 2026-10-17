"""跑 prism15 并显示错误细节。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
env = dict(os.environ)
env.update(dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
                TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="5",
                DISP_MM="0.05", CART_ELEM="prism15"))
env["PYTHONIOENCODING"] = "utf-8"
r = subprocess.run([str(ROOT / ".venv/Scripts/python.exe"), "-u",
                    str(ROOT / "scripts/ankle_fe/ankle_lig_feb.py")],
                   cwd=str(ROOT), env=env, capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
o = (r.stdout or "") + (r.stderr or "")
print("--- 关键行 ---")
for l in o.split("\n"):
    if any(k in l for k in ("rc=", "ERROR", "error", "Negative", "failed",
                            "FAILED", "unrecognized", "prism15", "接触面",
                            "tibiotalar", "NodeSet", "Surface")):
        print("  ", l.strip()[:130])
print("\n--- 尾部 600 字符 ---")
print(o[-600:])
