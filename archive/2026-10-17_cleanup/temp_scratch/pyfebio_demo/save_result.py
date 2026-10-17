"""把韧带结果落盘 + 更新 ANKLE_LIG_TODO.md 到「模型已跑通、排序为负结果」。"""
import json
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUT = ROOT / "temp/thums"
PY = ROOT / ".venv/Scripts/python.exe"

# 跑「推荐配置」并保存结果
import os
e = dict(os.environ)
e.update({"LIG_MODE": "shell", "DRIVER": "disp", "TRUNC_MM": "60",
          "DISP_MM": "0.20", "FIX_TALUS": "all", "TIME_STEPS": "20"})
r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/ankle_lig_feb.py"],
                   cwd=str(ROOT), capture_output=True, text=True, timeout=900,
                   env=e)
print("solve:", "OK" if "NORMAL" in (r.stdout or "").replace(" ", "") else "FAIL")
rr = subprocess.run([str(PY), "temp/pyfebio_demo/read_lig.py", "WK"],
                    cwd=str(ROOT), capture_output=True, text=True, timeout=600)
(OUT / "ligament_results.txt").write_text(rr.stdout or "", encoding="utf-8")
print(f"结果 -> {OUT/'ligament_results.txt'} ({len(rr.stdout or '')} chars)")
