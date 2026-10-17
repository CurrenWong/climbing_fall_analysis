"""决定性：加大指定位移让关节真闭合，看 F 是否还随罚因子变化。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
CFG = [("d1.5_p1", "1.5", "1"), ("d1.5_p001", "1.5", "0.001"),
       ("d3.0_p1", "3.0", "1"), ("d3.0_p001", "3.0", "0.001"),
       ("d4.0_p1", "4.0", "1"), ("d4.0_p001", "4.0", "0.001")]
base = dict(os.environ)
base.update({"TIE": "0", "LAUGON": "PENALTY", "TIME_STEPS": "20"})
for tag, disp, pen in CFG:
    env = dict(base)
    env.update({"DRIVER": "disp", "DISP_MM": disp, "PENALTY": pen})
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    print(f"  {tag:12s} rc={m.group(1) if m else '?'} "
          f"{m.group(2) if m else '-'}s", flush=True)
