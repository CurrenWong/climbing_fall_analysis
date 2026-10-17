"""罚因子收敛验证：PENALTY = 1/2/4/8，看 k 是否收敛到物理值 E·A/t ≈ 1059 N/mm。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
base = dict(os.environ)
base.update({"TIE": "0", "LAUGON": "PENALTY", "TIME_STEPS": "20",
             "DRIVER": "disp", "DISP_MM": "0.1472"})
for pen in ("1", "2", "4", "8", "16"):
    env = dict(base)
    env["PENALTY"] = pen
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"], cwd=str(ROOT),
                       env=env, capture_output=True, text=True, timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    neg = "Negative jacobian" in out
    conv = "N O R M A L" in out.replace("  ", " ")
    print(f"  PENALTY={pen:>4s}  rc={m.group(1) if m else '?'}  "
          f"{m.group(2) if m else '-':>8s}s  正常终止={'Y' if conv else 'n'}"
          f"{' (负 Jacobian)' if neg else ''}", flush=True)
