"""批量诊断：分别关接触 / 距骨全约束 / 加密时间步，看哪个能让模型收敛。"""
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
OUTD = ROOT / "temp/pyfebio_demo"
PY = ROOT / ".venv/Scripts/python.exe"

CASES = {
    "A_fixtalus_all": {"LIG_MODE": "shell", "FIX_TALUS": "all",
                       "NO_ROTLOCK": "1"},
    "B_nocontact": {"LIG_MODE": "shell", "NO_CONTACT": "1",
                    "NO_ROTLOCK": "1"},
    "C_steps50": {"LIG_MODE": "shell", "TIME_STEPS": "50",
                  "NO_ROTLOCK": "1"},
}

import os
for name, env in CASES.items():
    e = dict(os.environ)
    e.update(env)
    e.setdefault("FORCE_N", "600")
    r = subprocess.run([str(PY), "scripts/ankle_fe/ankle_lig_feb.py"],
                       cwd=str(ROOT), capture_output=True, text=True,
                       timeout=5400, env=e)
    o = (r.stdout or "") + (r.stderr or "")
    (OUTD / f"_diag_{name}.txt").write_text(o, encoding="utf-8")
    lines = o.split("\n")
    neg = next((l.strip() for l in lines if "negative jacobian" in l.lower()),
               None)
    conv = next((l.strip() for l in lines if "failed to converge at time"
                 in l.lower()), None)
    done = next((l.strip() for l in lines
                 if "Number of time steps completed" in l), None)
    maxit = next((l.strip() for l in lines
                  if "Max nr of iterations reached" in l), None)
    ok = "NORMAL TERMINATION" in o or "T E R M I N A T I O N" in o and \
         "N O R M A L" in o
    print(f"\n=== {name} {env} ===")
    print(f"   NORMAL={ok}  {done}")
    print(f"   {neg or '(无负Jacobian)'}")
    print(f"   {conv or '(未报收敛失败)'}")
    print(f"   {maxit or ''}")
