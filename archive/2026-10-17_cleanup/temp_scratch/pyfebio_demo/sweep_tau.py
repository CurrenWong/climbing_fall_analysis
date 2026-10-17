"""扫 Taubin 平滑强度，看负体积数怎么降。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
for tau in [0, 5, 10, 20, 40, 80]:
    e = dict(os.environ, CH_TAU=str(tau), CH_OUT=f"_probe_tau{tau}.npz")
    e["PATH"] = r"D:\Program\FEBioStudio\bin;" + e["PATH"]
    r = subprocess.run([str(ROOT / ".venv/Scripts/python.exe"), "-u",
                        str(ROOT / "scripts/ankle_fe/cart_hex.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=1800)
    out = r.stdout + r.stderr
    negs = [int(m) for m in re.findall(r"负体积\s+(\d+)", out)]
    gm = [float(m) for m in re.findall(r"γ中位=([\d.]+)", out)]
    tot = sum(negs)
    print(f"TAU={tau:3d}  负体积合计 {tot:6d} / 8952 ({100*tot/8952:5.1f}%)"
          f"  各处 {negs}  γ中位 {[f'{g:.3f}' for g in gm]}")
    if r.returncode != 0:
        print("   rc!=0:", out[-500:])
for tau in [0, 5, 10, 20, 40, 80]:
    p = ROOT / f"temp/thums/_probe_tau{tau}.npz"
    if p.exists():
        p.unlink()
