"""扫姿势角度量级：找几何旋转式的可用范围。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
print(f"{'case':<20} {'落位后间距':>10} {'t':>10}  结果", flush=True)
for ax, val in (("IE", "1"), ("IE", "2"), ("IE", "3"), ("IE", "5"),
                ("DF", "2"), ("DF", "5")):
    if (W / "ankle_lig_WK.feb").exists():
        (W / "ankle_lig_WK.feb").unlink()
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
             FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
             TAL_TOMLOCK="0", TAL_ROTLOCK="0", TIME_STEPS="20",
             TALUS_RIGID="0", LAUGON="PENALTY",
             CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
    e[f"POSE_{ax}"] = val
    p = subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    so = (p.stdout or "") + (p.stderr or "")
    if not (W / "ankle_lig_WK.feb").exists():
        print(f"{ax}{val:<18}  构建失败 rc={p.returncode}", flush=True)
        continue
    g = re.search(r"落位后 ([\d.]+) mm", so)
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{ax}{val:<18} {g.group(1) if g else '?':>10} "
          f"{mx if mx else '-':>10}  "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
