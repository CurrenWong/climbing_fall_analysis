"""POSE_IE=+10 的原始 FEBio 输出：失败在哪（t=0 还是加载中）。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
for tag, ov in (("tilt 10/0（已通）", dict(POSE_TILT="10")),
                ("IE +10（失败）", dict(POSE_IE="+10"))):
    if (W / "ankle_lig_WK.feb").exists():
        (W / "ankle_lig_WK.feb").unlink()
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
             FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
             TALUS_RIGID="1", TAL_ROTLOCK="0", TIME_STEPS="20",
             LAUGON="PENALTY", CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
    e.update(ov)
    subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                   capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    if not (W / "ankle_lig_WK.feb").exists():
        print(f"{tag}: 构建失败")
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    print(f"\n===== {tag} =====")
    for l in out.splitlines():
        ls = l.strip()
        if any(k in ls for k in ("WARNING", "ERROR", "negative", "failed",
                                 "beginning time step", "jacobian",
                                 "invalid facet", "converged")):
            print("  " + ls[:100])
