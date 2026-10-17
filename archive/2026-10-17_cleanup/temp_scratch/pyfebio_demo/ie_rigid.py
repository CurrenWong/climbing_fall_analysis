"""嫌疑：TALUS_RIGID 与足部旋转不兼容。测变形距骨 + 各种组合。"""
import os
import re
import subprocess
from pathlib import Path

R = Path(r"D:\Project\climbing_fall_analysis")
W = R / "temp" / "thums" / "febio_lig"
FEB = r"D:\Program\FEBioStudio\bin\febio4.exe"
PY = str(R / ".venv" / "Scripts" / "python.exe")
CASES = [
    ("刚体距骨 IE+10", dict(POSE_IE="+10", TALUS_RIGID="1")),
    ("变形距骨 IE+10", dict(POSE_IE="+10", TALUS_RIGID="0")),
    ("变形距骨 IE-10", dict(POSE_IE="-10", TALUS_RIGID="0")),
    ("变形距骨 DF+10", dict(POSE_IE="0", POSE_DF="+10", TALUS_RIGID="0")),
    ("变形距骨 无姿势", dict(TALUS_RIGID="0")),
]
for tag, ov in CASES:
    if (W / "ankle_lig_WK.feb").exists():
        (W / "ankle_lig_WK.feb").unlink()
    e = dict(os.environ)
    e.update(DRIVER="dispx", DISP_MM="-1.0", LOAD_SIGN="+1",
             FIX_TALUS="transverse", TRUNC_MM="60", AXIS_SRC="full",
             ALIGN_AXIS="1", CART_NEW="1", TIE="1", DET_MIN="0.0",
             TAL_ROTLOCK="0", TIME_STEPS="20",
             LAUGON="PENALTY", CART_DIR=str(R / "temp" / "thums" / "cart_m1"))
    e.update(ov)
    p = subprocess.run([PY, "-u", str(R / "scripts" / "ankle_fe" / "ankle_lig_feb.py")],
                       capture_output=True, text=True, env=e, cwd=str(R), timeout=900)
    if not (W / "ankle_lig_WK.feb").exists():
        print(f"{tag:<22} 构建失败 rc={p.returncode}", flush=True)
        continue
    r = subprocess.run([FEB, "-i", "ankle_lig_WK.feb", "-o", "s.xplt"],
                       capture_output=True, text=True, cwd=str(W), timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    fin = [float(x) for x in re.findall(r"converged at time : ([\d.]+)", out)]
    mx = max(fin) if fin else None
    print(f"{tag:<22} t={mx if mx else '-':<8} "
          f"{'NORMAL ✅' if 'N O R M A L' in out else 'FAIL ❌'}", flush=True)
