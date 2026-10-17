"""双靶验证 @+1.6mm（592.6 N）：读 软骨接触压力 + 接触面积，对吴恺 3.74 MPa / 335.5 mm²。"""
import io
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

for d in ("+1.6", "+1.0"):
    r = subprocess.run([str(PY), "-u", str(ACC), "DRIVER=dispx",
                        f"DISP_MM={d}", "DET_MIN=0.0", "TRUNC_MM=60",
                        "AXIS_SRC=full", "ALIGN_AXIS=1",
                        "FIX_TALUS=transverse", "CART_NEW=1", "TIE=1",
                        "CART_DIR=temp/thums/cart_m1", "POSE_IE=0"],
                       cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    tag = d.replace("+", "p").replace(".", "_")
    (ROOT / f"temp/pyfebio_demo/_dual_{tag}.log").write_text(out, encoding="utf-8")
    print("=" * 88)
    print(f"=== DISP_MM={d}")
    print("=" * 88)
    print(out.strip()[:3000])
