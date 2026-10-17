"""抓 FEBio 原始输出里的负 Jacobian 上下文（看是哪些单元/材料）。"""
import io
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
SCRIPT = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
env = os.environ.copy()
env.update({"DRIVER": "dispx", "DISP_MM": "+2.0", "DET_MIN": "0.0",
            "TRUNC_MM": "60", "AXIS_SRC": "full", "ALIGN_AXIS": "1",
            "FIX_TALUS": "transverse", "CART_NEW": "1", "TIE": "1",
            "CART_DIR": str(ROOT / "temp/thums/cart_m1"), "POSE_IE": "0",
            "PYTHONIOENCODING": "utf-8"})
r = subprocess.run([str(PY), "-u", str(SCRIPT)], cwd=str(ROOT), env=env,
                   capture_output=True, text=True, encoding="utf-8",
                   errors="replace", timeout=3600)
out = (r.stdout or "") + (r.stderr or "")
(ROOT / "temp/pyfebio_demo/_raw20.log").write_text(out, encoding="utf-8")
print(f"原始输出 {len(out)} 字符\n" + "=" * 88)
for m in re.finditer(r"[^\n]*(?:jacobian|Jacobian|ELEMENT|Element \d|negative)[^\n]*",
                     out):
    s = m.group(0).strip()
    if s:
        print("  " + s[:160])
print("=" * 88)
i = out.lower().find("failed to converge")
if i > 0:
    print("失败上下文：")
    print(out[max(0, i - 900):i + 300])
else:
    print("尾部 1200 字符：")
    print(out[-1200:])
