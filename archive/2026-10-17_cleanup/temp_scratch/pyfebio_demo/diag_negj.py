"""诊断 +2.0 压缩失败的 12 个负 Jacobian 是哪些单元（是不是软骨）。"""
import io
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

r = subprocess.run([str(PY), "-u", str(ACC), "DRIVER=dispx", "DISP_MM=+2.0",
                    "DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1",
                    "FIX_TALUS=transverse", "CART_NEW=1", "TIE=1",
                    "CART_DIR=temp/thums/cart_m1", "POSE_IE=0"],
                   cwd=str(ROOT), capture_output=True, text=True,
                   encoding="utf-8", errors="replace", timeout=3600)
out = (r.stdout or "") + (r.stderr or "")
(ROOT / "temp/pyfebio_demo/_d20.log").write_text(out, encoding="utf-8")
print(f"总输出 {len(out)} 字符")
print("=" * 86)
print("含 negative jacobian / 单元 / mat 的行")
for i, l in enumerate(out.splitlines()):
    if re.search(r"negative jacobian|Element|element|mat\d|非正|Jacobian", l, re.I):
        print(f"  {i+1:>5}: {l.strip()[:150]}")
print("=" * 86)
print("失败点上下文：")
m = re.search(r"[^\n]*failed to converge[^\n]*", out)
if m:
    s = max(0, m.start() - 700)
    print(out[s:m.end() + 200])
else:
    print("（未找到 failed to converge）")
