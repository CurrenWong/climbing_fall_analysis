"""干净的载荷路径分解：FIX_TALUS=all 保证"无韧带"也不奇异。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

BASE = ["DRIVER=dispx", "DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full",
        "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1",
        "CART_DIR=temp/thums/cart_m1"]
CASES = [
    ("距骨固定·压缩·韧带开·接触开", ["FIX_TALUS=all", "DISP_MM=+1.0"]),
    ("距骨固定·压缩·**韧带关**·接触开", ["FIX_TALUS=all", "DISP_MM=+1.0", "LIG_MODE=none"]),
    ("距骨固定·压缩·韧带开·**接触关**", ["FIX_TALUS=all", "DISP_MM=+1.0", "NO_CONTACT=1"]),
    ("距骨固定·压缩·**两者都关**", ["FIX_TALUS=all", "DISP_MM=+1.0", "LIG_MODE=none", "NO_CONTACT=1"]),
    ("距骨横约束·压缩·韧带开·接触开", ["FIX_TALUS=transverse", "DISP_MM=+1.0"]),
    ("距骨横约束·压缩·**韧带关**", ["FIX_TALUS=transverse", "DISP_MM=+1.0", "LIG_MODE=none"]),
]
print(f"{'用例':<34} {'rc':>3} {'NORMAL':>7} {'两顶面合力':>12} {'|u|max':>12}")
print("─" * 82)
for name, extra in CASES:
    r = subprocess.run([str(PY), "-u", str(ACC)] + BASE + extra,
                       cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_dp_{re.sub(r'[^A-Za-z0-9]+','_',name)}.log").write_text(
        out, encoding="utf-8")
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    m = re.search(r"两顶面合力 = \[[^\]]*\] N\s+\|F\|=([\d.eE+-]+)", out)
    u = re.search(r"\|u\|max = ([\d.]+) mm", out)
    print(f"{name:<34} {r.returncode:>3} {nm:>7} "
          f"{(m.group(1) if m else '—'):>12} {(u.group(1) if u else '—'):>12}")
print("\n★ 『接触关』与『全开』之差 = **接触贡献**；")
print("  『韧带关』与『全开』之差 = **韧带贡献**。")
