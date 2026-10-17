"""定性：这个配置下"关节力"到底走哪条路径？（接触 vs 韧带 vs 舟骨固定）"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

BASE = ["DRIVER=dispx", "DISP_MM=-1.0", "LOAD_SIGN=+1", "DET_MIN=0.0",
        "TRUNC_MM=60", "AXIS_SRC=full", "ALIGN_AXIS=1", "FIX_TALUS=transverse",
        "CART_NEW=1", "TIE=1", "CART_DIR=temp/thums/cart_m1"]

CASES = [
    ("基线（韧带开·接触开·舟骨固定）", []),
    ("无韧带", ["LIG_MODE=none"]),
    ("无接触", ["NO_CONTACT=1"]),
    ("无韧带+无接触", ["LIG_MODE=none", "NO_CONTACT=1"]),
    ("舟骨不固定", ["FIX_NAVICULAR=0"]),
    ("舟骨不固定+无韧带", ["FIX_NAVICULAR=0", "LIG_MODE=none"]),
    ("距骨全固定", ["FIX_TALUS=all"]),
    ("距骨全固定+无韧带", ["FIX_TALUS=all", "LIG_MODE=none"]),
]

print(f"{'用例':<28} {'rc':>3} {'NORMAL':>7} {'两顶面合力(N)':>14} "
      f"{'|u|max':>9}  备注")
print("─" * 96)
for name, extra in CASES:
    r = subprocess.run([str(PY), "-u", str(ACC)] + BASE + extra,
                       cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_path_{re.sub(r'[^A-Za-z0-9]+','_',name)}.log").write_text(
        out, encoding="utf-8")
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    m = re.search(r"两顶面合力 = \[[^\]]*\] N\s+\|F\|=([\d.eE+-]+)", out)
    u = re.search(r"\|u\|max = ([\d.]+) mm", out)
    note = ""
    if nm == "ERROR":
        mt = re.search(r"failed to converge at time : ([\d.]+)", out)
        note = f"失败 t={mt.group(1)}" if mt else ""
    print(f"{name:<28} {r.returncode:>3} {nm:>7} "
          f"{(m.group(1) if m else '—'):>14} {(u.group(1) if u else '—'):>9}  {note}")
print("\n★ 判读：若『无韧带』≈0 而『无接触』≈基线 ⇒ 载荷不走接触；")
print("  若『舟骨不固定+无韧带』也≈0 ⇒ 舟骨固定是足部唯一的锚。")
