"""真·压缩曲线（DISP_MM>0），并分解载荷路径（接触 vs 韧带）。"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

BASE = ["DRIVER=dispx", "DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full",
        "ALIGN_AXIS=1", "FIX_TALUS=transverse", "CART_NEW=1", "TIE=1",
        "CART_DIR=temp/thums/cart_m1"]

SERIES = [
    ("全开（接触+韧带）", []),
    ("无韧带（只看接触）", ["LIG_MODE=none"]),
]

for lbl, extra in SERIES:
    print(f"\n{'='*82}\n=== {lbl}\n{'='*82}")
    print(f"{'DISP_MM':>8} {'rc':>3} {'NORMAL':>7} {'两顶面合力(N)':>14} "
          f"{'|u|max':>9}  备注")
    print("─" * 74)
    for d in ("+0.05", "+0.1", "+0.25", "+0.5", "+1.0", "+1.6", "+2.5"):
        r = subprocess.run([str(PY), "-u", str(ACC)] + BASE + extra
                           + [f"DISP_MM={d}"], cwd=str(ROOT),
                           capture_output=True, text=True, encoding="utf-8",
                           errors="replace", timeout=3600)
        out = (r.stdout or "") + (r.stderr or "")
        (ROOT / f"temp/pyfebio_demo/_cmp_{re.sub(r'[^A-Za-z0-9]+','_',lbl)}_{d}.log"
         ).write_text(out, encoding="utf-8")
        nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
        m = re.search(r"两顶面合力 = \[[^\]]*\] N\s+\|F\|=([\d.eE+-]+)", out)
        u = re.search(r"\|u\|max = ([\d.]+) mm", out)
        note = ""
        if nm == "ERROR":
            mt = re.search(r"failed to converge at time : ([\d.]+)", out)
            note = f"失败 t={mt.group(1)}" if mt else ""
        print(f"{d:>8} {r.returncode:>3} {nm:>7} "
              f"{(m.group(1) if m else '—'):>14} {(u.group(1) if u else '—'):>9}  {note}")
print("\n★ 若『无韧带』与『全开』同量级 ⇒ 压缩时确实由接触承力（与拉伸情形相反）。")
