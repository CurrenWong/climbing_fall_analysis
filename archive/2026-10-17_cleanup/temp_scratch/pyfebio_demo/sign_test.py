"""决定性符号实验：哪个 DISP_MM/LOAD_SIGN 才是真正"压实"？

关掉韧带（LIG_MODE=none）⇒ 测到的力**只可能来自接触**。
压实 ⇒ 接触承力（应给大而正的阻力）；拉开 ⇒ 接触断开（应≈0）。
"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

BASE = ["DRIVER=dispx", "DET_MIN=0.0", "TRUNC_MM=60", "AXIS_SRC=full",
        "ALIGN_AXIS=1", "CART_NEW=1", "TIE=1",
        "CART_DIR=temp/thums/cart_m1", "LIG_MODE=none"]

CASES = []
for fix in ("transverse", "all"):
    for sign in ("+1", "-1"):
        for d in ("-1.0", "+1.0"):
            CASES.append((f"FIX={fix:10s} SIGN={sign:>2} DISP={d:>5}",
                          [f"FIX_TALUS={fix}", f"LOAD_SIGN={sign}",
                           f"DISP_MM={d}"]))

print(f"{'用例':<34} {'rc':>3} {'NORMAL':>7} {'两顶面合力(N)':>14} {'|u|max':>9}")
print("─" * 84)
for name, extra in CASES:
    r = subprocess.run([str(PY), "-u", str(ACC)] + BASE + extra,
                       cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    (ROOT / f"temp/pyfebio_demo/_sc_{re.sub(r'[^A-Za-z0-9]+','_',name)}.log").write_text(
        out, encoding="utf-8")
    nm = "NORMAL" if "N O R M A L   T E R M I N A T I O N" in out else "ERROR"
    m = re.search(r"两顶面合力 = \[[^\]]*\] N\s+\|F\|=([\d.eE+-]+)", out)
    u = re.search(r"\|u\|max = ([\d.]+) mm", out)
    print(f"{name:<34} {r.returncode:>3} {nm:>7} "
          f"{(m.group(1) if m else '—'):>14} {(u.group(1) if u else '—'):>9}")
    # 顺带把 BC 里写的 z 位移值抓出来
    feb = ROOT / "temp/thums/febio_lig/ankle_lig_WK.feb"
    if feb.exists():
        t = feb.read_text(encoding="latin-1", errors="replace")
        mm = re.search(r'<bc name="tibia_top_z".*?<value[^>]*>([-\d.]+)</value>',
                       t, re.S)
        if mm:
            print(f"{'':<34}   ↳ .feb 写入 tibia_top_z 值 = {mm.group(1)}")
print("\n★ 判读：给出**大阻力**的那个方向才是压实；≈0 的是拉开。")
