"""对齐 Anderson 文献参数的受控对比。

基准（MATSET=WK，软骨 E=10/ν=0.40，厚 1.5mm）：
    tet4 259.3 N / penta6 116.9 N
本脚本：CART_E=12, CART_NU=0.42（Anderson 2006/2007 的软骨参数）
"""
import os
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
PY = ROOT / ".venv/Scripts/python.exe"

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20", DISP_MM="0.1472")

VARIANTS = [
    ("基准 WK 10/0.40", {}),
    ("Anderson 12/0.42", {"CART_E": "12.0", "CART_NU": "0.42"}),
    ("ν=0.42（E 不变）", {"CART_NU": "0.42"}),
    ("E=12（ν 不变）", {"CART_E": "12.0"}),
]


def run(elem, extra):
    env = dict(os.environ)
    env.update(BASE)
    env.update(extra)
    env["CART_ELEM"] = elem
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    t = re.search(r"用时 ([\d.]+)s", o)
    # 支反力
    F = None
    for pat in (r"支反力\s*[:=]?\s*([\d.]+)\s*N", r"轴向支反力\s*[:=]?\s*([\d.]+)"):
        mm = re.search(pat, o)
        if mm:
            F = float(mm.group(1)); break
    carts = re.findall(r"(tibiotalar_\w+)\s+(hex8|penta6|penta15|tet4)\s+保留\s+(\d+)\s+剔除\s+(\d+)", o)
    return rc, (t.group(1) if t else "?"), F, carts, o


print("软骨材料/厚度对齐 Anderson 的受控对比（DISP=0.1472, TRUNC=60, FIX_TALUS=all）")
print("=" * 84)
for label, extra in VARIANTS:
    for elem in ("tet4", "prism"):
        rc, t, F, carts, o = run(elem, extra)
        kept = "+".join(c[2] for c in carts) if carts else "-"
        drop = "+".join(c[3] for c in carts) if carts else "-"
        print(f"  {label:20s} {elem:6s} rc={rc} {t:>6}s  "
              f"软骨保留 {kept:>6s} 剔除 {drop:>5s}  支反力 {F}")
    print()
