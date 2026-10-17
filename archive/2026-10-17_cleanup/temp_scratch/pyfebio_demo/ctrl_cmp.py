"""受控对比：同一份『法向偏移』软骨下 hex8 vs penta6 vs tet4。
用 ankle_lig_feb.py（不是 build_feb.py！），env 与 pconv.py 基准一致。
"""
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
PY = ROOT / ".venv/Scripts/python.exe"
LOG = ROOT / "temp/pyfebio_demo/_ctrl.txt"

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20", DISP_MM="0.1472",
            CART_OFFSET="normal")


def run(elem):
    env = dict(os.environ)
    env.update(BASE)
    env["CART_ELEM"] = elem
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    log = LOG.with_name(f"_ctrl_{elem}.txt")
    log.write_text(o, encoding="utf-8")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    t = re.search(r"用时 ([\d.]+)s", o)
    types = re.search(r"按类型统计: (\{[^}]*\})", o)
    err = ""
    if "Negative jacobian" in o:
        err = "负Jacobian"
    elif "failed to converge" in o:
        err = "不收敛"
    carts = re.findall(r"(tibiotalar_\w+)\s+(hex8|penta6|penta15|tet4)\s+保留\s+(\d+)\s+剔除\s+(\d+)", o)
    surfs = re.findall(r"接触面 (\w+): (\d+) 面", o)
    return rc, (t.group(1) if t else "?"), (types.group(1) if types else "?"), err, carts, surfs


print("受控对比（同一份 CART_OFFSET=normal 软骨，env 与基准一致）")
print(f"{'CART_ELEM':>10} {'rc':>3} {'用时':>8}  软骨单元(保留/剔除)")
print("-" * 78)
for elem in ("tet4", "prism", "hex8"):
    rc, t, types, err, carts, surfs = run(elem)
    print(f"{elem:>10} {rc:>3} {t:>7}s  {err}")
    if types:
        print(f"           全模型类型: {types}")
    for c in carts:
        print(f"           {c[0]:20s} {c[1]:7s} 保留 {c[2]:>5s}  剔除 {c[3]:>4s}")
    if surfs:
        print(f"           接触面: {surfs}")
    print()
print(f"分日志: {LOG.parent}/_ctrl_<elem>.txt")
