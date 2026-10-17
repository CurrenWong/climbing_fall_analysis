"""二分：关掉接触，看 tet4 与 penta6 是否还差 2.2 倍。"""
import os
import re
import sys
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
sys.path.insert(0, str(Path(__file__).parent))
import fepatch


def run(elem, no_contact, lig="off"):
    env = dict(os.environ)
    env.update(dict(MATSET="WK", LIG_MODE=lig, DRIVER="disp",
                    TOPFIX="lateral", TRUNC_MM="60", FIX_TALUS="all",
                    TIME_STEPS="20", DISP_MM="0.1472", CART_ELEM=elem,
                    NO_CONTACT="1" if no_contact else "0"))
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    if rc != 0:
        return rc, float("nan"), float("nan"), o
    d = fepatch.load_last("_bc2.h5")
    return rc, d["F"], d["umax"], o


print("LIG_MODE=off（排除韧带），只比实体网格")
print(f"{'单元':>8} {'接触':>6} {'rc':>3} {'F(N)':>10} {'|u|max':>9}")
print("-" * 46)
for nc in (True, False):
    vals = {}
    for elem in ("tet4", "prism"):
        rc, F, um, o = run(elem, nc)
        vals[elem] = F
        print(f"{elem:>8} {'关' if nc else '开':>6} {rc:>3} {F:10.1f} {um:9.4f}")
    a, b = vals["tet4"], vals["prism"]
    if b == b and b != 0:
        print(f"         -> tet4/penta6 = {a / b:.3f}   "
              f"({'一致(<15%)' if abs(a / b - 1) < 0.15 else '仍有差异'})")
    print()
