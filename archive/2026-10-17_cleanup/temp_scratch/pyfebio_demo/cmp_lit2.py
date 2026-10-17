"""对齐 Anderson 参数的受控对比（含从 .xplt 读支反力）。"""
import os
import re
import textwrap
import inspect
import subprocess
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
PY = ROOT / ".venv/Scripts/python.exe"
W = ROOT / "temp/thums/febio_lig"
AX = np.array([0.8153, 0.135, 0.5631])

import pyfebio.xplt as X
_src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    _src = _src.replace(_o, _n)
_ns = X.__dict__
exec(compile(_src, "<p>", "exec"), _ns)
X.parse_mesh = _ns["parse_mesh"]
_ps = textwrap.dedent(inspect.getsource(X.parse_state))
for _a, _b in (
    ('set_name = mesh_dict[set_lut[key]][set_id + 1]',
     'set_name = mesh_dict[set_lut[key]].get(set_id + 1) or '
     'mesh_dict[set_lut[key]].get(int(set_id) + 1) or f"set{set_id}"'),
    ('set_name = mesh_dict[set_lut[key]][set_id]',
     'set_name = mesh_dict[set_lut[key]].get(set_id) or '
     'mesh_dict[set_lut[key]].get(int(set_id)) or f"set{set_id}"'),
):
    _ps = _ps.replace(_a, _b)
_ns2 = X.__dict__
exec(compile(_ps, "<p2>", "exec"), _ns2)
X.parse_state = _ns2["parse_state"]

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20", DISP_MM="0.1472")


def run(elem, extra):
    env = dict(os.environ)
    env.update(BASE); env.update(extra); env["CART_ELEM"] = elem
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    if rc != 0:
        return rc, None, None, None
    src = sorted(W.glob("ankle_lig_WK*.xplt"), key=lambda p: -p.stat().st_mtime)[0]
    dst = W / "_lit.h5"
    X.to_hdf5(str(src), str(dst))
    with h5py.File(dst, "r") as h:
        last = max(int(k) for k in h["states"])
        R = np.asarray(h[f"states/{last}/node_data/reaction forces/all"][()], float)
        Ud = np.asarray(h[f"states/{last}/node_data/displacement/all"][()], float)
        F = float((R.reshape(-1, 3) @ AX).sum())
        um = float(np.linalg.norm(Ud.reshape(-1, 3), axis=1).max())
    carts = re.findall(r"(tibiotalar_\w+)\s+(hex8|penta6|penta15|tet4)\s+保留\s+(\d+)\s+剔除\s+(\d+)", o)
    return rc, F, um, carts


VARIANTS = [("基准 WK 10/0.40", {}),
            ("Anderson 12/0.42", {"CART_E": "12.0", "CART_NU": "0.42"}),
            ("仅 ν=0.42", {"CART_NU": "0.42"}),
            ("仅 E=12", {"CART_E": "12.0"})]

print("软骨材料参数受控对比（DISP=0.1472, TRUNC=60, FIX_TALUS=all, TIME_STEPS=20）")
print("=" * 86)
print(f"{'参数':>20s} {'单元':>7s} {'rc':>3s} {'F(N)':>9s} {'|u|max':>8s}  软骨保留")
print("-" * 86)
base_v = {}
for label, extra in VARIANTS:
    for elem in ("tet4", "prism"):
        rc, F, um, carts = run(elem, extra)
        kept = "+".join(c[2] for c in carts) if carts else "-"
        sF = f"{F:9.1f}" if F is not None else f"{'-':>9s}"
        su = f"{um:8.3f}" if um is not None else f"{'-':>8s}"
        print(f"{label:>20s} {elem:>7s} {rc:>3d} {sF} {su}  {kept}")
        if label == "基准 WK 10/0.40" and F:
            base_v[elem] = F
    print()
