"""ν 扫描：判定 2.2 倍刚度差是否来自线性 tet4 的体积锁死。

锁死判据：差效应应随 ν→0.5 放大（K = E/[3(1-2ν)] → ∞）。
"""
import os
import re
import sys
import textwrap
import inspect
import subprocess
import h5py
import numpy as np
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio_lig"
PY = ROOT / ".venv/Scripts/python.exe"
MODEL = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
AX = np.array([0.8153, 0.135, 0.5631])

BASE = dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp", TOPFIX="lateral",
            TRUNC_MM="60", FIX_TALUS="all", TIME_STEPS="20", DISP_MM="0.1472")


def run(cart_elem, nu):
    env = dict(os.environ)
    env.update({k: str(v) for k, v in BASE.items()})
    env["CART_ELEM"] = cart_elem
    env["CART_NU"] = str(nu)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时 ([\d.]+)s", out)
    rc = int(m.group(1)) if m else None
    note = ""
    if rc != 0:
        if "Negative jacobian" in out:
            note = "负Jacobian"
        elif "failed to converge" in out:
            note = "不收敛"
        else:
            note = "其他"
    return rc, note, out


def read_last():
    import pyfebio.xplt as X
    global _patched
    if not _patched:
        _src = textwrap.dedent(inspect.getsource(X.parse_mesh))
        for _o, _n in (('domain["name"]',
                        'domain.get("name", str(domain["id"][0]))'),
                       ('surface["name"]',
                        'surface.get("name", str(surface["id"][0]))'),
                       ("{domain['name']}",
                        "{domain.get('name', str(domain['id'][0]))}"),
                       ("{surface['name']}",
                        "{surface.get('name', str(surface['id'][0]))}"),
                       ('if nodeset["name"] == "":',
                        'if nodeset.get("name", "") == "":'),
                       ("{nodeset['name']}",
                        "{nodeset.get('name', str(nodeset['id'][0]))}"),
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
        _patched = True
    src = sorted(W.glob("ankle_lig_WK*.xplt"),
                 key=lambda p: -p.stat().st_mtime)[0]
    dst = W / "_sw.h5"
    X.to_hdf5(str(src), str(dst))
    with h5py.File(dst, "r") as h:
        st = sorted(int(k) for k in h["states"])
        R = np.asarray(
            h[f"states/{max(st)}/node_data/reaction forces/all"][()], float)
        D = np.asarray(h[f"states/{max(st)}/node_data/displacement/all"][()],
                       float)
    return float((R.reshape(-1, 3) @ AX).sum()), float(
        np.linalg.norm(D.reshape(-1, 3), axis=1).max())


_patched = False

print(f"{'ν':>6} | {'tet4 F(N)':>10} {'备注':>8} | {'penta6 F(N)':>11} "
      f"{'备注':>8} | {'比值 tet/prism':>13}")
print("-" * 84)
rows = []
for nu in (0.0, 0.2, 0.3, 0.4, 0.45, 0.49):
    r1, n1, o1 = run("tet4", nu)
    F1 = read_last()[0] if r1 == 0 else float("nan")
    r2, n2, o2 = run("prism", nu)
    F2 = read_last()[0] if r2 == 0 else float("nan")
    ratio = F1 / F2 if (F2 and F2 == F2 and F2 != 0) else float("nan")
    print(f"{nu:6.2f} | {F1:10.1f} {n1:>8} | {F2:11.1f} {n2:>8} | "
          f"{ratio:13.2f}")
    rows.append((nu, F1, F2, ratio))

print()
print("判读：若 tet4 的体积锁死是主因，则 ν=0 时比值应≈1，ν→0.5 时显著 >1。")
for nu, F1, F2, rt in rows:
    if rt == rt:
        print(f"  ν={nu:.2f}  tet4/penta6 = {rt:.2f}")
