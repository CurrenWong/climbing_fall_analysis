"""AUGLAG 的 F 是否与初始罚因子无关 ⇒ AUGLAG 是否给出可信解。"""
import inspect
import os
import re
import subprocess
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
PY = ROOT / ".venv/Scripts/python.exe"
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)
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

base = dict(os.environ)
base.update({"TIE": "0", "LAUGON": "AUGLAG", "TIME_STEPS": "20",
             "DRIVER": "disp", "DISP_MM": "0.1472"})
rows = []
for pen in ("0.001", "0.01", "0.1", "1.0"):
    env = dict(base)
    env["PENALTY"] = pen
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    rc = int(m.group(1)) if m else 99
    F = um = float("nan")
    nst = 0
    if rc == 0:
        xp = W / "ankle_tibiotalar.xplt"
        X.to_hdf5(str(xp), str(W / "_aug.h5"))
        with h5py.File(W / "_aug.h5", "r") as h:
            ns = sorted((int(q) for q in h["states"]), key=int)
            nst = len(ns)
            u1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
            R = np.asarray(h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
        F = float((R @ AX).sum())
        um = float(np.linalg.norm(u1, axis=1).max())
    rows.append((pen, rc, nst, F, um))
    print(f"  AUGLAG 初始penalty={pen:>6s}  rc={rc}  帧={nst:>3d}  "
          f"F={F:8.2f} N  |u|max={um:.4f}", flush=True)

print("\n" + "=" * 70)
print("AUGLAG 对初始罚因子的敏感性")
print("=" * 70)
print(f"  {'初始penalty':>12} {'rc':>3} {'帧':>4} {'F(N)':>10} {'|u|max':>9}")
for r in rows:
    print(f"  {r[0]:>12} {r[1]:>3} {r[2]:>4} {r[3]:>10.2f} {r[4]:>9.4f}")
fs = [r[3] for r in rows if r[1] == 0]
if len(fs) > 1:
    print(f"\n  F 极差 = {max(fs)-min(fs):.2f} N  "
          f"（相对 {(max(fs)-min(fs))/np.mean(fs)*100:.1f}%）")
