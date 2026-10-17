"""关节接触罚因子收敛扫描（baseline + m2）：F 是否随 PENALTY 收敛到定值。"""
import inspect
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"
PY = ROOT / ".venv/Scripts/python.exe"

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
AX = np.array([0.8153, 0.135, 0.5631])
AX = AX / np.linalg.norm(AX)

CFG = []
for p in ("0.0001", "0.001", "0.01", "0.1", "1.0"):
    CFG.append((f"base_pen{p}", dict(TIE="0", LAUGON="PENALTY", PENALTY=p), "base"))
for p in ("0.0001", "0.001", "0.01", "0.1", "1.0"):
    CFG.append((f"m2_pen{p}", dict(CART_MESH="cartilage_mesh_patch_m2.npz",
                                   TIE="1", LAUGON="PENALTY", PENALTY=p), "m2"))

base = dict(os.environ)
base.update({"DRIVER": "disp", "DISP_MM": "0.1472", "TIME_STEPS": "20"})
rows = []
for tag, cfg, kind in CFG:
    env = dict(base)
    env.update(cfg)
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    rc, secs = (int(m.group(1)), float(m.group(2))) if m else (99, 0.0)
    F = um = uax = float("nan")
    nst = 0
    if rc == 0:
        xp = W / "ankle_tibiotalar.xplt"
        if xp.exists():
            X.to_hdf5(str(xp), str(W / f"_s_{tag}.h5"))
            with h5py.File(W / f"_s_{tag}.h5", "r") as h:
                ns = sorted((int(q) for q in h["states"]), key=int)
                nst = len(ns)
                U = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
                R = np.asarray(
                    h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
            F = float((R @ AX).sum())
            um = float(np.linalg.norm(U, axis=1).max())
            uax = float(np.abs(U @ AX).max())
    rows.append((tag, kind, cfg.get("PENALTY", "-"), rc, secs, nst, F, um, uax))
    print(f"  {tag:16s} rc={rc} {secs:6.1f}s  帧={nst:>3d}  F={F:8.2f} N  "
          f"|u|max={um:.4f}", flush=True)

print("\n" + "=" * 84)
print("关节接触罚因子收敛扫描")
print("=" * 84)
print(f"  {'tag':16s} {'kind':>5} {'PENALTY':>9} {'rc':>3} {'s':>7} {'帧':>4} "
      f"{'F(N)':>9} {'|u|max':>9} {'|u·ax|max':>10}")
for r in rows:
    print(f"  {r[0]:16s} {r[1]:>5} {r[2]:>9} {r[3]:>3} {r[4]:>7.1f} {r[5]:>4} "
          f"{r[6]:>9.2f} {r[7]:>9.4f} {r[8]:>10.4f}")
