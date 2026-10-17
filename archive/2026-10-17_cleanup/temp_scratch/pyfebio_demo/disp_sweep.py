"""位移扫描（PENALTY=0.001，关节低罚）：F 何时脱离线性 ⇒ 关节何时真闭合。"""
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

import pttri

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

DISPS = ["0.1472", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"]
base = dict(os.environ)
base.update({"TIE": "0", "LAUGON": "PENALTY", "TIME_STEPS": "20",
             "PENALTY": "0.001", "DRIVER": "disp"})
rows = []
for d in DISPS:
    env = dict(base)
    env["DISP_MM"] = d
    r = subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"],
                       cwd=str(ROOT), env=env, capture_output=True, text=True,
                       timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)\s+用时\s+([\d.]+)s", out)
    rc = int(m.group(1)) if m else 99
    F = um = gmin = gmed = float("nan")
    if rc == 0:
        xp = W / "ankle_tibiotalar.xplt"
        fp = W / "ankle_tibiotalar.feb"
        shutil.copy2(xp, W / "_d.xplt")
        shutil.copy2(fp, W / "_d.feb")
        txt = fp.read_text(encoding="latin-1")
        XYZ = {}
        for mm in re.finditer(
                r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)',
                txt):
            XYZ[int(mm.group(1))] = np.array([float(mm.group(2)), float(mm.group(3)),
                                              float(mm.group(4))])
        sf = {}
        for mm in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
            sf[mm.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                               re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</",
                                          mm.group(2))]
        X.to_hdf5(str(xp), str(W / "_d.h5"))
        with h5py.File(W / "_d.h5", "r") as h:
            ns = sorted((int(q) for q in h["states"]), key=int)
            u1 = np.asarray(h[f"states/{ns[-1]}/node_data/displacement/all"][()], float)
            R = np.asarray(h[f"states/{ns[-1]}/node_data/reaction forces/all"][()], float)
        NX = max(XYZ) + 1
        U1 = np.zeros((NX, 3))
        U1[:len(u1)] = u1
        pn = sorted({x for f in sf["surf_tibiotalar_tibia"] for x in f})
        tris = np.array([np.array([XYZ[i] + U1[i - 1] for i in f])
                         for f in sf["surf_tibiotalar_talus"]])
        P = np.array([XYZ[i] + U1[i - 1] for i in pn])
        dvec = pttri.dist_pt_tris(P, tris)
        F = float((R @ AX).sum())
        um = float(np.linalg.norm(U1, axis=1).max())
        gmin = float(dvec.min())
        gmed = float(np.median(dvec))
    rows.append((d, rc, F, um, gmin, gmed))
    print(f"  DISP={d:>7s}  rc={rc}  F={F:9.2f} N  |u|max={um:8.4f}  "
          f"最小间隙={gmin:8.4f}  中位={gmed:8.4f}", flush=True)

print("\n" + "=" * 84)
print("位移扫描（关节 PENALTY=0.001，无 tie）")
print("=" * 84)
print(f"  {'DISP':>7} {'rc':>3} {'F(N)':>10} {'|u|max':>9} {'最小间隙':>10} "
      f"{'中位间隙':>10} {'F/DISP':>9}")
for r in rows:
    print(f"  {r[0]:>7} {r[1]:>3} {r[2]:>10.2f} {r[3]:>9.4f} {r[4]:>10.4f} "
          f"{r[5]:>10.4f} {r[2]/float(r[0]):>9.2f}")
