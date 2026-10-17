"""① tie 罚因子扫描 ② 够硬罚因子下的收敛研究（自包含）。"""
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
THUMS = ROOT / "temp/thums"
WORK = THUMS / "febio"
PY = str(ROOT / ".venv/Scripts/python.exe")
BONE = ("tibia_cort", "tibia_end_spon", "tibia_marrow", "talus_cort",
        "talus_spon", "calcaneus_cort", "calcaneus_spon")

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


def post():
    X.to_hdf5(str(WORK / "ankle_tibiotalar.xplt"), str(WORK / "_pc.h5"))
    r = {}
    with h5py.File(WORK / "_pc.h5", "r") as h:
        S = max(int(k.split("/")[-1]) for k in h["states"])
        b = f"/states/{S}/node_data"
        u = np.asarray(h[b + "/displacement/all"][()], float)
        rf = np.asarray(h[b + "/reaction forces/all"][()], float)
        r["umax"] = float(np.linalg.norm(u, axis=1).max())
        r["F"] = float(np.linalg.norm(rf.sum(axis=0)))
        found = {}

        def walk(gr, pre=""):
            for k in gr:
                if isinstance(gr[k], h5py.Group):
                    walk(gr[k], pre + "/" + k)
                elif "stress" in pre:
                    found[pre + "/" + k] = np.asarray(gr[k][()], float)
        walk(h[f"/states/{S}/element_data"])
        vs = [np.linalg.norm(v[:, :6], axis=1) for k, v in found.items()
              if "stress" in k and not any(x in k for x in BONE)
              and v.ndim == 2 and v.shape[1] >= 6]
        if vs:
            a = np.concatenate(vs)
            r["p99"] = float(np.percentile(a, 99))
    return r


env0 = dict(os.environ)
env0["PATH"] = r"D:\Program\FEBioStudio\bin;" + env0["PATH"]
env0.update(TIE="1", DRIVER="disp", DISP_MM="0.1472", TIME_STEPS="20",
            LAUGON="AUGLAG", CP_DOMAIN="bary", CP_DIAG="02")


def run(mesh, pen):
    e = dict(env0, CART_MESH=mesh, TIE_PENALTY=str(pen))
    r = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/build_feb.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=7200)
    o = r.stdout + r.stderr
    rc = re.search(r"FEBio rc=(\d+)", o)
    tm = re.search(r"用时 ([\d.]+)s", o)
    if not (rc and rc.group(1) == "0"):
        return None, f"FAIL rc={rc.group(1) if rc else '?'}"
    return post(), (tm.group(1) + "s" if tm else "?")


print("=" * 76)
print("① tie 罚因子扫描 @ m=1")
print("=" * 76)
print(f"{'pen':>5s} {'耗时':>8s} {'|u|max':>9s} {'|F|(N)':>9s} {'软骨p99':>9s}")
for pen in (1, 4, 16, 64):
    p, info = run("cartilage_mesh_patch_m1.npz", pen)
    print(f"{pen:>5d} {info:>8s} " + (f"{p['umax']:>9.4f} {p['F']:>9.3f} "
                                      f"{p.get('p99', 0):>9.4f}" if p else ""))

print()
print("=" * 76)
print("② 收敛研究 @ TIE_PENALTY=16")
print("=" * 76)
print(f"{'m':>2s} {'耗时':>8s} {'|u|max':>9s} {'|F|(N)':>9s} {'软骨p99':>9s}")
for m in (1, 2, 3, 4):
    mesh = f"cartilage_mesh_patch_m{m}.npz"
    if not (THUMS / mesh).exists():
        e = dict(env0, CP_M=str(m), CP_OUT=mesh)
        subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/cart_patch.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    p, info = run(mesh, 16)
    print(f"{m:>2d} {info:>8s} " + (f"{p['umax']:>9.4f} {p['F']:>9.3f} "
                                    f"{p.get('p99', 0):>9.4f}" if p else ""))
