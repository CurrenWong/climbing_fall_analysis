"""收敛研究（AUGLAG）—— 修正后处理：取**最后** state，按域名分离软骨应力。"""
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
    X.to_hdf5(str(WORK / "ankle_tibiotalar.xplt"), str(WORK / "_cv.h5"))
    r = {}
    with h5py.File(WORK / "_cv.h5", "r") as h:
        states = sorted(int(k.split("/")[-1]) for k in h["states"])
        S = states[-1]
        base = f"/states/{S}/node_data"
        u = np.asarray(h[base + "/displacement/all"][()], float)
        rf = np.asarray(h[base + "/reaction forces/all"][()], float)
        r["nstate"] = len(states)
        r["umax"] = float(np.linalg.norm(u, axis=1).max())
        r["F"] = float(np.linalg.norm(rf.sum(axis=0)))
        r["Fvec"] = rf.sum(axis=0)
        # 软骨应力：域名不含骨名
        ebase = f"/states/{S}/element_data"
        found = {}

        def walk(gr, pre=""):
            for k in gr:
                if isinstance(gr[k], h5py.Group):
                    walk(gr[k], pre + "/" + k)
                elif "stress" in pre:
                    found[pre + "/" + k] = np.asarray(gr[k][()], float)
        walk(h[ebase])
        cart = {k: v for k, v in found.items()
                if "stress" in k and not any(b in k for b in BONE)}
        allv = []
        for k, v in cart.items():
            # von Mises 用前 6 个分量近似：这里取 |σ| 的 p99（分量范数）
            if v.ndim == 2 and v.shape[1] >= 6:
                allv.append(np.linalg.norm(v[:, :6], axis=1))
        if allv:
            a = np.concatenate(allv)
            r["cart_p99"] = float(np.percentile(a, 99))
            r["cart_max"] = float(a.max())
            r["cart_n"] = int(a.size)
            r["cart_keys"] = list(cart)
    return r


env0 = dict(os.environ)
env0["PATH"] = r"D:\Program\FEBioStudio\bin;" + env0["PATH"]
env0.update(CP_DOMAIN="bary", CP_DIAG="02", TIE="1", DRIVER="disp",
            DISP_MM="0.1472", TIME_STEPS="20", LAUGON="AUGLAG")

print(f"{'m':>2s} {'底三角':>7s} {'软骨tet':>8s} {'接触面(tib/tal)':>16s} "
      f"{'rc':>3s} {'耗时':>8s} {'|u|max':>9s} {'|F|(N)':>9s} {'软骨p99':>9s}")
for m in (1, 2, 3):
    out = f"cartilage_mesh_patch_m{m}.npz"
    e = dict(env0, CP_M=str(m), CP_OUT=out)
    r = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/cart_patch.py")],
                       cwd=ROOT, env=e, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    nc = [int(x) for x in re.findall(r"面\s+(\d+)\(", r.stdout)]
    e2 = dict(env0, CART_MESH=out)
    r2 = subprocess.run([PY, "-u", str(ROOT / "scripts/ankle_fe/build_feb.py")],
                        cwd=ROOT, env=e2, capture_output=True, text=True,
                        encoding="utf-8", errors="replace", timeout=7200)
    o = r2.stdout + r2.stderr
    cf = re.findall(r"外层接触面 (\d+)", o)
    rc = re.search(r"FEBio rc=(\d+)", o)
    tm = re.search(r"用时 ([\d.]+)s", o)
    nt = re.findall(r"保留\s+([\d,]+)\s+剔除", o)
    cart_n = sum(int(x.replace(",", "")) for x in nt[-2:]) if nt else 0
    if rc and rc.group(1) == "0":
        p = post()
        print(f"{m:>2d} {sum(nc):>7d} {cart_n:>8d} {str(cf):>16s} "
              f"{rc.group(1):>3s} {tm.group(1)+'s':>8s} {p['umax']:>9.4f} "
              f"{p['F']:>9.3f} {p.get('cart_p99', 0):>9.4f}")
    else:
        print(f"{m:>2d} {sum(nc):>7d} {cart_n:>8d} {str(cf):>16s} "
              f"{rc.group(1) if rc else '?':>3s} FAIL")
