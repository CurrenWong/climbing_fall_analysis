"""p 收敛：penta6 (线性) vs penta15 (二次) vs tet4 (线性，作参照)
同网格、同位移，比支反力。
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


def run(elem, disp=0.1472, steps=20, extra=None):
    env = dict(os.environ)
    env.update(dict(MATSET="WK", LIG_MODE="shell", DRIVER="disp",
                    TOPFIX="lateral", TRUNC_MM="60", FIX_TALUS="all",
                    TIME_STEPS=str(steps), DISP_MM=str(disp),
                    CART_ELEM=elem))
    if extra:
        env.update({k: str(v) for k, v in extra.items()})
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(MODEL)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=3600)
    o = (r.stdout or "") + (r.stderr or "")
    m = re.search(r"FEBio rc=(\d+)", o)
    rc = int(m.group(1)) if m else None
    info = {}
    for pat, key in ((r"(tibiotalar_\w+)\s+(penta15|penta6|tet4)\s+保留\s+(\d+)\s+剔除\s+(\d+)", "cart"),
                     (r"接触面 (\w+): (\d+) 面", "surf"),
                     (r"\[prism15\][^\n]*", "p15")):
        info[key] = re.findall(pat, o)
    if "failed to converge" in o:
        info["err"] = "不收敛"
    elif "Negative jacobian" in o:
        info["err"] = "负Jacobian"
    return rc, info


def readF():
    src = sorted(W.glob("ankle_lig_WK*.xplt"),
                 key=lambda p: -p.stat().st_mtime)[0]
    dst = W / "_pc.h5"
    X.to_hdf5(str(src), str(dst))
    with h5py.File(dst, "r") as h:
        st = sorted(int(k) for k in h["states"])
        last = max(st)
        R = np.asarray(
            h[f"states/{last}/node_data/reaction forces/all"][()], float)
        D = np.asarray(h[f"states/{last}/node_data/displacement/all"][()],
                       float)
        # 软骨单元数与节点总数
        ed = [k for k in h[f"states/{last}/element_data/stress"]]
    return (float((R.reshape(-1, 3) @ AX).sum()),
            float(np.linalg.norm(D.reshape(-1, 3), axis=1).max()),
            len(st))


print("同网格、同位移 0.1472 mm 的 p 收敛对比")
print(f"{'软骨单元':>10} {'rc':>3} {'F(N)':>9} {'|u|max':>8} {'状态数':>6}  "
      f"{'备注':>10}")
print("-" * 68)
base = None
for elem in ("tet4", "prism", "prism15"):
    rc, info = run(elem)
    if rc == 0:
        F, um, ns = readF()
        note = ""
        if base is None:
            base = F
            note = "(基准)"
        else:
            note = f"×{F / base:.2f}"
        print(f"{elem:>10} {rc:>3} {F:9.1f} {um:8.3f} {ns:6d}  {note:>10}")
    else:
        print(f"{elem:>10} {rc:>3} {'-':>9} {'-':>8} {'-':>6}  "
              f"{info.get('err', '失败'):>10}")
    print(f"           {info.get('cart')} {info.get('p15')}")
