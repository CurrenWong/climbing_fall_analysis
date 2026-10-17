"""判定实验：把加载截面从顶端移到关节上方 H mm，看关节是否闭合、F 是否摆脱罚因子。"""
import inspect
import os
import re
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
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
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


def get_geom(txt):
    XYZ = {}
    for m in re.finditer(
            r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)', txt):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    sf = {}
    for m in re.finditer(r'<Surface name="([^"]+)"[^>]*>(.*?)</Surface>', txt, re.S):
        sf[m.group(1)] = [tuple(int(y) for y in s.split(",")) for s in
                          re.findall(r"<[a-z0-9]+ id=[^>]*>\s*([\d,]+)\s*</", m.group(2))]
    els = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        els[m.group(2)] = {int(y) for s in
                           re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3))
                           for y in s.split(",")}
    return XYZ, sf, els


def run_case(H, pen, tag):
    """H=加载截面距关节的高度(mm)；pen=关节 PENALTY。返回 (rc,F,|u|max,min gap,median gap)."""
    env = dict(os.environ)
    env.update({"TIE": "0", "LAUGON": "PENALTY", "TIME_STEPS": "20",
                "DRIVER": "disp", "DISP_MM": "0.1472", "PENALTY": str(pen)})
    subprocess.run([str(PY), "-u", "scripts/ankle_fe/build_feb.py"], cwd=str(ROOT),
                   env=env, capture_output=True, text=True, timeout=1800)
    txt = (W / "ankle_tibiotalar.feb").read_text(encoding="latin-1")
    XYZ, sf, els = get_geom(txt)
    if H is not None:
        # 关节面投影
        jp = np.array([XYZ[i] @ AX for i in
                       sorted({x for f in sf["surf_tibiotalar_tibia"] for x in f})])
        j0 = jp.min()
        band = sorted(i for i in els["tibia_cort__hex8"]
                      if j0 + H - 2.0 <= XYZ[i] @ AX <= j0 + H + 2.0)
        if not band:
            print(f"  H={H} 找不到带状节点")
            return None
        old = re.search(r'<NodeSet name="tibia_top_ns">([^<]*)</NodeSet>', txt).group(1)
        txt2 = txt.replace(f'<NodeSet name="tibia_top_ns">{old}</NodeSet>',
                           f'<NodeSet name="tibia_top_ns">'
                           f'{",".join(str(x) for x in band)}</NodeSet>')
        txt = txt2
    feb = W / f"_h{H}_{tag}.feb"
    feb.write_text(txt, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(feb)], cwd=str(W),
                       capture_output=True, text=True, timeout=1800)
    so = (r.stdout or "") + (r.stderr or "")
    okc = "N O R M A L" in so.replace("  ", " ")
    xp = feb.with_suffix(".xplt")
    if not okc or not xp.exists():
        print(f"    [debug H={H} pen={pen}] rc={r.returncode} okc={okc} "
              f"xplt={xp.exists()}")
        print("    " + "\n    ".join(so.strip().splitlines()[-12:]))
        return ("FAIL", float("nan"), float("nan"), float("nan"), float("nan"))
    X.to_hdf5(str(xp), str(W / f"_h{H}_{tag}.h5"))
    with h5py.File(W / f"_h{H}_{tag}.h5", "r") as h:
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
    dv = pttri.dist_pt_tris(P, tris)
    return ("OK", float((R @ AX).sum()), float(np.linalg.norm(U1, axis=1).max()),
            float(dv.min()), float(np.median(dv)))


print("=" * 84)
print("加载截面位置 × 关节罚因子")
print("=" * 84)
print(f"  {'H(mm)':>7} {'PENALTY':>9} {'rc':>5} {'F(N)':>10} {'|u|max':>9} "
      f"{'最小间隙':>10} {'中位间隙':>10}")
res = {}
for H in (None, 60.0, 30.0):
    for pen in (1.0, 0.001):
        out = run_case(H, pen, f"p{pen}")
        if out:
            res[(H, pen)] = out
            print(f"  {str(H) if H else 'top':>7} {pen:>9} {out[0]:>5} {out[1]:>10.2f} "
                  f"{out[2]:>9.4f} {out[3]:>10.4f} {out[4]:>10.4f}", flush=True)
print("\n  判据：若同一 H 下 F(pen=1) ≈ F(pen=0.001) ⇒ 关节真闭合、F 可信")
