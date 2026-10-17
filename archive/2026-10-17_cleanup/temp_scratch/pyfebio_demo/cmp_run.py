"""通用读数：对任意 .feb/.xplt 对，报增量放大比 + |u|max 部件归属 + 支反力。"""
import inspect
import re
import sys
import textwrap
from pathlib import Path

import h5py
import numpy as np
import pyfebio.xplt as X

ROOT = Path(r"D:\Project\climbing_fall_analysis")
W = ROOT / "temp/thums/febio"

for _o, _n in (('domain["name"]', 'domain.get("name", str(domain["id"][0]))'),
               ('surface["name"]', 'surface.get("name", str(surface["id"][0]))'),
               ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
               ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
               ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
               ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
               ("{elementset['name']}",
                "{elementset.get('name', str(elementset['id'][0]))}")):
    _src = _o
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
GIVEN = 0.1472


def report(tag, feb_p, xplt_p):
    txt = Path(feb_p).read_text(encoding="latin-1")
    part_of = {}
    for m in re.finditer(
            r'<Elements\b[^>]*type="(tet4|hex8|penta6)"[^>]*name="([^"]+)"[^>]*>(.*?)</Elements>',
            txt, re.S):
        for s in re.findall(r"<elem\b[^>]*>\s*([\d,]+)\s*</elem>", m.group(3)):
            for x in s.split(","):
                part_of[int(x)] = m.group(2)
    for n in ("tiesurf_tibiotalar_tibia", "tiesurf_tibiotalar_talus"):
        mm = re.search(rf'<Surface name="{n}"[^>]*>(.*?)</Surface>', txt, re.S)
        if mm:
            for s in re.findall(r"<tri3\b[^>]*>\s*([\d,]+)\s*</", mm.group(1)):
                for x in s.split(","):
                    part_of.setdefault(int(x), "cart(" + n.replace("tiesurf_", "") + ")")

    dst = W / f"_{tag}.h5"
    X.to_hdf5(str(xplt_p), str(dst))
    with h5py.File(dst, "r") as h:
        ns = sorted((int(k) for k in h["states"]), key=int)
        U = {s: np.asarray(h[f"states/{s}/node_data/displacement/all"][()], float)
             for s in ns}
        R = {s: np.asarray(h[f"states/{s}/node_data/reaction forces/all"][()], float)
             for s in ns}
    n = len(ns) - 1
    print("\n" + "#" * 76)
    print(f"# {tag}   （{Path(feb_p).name}）")
    print("#" * 76)
    print(f"  {'帧':>4} {'施加':>9} {'|u|max':>10} {'放大×':>7} {'F(N)':>9}")
    for s in ns:
        t = s / n
        um = np.linalg.norm(U[s], axis=1).max()
        print(f"  {s:>4d} {t*GIVEN:>9.5f} {um:>10.4f} "
              f"{(um/(t*GIVEN) if t else float('nan')):>7.1f} "
              f"{(R[s] @ AX).sum():>9.2f}")
    f = ns[-1]
    mag = np.linalg.norm(U[f], axis=1)
    bp = {}
    for i, m_ in enumerate(mag):
        bp.setdefault(part_of.get(i + 1, "?"), []).append(m_)
    print("\n  末帧各部件 |u|（按 max 降序）：")
    for nm, v in sorted(bp.items(), key=lambda kv: -max(kv[1])):
        v = np.array(v)
        print(f"    {nm:28s} n={len(v):>5d} max={v.max():>7.4f} "
              f"中位={np.median(v):>7.4f} 中位/max={np.median(v)/v.max():>5.2f}")
    print(f"\n  末帧沿轴支反力 = {(R[f] @ AX).sum():.2f} N")


if __name__ == "__main__":
    report("baseline", W / "_base.feb", W / "_base.xplt")
    if len(sys.argv) > 1:
        report(sys.argv[1], W / "ankle_tibiotalar.feb",
               sorted(W.glob("ankle_tibiotalar.xplt"), key=lambda p: -p.stat().st_mtime)[0])
