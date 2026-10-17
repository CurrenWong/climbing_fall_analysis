"""二分（修正版）：用「指定位移」驱动（收敛判据正常），隔离 Surface / surface_load。
"""
import subprocess
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = P / "temp/thums/febio/bisect2"
W.mkdir(parents=True, exist_ok=True)

db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
FACES = {4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
         8: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
         6: [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]}


def dedupe(row):
    out = []
    for x in row:
        xi = int(x)
        if xi not in (0, -1) and xi not in out:
            out.append(xi)
    return out


tib_hex = [(8, dedupe(r)) for r in db["81000700_solid"] if len(dedupe(r)) == 8]


def surface_of(els, flip=False):
    from collections import Counter
    cnt, info = Counter(), {}
    for k, conn in els:
        Pp = xyz[[idx[v] for v in conn]]
        cen = Pp.mean(0)
        for f in FACES[k]:
            fv = [conn[z] for z in f]
            pts = xyz[[idx[v] for v in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                fv = fv[::-1]
            if flip:
                fv = fv[::-1]
            key = tuple(sorted(fv))
            info[key] = tuple(fv)
            cnt[key] += 1
    return [info[k] for k, c in cnt.items() if c == 1]


S_OUT = surface_of(tib_hex, flip=False)   # 朝外
S_IN = surface_of(tib_hex, flip=True)     # 朝内（故意反向）

MESH = {}


def build(fname, with_surface, surface_faces, with_load, title):
    used = set()
    for _, c in tib_hex:
        used.update(c)
    if with_surface:
        for f in surface_faces:
            used.update(f)
    ids = sorted(used)
    remap = {v: i + 1 for i, v in enumerate(ids)}
    allc = np.array([xyz[idx[v]] for v in ids])
    a = int(np.argmax(allc.max(0) - allc.min(0)))
    lo, hi = allc[:, a].min(), allc[:, a].max()
    bot = [remap[v] for v, p in zip(ids, allc) if p[a] < lo + 0.1 * (hi - lo)]
    top = [remap[v] for v, p in zip(ids, allc) if p[a] > hi - 0.1 * (hi - lo)]

    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">',
         '  <Module type="solid"/>',
         '  <Globals><Constants><T>298.0</T><R>8.314e-06</R><Fc>9.6485e-05</Fc></Constants></Globals>',
         '  <Material><material name="m" id="1" type="isotropic elastic">'
         '<density>2e-09</density><E>18000</E><v>0.3</v></material></Material>',
         '  <Mesh>', '    <Nodes name="nodes">']
    for v in ids:
        p = xyz[idx[v]]
        L.append(f'      <node id="{remap[v]}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
    L.append('    </Nodes>')
    L.append('    <Elements type="hex8" name="t__hex8">')
    for i, (_, c) in enumerate(tib_hex, 1):
        L.append(f'      <elem id="{i}">{",".join(str(remap[x]) for x in c)}</elem>')
    L.append('    </Elements>')
    L.append(f'    <NodeSet name="bot">{",".join(map(str, sorted(bot)))}</NodeSet>')
    L.append(f'    <NodeSet name="top">{",".join(map(str, sorted(top)))}</NodeSet>')
    if with_surface:
        L.append('    <Surface name="s">')
        for i, f in enumerate(surface_faces, 1):
            kind = "quad4" if len(f) == 4 else "tri3"
            L.append(f'      <{kind} id="{i}">{",".join(str(remap[x]) for x in f)}</{kind}>')
        L.append('    </Surface>')
    L.append('  </Mesh>')
    L.append('  <MeshDomains><SolidDomain name="t__hex8" mat="m"/></MeshDomains>')
    if with_load:
        L.append('  <Loads>')
        L.append('    <surface_load type="traction" surface="s">'
                 '<traction>0.001,0,0</traction></surface_load>')
        L.append('  </Loads>')
    L.append('  <Boundary>')
    L.append('    <bc type="zero displacement" node_set="bot">'
             '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc>')
    L.append('    <bc type="prescribed displacement" node_set="top">'
             f'<dof>{"xyz"[a]}</dof><value>-0.05</value><relative>0</relative></bc>')
    L.append('  </Boundary>')
    L.append('  <Control><analysis>STATIC</analysis><time_steps>2</time_steps>'
             '<step_size>0.5</step_size><solver type="solid">'
             '<linear_solver type="pardiso"/></solver></Control>')
    L.append('</febio_spec>')
    f = W / fname
    f.write_text("\n".join(L), encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=1200)
    out = (r.stdout or "") + (r.stderr or "")
    err = ""
    if "ERROR" in out:
        i = out.find("ERROR")
        err = " | ".join(x.strip(" *") for x in out[i:i + 200].split("\n")[1:3] if x.strip())
    print(f"  {'OK  ' if r.returncode == 0 else 'FAIL'} {title:44s} rc={r.returncode}  {err[:90]}")
    return r.returncode


print("=== 修正版二分（指定位移驱动）===")
build("a_nosurf.feb", False, None, False, "A hex8 + 指定位移（无 Surface）")
build("b_surf_out.feb", True, S_OUT, False, "B A + Surface（法向朝外，正确）")
build("c_surf_in.feb", True, S_IN, False, "C A + Surface（法向朝内，故意反向）")
build("d_surf_load.feb", True, S_OUT, True, "D B + surface_load")
