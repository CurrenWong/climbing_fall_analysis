"""修正版重构验证：正解梯度 + tet4 支持 + 排序分布对比。"""
import re
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])
DOMS = {}
for m in re.finditer(r'<Elements(\s[^>]*?)>(.*?)</Elements>', txt, re.S):
    a, body = m.group(1), m.group(2)
    k = re.search(r'type="(\w+)"', a)
    n = re.search(r'name="([^"]*)"', a)
    if not (k and n):
        continue
    DOMS[n.group(1)] = (k.group(1), [[int(x) for x in c.split(",")] for c in
                                     re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*'
                                                r'</elem>', body)])

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    nid = np.asarray(nd["id"], int)
    X = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    sids = sorted(f["states"].keys(), key=int)
    D = np.asarray(f["states"][sids[-1]]["node_data"]["displacement"]["all"][:],
                   float)
    SG = f["states"][sids[-1]]["element_data"]["stress"]
    STORED = {k: np.asarray(SG[k][:], float) for k in SG.keys()}
pos = {int(v): i for i, v in enumerate(nid)}


def dN_hex(xi, et, ze):
    sg = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                   [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)
    dN = np.zeros((8, 3))
    for i in range(8):
        a, b, c = sg[i]
        dN[i] = [a*(1+b*et)*(1+c*ze)/8, b*(1+a*xi)*(1+c*ze)/8,
                 c*(1+a*xi)*(1+b*et)/8]
    return dN


def dN_wedge(r, s, t):
    return np.array([-(1-t), -(1-t), -(1-r-s), (1-t), 0, -r, 0, (1-t), -s,
                     -t, -t, (1-r-s), t, 0, r, 0, t, s]).reshape(6, 3)


dN_TET = np.array([[-1., -1., -1.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
GH = [tuple(g) for g in np.array(
    [[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],
    float)*(1/np.sqrt(3))]
GW = [(1/6, 1/6, .5), (2/3, 1/6, .5), (1/6, 2/3, .5),
      (1/6, 1/6, -.5), (2/3, 1/6, -.5), (1/6, 2/3, -.5)]


def vm(s):
    xx, yy, zz, xy, yz, xz = (s[..., i] for i in range(6))
    return np.sqrt(0.5*((xx-yy)**2+(yy-zz)**2+(zz-xx)**2)
                   + 3.0*(xy**2+yz**2+xz**2))


def recon(kind, conn, E, nu):
    lam, mu = E*nu/((1+nu)*(1-2*nu)), E/(2*(1+nu))
    out = []
    for c in conn:
        try:
            P0 = np.array([XYZ[n] for n in c], float)
            P1 = np.array([X[pos[n]] + D[pos[n]] for n in c], float)
        except KeyError:
            continue
        pts = {"hex8": GH, "penta6": GW, "tet4": [None]}[kind]
        vs = []
        for p in pts:
            dN = (dN_TET if kind == "tet4" else
                  (dN_hex(*p) if kind == "hex8" else dN_wedge(*p)))
            J = P0.T @ dN
            M = (P1 - P0).T @ dN
            g = M @ np.linalg.inv(J)                     # ★ 正解（畸变单元对照过）
            e = 0.5*(g + g.T)
            s = lam*np.trace(e)*np.eye(3) + 2*mu*e
            vs.append([s[0, 0], s[1, 1], s[2, 2], s[0, 1], s[1, 2], s[0, 2]])
        if vs:
            out.append(np.nanmax(vm(np.array(vs))))
    return np.array(out)


print("=" * 108)
print("排序分布对比（与顺序无关）—— 修正梯度后")
print("=" * 108)
print(f"{'域':<26}{'n':>5}{'存储 p50':>10}{'重构 p50':>10}{'存储 p99':>10}"
      f"{'重构 p99':>10}{'存储 max':>10}{'重构 max':>10}{'p99比':>8}")
for nm, (kind, conn) in sorted(DOMS.items()):
    if nm not in STORED:
        continue
    low = nm.lower()
    if any(k in low for k in ("cort", "spon", "marrow")):
        E, nu = 7300.0, 0.3
    elif "tibiotalar" in low:
        E, nu = 10.0, 0.4
    elif low.startswith("lig_"):
        continue
    else:
        continue
    if kind not in ("hex8", "penta6", "tet4"):
        continue
    rec = recon(kind, conn, E, nu)
    if not len(rec):
        continue
    st = vm(STORED[nm])
    n = min(len(rec), len(st))
    a, b = np.sort(st)[-n:], np.sort(rec)[-n:]
    print(f"{nm:<26}{n:>5}{np.median(st):>10.3f}{np.median(rec):>10.3f}"
          f"{np.percentile(st,99):>10.3f}{np.percentile(rec,99):>10.3f}"
          f"{st.max():>10.3f}{rec.max():>10.3f}"
          f"{np.percentile(rec,99)/max(np.percentile(st,99),1e-9):>8.2f}")
