"""独立验证：从位移经"应变→本构"重构骨应力，与 h5 的 stress 输出对比。
⇒ 验证 P1 骨折裕度表的输入是否可信。
"""
import re
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"
txt = FEB.read_text(encoding="latin-1", errors="replace")

# ── 材料（从 .feb 读 E/ν）──
MAT = {}
for m in re.finditer(r'<material id="(\d+)" name="([^"]+)" type="([^"]+)">(.*?)'
                     r'</material>', txt, re.S):
    mid, nm, ty, body = int(m.group(1)), m.group(2), m.group(3), m.group(4)
    E = re.search(r'<E>([-\d.eE+]+)</E>', body)
    V = re.search(r'<v>([-\d.eE+]+)</v>', body)
    MAT[mid] = (nm, ty, float(E.group(1)) if E else None,
                float(V.group(1)) if V else None)
print("=== 材料 ===")
for k, v in sorted(MAT.items()):
    if v[3] is not None:
        print(f"  id={k:<3} {v[0]:<22} E={v[2]:<10.1f} ν={v[3]}")

# ── 域 → 材料 id + 单元（属性顺序不定 ⇒ 逐个解析）──
DOMS = {}
for m in re.finditer(r'<Elements(\s[^>]*?)>(.*?)</Elements>', txt, re.S):
    attrs, body = m.group(1), m.group(2)
    kind = re.search(r'type="(\w+)"', attrs)
    mid = re.search(r'mat="(\d+)"', attrs)
    nm = re.search(r'name="([^"]*)"', attrs)
    if not (kind and nm):
        continue
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm.group(1)] = (kind.group(1), int(mid.group(1)) if mid else None, conn)
print(f"\n域数 {len(DOMS)}: {sorted(DOMS)[:6]} ...")

XYZ = {}
for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                     r'\s*([-\d.eE+]+)', txt):
    XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                     float(m.group(4))])

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


def dN_tet():
    return np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)


def dN_wedge(r, s, t):
    return np.array([-(1-t), -(1-t), -(1-r-s), (1-t), 0, -r, 0, (1-t), -s,
                     -t, -t, (1-r-s), t, 0, r, 0, t, s]).reshape(6, 3)


GAUSS_HEX = [tuple(g) for g in np.array(
    [[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],
    float) * (1/np.sqrt(3))]
GAUSS_WEDGE = [(1/6, 1/6, .5), (2/3, 1/6, .5), (1/6, 2/3, .5),
               (1/6, 1/6, -.5), (2/3, 1/6, -.5), (1/6, 2/3, -.5)]


def vm(s):
    xx, yy, zz, xy, yz, xz = s[..., 0], s[..., 1], s[..., 2], s[..., 3], \
        s[..., 4], s[..., 5]
    return np.sqrt(0.5*((xx-yy)**2 + (yy-zz)**2 + (zz-xx)**2)
                   + 3.0*(xy**2 + yz**2 + xz**2))


print("\n" + "=" * 100)
print("逐域：h5 存储应力 vs 从位移重构（本构 σ = λtr(ε)I + 2μ ε）")
print("=" * 100)
print(f"{'域':<32}{'单元':>6}{'存储 vm_max':>13}{'重构 vm_max':>13}"
      f"{'比值':>8}{'相关系数':>10}")
for nm, (kind, mid, conn) in sorted(DOMS.items()):
    if nm not in STORED:
        continue
    # .feb 的 <Elements> 无 mat= ⇒ 按域名判材料（实测所有骨材料均 E=7300/ν=0.3）
    low = nm.lower()
    if any(k in low for k in ("cort", "spon", "marrow")):
        E, nu = 7300.0, 0.3
    elif "cartilage" in low or "tibiotalar" in low:
        E, nu = 10.0, 0.4
    elif low.startswith("lig_") or "ligament" in low:
        E, nu = 260.0, 0.4
    else:
        continue
    lam = E * nu / ((1+nu)*(1-2*nu))
    mu = E / (2*(1+nu))
    stored = STORED[nm]
    rec = []
    for c in conn:
        try:
            P0 = np.array([XYZ[n] for n in c], float)
            P1 = np.array([X[pos[n]] + D[pos[n]] for n in c], float)
        except KeyError:
            continue
        pts = {"hex8": GAUSS_HEX, "penta6": GAUSS_WEDGE}.get(kind)
        if pts is None:
            continue
        vs = []
        for p in pts:
            dN = (dN_hex(*p) if kind == "hex8" else dN_wedge(*p))
            J = P0.T @ dN
            try:
                g = np.linalg.solve(J, (P1 - P0).T @ dN).T   # (3, npd)
            except np.linalg.LinAlgError:
                vs.append(np.nan); continue
            eps = 0.5 * (g + g.T)
            sig = lam * np.trace(eps) * np.eye(3) + 2*mu*eps
            # Voigt [xx,yy,zz,xy,yz,xz]
            vs.append([sig[0,0], sig[1,1], sig[2,2], sig[0,1], sig[1,2], sig[0,2]])
        if vs:
            rec.append(np.nanmax(vm(np.array(vs))))
    if not rec:
        continue
    rec = np.array(rec)
    n = min(len(rec), len(stored))
    a, b = vm(stored)[:n], rec[:n]
    r = np.corrcoef(a, b)[0, 1] if n > 2 and a.std() > 0 and b.std() > 0 else np.nan
    ratio = rec.max()/vm(stored).max() if vm(stored).max() > 0 else np.nan
    print(f"{nm:<32}{n:>6}{vm(stored).max():>13.3f}{rec.max():>13.3f}"
          f"{ratio:>8.2f}{r:>10.3f}")
