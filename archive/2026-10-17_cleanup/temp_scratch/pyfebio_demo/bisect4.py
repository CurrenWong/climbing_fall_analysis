"""最小复现：只放软骨 tet4，看 FEBio 是否在 mesh init 就报负 Jacobian。
同时逐单元报告几何质量，找出坏单元。
"""
import subprocess
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = P / "temp/thums/febio/bisect4"
W.mkdir(parents=True, exist_ok=True)

db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(P / "temp/thums/cartilage_mesh.npz")
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
CARTS = ["tibiotalar_tibia", "tibiotalar_talus"]

# 统一编号：骨节点 + 各关节 offset
bone_ids = sorted({int(x) for j in CARTS for row in cg[f"{j}_elems"]
                   for x in row if int(x) > 0})
remap = {v: i + 1 for i, v in enumerate(bone_ids)}
coords = [tuple(xyz[idx[v]]) for v in bone_ids]
off_of = {}
for j in CARTS:
    for k, p in enumerate(cg[f"{j}_offsets"]):
        off_of[f"off:{j}:{k}"] = len(coords) + 1
        coords.append(tuple(p))
coords = np.asarray(coords, float)
NID = lambda x: remap[x] if not isinstance(x, str) else off_of[x]

tets = {}     # j -> [(tet_nodeids, parent_cell_index)]
for j in CARTS:
    tets[j] = []
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        full = vals[:m] + [f"off:{j}:{-x - 1}" for x in vals[m:]]
        for t in (HEX_T if m == 4 else WEDGE_T):
            tets[j].append([full[z] for z in t])

print(f"软骨 tet: { {k: len(v) for k, v in tets.items()} }  "
      f"节点 {len(coords):,}")

# ── 逐 tet 几何质量 ──
print("\n=== 软骨 tet 质量（按 |6V| 与形状）===")
for j, ts in tets.items():
    vols, mind, dup = [], [], 0
    for t in ts:
        ids = [NID(x) for x in t]
        if len(set(ids)) < 4:
            dup += 1
            continue
        Pp = coords[[i - 1 for i in ids]]
        v = float(np.dot(np.cross(Pp[1] - Pp[0], Pp[2] - Pp[0]), Pp[3] - Pp[0])) / 6.0
        vols.append(v)
        d = np.linalg.norm(Pp[:, None, :] - Pp[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        mind.append(d.min())
    vols = np.array(vols)
    mind = np.array(mind)
    print(f"  {j}: n={len(vols)}  体积 min={vols.min():.6f} p1={np.percentile(vols,1):.6f} "
          f"中位={np.median(vols):.6f}  **负体积 {int((vols<0).sum())}**  "
          f"|V|<1e-6 的 {int((np.abs(vols)<1e-6).sum())}  重复节点 {dup}")
    print(f"      最小节点间距 min={mind.min():.6f} p1={np.percentile(mind,1):.6f} "
          f"中位={np.median(mind):.6f}  间距<1e-3 的 {int((mind<1e-3).sum())}")


def build(fname, which, title):
    use = [j for j in which]
    used = set()
    for j in use:
        for t in tets[j]:
            used.update(NID(x) for x in t)
    ids = sorted(used)
    r2 = {v: i + 1 for i, v in enumerate(ids)}
    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">',
         '  <Module type="solid"/>',
         '  <Globals><Constants><T>298.0</T><R>8.314e-06</R><Fc>9.6485e-05</Fc></Constants></Globals>',
         '  <Material><material name="m" id="1" type="isotropic elastic">'
         '<density>1.1e-09</density><E>260</E><v>0.4</v></material></Material>',
         '  <Mesh>', '    <Nodes name="all">']
    for v in ids:
        p = coords[v - 1]
        L.append(f'      <node id="{r2[v]}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
    L.append('    </Nodes>')
    for j in use:
        L.append(f'    <Elements type="tet4" name="{j}__tet4">')
        for i, t in enumerate(tets[j], 1):
            L.append(f'      <elem id="{i}">{",".join(str(r2[NID(x)]) for x in t)}</elem>')
        L.append('    </Elements>')
    p0 = np.array([coords[v - 1] for v in ids])
    a = int(np.argmax(p0.max(0) - p0.min(0)))
    lo, hi = p0[:, a].min(), p0[:, a].max()
    bot = [r2[v] for v in ids if coords[v - 1][a] < lo + 0.1 * (hi - lo)]
    L.append(f'    <NodeSet name="bot">{",".join(map(str, sorted(bot)))}</NodeSet>')
    L.append('  </Mesh>')
    L.append('  <MeshDomains>')
    for j in use:
        L.append(f'    <SolidDomain name="{j}__tet4" mat="m"/>')
    L.append('  </MeshDomains>')
    L.append('  <Boundary><bc type="zero displacement" node_set="bot">'
             '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc></Boundary>')
    L.append('  <Control><analysis>STATIC</analysis><time_steps>1</time_steps>'
             '<step_size>1.0</step_size><solver type="solid">'
             '<linear_solver type="pardiso"/></solver></Control>')
    L.append('</febio_spec>')
    f = W / fname
    f.write_text("\n".join(L), encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=900)
    out = (r.stdout or "") + (r.stderr or "")
    err = ""
    if "ERROR" in out:
        i = out.find("ERROR")
        err = " | ".join(x.strip(" *") for x in out[i:i + 220].split("\n")[1:3] if x.strip())
    print(f"  {'OK  ' if r.returncode == 0 else 'FAIL'} {title:40s} rc={r.returncode}  {err[:95]}")
    return r.returncode


print("\n=== 只放软骨 tet4，逐一测试 ===")
build("c1_tib_cart.feb", [CARTS[0]], "C1 仅 tibiotalar_tibia 软骨")
build("c2_tal_cart.feb", [CARTS[1]], "C2 仅 tibiotalar_talus 软骨")
build("c3_both.feb", CARTS, "C3 两侧软骨")
