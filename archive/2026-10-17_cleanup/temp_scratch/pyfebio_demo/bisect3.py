"""验证假设：全量模型失败是因为「胫骨只靠接触约束」→ 接触未咬合时系统奇异。

若换成「指定位移」驱动（位移有界，不会因奇异而爆掉）就能跑通，
则证明网格本身是好的，问题在接触咬合。
"""
import subprocess
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = P / "temp/thums/febio/bisect3"
W.mkdir(parents=True, exist_ok=True)

db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(P / "temp/thums/cartilage_mesh.npz")

FACES = {4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
         8: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
         6: [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]}
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]
WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]

BONES = {"tib_c": "81000700_solid", "tib_s": "81000600_solid",
         "tib_m": "81000601_solid", "tal_c": "81001100_solid",
         "tal_s": "81001000_solid"}
CARTS = ["tibiotalar_tibia", "tibiotalar_talus"]


def dedupe(row):
    out = []
    for x in row:
        xi = int(x)
        if xi not in (0, -1) and xi not in out:
            out.append(xi)
    return out


# ── 收集单元 ──
blocks = {}          # name -> [(etype, conn_of_nodeids)]
used_bone = set()
for nm, key in BONES.items():
    for row in db[key]:
        c = dedupe(row)
        if len(c) == 8:
            blocks.setdefault(nm, []).append(("hex8", c))
        elif len(c) == 4:
            blocks.setdefault(nm, []).append(("tet4", c))
        elif len(c) == 6:
            for t in WEDGE_T:
                blocks.setdefault(nm, []).append(("tet4", [c[z] for z in t]))
        used_bone.update(c)

cart = {}
for j in CARTS:
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        full = vals[:m] + [f"off:{j}:{-x - 1}" for x in vals[m:]]
        for t in (HEX_T if m == 4 else WEDGE_T):
            cart.setdefault(j, []).append([full[z] for z in t])
        used_bone.update(vals[:m])
print(f"骨单元 block: { {k: len(v) for k, v in blocks.items()} }")
print(f"软骨 tet 数: { {k: len(v) for k, v in cart.items()} }")

# ── 统一编号 ──
remap = {v: i + 1 for i, v in enumerate(sorted(used_bone))}
coords = [tuple(xyz[idx[v]]) for v in sorted(used_bone)]
off_of = {}
for j in CARTS:
    for k, p in enumerate(cg[f"{j}_offsets"]):
        off_of[f"off:{j}:{k}"] = len(coords) + 1
        coords.append(tuple(p))
coords = np.asarray(coords, float)
NID = lambda x: remap[x] if not isinstance(x, str) else off_of[x]
print(f"合并节点 {len(coords):,}")

# ── 外表面（bone）与软骨外层 ──
def bone_surface(nm):
    from collections import Counter
    cnt, info = Counter(), {}
    for et, conn in blocks[nm]:
        Pp = coords[[NID(v) - 1 for v in conn]]
        cen = Pp.mean(0)
        for f in FACES[len(conn)]:
            fv = [conn[z] for z in f]
            pts = coords[[NID(v) - 1 for v in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                fv = fv[::-1]
            key = tuple(sorted(fv))
            info[key] = tuple(fv)
            cnt[key] += 1
    return [info[k] for k, c in cnt.items() if c == 1]


def cart_outer(j):
    out = []
    for row in cg[f"{j}_elems"]:
        vals = [int(x) for x in row if int(x) != 0]
        if not vals or vals[0] <= 0:
            continue
        m = len(vals) // 2
        out.append([off_of[f"off:{j}:{-x - 1}"] for x in vals[m:]])
    return out


tib_surf = bone_surface("tib_c")
tal_surf = bone_surface("tal_c")
tb = coords[[NID(v) - 1 for v in sorted(remap)]]
axis = (coords[[NID(v) - 1 for v in sorted(remap)]].mean(0))
T = coords[[NID(v) - 1 for v in tib_surf[0]]].mean(0)
Ta = coords[[NID(v) - 1 for v in tal_surf[0]]].mean(0)
ax = (T - Ta)
ax = ax / np.linalg.norm(ax)
pr = coords @ ax
hi, lo = pr.max(), pr.min()
tib_top = [f for f in tib_surf if pr[[NID(v) - 1 for v in f]].mean() > hi - 2.5]
tal_bot = [f for f in tal_surf if pr[[NID(v) - 1 for v in f]].mean() < lo + 2.5]
print(f"加载面 {len(tib_top)} 面；固定面 {len(tal_bot)} 面")

OUTER = {j: cart_outer(j) for j in CARTS}
for j, v in OUTER.items():
    print(f"  {j} 接触面 {len(v)} 面")


def build(fname, use_contact, driver, mag, title):
    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">',
         '  <Module type="solid"/>',
         '  <Globals><Constants><T>298.0</T><R>8.314e-06</R><Fc>9.6485e-05</Fc></Constants></Globals>',
         '  <Material>',
         '    <material name="m_bone" id="1" type="isotropic elastic">'
         '<density>2e-09</density><E>15000</E><v>0.3</v></material>',
         '    <material name="m_cart" id="2" type="isotropic elastic">'
         '<density>1.1e-09</density><E>260</E><v>0.4</v></material>',
         '  </Material>', '  <Mesh>', '    <Nodes name="all">']
    for i, p in enumerate(coords, 1):
        L.append(f'      <node id="{i}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
    L.append('    </Nodes>')
    for nm, els in blocks.items():
        for et in ("hex8", "tet4"):
            sub = [c for e, c in els if e == et]
            if not sub:
                continue
            L.append(f'    <Elements type="{et}" name="{nm}__{et}">')
            for i, c in enumerate(sub, 1):
                L.append(f'      <elem id="{i}">{",".join(str(NID(x)) for x in c)}</elem>')
            L.append('    </Elements>')
    for j, tets in cart.items():
        L.append(f'    <Elements type="tet4" name="{j}__tet4">')
        for i, c in enumerate(tets, 1):
            L.append(f'      <elem id="{i}">{",".join(str(NID(x)) for x in c)}</elem>')
        L.append('    </Elements>')
    L.append(f'    <NodeSet name="tib_top_ns">'
             f'{",".join(str(NID(v)) for v in sorted({x for f in tib_top for x in f}))}'
             f'</NodeSet>')
    L.append(f'    <NodeSet name="tal_bot_ns">'
             f'{",".join(str(NID(v)) for v in sorted({x for f in tal_bot for x in f}))}'
             f'</NodeSet>')
    for nm, fs in (("srfc_tib", OUTER[CARTS[0]]), ("srfc_tal", OUTER[CARTS[1]])):
        L.append(f'    <Surface name="{nm}">')
        for i, f in enumerate(fs, 1):
            kind = "tri3" if len(f) == 3 else "quad4"
            L.append(f'      <{kind} id="{i}">{",".join(str(x) for x in f)}</{kind}>')
        L.append('    </Surface>')
    if use_contact:
        L.append('    <SurfacePair name="pair">'
                 '<primary>srfc_tib</primary><secondary>srfc_tal</secondary></SurfacePair>')
    L.append('  </Mesh>')
    L.append('  <MeshDomains>')
    for nm, els in blocks.items():
        for et in sorted({e for e, _ in els}):
            L.append(f'    <SolidDomain name="{nm}__{et}" mat="m_bone"/>')
    for j in cart:
        L.append(f'    <SolidDomain name="{j}__tet4" mat="m_cart"/>')
    L.append('  </MeshDomains>')
    L.append('  <Boundary>')
    L.append('    <bc type="zero displacement" node_set="tal_bot_ns">'
             '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc>')
    if driver == "disp":
        # 沿「距骨→胫骨」轴分量最大的那个 dof 施加指定位移（bisect2 已验证此写法可跑）
        k = int(np.argmax(np.abs(ax)))
        sgn = -1.0 if ax[k] > 0 else 1.0
        L.append('    <bc type="prescribed displacement" node_set="tib_top_ns">'
                 f'<dof>{"xyz"[k]}</dof><value lc="1">{sgn*mag:.6f}</value>'
                 '<relative>0</relative></bc>')
    L.append('  </Boundary>')
    if use_contact:
        L.append('  <Contact>')
        L.append('    <contact name="ct" surface_pair="pair" type="sliding-elastic">'
                 '<laugon>PENALTY</laugon><penalty>1.0</penalty>'
                 '<auto_penalty>1</auto_penalty><tolerance>0.1</tolerance></contact>')
        L.append('  </Contact>')
    L.append('  <LoadData><load_controller id="1" type="loadcurve">'
             '<interpolate>LINEAR</interpolate><extend>CONSTANT</extend>'
             '<points><pt>0.0,0.0</pt><pt>1.0,1.0</pt></points>'
             '</load_controller></LoadData>')
    L.append('  <Control><analysis>STATIC</analysis><time_steps>10</time_steps>'
             '<step_size>0.1</step_size><solver type="solid">'
             '<linear_solver type="pardiso"/></solver></Control>')
    L.append('</febio_spec>')
    f = W / fname
    f.write_text("\n".join(L), encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                       capture_output=True, text=True, timeout=1800)
    out = (r.stdout or "") + (r.stderr or "")
    err = ""
    if "ERROR" in out:
        i = out.find("ERROR")
        err = " | ".join(x.strip(" *") for x in out[i:i + 200].split("\n")[1:3] if x.strip())
    print(f"  {'OK  ' if r.returncode == 0 else 'FAIL'} {title:46s} rc={r.returncode}  {err[:88]}")
    return r.returncode


print("\n=== 全量模型：指定位移 vs 载荷，用/不用接触 ===")
build("e1_disp_contact.feb", True, "disp", 0.02, "E1 指定位移 + 接触")
build("e2_disp_nocontact.feb", False, "disp", 0.02, "E2 指定位移 + 无接触")
