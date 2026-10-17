"""二分定位：从「只有 tibia_cort hex8」逐步加到全量，隔离负 Jacobian 的引入点。

为排除"收敛困难"的干扰，所有变体都只用「固定底面 + 无载荷」，
这样 FEBio 只做一次装配+线性求解 —— 若网格合法就必然 rc=0。
"""
import subprocess
from pathlib import Path
import numpy as np

P = Path(r"D:\Project\climbing_fall_analysis")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = P / "temp/thums/febio/bisect"
W.mkdir(parents=True, exist_ok=True)

db = np.load(P / "temp/thums/thums_lowerlimb.npz")
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}
cg = np.load(P / "temp/thums/cartilage_mesh.npz")

FACES = {4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
         8: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)],
         6: [(0, 2, 1), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]}
FIELD = {4: "tet4", 8: "hex8", 6: "penta6"}


def dedupe(row):
    out = []
    for x in row:
        xi = int(x)
        if xi not in (0, -1) and xi not in out:
            out.append(xi)
    return out


def bone_elems(key, allow):
    """allow: 允许的节点数集合，如 {8} / {8,4}"""
    out = []
    for row in db[key]:
        c = dedupe(row)
        if len(c) in allow:
            out.append((FIELD[len(c)], c))
    return out


def build(fname, blocks, surfaces, loads=None, title=""):
    """blocks: [(name, [(etype, conn)，...])]  conn 为 node id
       surfaces: [(name, [(kind, conn)])]"""
    used = set()
    for _, els in blocks:
        for _, c in els:
            used.update(c)
    for _, fs in surfaces:
        for _, c in fs:
            used.update(c)
    ids = sorted(used)
    remap = {v: i + 1 for i, v in enumerate(ids)}
    L = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<febio_spec version="4.0">',
         '  <Module type="solid"/>',
         '  <Globals><Constants><T>298.0</T><R>8.314e-06</R><Fc>9.6485e-05</Fc></Constants></Globals>',
         '  <Material><material name="m" id="1" type="isotropic elastic">'
         '<density>2e-09</density><E>18000</E><v>0.3</v></material></Material>',
         '  <Mesh>', '    <Nodes name="nodes">']
    for v in ids:
        p = xyz[idx[v]] if isinstance(v, int) else None
        L.append(f'      <node id="{remap[v]}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
    L.append('    </Nodes>')
    for nm, els in blocks:
        byt = {}
        for et, c in els:
            byt.setdefault(et, []).append(c)
        for et, lst in byt.items():
            L.append(f'    <Elements type="{et}" name="{nm}__{et}">')
            for i, c in enumerate(lst, 1):
                L.append(f'      <elem id="{i}">{",".join(str(remap[x]) for x in c)}</elem>')
            L.append('    </Elements>')
    # 固定面：用第一个 block 的最低点
    allc = np.array([xyz[idx[v]] for v in ids])
    a = int(np.argmax(allc.max(0) - allc.min(0)))
    lo = allc[:, a].min()
    bot = [remap[v] for v, p in zip(ids, allc)
           if p[a] < lo + 0.15 * (allc[:, a].max() - lo)]
    L.append(f'    <NodeSet name="bot">{",".join(map(str, sorted(bot)))}</NodeSet>')
    for nm, fs in surfaces:
        L.append(f'    <Surface name="{nm}">')
        for kind, c in fs:
            L.append(f'      <{kind} id="1">{",".join(str(remap[x]) for x in c)}</{kind}>')
        L.append('    </Surface>')
    L.append('  </Mesh>')
    L.append('  <MeshDomains>')
    for nm, els in blocks:
        byt = sorted({et for et, _ in els})
        for et in byt:
            L.append(f'    <SolidDomain name="{nm}__{et}" mat="m"/>')
    L.append('  </MeshDomains>')
    if loads:
        L.append('  <Loads>')
        L.append(f'    <surface_load type="traction" surface="{loads}">'
                 '<traction>0,0,-1</traction></surface_load>')
        L.append('  </Loads>')
    L.append('  <Boundary><bc type="zero displacement" node_set="bot">'
             '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc></Boundary>')
    L.append('  <Control><analysis>STATIC</analysis><time_steps>1</time_steps>'
             '<step_size>1.0</step_size><solver type="solid">'
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
    print(f"  {'OK  ' if r.returncode == 0 else 'FAIL'} {title:46s} rc={r.returncode}  {err[:95]}")
    return r.returncode


# 表面集合（用 tibia_top 那批面）
tib_hex = bone_elems("81000700_solid", {8})
tib_all = bone_elems("81000700_solid", {8, 6})
spon = bone_elems("81000600_solid", {4})
marrow = bone_elems("81000601_solid", {4})
tal_h = bone_elems("81001100_solid", {8})
tal_s = bone_elems("81001000_solid", {4})


def surface_of(els):
    """从 hex8 单元集提取外表面（带朝外定向）"""
    from collections import Counter
    cnt, info = Counter(), {}
    for et, conn in els:
        P_ = xyz[[idx[v] for v in conn]]
        cen = P_.mean(0)
        for f in FACES[len(conn)]:
            fv = [conn[z] for z in f]
            pts = xyz[[idx[v] for v in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                n = -n
                fv = fv[::-1]
            k = tuple(sorted(fv))
            info[k] = tuple(fv)
            cnt[k] += 1
    outer = [info[k] for k, c in cnt.items() if c == 1]
    return [("quad4", f) for f in outer if len(f) == 4] + \
           [("tri3", f) for f in outer if len(f) == 3]


S = surface_of(tib_hex)
print(f"tibia_cort 外表面 {len(S)} 面（quad4 {sum(1 for k,_ in S if k=='quad4')}）\n")
print("=== 二分结果（固定底面 + 无载荷，只看网格是否合法）===")
build("v1_hex8_only.feb", [("t", tib_hex)], [], None, "V1 只有 tibia_cort hex8 (2415)")
build("v2_hex8_surface.feb", [("t", tib_hex)], [("s", S)], None, "V2 V1 + Surface(quad4/tri3)")
build("v3_hex8_surf_load.feb", [("t", tib_hex)], [("s", S)], "s",
      "V3 V2 + surface_load(traction)")
build("v4_hex8_tet4.feb", [("t", tib_all)], [], None, "V4 tibia_cort hex8+tet4")
build("v5_tibia_all.feb", [("tib_c", tib_hex), ("tib_s", spon), ("tib_m", marrow)],
      [], None, "V5 胫骨全部(cort+spon+marrow)")
build("v6_full_nocart.feb",
      [("tib_c", tib_hex), ("tib_s", spon), ("tib_m", marrow),
       ("tal_c", tal_h), ("tal_s", tal_s)], [], None, "V6 胫骨+距骨(无软骨)")
