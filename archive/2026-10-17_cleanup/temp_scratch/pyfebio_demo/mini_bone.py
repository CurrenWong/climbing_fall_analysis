"""最小复现：把 THUMS tibia_cort 的前 N 个 hex8 原样放进一个 .feb，看 FEBio 收不收。
用来判定「我的 det(J) 公式是否等价于 FEBio 的判据」。
"""
import subprocess
from pathlib import Path
import numpy as np

NPZ = Path(r"D:\Project\climbing_fall_analysis\temp\thums\thums_lowerlimb.npz")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio\probe")
W.mkdir(parents=True, exist_ok=True)

db = np.load(NPZ)
nid, xyz = db["node_ids"], db["node_xyz"]
idx = {int(v): i for i, v in enumerate(nid)}


def dedupe(row):
    out = []
    for x in row:
        xi = int(x)
        if xi != -1 and xi not in out:
            out.append(xi)
    return out


DOMAIN = "81000700_solid"       # tibia_cort，纯 hex8
N = 200

elems = []
used = set()
for row in db[DOMAIN]:
    c = dedupe(row)
    if len(c) != 8:
        continue
    elems.append(c)
    used.update(c)
    if len(elems) >= N:
        break

ids = sorted(used)
remap = {v: i + 1 for i, v in enumerate(ids)}
print(f"取 {len(elems)} 个 hex8，{len(ids)} 个节点")

# 底面 = 沿 y 最小的 30% 节点（随便选个固定面）
P = np.array([xyz[idx[v]] for v in ids])
axis = 1
lo = P[:, axis].min()
bot = [remap[v] for v, p in zip(ids, P) if p[axis] < lo + 0.15 * (P[:, axis].max() - lo)]
top = [remap[v] for v, p in zip(ids, P) if p[axis] > P[:, axis].max() - 0.15 * (P[:, axis].max() - lo)]
print(f"固定 {len(bot)} 节点, 加载 {len(top)} 节点")

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
L.append('    <Elements type="hex8" name="p1">')
for i, c in enumerate(elems, 1):
    L.append(f'      <elem id="{i}">{",".join(str(remap[x]) for x in c)}</elem>')
L.append('    </Elements>')
L.append(f'    <NodeSet name="bot">{",".join(map(str, sorted(bot)))}</NodeSet>')
L.append(f'    <NodeSet name="top">{",".join(map(str, sorted(top)))}</NodeSet>')
L.append('  </Mesh>')
L.append('  <MeshDomains><SolidDomain name="p1" mat="m"/></MeshDomains>')
L.append('  <Boundary><bc type="zero displacement" node_set="bot">'
         '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc></Boundary>')
L.append('  <Control><analysis>STATIC</analysis><time_steps>2</time_steps>'
         '<step_size>0.5</step_size><solver type="solid">'
         '<linear_solver type="pardiso"/></solver></Control>')
L.append('</febio_spec>')

f = W / "mini_bone.feb"
f.write_text("\n".join(L), encoding="latin-1")
r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(W),
                   capture_output=True, text=True, timeout=900)
out = (r.stdout or "") + (r.stderr or "")
print(f"\nFEBio rc={r.returncode}")
print("--- 关键行 ---")
for l in out.split("\n"):
    if any(k in l for k in ("ERROR", "jacobian", "SUCCESS", "NORMAL", "WARNING",
                            "iter", "Time", "failed", "FAILED", "isolated")):
        print("  " + l.strip())
