"""把软骨接触面 + Contact 段搬进 ankle_lig_feb.py。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① OFFSET_IDS
a1 = "XYZ = np.asarray(coords, float)"
assert s.count(a1) == 1
s = s.replace(a1, a1 + """
# 软骨外表面（offset）节点的 FE 编号集合 —— 用于识别接触面
OFFSET_IDS = {node_id_of[f"off:{j}:{k}"]
              for j in CARTS for k in range(len(cg[f"{j}_offsets"]))}
print(f"其中 offset（软骨外表面）节点 {len(OFFSET_IDS):,}")""")

# ② 接触面抽取 + facet_ok，插在 Jacobian 修正之后、写 Mesh 之前
a2 = "w('  <Mesh>')"
assert s.count(a2) == 1
port = """# ── 软骨接触面（★ 必须取 tet 的**真实外表面**，不能用 npz 原始四边形）──
# 这两条都是实测换来的：
#   ① 接触面与网格实际面不重合 ⇒ 接触**完全不检测**，但 rc=0、无警告、只是应力全 0
#   ② 反向：面必须朝外定向（法向背离对顶点）
TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
MIN_FACET_AREA = 1e-8


def facet_ok(f):
    if len(set(f)) != len(f):
        return False
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return np.linalg.norm(n) / 2.0 >= MIN_FACET_AREA


def cart_outer(j):
    out, seen = [], set()
    for ids in blocks[(j, "tet4")]:
        for f in TET_FACES:
            fv = [ids[z] for z in f]
            if not all(x in OFFSET_IDS for x in fv):
                continue
            key = tuple(sorted(fv))
            if key in seen or not facet_ok(fv):
                continue
            seen.add(key)
            opp = ids[[z for z in range(4) if z not in f][0]]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(3):
                n += np.cross(pts[i], pts[(i + 1) % 3])
            if np.dot(n, pts.mean(0) - XYZ[opp - 1]) < 0:
                fv = fv[::-1]
            out.append(fv)
    return out


contact_surfs = {j: cart_outer(j) for j in CARTS}
for j, v in contact_surfs.items():
    print(f"  接触面 {j}: {len(v)} 面")

w('  <Mesh>')"""
s = s.replace(a2, port)

# ③ Mesh 里加接触 Surface + SurfacePair
a3 = "w('  </Mesh>')"
assert s.count(a3) == 1
s = s.replace(a3, """for j, sel in contact_surfs.items():
    w(f'    <Surface name="surf_{j}">')
    for n, f in enumerate(sel, 1):
        w(f'      <tri3 id="{n}">{",".join(str(x) for x in f)}</tri3>')
    w('    </Surface>')
w('    <SurfacePair name="tibiotalar_pair">')
w('      <primary>surf_tibiotalar_tibia</primary>')
w('      <secondary>surf_tibiotalar_talus</secondary>')
w('    </SurfacePair>')
w('  </Mesh>')""")

# ④ Contact 段（段序：Loads -> Boundary -> Contact）
a4 = "w('  <Output>')"
assert s.count(a4) == 1
s = s.replace(a4, """w('  <Contact>')
w('    <contact name="tibiotalar" surface_pair="tibiotalar_pair" '
  'type="sliding-elastic">')
w('      <laugon>PENALTY</laugon>')
w('      <penalty>1</penalty>')
w('      <auto_penalty>1</auto_penalty>')
w('      <two_pass>0</two_pass>')
w('      <node_reloc>0</node_reloc>')
w('      <symmetric_stiffness>0</symmetric_stiffness>')
w('      <tolerance>0.02</tolerance>')
w('    </contact>')
w('  </Contact>')
w('  <Output>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
