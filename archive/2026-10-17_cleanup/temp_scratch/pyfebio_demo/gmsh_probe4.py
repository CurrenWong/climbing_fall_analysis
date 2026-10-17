"""用 gmsh 的 mesh.createTopology/createGeometry 把 Discrete surface 转成实体后四面体化。"""

from pathlib import Path

import gmsh

STL = Path(__file__).resolve().parents[2] / "temp" / "ankle_geom" / "talus_r.stl"
CL = 0.006

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("t")
gmsh.merge(str(STL))
gmsh.model.occ.synchronize()
print(f"merge 后 ents: {  {d: len(gmsh.model.getEntities(d)) for d in (0,1,2,3)} }")
print(f"dim2 type: {gmsh.model.getType(2, gmsh.model.getEntities(2)[0][1])}")

# 1) createTopology: 建立几何拓扑 -> 可能产生 dim3
print("\n[1] mesh.createTopology()")
try:
    gmsh.model.mesh.createTopology()
    print(f"  ents: {  {d: len(gmsh.model.getEntities(d)) for d in (0,1,2,3)} }")
except Exception as e:
    print("  FAIL:", e)

# 2) createGeometry: 从网格创建几何实体
print("\n[2] mesh.createGeometry(dim=3)")
try:
    gmsh.model.mesh.createGeometry(3)
    print(f"  ents: {  {d: len(gmsh.model.getEntities(d)) for d in (0,1,2,3)} }")
    for d in (2, 3):
        for dim, tag in gmsh.model.getEntities(d)[:3]:
            print(f"    dim={d} tag={tag} type={gmsh.model.getType(d, tag)}")
except Exception as e:
    print("  FAIL:", e)

# 3) 再加 repair / removeAllDuplicates
print("\n[3] + heal/healShapes 组合")
for tag, fn in [
    ("removeAllDuplicates", lambda: gmsh.model.occ.removeAllDuplicates()),
    ("healShapes(all)", lambda: gmsh.model.occ.healShapes(gmsh.model.getEntities(), 1e-8, True)),
    ("removeDuplicateNodes", lambda: gmsh.model.mesh.removeDuplicateNodes()),
]:
    try:
        fn()
        gmsh.model.occ.synchronize()
        print(f"  {tag:<24} ents={ {d: len(gmsh.model.getEntities(d)) for d in (0,1,2,3)} }")
    except Exception as e:
        print(f"  {tag:<24} FAIL {str(e)[:60]}")

# 4) 最终尝试四面体化
print("\n[4] generate(3)")
try:
    gmsh.option.setNumber("Mesh.MeshSizeMax", CL)
    gmsh.option.setNumber("Mesh.MeshSizeMin", CL * 0.4)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)
    etypes, tags, _ = gmsh.model.mesh.getElements(3)
    nn, _, _ = gmsh.model.mesh.getNodes()
    n3 = 0
    for et, tg in zip(etypes, tags):
        if et == 4:
            n3 += len(tg)
    print(f"  nodes={len(nn):,}  3D types={etypes}  tets={n3:,}")
except Exception as e:
    print("  FAIL:", e)

gmsh.finalize()
