"""gmsh 4.15 导入 STL 的正确姿势: 先把 STL 写进 model, 再取 native pointer 重建实体。"""

from pathlib import Path

import gmsh

STL = Path(__file__).resolve().parents[2] / "temp" / "ankle_geom" / "talus_r.stl"
CL = 0.006

print(f"STL {STL.name}\n")

# ---- A: importShapesNativePointer ----
print("[A] importShapesNativePointer")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("a")
    # 1) 把 STL 作为网格读进当前 model
    gmsh.open(str(STL))
    gmsh.model.add("a")  # 切回干净 model? 不, open 已换 model
    # 2) 取当前 mesh 的 native pointer
    ptr = gmsh.model.mesh.getNativePointer()
    print("  mesh native ptr:", type(ptr), repr(ptr)[:80])
    # 3) 导成 OCC 实体
    ents = gmsh.model.occ.importShapesNativePointer(ptr)
    gmsh.model.occ.synchronize()
    dims = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
    print("  imported entities:", dims)
    gmsh.finalize()
except Exception as e:
    print("  FAIL:", e)
    try:
        gmsh.finalize()
    except Exception:
        pass

# ---- B: getNativePointer 取 STL 字符串, 再 importShapesNativePointer ----
print("\n[B] getNativePointer 类别探测")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("b")
    gmsh.merge(str(STL))  # merge 会把 STL 变成 dim2 实体
    gmsh.model.occ.synchronize()
    dims = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
    print("  after merge:", dims)
    for d in (2, 3):
        for dim, tag in gmsh.model.getEntities(d)[:3]:
            try:
                cls = gmsh.model.getType(d, tag)
                print(f"    dim={d} tag={tag} type={cls}")
            except Exception:
                pass
    gmsh.finalize()
except Exception as e:
    print("  FAIL:", e)
    try:
        gmsh.finalize()
    except Exception:
        pass

# ---- C: 直接用 gmsh.merge 出的 dim2 实体 + healShapes 造实体 ----
print("\n[C] merge -> healShapes -> 实体 -> tet")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("c")
    gmsh.model.occ.importShapes(str(STL), format="stl", highestDimOnly=False)
    gmsh.finalize()
except Exception as e:
    print("  importShapes FAIL:", e)
    try:
        gmsh.finalize()
    except Exception:
        pass

# ---- D: 用 gmsh.merge 然后 Mesh.Algorithm3D 对闭合壳直接 generate(3) ----
print("\n[D] gmsh.merge + 直接 generate(3) (依赖 STL 闭合)")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("d")
    gmsh.merge(str(STL))
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", CL)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", CL * 0.4)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    # 关键: 有时需要在 surface 上加 volume。尝试 MeshVolume
    try:
        gmsh.model.mesh.generate(3)
    except Exception as e:
        print("  generate(3) 直接失败:", e)
    etypes3, _, _ = gmsh.model.mesh.getElements(3)
    nn, _, _ = gmsh.model.mesh.getNodes()
    print(f"  nodes={len(nn)}  3D elements types={etypes3}")
    gmsh.finalize()
except Exception as e:
    print("  FAIL:", e)
    try:
        gmsh.finalize()
    except Exception:
        pass
