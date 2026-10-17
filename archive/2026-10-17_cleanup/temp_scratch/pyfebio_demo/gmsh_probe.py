"""实测: gmsh 4.15 (pip wheel) 到底支持哪种 STL -> tet 路径。"""

from pathlib import Path

import gmsh

STL = Path(__file__).resolve().parents[2] / "temp" / "ankle_geom" / "talus_r.stl"
print(f"STL: {STL}  ({STL.stat().st_size:,} bytes)\n")


def report(tag):
    """打印当前 model 的实体和网格统计。"""
    ents = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
    etypes, _, _ = gmsh.model.mesh.getElements(3)
    tri_types, _, _ = gmsh.model.mesh.getElements(2)
    nn, _, _ = gmsh.model.mesh.getNodes()
    print(f"  {tag:<44} ents(0/1/2/3)={ents[0]}/{ents[1]}/{ents[2]}/{ents[3]}"
          f"  nodes={len(nn)}  3D={etypes}  2D={tri_types}")


# ---- 方案 1: gmsh.open(stl) ----
print("[1] gmsh.open(stl)")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(STL))
    report("gmsh.open")
    print("  -> getEntities(3):", gmsh.model.getEntities(3))
except Exception as e:
    print("  FAIL:", e)
finally:
    try:
        gmsh.finalize()
    except Exception:
        pass

# ---- 方案 2: occ.importShapes 各种 format ----
print("\n[2] occ.importShapes 各种 format")
for fmt in ["stl", "STL", "STL14", "stl14", "auto", ""]:
    try:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("t")
        gmsh.model.occ.importShapes(str(STL), format=fmt)
        gmsh.model.occ.synchronize()
        ents = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
        print(f"  format={fmt!r:<8} OK  ents={ents}")
        gmsh.finalize()
        break
    except Exception as e:
        msg = str(e)[:60]
        print(f"  format={fmt!r:<8} FAIL {msg}")
        try:
            gmsh.finalize()
        except Exception:
            pass

# ---- 方案 3: mesh.importStl + option ----
print("\n[3] mesh.importStl 各种 option 组合")
for opt_val, fname_fmt in [(0, None), (1, str(STL))]:
    try:
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("t")
        if fname_fmt:
            gmsh.option.setNumber("MeshSTLInFile", 1)
            # 看看有没有别的 option 指定文件名
            for o in ("Mesh.ImportFile", "General.FileName", "Mesh.FileName"):
                try:
                    gmsh.option.setNumber(o, 0)
                except Exception:
                    pass
        gmsh.model.mesh.importStl()
        ents = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
        nn, _, _ = gmsh.model.mesh.getNodes()
        print(f"  MeshSTLInFile={opt_val} OK  ents={ents} nodes={len(nn)}")
        gmsh.finalize()
    except Exception as e:
        print(f"  MeshSTLInFile={opt_val} FAIL {str(e)[:70]}")
        try:
            gmsh.finalize()
        except Exception:
            pass

# ---- 方案 4: gmsh.open 后 reclassify/heal ----
print("\n[4] gmsh.open + occ.reclassify/healShapes")
try:
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(STL))
    # 把 mesh 转成几何实体再重建
    for fn, tag in [
        (lambda: gmsh.model.occ.reclassify(1, 1), "reclassify(1,1)"),
        (lambda: gmsh.model.occ.healShapes(gmsh.model.getEntities(2), 1e-8, True), "healShapes(dim2)"),
    ]:
        try:
            fn()
            gmsh.model.occ.synchronize()
            ents = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
            print(f"  {tag:<24} OK ents={ents}")
        except Exception as e:
            print(f"  {tag:<24} FAIL {str(e)[:60]}")
    gmsh.finalize()
except Exception as e:
    print("  FAIL:", e)
    try:
        gmsh.finalize()
    except Exception:
        pass

# ---- 方案 5: 看 gmsh 有没有 stl 相关的 heal/rebuild option ----
print("\n[5] 关键 option 当前值")
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
for o in ["MeshSTLInFile", "MeshBinary", "Geometry.OCCImportLabels"]:
    try:
        print(f"  {o:<28} = {gmsh.option.getNumber(o)}")
    except Exception as e:
        print(f"  {o:<28} <无>")
gmsh.finalize()
