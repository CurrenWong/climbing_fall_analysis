"""找到把 gmsh 'Discrete surface' 转成 3D 实体的方法。"""

from pathlib import Path

import gmsh

STL = Path(__file__).resolve().parents[2] / "temp" / "ankle_geom" / "talus_r.stl"
CL = 0.006


def try_tet(tag, setup):
    """setup(gmsh) 之后尝试四面体网格化, 报告结果。"""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.model.add("m")
        setup()
        gmsh.model.occ.synchronize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", CL)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", CL * 0.4)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        et3, tags3, _ = gmsh.model.mesh.getElements(3)
        nn, _, _ = gmsh.model.mesh.getNodes()
        if et3 and 4 in et3:
            cnt = len([t for et, tg in zip(et3, tags3) if et == 4 for _ in tg])
            print(f"  [OK] {tag:<40} tets={cnt:,}  nodes={len(nn):,}")
        else:
            dims = {d: len(gmsh.model.getEntities(d)) for d in (0, 1, 2, 3)}
            print(f"  [--] {tag:<40} 无四面体  ents={dims}")
    except Exception as e:
        print(f"  [XX] {tag:<40} {str(e)[:70]}")
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


print("尝试把 Discrete surface 转实体:\n")

# 1) removeAllDuplicates: 把离散几何的重复点合并, 有时能触发实体识别
try_tet(
    "merge + removeAllDuplicates",
    lambda: (
        gmsh.merge(str(STL)),
        gmsh.model.occ.removeAllDuplicates(),
    ),
)

# 2) healShapes 在 Discrete 上
def s2():
    gmsh.merge(str(STL))
    gmsh.model.occ.healShapes([(2, t) for d, t in gmsh.model.getEntities(2)], 1e-6, True)


try_tet("merge + healShapes", s2)

# 3) geo 模块的 addVolumeFromEntities? 查 geo 有无
print("\ngeo 模块里 volume 相关:")
print("  ", [a for a in dir(gmsh.model.geo) if "olume" in a])

# 4) 关键: gmsh 有 "Mesh.MeshOnly" 与 "GeneralizeToSolid"?
for opt in ["Mesh.MeshOnly", "Mesh.SaveAll", "Geometry.OCCBooleanPreserveNumbering"]:
    gmsh.initialize()
    try:
        v = gmsh.option.getNumber(opt)
        print(f"  option {opt} = {v}")
    except Exception:
        print(f"  option {opt} = <无>")
    gmsh.finalize()

# 5) 最强候选: gmsh.model.occ.fragment / sewing
def s5():
    gmsh.merge(str(STL))
    # sewing 把缝合的面片变实体
    gmsh.model.occ.fragment(
        [(2, t) for d, t in gmsh.model.getEntities(2)],
        [],
    )


try_tet("merge + fragment", s5)

# 6) kernel 层
print("\nkernel 模块接口:")
gmsh.initialize()
print("  ", [a for a in dir(gmsh.model.getCurrentEntityName) if "ol" in a.lower()][:5] if hasattr(gmsh.model, "getCurrentEntityName") else "n/a")
gmsh.finalize()

# 7) 看 gmsh.kernel 是否暴露 STL->Brep
gmsh.initialize()
print("\ngmsh.kernel attrs (含 brep/stl/heal):")
print("  ", [a for a in dir(gmsh.kernel) if any(k in a.lower() for k in ("stl", "brep", "heal", "shell", "solid"))])
gmsh.finalize()
