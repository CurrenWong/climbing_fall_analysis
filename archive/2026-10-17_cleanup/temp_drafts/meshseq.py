"""生成**同一表面 × 同一生成器 × 多档尺寸**的 tet 网格序列（干净收敛序列）。

输入表面 = `parts_gmsh/calcaneus_r/calcaneus_r_surface_used.stl`（636 点，**anatframe**，
即现有细化网格 `parts_gmsh/calcaneus_r/calcaneus_r_mesh.npz` 实际用的那个面）。
gmsh 调用与 `probe_gmsh_remesh.py` **逐行同法**（Algorithm 6 / Algorithm3D 1 / Min 0.5cl）。

输出 `temp/opensim_fe/meshseq/calc_cl<cl>.npz`（键：nodes_xyz / tets）
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

SDK = Path(r"D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk")
LIB = SDK / "lib"
os.environ["PATH"] = f"{LIB};{SDK / 'bin'};{os.environ['PATH']}"
sys.path.insert(0, str(LIB))

import numpy as np  # noqa: E402
from climbing.coupling.meshing import _scaled_jacobian, _tet_volumes  # noqa: E402

SURF = ROOT / "temp" / "opensim_fe" / "parts_gmsh" / "calcaneus_r" / "calcaneus_r_surface_used.stl"
OUT = ROOT / "temp" / "opensim_fe" / "meshseq"
OUT.mkdir(parents=True, exist_ok=True)
CLS = [4.0, 3.0, 2.0, 1.5]


def build(cl):
    import gmsh
    gmsh.initialize(str(LIB))
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.option.setNumber("General.Verbosity", 0)
        gmsh.merge(str(SURF))
        gmsh.option.setNumber("Geometry.Tolerance", 1e-6)
        gmsh.model.mesh.removeDuplicateNodes()
        gmsh.model.mesh.classifySurfaces(0.001)
        gmsh.model.mesh.createGeometry()
        gmsh.model.geo.synchronize()
        surfs = [t for _, t in gmsh.model.getEntities(2)]
        gmsh.model.geo.addVolume([gmsh.model.geo.addSurfaceLoop(surfs)])
        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.MeshSizeMin", cl * 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", cl)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.model.mesh.generate(2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        tags, coords, _ = gmsh.model.mesh.getNodes()
        tag_arr = np.array(tags, dtype=np.int64)
        nodes = np.array(coords, dtype=np.float64).reshape(-1, 3)
        _, nidx = gmsh.model.mesh.getElementsByType(4)
        lut = np.zeros(int(tag_arr.max()) + 1, dtype=np.int64)
        lut[tag_arr] = np.arange(len(tag_arr))
        tets = lut[np.array(nidx, dtype=np.int64)].reshape(-1, 4)
        return nodes, tets
    finally:
        gmsh.finalize()


print("同一表面 × 同一生成器 × 多档 cl（anatframe）\n")
print(f"  {'cl(mm)':>7} {'nodes':>7} {'tets':>8} {'minJ':>10} {'fracJ<0.3':>10} "
      f"{'vol(mm³)':>10} {'边长中位':>9}")
rows = []
for cl in CLS:
    try:
        nodes, tets = build(cl)
    except Exception as ex:                            # noqa: BLE001
        print(f"  {cl:>7.1f} FAIL {str(ex)[:80]}", flush=True)
        continue
    jac = _scaled_jacobian(nodes, tets)
    vol = float(_tet_volumes(nodes, tets).sum())
    E = [np.linalg.norm(nodes[tets[:, i]] - nodes[tets[:, j]], axis=1)
         for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]
    med_edge = float(np.median(np.array(E)))
    np.savez(OUT / f"calc_cl{cl:g}.npz", nodes_xyz=nodes, tets=tets)
    print(f"  {cl:>7.1f} {len(nodes):>7d} {len(tets):>8d} {jac.min():>10.3e} "
          f"{float((jac<0.3).mean()):>10.4f} {vol:>10.1f} {med_edge:>9.3f}", flush=True)
    rows.append({"cl_mm": cl, "nodes": int(len(nodes)), "tets": int(len(tets)),
                 "minJ": float(jac.min()), "vol_mm3": vol, "median_edge_mm": med_edge})
import json  # noqa: E402
(OUT / "meshseq_summary.json").write_text(
    json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
print("\n->", OUT)
