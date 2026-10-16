"""★ S10 收敛研究的**核心读数**：细化 vs 粗化，同口径对比 σ_max。

口径（与粗网格一致）：anat 坐标系 / F=22723.586 N（=h 2 m）/
跖面 `fixed`（三向零位移）/ 皮质域 gauge_max（ℓ=4 mm）+ 裸 max。
粗网格基准（本链条重跑）= 161.402 MPa（`v1_ucurve` 的 fixed h=2）。
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for _p in (str(ROOT / "src"), str(ROOT / "scripts" / "opensim_fe")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import numpy as np  # noqa: E402
import thums_feb as T  # noqa: E402
from s2_thums_sweep import _read_domain_vm  # noqa: E402

REF = ROOT / "temp" / "opensim_fe" / "parts_gmsh" / "calcaneus_r" / "calcaneus_r_mesh.npz"
RUN = ROOT / "temp" / "opensim_fe" / "s10_retry" / "shell1.4_s1_n1" / "shell1.4_s1_n1"
REG = 4.0
d = np.load(REF)
X = np.asarray(d["nodes_xyz"], float)
TET = np.asarray(d["tets"], int)
# 复算分层（与构造器同法）
elements = [TET[i] for i in range(len(TET))]
polys, owners = T._mixed_boundary_faces(elements)
elem_c = np.asarray([X[e].mean(0) for e in elements], float)
bn = np.unique(np.concatenate([np.asarray(p, int) for p in polys]))
try:
    from scipy.spatial import cKDTree
    dmin, _ = cKDTree(X[bn]).query(elem_c)
except Exception:                                   # noqa: BLE001
    dmin = np.linalg.norm(elem_c[:, None, :] - X[bn][None, :, :], axis=2).min(1)
cort = np.flatnonzero(dmin <= 1.40)
spon = np.flatnonzero(dmin > 1.40)
cc = TET[cort]
nodes = X
c0 = nodes[cc].mean(1)
v = np.abs(np.einsum('ij,ij->i',
                     np.cross(nodes[cc[:, 1]] - nodes[cc[:, 0]],
                              nodes[cc[:, 2]] - nodes[cc[:, 0]]),
                     nodes[cc[:, 3]] - nodes[cc[:, 0]])) / 6.0
print(f"皮质 tet {len(cort)}  松质 tet {len(spon)}")
print(f"皮质体积 ΣV = {v.sum():.2f} mm³  单元体积 中位 {np.median(v):.4f}")

vm = _read_domain_vm(Path(str(RUN) + ".xplt"), "calcaneus", Path(str(RUN) + ".hdf5"))
vm = np.asarray(vm, float).ravel()
print(f"\nvm 数组长度 {len(vm)}（应 = {len(cort)}）")
res = {}
if len(vm) == len(cort):
    g = T.gauge_von_mises_mixed(vm, c0, v, REG)
    res = {"refined_gauge_max_mpa": float(g.max()), "refined_max_mpa": float(vm.max()),
           "refined_mean_mpa": float(vm.mean()), "refined_p95_mpa": float(np.percentile(vm, 95)),
           "n_cort_tet": int(len(cort))}
else:
    res = {"refined_max_mpa": float(vm.max()),
           "refined_mean_mpa": float(vm.mean()),
           "refined_p95_mpa": float(np.percentile(vm, 95)),
           "n_cort_tet": int(len(cort)), "n_vm": int(len(vm)),
           "note": "长度不匹配 ⇒ 只给裸 max/mean/p95"}
COARSE = {"coarse_gauge_max_mpa": 161.402, "coarse_max_mpa": 172.167, "note": "粗 1544/634hex8"}
print(json.dumps({**COARSE, **res}, ensure_ascii=False, indent=2))
out = ROOT / "results" / "opensim_fe" / "s10_convergence_readout.json"
out.write_text(json.dumps({"mode": "fixed @ F=22723.586 N (h=2 m), anat",
                           "reg_mm": REG, **COARSE, **res}, ensure_ascii=False, indent=2),
               encoding="utf-8")
print("\n->", out)
