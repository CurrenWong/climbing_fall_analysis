# -*- coding: utf-8 -*-
"""F-3′ · 从**已有 h5 的已对齐几何**量出股骨的 `A` / `W` / `e`，给出可证伪的首次跨越预测。

**零新增 FE 解算**：几何取自 `temp/fe_out/fe_bc_stab/femur_r/A2_axial_lat_90_h50/femur_r.hdf5`
（`align_pca=True` ⇒ 载荷轴 = z；几何与坠落高度无关）。

**定义（事先写明）**
* `A`   = 骨干皮质**载流截面** = `Σ V_elem(骨干皮质) / L_shaft`（M1 同法：几何量，非全长增厚）
* `R`   = 骨干皮质节点到骨干轴的**最大垂距**（外半径）
* `W`   = **环形等效**截面模量：由 `A` 与 `R` 反解内径 `d`，`I = π(D⁴−d⁴)/64`，`W = I/(D/2)`
* **`e`** = **载荷带形心到骨干轴的垂距**（`top_frac=0.30`，沿 z 的顶 30%）
  —— **这是模型自身的几何偏心，不是查来的文献值**
* 预测：`σ_shaft(h) ≈ F(h)/A + F(h)·e/W`，`F(h) = loads_n["hip_r"](h) × load_fraction`

**纪律**：本件**不改任何判定**、**不调任何参数**；若预测与实测不符，**如实报不符**。
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import h5py  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
H5 = ROOT / "temp" / "opensim_fe" / "fe_bc_stab" / "femur_r" / "A2_axial_lat_90_h50" / "femur_r.hdf5"
TOP_FRAC = 0.30


def signed_volumes(coords: np.ndarray, conn: np.ndarray) -> np.ndarray:
    """逐单元有符号体积（hex8 / tet4 通用，1-based 局部索引）。"""
    idx = conn[:, 1:] - 1
    npe = idx.shape[1]
    p = coords[idx]                       # (n, npe, 3)
    if npe == 8:                          # hex8：拆 6 个 tet（与 _hex8_signed_volumes 同族）
        TETS = [(0, 1, 2, 5), (0, 2, 3, 7), (0, 5, 7, 4),
                (5, 6, 7, 2), (0, 5, 2, 7), (0, 7, 2, 3)]
        v = np.zeros(len(idx), dtype=float)
        for a, b, c, d in TETS:
            A4 = p[:, a], p[:, b], p[:, c], p[:, d]
            v += np.abs(np.einsum("ij,ij->i", np.cross(A4[1] - A4[0], A4[2] - A4[0]),
                                  A4[3] - A4[0]) / 6.0)
        return v / 6.0                    # 6 个 tet 平均 ⇒ 与单元体积同量级
    if npe == 4:
        return np.abs(np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                                p[:, 3] - p[:, 0]) / 6.0)
    raise SystemExit(f"未支持的单元类型 npe={npe}")


def main() -> int:
    fj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    pr = next(p for p in fj["parts"] if p["bone"] == "femur_r")
    joint, lf = pr["joint"], float(pr.get("load_fraction", 1.0))
    F_all = np.asarray(fj["loads_n"][joint], float) * lf
    heights = np.asarray(fj["heights_m"], float)
    reg = next(r for r in pa.REGIONS if r[1] == "femur_r")
    sig_neck, sig_shaft = pa._sc(reg[3][0][1]), pa._sc(reg[3][1][1])
    print(f"[f3p] joint={joint} load_fraction={lf}  "
          f"F(1m)={F_all[0]/1e3:.2f} kN  F(50m)={F_all[-1]/1e3:.2f} kN")
    print(f"[f3p] 阈值: 股骨颈={sig_neck} 骨干={sig_shaft} MPa")

    with h5py.File(H5, "r") as h:
        doms = list(h["meshes/0/domains"].keys())
        print(f"[f3p] domains={doms}")
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        cort = next((d for d in doms if "cort" in d.lower()), doms[0])
        conn = h["meshes/0/domains"][cort][:]
        sid = sorted(h["states"].keys(), key=int)[-1]
        stress = np.asarray(h[f"states/{sid}/element_data/stress/{cort}"][:], float)

    vol = signed_volumes(coords, conn)
    cent = coords[conn[:, 1:] - 1].mean(axis=1)
    part = pa.partition_ends_shaft(coords, cent)
    sm, em = part["shaft_mask"], part["ends_mask"]
    print(f"[f3p] {cort}: {conn.shape[0]} 单元  shaft={int(sm.sum())}  ends={int(em.sum())}")

    # --- 骨干轴（PCA on shaft centroids） ---
    C = cent[sm]
    mu = C.mean(axis=0)
    _, _, Vt = np.linalg.svd(C - mu, full_matrices=False)
    axis = Vt[0] / np.linalg.norm(Vt[0])
    axis *= np.sign(np.dot(axis, np.array([0.0, 0.0, 1.0]))) or 1.0
    print(f"[f3p] 骨干轴 = {np.round(axis, 5)}   （z 轴夹角 "
          f"{np.degrees(np.arccos(abs(axis[2]))):.2f}°）")

    def perp_dist(P):
        d = P - mu
        return np.linalg.norm(d - np.outer(d @ axis, axis), axis=1)

    # --- A：骨干皮质载流截面 ---
    L_shaft = float(np.ptp((C - mu) @ axis))
    A = float(vol[sm].sum() / L_shaft)
    print(f"[f3p] 骨干长 L={L_shaft:.2f} mm   A = ΣV/L = {A:.1f} mm²")

    # --- R 与环形等效 W ---
    snodes = np.unique(conn[sm][:, 1:] - 1)
    R = float(perp_dist(coords[snodes]).max())
    D = 2 * R
    inner2 = D ** 2 - 4 * A / np.pi
    d_in = float(np.sqrt(inner2)) if inner2 > 0 else 0.0
    I = np.pi * (D ** 4 - d_in ** 4) / 64.0
    W = I / (D / 2.0)
    print(f"[f3p] 外半径 R={R:.2f} mm (D={D:.2f})  反解内径 d={d_in:.2f} mm  "
          f"⇒ I={I:.0f} mm⁴  W={W:.0f} mm³")

    # --- e：载荷带形心到骨干轴的垂距 ---
    z = coords @ axis
    z0, z1 = float(z.min()), float(z.max())
    thr = z0 + (1.0 - TOP_FRAC) * (z1 - z0)
    top = coords[z >= thr]
    e = float(perp_dist(top).mean())
    e_c = float(np.linalg.norm((top.mean(axis=0) - mu)
                               - np.dot(top.mean(axis=0) - mu, axis) * axis))
    print(f"[f3p] top 带（顶 {TOP_FRAC:.0%}）: {len(top)} 节点  "
          f"**e(形心垂距) = {e_c:.2f} mm**  (带内平均垂距 {e:.2f} mm)")

    # --- 预测 ---
    sig_ax = F_all / A
    sig_bd = F_all * e_c / W
    sig_tot = sig_ax + sig_bd
    print(f"\n[f3p] 预测 σ_shaft(h) = F/A + F·e/W  （e={e_c:.2f} mm, W={W:.0f} mm³, A={A:.1f} mm²）")
    print(f"{'h':>3} {'F_kN':>7} {'F/A':>8} {'F·e/W':>8} {'合计':>8} {'/220':>6}  骨折")
    rows = []
    for i in (0, 4, 9, 19, 29, 39, 49):
        r = sig_tot[i] / sig_shaft
        print(f"{i+1:>3} {F_all[i]/1e3:>7.2f} {sig_ax[i]:>8.1f} {sig_bd[i]:>8.1f} "
              f"{sig_tot[i]:>8.1f} {r:>6.2f}  {'1' if r >= 1 else '0'}")
        rows.append([i + 1, F_all[i], sig_ax[i], sig_bd[i], sig_tot[i], r])
    hit = np.where(sig_tot >= sig_shaft)[0]
    first = int(hit[0]) + 1 if hit.size else None
    print(f"\n[f3p] 预测首次跨越 = {first}（{heights[first-1] if first else '—'} m），"
          f"跨 阈数 = {hit.size}/50")
    print(f"[f3p] 对照实测（稳定化 + median 口径）：0/50 ⇒ all_zero")

    out = RES / "F3P_FEMUR_GEOM_PREDICTION.json"
    out.write_text(json.dumps({
        "piece": "F-3′ 股骨几何 + 载荷偏心实测与预测（零 FE）",
        "h5": str(H5.relative_to(ROOT)).replace("\\", "/"),
        "domain": cort, "joint": joint, "load_fraction": lf,
        "sigma_c_neck_mpa": sig_neck, "sigma_c_shaft_mpa": sig_shaft,
        "measured": {"A_mm2": A, "L_shaft_mm": L_shaft, "R_mm": R, "D_mm": D,
                     "d_inner_mm": d_in, "I_mm4": float(I), "W_mm3": float(W),
                     "e_centroid_mm": e_c, "e_band_mean_mm": e,
                     "n_top_nodes": int(len(top)), "axis": axis.tolist()},
        "prediction": {"rows": rows, "first_cross_h": first,
                       "n_cross": int(hit.size)},
        "observed": {"n_cross_median_caliber": 0, "behav": "all_zero"},
        "caveats": [
            "W 用环形等效（由 A 与外半径 R 反解内径），非真实截面二阶矩",
            "e 定义为 top 带形心到骨干轴的垂距；真实解剖偏心是股骨头中心到骨干轴",
            "未改动任何判定；预测未用任何可调参数",
        ],
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[f3p] 写出 {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())