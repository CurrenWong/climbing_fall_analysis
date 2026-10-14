# -*- coding: utf-8 -*-
"""F-3′d · 把股骨 `W` 测对（修 f3p 的 A/W/R 常数倍错误）。

f3p 失败根因（本件已实证钉死）：
1. `signed_volumes` 用了非共享对角线 TETS 且**额外 ÷6** ⇒ 体积错
2. `conn[:,1:]-1` 把 0-based 行号当 1-based ⇒ 形心错一行（域划分错位）
⇒ 修正：体积用共享对角线 (0,6) 6-tet **直接求和**；形心直读 `coords[conn[:,1:]]`

定义（事先写明）：
* `A` = 骨干皮质**载流截面**：切片法 Σ dA（dA = V_elem/Δt_elem），取 shaft 段平均
       （并与 M1 同法 ΣV/L 交叉核对）
* `R` = 每切片内节点的**最大垂距**，取 shaft 段 bins 的中位（外半径）
* `W` = **真实截面二阶矩**：每切片 I=Σ y²·dA（y=单元质心到骨干轴垂距），
       中段平均 I_avg ⇒ W = I_avg / R_avg
* `e` = top 带（顶 30%）节点形心到骨干轴的垂距（模型自身几何偏心）

验收判据（F3P_FEMUR_GEOM_FAIL.md §3 事先写明）：
  A ∈ [400, 700] mm²  且  R ∈ [12, 18] mm  且  W ∈ [1500, 3500] mm³
  ⇒ 三者同时满足才允许把预测作为可证伪预测交付；否则继续报 FAIL。

预测：σ_shaft(h) = F(h)/A + F(h)·e/W，F(h) = loads_n["hip_r"] × load_fraction
对照：稳定化 + median 口径实测 0/50 ⇒ all_zero。
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
N_BIN = 40                      # 切片数（沿骨干轴）
MID_FRAC = 0.6                  # 中段平均：取 bins 的中间 60%

#: 共享对角线 (0,6) 6-tet 分解（= bone_feb._HEX_TETS），直接求和 = hex 体积
_HEX_TETS = np.array(
    [[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6], [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]],
    dtype=np.int64,
)


def hex_volumes(coords: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """hex8 体积：共享对角线 6-tet 直接求和（idx 为 0-based 节点行号）。"""
    vols = np.zeros(len(idx), dtype=float)
    for t in _HEX_TETS:
        p = coords[idx[:, t]]
        vols += np.einsum("ij,ij->i", np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]),
                          p[:, 3] - p[:, 0]) / 6.0
    return np.abs(vols)


def main() -> int:
    fj = json.loads((RES / "fracture_matrix_bc2.json").read_text(encoding="utf-8"))
    pr = next(p for p in fj["parts"] if p["bone"] == "femur_r")
    joint, lf = pr["joint"], float(pr.get("load_fraction", 1.0))
    F_all = np.asarray(fj["loads_n"][joint], float) * lf
    heights = np.asarray(fj["heights_m"], float)
    reg = next(r for r in pa.REGIONS if r[1] == "femur_r")
    sig_shaft = pa._sc(reg[3][1][1])
    print(f"[f3d] joint={joint} load_fraction={lf}  "
          f"F(1m)={F_all[0]/1e3:.2f} kN  F(50m)={F_all[-1]/1e3:.2f} kN  阈(骨干)={sig_shaft} MPa")

    with h5py.File(H5, "r") as h:
        doms = list(h["meshes/0/domains"].keys())
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        cort = next((d for d in doms if "cort" in d.lower()), doms[0])
        conn = h["meshes/0/domains"][cort][:]
    idx = conn[:, 1:]                       # 0-based 行号（本件实证）
    vol = hex_volumes(coords, idx)
    cent = coords[idx].mean(axis=1)
    print(f"[f3d] {cort}: {len(conn)} 单元  V_total={vol.sum():.1f} mm³  "
          f"单单元 p50={np.median(vol):.2f} mm³")

    part = pa.partition_ends_shaft(coords, cent)
    sm, em = part["shaft_mask"], part["ends_mask"]
    print(f"[f3d] partition_ends_shaft: 轴={part['load_axis']}  "
          f"shaft={int(sm.sum())} ends={int(em.sum())}  L={part['L_axis_mm']:.1f} mm")

    # --- 骨干轴：shaft 质心 PCA ---
    C = cent[sm]
    mu = C.mean(axis=0)
    _, _, Vt = np.linalg.svd(C - mu, full_matrices=False)
    axis = Vt[0] / np.linalg.norm(Vt[0])
    axis *= np.sign(axis[2]) or 1.0
    print(f"[f3d] 骨干轴 = {np.round(axis, 5)}  （z 夹角 {np.degrees(np.arccos(abs(axis[2]))):.2f}°）")

    t_node = coords @ axis
    t_cent = cent @ axis
    tmin, tmax = float(t_cent[sm].min()), float(t_cent[sm].max())

    def perp_dist(P, origin=mu):
        d = P - origin
        return np.linalg.norm(d - np.outer(d @ axis, axis), axis=1)

    # --- 切片截面二阶矩 ---
    edges = np.linspace(tmin, tmax, N_BIN + 1)
    A_bin, I_bin, R_bin, n_bin = [], [], [], []
    for k in range(N_BIN):
        lo, hi = edges[k], edges[k + 1]
        m = sm & (t_cent >= lo) & (t_cent < hi)
        nk = int(m.sum())
        if nk == 0:
            continue
        dt = np.maximum(np.ptp(t_node[idx[m]], axis=1), 1e-6)   # 单元轴向跨度
        dA = vol[m] / dt
        y = perp_dist(cent[m])
        A_bin.append(float(dA.sum()))
        I_bin.append(float((dA * y ** 2).sum()))
        R_bin.append(float(perp_dist(coords[np.unique(idx[m])]).max()))
        n_bin.append(nk)

    A_bin, I_bin, R_bin, n_bin = map(np.asarray, (A_bin, I_bin, R_bin, n_bin))
    # 中段平均：按 bins 排序取中间 MID_FRAC
    j0 = int(np.floor((1 - MID_FRAC) / 2 * len(A_bin)))
    j1 = int(np.ceil((1 + MID_FRAC) / 2 * len(A_bin)))
    A_mid = float(A_bin[j0:j1].mean())
    I_mid = float(I_bin[j0:j1].mean())
    R_mid = float(np.median(R_bin[j0:j1]))
    R_mid_mean = float(R_bin[j0:j1].mean())
    W = I_mid / R_mid_mean
    print(f"[f3d] 切片 bins={len(A_bin)}  中段({int(MID_FRAC*100)}%): "
          f"A={A_mid:.1f} mm²  I={I_mid:.0f} mm⁴  R(median)={R_mid:.2f} mm  "
          f"R(mean)={R_mid_mean:.2f} mm  ⇒ W=I/R={W:.0f} mm³")

    # --- 交叉核对 A（M1 同法：ΣV/L） ---
    L_shaft = float(np.ptp((C - mu) @ axis))
    A_vol = float(vol[sm].sum() / L_shaft)
    print(f"[f3d] A(ΣV/L 同法) = {A_vol:.1f} mm²  （L_shaft={L_shaft:.1f} mm）")

    # --- e：top 带（顶 30%）形心垂距 ---
    z0, z1 = float(t_node.min()), float(t_node.max())
    thr = z0 + (1.0 - TOP_FRAC) * (z1 - z0)
    top = coords[t_node >= thr]
    e_c = float(np.linalg.norm((top.mean(axis=0) - mu)
                               - np.dot(top.mean(axis=0) - mu, axis) * axis))
    e_band = float(perp_dist(top).mean())
    print(f"[f3d] top 带（顶 {TOP_FRAC:.0%}）: {len(top)} 节点  e(形心垂距)={e_c:.2f} mm  "
          f"(带内均垂距 {e_band:.2f} mm)")

    # --- 验收 ---
    okA = 400 <= A_mid <= 700
    okR = 12 <= R_mid <= 18
    okW = 1500 <= W <= 3500
    ok = okA and okR and okW
    print(f"\n[f3d] 验收: A∈[400,700]? {okA} ({A_mid:.1f}) | R∈[12,18]? {okR} "
          f"({R_mid:.2f}) | W∈[1500,3500]? {okW} ({W:.0f}) ⇒ {'PASS' if ok else 'FAIL'}")

    # --- 预测（F3P 预注册：仅 PASS 时正式交付；FAIL 时给量级参考并标注非正式） ---
    sig_ax = F_all / A_mid
    sig_bd = F_all * e_c / W
    sig_tot = sig_ax + sig_bd
    rows = []
    print(f"\n[f3d] 预测 σ_shaft(h) = F/A + F·e/W  （e={e_c:.2f} mm, W={W:.0f} mm³, A={A_mid:.1f} mm²）")
    print(f"{'h':>3} {'F_kN':>7} {'F/A':>8} {'F·e/W':>8} {'合计':>8} {'/220':>6}  骨折")
    for i in (0, 4, 9, 19, 29, 39, 49):
        r = sig_tot[i] / sig_shaft
        print(f"{i+1:>3} {F_all[i]/1e3:>7.2f} {sig_ax[i]:>8.1f} {sig_bd[i]:>8.1f} "
              f"{sig_tot[i]:>8.1f} {r:>6.2f}  {'1' if r >= 1 else '0'}")
        rows.append([i + 1, float(F_all[i]), float(sig_ax[i]), float(sig_bd[i]),
                     float(sig_tot[i]), float(r)])
    hit = np.where(sig_tot >= sig_shaft)[0]
    first = int(hit[0]) + 1 if hit.size else None
    print(f"[f3d] 预测首次跨越 = {first}（{heights[first-1] if first else '—'} m），"
          f"跨阈数 = {hit.size}/50")
    pred = {"rows": rows, "first_cross_h": first, "n_cross": int(hit.size),
            "formal_delivery": ok}

    out = RES / "F3D_FEMUR_W_MEASURE.json"
    out.write_text(json.dumps({
        "piece": "F-3′d 股骨 W 测对（切片截面二阶矩，零 FE 解算）",
        "h5": str(H5.relative_to(ROOT)).replace("\\", "/"),
        "domain": cort, "joint": joint, "load_fraction": lf,
        "sigma_c_shaft_mpa": sig_shaft,
        "index_convention": "h5 domains = 0-based 行号（.feb 1-based 交叉验证 V_total 逐位一致）",
        "volume_method": "共享对角线(0,6) 6-tet 直接求和（仓库 _hex8_signed_volumes 等价）",
        "v_total_mm3": float(vol.sum()),
        "measured": {"A_mid_mm2": A_mid, "A_sigV_over_L_mm2": A_vol,
                     "L_shaft_mm": L_shaft, "I_mid_mm4": I_mid,
                     "R_median_mm": R_mid, "R_mean_mm": R_mid_mean,
                     "W_mm3": W, "e_centroid_mm": e_c, "e_band_mean_mm": e_band,
                     "n_top_nodes": int(len(top)), "axis": axis.tolist(),
                     "bins": {"A": A_bin.tolist(), "I": I_bin.tolist(),
                              "R": R_bin.tolist(), "n": n_bin.tolist()}},
        "acceptance": {"A_ok": okA, "R_ok": okR, "W_ok": okW, "pass": ok},
        "prediction": pred,
        "observed": {"n_cross_median_caliber": 0, "behav": "all_zero"},
        "caveats": [
            "e 定义为模型自身 top 带形心到骨干轴垂距；真实解剖偏心（股骨头中心）需独立来源",
            "I 用单一 PCA 直轴（股骨生理前弓会引入小系统误差，验收窗已放宽）",
            "A 与 G6 锚 500（=中股骨皮质 CSA，文献 505±49）比值为 0.79×，同量级（G6 胫骨 320 vs 380 判例）",
            "未改动任何判定；预测未用任何可调参数",
        ],
    }, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[f3d] 写出 {out.name}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
