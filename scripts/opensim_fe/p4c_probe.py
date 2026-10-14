# -*- coding: utf-8 -*-
"""P4-c 对照探针：解除/缩小 lateral_support 验证弯曲路径。

预注册：P4C_PROBE_TODO.md（2026-10-14，冻结）。
复用 p4_pb_probe.load_ctx / head_joint_nodes（p4_pb_probe.py 零改动）。
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402

import bone_feb as bf  # noqa: E402
import paper_alignment as pa  # noqa: E402
from p4_pb_probe import load_ctx, head_joint_nodes  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
BASE_DIR = ROOT / "temp" / "opensim_fe" / "p4c_probe" / "femur_r"
P4A_DIR = ROOT / "temp" / "opensim_fe" / "p4_pb_probe" / "femur_r"

HEIGHTS = [10, 20, 30, 50]           # 预注册 §1.2
SIGMA_C = 220.0
VARIANTS = {
    "p4c_narrow": {"lateral_support": [(0.70, 1.00)], "tag": "约束仅固定端带[0.70,1.00]"},
    "p4c_none": {"lateral_support": None, "tag": "完全解除"},
}


def statistic_p4c(h5: Path, dom: str, load_n: float) -> dict:
    """shaft 域观测：vm median（判定口径）+ 弯曲指标 + 位移哨兵（预注册 §1.3）。

    口径全部复用 P4-a' 诊断脚本（temp/p4_bend_check.py / p4_diag_check.py）：
    - szz_max > 0 = 拉伸侧出现 = 弯曲存在
    - r2_linear = σ_zz 对单元横向坐标 (x,y) 的线性回归 r²（基线 0.291 vs P4-a' 0.054）
    - a_eff = F / |σ_zz median|（≈408 = 纯轴向）
    - dis_z_max = 节点 z 位移最大幅值（屈曲哨兵）
    """
    import h5py
    with h5py.File(h5, "r") as h:
        raw = h["meshes/0/nodes"][:]
        coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
        conn = h["meshes/0/domains"][dom][:]
        sids = sorted(h["states"].keys(), key=int)
        s = h[f"states/{sids[-1]}"]
        stress = s["element_data/stress"][dom][:].astype(float)
        dis = None
        if "node_data" in s and "displacement" in s["node_data"]:
            dd = s["node_data/displacement"]
            for k in dd.keys():
                dis = dd[k][:].astype(float)
                break
    idx = conn[:, 1:] - 1
    cent = coords[idx].mean(axis=1)
    p = pa.partition_ends_shaft(coords, cent)
    m = p["shaft_mask"]
    vm = pa.vm_voigt(stress)
    szz = stress[:, 2]
    cc = cent[m]
    s_sh = szz[m]
    cg = coords.mean(0)
    lat = cc - cg
    A = np.stack([np.ones(len(s_sh)), lat[:, 0], lat[:, 1]], 1)
    beta, *_ = np.linalg.lstsq(A, s_sh, rcond=None)
    resid = s_sh - A @ beta
    r2 = 1 - resid.var() / s_sh.var() if s_sh.var() > 0 else 0.0
    med = float(np.median(s_sh))
    out = {
        "ends_stat_mpa": float(np.median(vm[p["ends_mask"]])),
        "shaft_stat_mpa": float(np.median(vm[m])),
        "vm_max_mpa": float(vm.max()),
        "szz_median_mpa": med,
        "szz_min_mpa": float(s_sh.min()),
        "szz_max_mpa": float(s_sh.max()),
        "r2_linear": float(r2),
        "a_eff_mm2": float(load_n / abs(med)) if med != 0 else float("nan"),
        "n_shaft_elem": int(m.sum()),
    }
    if dis is not None:
        out["dis_z_max_mm"] = float(np.abs(dis[:, 2]).max())
        out["dis_mag_max_mm"] = float(np.linalg.norm(dis, axis=1).max())
    else:
        out["dis_z_max_mm"] = None
        out["dis_mag_max_mm"] = None
    return out


def main() -> int:
    t_all = time.time()
    part, mesh, solid_names, shell_specs, dom, joint, loads, heights = load_ctx()
    print(f"[p4c] {part.bone}: mesh {mesh['n_nodes']} nodes, dom={dom}, joint={joint}",
          flush=True)
    hj = head_joint_nodes(mesh)
    print(f"[p4c] 载荷=头端球冠 n={hj['n_head_nodes']} e_band={hj['band_eccentricity_mm']:.2f}mm, "
          f"固定=髁端 n={hj['n_fix_nodes']}", flush=True)

    # 正对照基准：P4-a' 已收口（不重跑），读 P4_PB_PROBE.json 同高度值
    p4a = json.loads((RES / "P4_PB_PROBE.json").read_text(encoding="utf-8"))
    p4a_by_h = {r["h_no"]: r for r in p4a["rows"] if r["case"] == "P4a_head_load"}
    p4a_h5 = P4A_DIR / "head_load_h50" / "femur_r.hdf5"
    p4a_dis = None
    if p4a_h5.is_file():
        st = statistic_p4c(p4a_h5, dom, float(p4a_by_h[50]["load_n"]))
        p4a_dis = st.get("dis_z_max_mm")

    rows: list[dict] = []
    for vname, vcfg in VARIANTS.items():
        for hno in HEIGHTS:
            i = hno - 1
            load = float(loads[i])
            d = BASE_DIR / f"{vname}_h{hno:02d}"
            d.mkdir(parents=True, exist_ok=True)
            feb = d / "femur_r.feb"
            h5 = feb.with_suffix(".hdf5")
            rec = {"case": vname, "variant": vcfg["tag"], "h_no": hno,
                   "h_m": float(heights[i]), "load_n": load,
                   "lateral_support": vcfg["lateral_support"]}
            if not h5.is_file():
                t = time.time()
                try:
                    bf.build_bone_feb(
                        part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30,
                        axis=None, with_shell=True, time_steps=1, mesh=mesh,
                        align_pca=True,
                        lateral_support=vcfg["lateral_support"],
                        bottom_fix_axial_only=True,
                        load_nodes_override=np.asarray(hj["head_nodes_1based"], dtype=np.int64),
                        fix_nodes_override=np.asarray(hj["fix_nodes_1based"], dtype=np.int64),
                        load_dir_vec=np.asarray(hj["load_dir_vec"], dtype=float),
                    )
                    bf.run_bone_feb(feb, workdir=feb.parent, timeout=1800)
                    if feb.with_suffix(".xplt").is_file():
                        bf.post_process(feb.with_suffix(".xplt"), h5,
                                        solid_names=solid_names,
                                        shell_specs=shell_specs, label=part.bone)
                    rec["converged"] = True
                except Exception as exc:  # noqa: BLE001
                    rec["converged"] = False
                    rec["error"] = f"{type(exc).__name__}: {exc}"
                rec["seconds"] = round(time.time() - t, 1)
            else:
                rec["converged"] = True
                rec["seconds"] = 0.0
            if h5.is_file():
                rec.update(statistic_p4c(h5, dom, load))
            rows.append(rec)
            print(f"[p4c] {vname} h={hno:2d} ({load:9.1f}N) "
                  f"{'OK ' if rec.get('converged') else 'FAIL'} "
                  + (f"shaft={rec['shaft_stat_mpa']:8.2f} szzmax={rec.get('szz_max_mpa', float('nan')):7.1f} "
                     f"r2={rec.get('r2_linear', float('nan')):.3f} "
                     f"dz={rec.get('dis_z_max_mm')}  [{rec.get('seconds')}s]"
                     if "shaft_stat_mpa" in rec else f"({rec.get('error', '')[:70]})"),
                  flush=True)

    # ---- 正对照对比（P4-a'，仅记录） ----
    print(f"\n[ref] P4-a' h50: shaft=83.3 szzmax=0.0 r2=0.054 "
          f"dis_z_max={p4a_dis if p4a_dis is not None else 'n/a'}", flush=True)

    # ---- 判定（预注册 §1.4，跑前写死，跑后如实） ----
    verdict, attribution, sentinel = judge(rows, p4a_dis)
    print("\n=== 判定 ===")
    print(f"  主判据：{verdict}")
    print(f"  归因：{attribution}")
    if sentinel:
        print(f"  屈曲哨兵：{sentinel}")

    out_json = RES / "P4C_PROBE.json"
    out_json.write_text(json.dumps({
        "probe": "P4-c lateral_support 对照（解除/缩小 ⇒ 弯曲能否进入 shaft）",
        "bone": part.bone, "domain": dom, "joint": joint, "sigma_c_mpa": SIGMA_C,
        "head_geom": {k: v for k, v in hj.items() if k not in
                      ("head_nodes_1based", "fix_nodes_1based")},
        "p4a_reference": {str(k): {"shaft_stat_mpa": r["shaft_stat_mpa"],
                                   "szz_max_mpa": r.get("szz_max_mpa")}
                          for k, r in p4a_by_h.items()},
        "criterion": "主判据：任一变体在最高收敛高度 szz_max>0 且 r2≥0.2 ⇒ 归因成立（证真）；"
                     "两变体均否 ⇒ 证伪。屈曲哨兵：发散/位移剧增/σ-F 非线性。",
        "verdict": verdict, "attribution": attribution, "sentinel": sentinel,
        "rows": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[p4c] 总耗时 {time.time()-t_all:.0f}s；写出 {out_json}")
    return 0


def judge(rows: list[dict], p4a_dis50: float | None) -> tuple[str, str, str]:
    """预注册 §1.4 判据（跑前写死）。返回 (verdict, attribution, sentinel)。"""
    sentinels = []
    for vname in VARIANTS:
        vrows = [r for r in rows if r["case"] == vname]
        conv = [r for r in vrows if r.get("converged") and "shaft_stat_mpa" in r]
        failed = [r for r in vrows if not r.get("converged")]
        if failed:
            sentinels.append(f"{vname}: 发散 @h={[r['h_no'] for r in failed]}")
        if not conv:
            continue
        sig_f = [r["shaft_stat_mpa"] / r["load_n"] for r in conv]
        lin_dev = (max(sig_f) - min(sig_f)) / min(sig_f) if sig_f else 0.0
        if len(sig_f) >= 2 and lin_dev > 0.20:
            sentinels.append(f"{vname}: σ/F 偏离线性 {lin_dev*100:.0f}%（屈曲非线性信号）")
        # 位移哨兵：h50 位移 > 3× P4-a' h50（仅高载荷屈曲风险最高）
        h50 = next((r for r in conv if r["h_no"] == 50), None)
        if h50 and h50.get("dis_z_max_mm") and p4a_dis50:
            ratio = h50["dis_z_max_mm"] / p4a_dis50
            if ratio > 3.0:
                sentinels.append(f"{vname}: h50 位移 {h50['dis_z_max_mm']:.3f}mm "
                                 f"= {ratio:.1f}× P4-a'（{p4a_dis50:.3f}mm）剧增")
    sentinel = "；".join(sentinels) if sentinels else ""

    # 主判据：任一变体在其最高收敛高度 证真条件
    proven = []
    for vname in VARIANTS:
        conv = [r for r in rows if r["case"] == vname and r.get("converged")
                and "szz_max_mpa" in r]
        if not conv:
            continue
        hmax_c = max(r["h_no"] for r in conv)
        hc = next(r for r in conv if r["h_no"] == hmax_c)
        if hc["szz_max_mpa"] > 0 and hc["r2_linear"] >= 0.2:
            proven.append((vname, hmax_c, hc["szz_max_mpa"], hc["r2_linear"],
                           hc["shaft_stat_mpa"]))
    if proven:
        (vname, hmax_c, szzm, r2, shaft) = proven[0]
        verdict = (f"证真（主判据）：{vname} 在 h={hmax_c} 弯曲恢复 "
                   f"（szz_max={szzm:.1f}>0、r2={r2:.3f}≥0.2、σ_shaft={shaft:.1f}）"
                   f"⇒ lateral_support 结构性抑制弯曲路径归因成立")
        attribution = ("解除/缩小 lateral_support 后弯曲进入 shaft：拉伸侧出现 + σ_zz 对横向坐标"
                       "线性相关恢复 ⇒ P4-a' 归因（防屈曲 BC 抑制弯曲）获对照证实")
    else:
        verdict = ("证伪（主判据）：两变体最高收敛高度均 szz_max≤0 或 r2<0.2 "
                   "⇒ 解除 lateral_support 后弯曲仍缺失，归因错误（另有原因）")
        attribution = ("lateral_support 不是弯曲缺失的原因；弯曲缺失由其他因素（载荷路径/几何）主导")
    return verdict, attribution, sentinel


if __name__ == "__main__":
    raise SystemExit(main())


