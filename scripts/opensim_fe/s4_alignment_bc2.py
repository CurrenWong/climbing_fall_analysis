"""s4_alignment_bc2.py -- B1.4: S4 stats + A1-A4 verdicts on the BC2 matrix.

B1.1 decision (6): Route 甲 (1D risk) is run through the same pipeline as a
code-path cross-reference; OFFICIAL A1-A4 verdicts come from Route 乙 = the
B1.3 BC2 matrix (`fracture_matrix_bc2.json`, primary stat = median).
A3/A4 accept ONLY route 乙 (route 甲 has complete separation -> no OR).

Routes
------
甲-1D  : exactly paper_alignment.build()'s Route-1D (risk_1d_bending.json
        sigma_total(h) vs weak-subregion sigma_c) -- inputs unchanged.
乙-BC2 : sub-region remap (REGIONS / weak sigma_c from [Y25] Table 1) on the
        BC2 products:
          * long bones (tibia/fibula/femur): ends/shaft partition read from
            temp/opensim_fe/fracture_matrix_bc2/<b>/<b>.hdf5, sub-region
            statistic = SUB_STAT (primary median -- B1.1 decision 3; p95 run
            as sensitivity), sigma(h) = stat * load(h)/fe_solve_n;
          * other parts: whole-domain sigma_vm_mpa[SUB_STAT][bone] from the
            bc2 JSON (already load-scaled, strict linear elasticity).

Verdicts (docs/methods/OpenSim_FE复现方案.md section 6, verbatim)
------------------------------------------------------------------
  A1  cluster/phase boundaries vs 7 m / 16 m, |delta| < 2 m both
  A2  stepwise targets {足, 胫骨两端, 股骨颈, 脊柱, 颈椎} >= 4/5 hit
      (hit = logistic status complete_separation; 脊柱 hit = 腰椎 AND 胸椎
       both hit (strict); loose OR-count also reported)
  A3  progressive ORs {腓骨 1.682, 骨盆 1.236, 颅骨 1.576}: >= 2/3 within +-30%
  A4  dual-path PCA separable: covariance of per-part stress-vs-height,
      top-2 PCs; group1 {胫,股,脊柱(腰/胸/颈)} vs group2 {腓,骨盆,颅} must be
      linearly separable in (PC1, PC2) loading space (LP feasibility).
      足部 belongs to neither §六 group -> unlabeled, excluded from the test.

Outputs (new files only):
  results/opensim_fe/s4_alignment_bc2.json
  results/opensim_fe/S4_ALIGNMENT_BC2_REPORT.md

Run (repo root, PYTHONPATH=src):
  & .venv\\Scripts\\python.exe scripts\\opensim_fe\\s4_alignment_bc2.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import paper_alignment as pa  # noqa: E402  (reuse: REGIONS, logistic_fit, ...)

RES = ROOT / "results" / "opensim_fe"
BC2_JSON = RES / "fracture_matrix_bc2.json"
BC2_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
BC2_REPORT = RES / "S4_ALIGNMENT_BC2_REPORT.md"
BC2_OUT = RES / "s4_alignment_bc2.json"
B1_LOG = RES / "B1_LOG.md"

# §六 A3 OR targets (per-metre OR, [Y25])
A3_TARGETS = {"腓骨": 1.682, "骨盆": 1.236, "颅骨": 1.576}
# §六 A2 stepwise targets (5 items); 脊柱 = 腰椎+胸椎 (strict: both)
A2_TARGETS_STRICT = [["足部"], ["胫骨"], ["股骨"], ["腰椎", "胸椎"], ["颈椎"]]
# §六 A4 groups
A4_G1 = ["胫骨", "股骨", "腰椎", "胸椎", "颈椎"]      # {胫/股/脊柱}
A4_G2 = ["腓骨", "骨盆", "颅骨"]                      # {腓/骨盆/颅}
A1_TARGETS = (7.0, 16.0)
A1_TOL = 2.0


def log_line(msg: str) -> None:
    with B1_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{time.strftime('%Y-%m-%d %H:%M')}] {msg}\n")
    print(f"[B1_LOG] {msg}", flush=True)


def _first(heights, mask):
    return pa._first_fracture(heights, mask)


def three_phases_seq(labels, matrix: np.ndarray) -> list[dict]:
    """3-phase reconstruction from ordered clusters -- SEQUENTIAL split.

    Spec correction (documented): pa.three_phases interleaves
    (``order[i::n_grp]``), which for 5 non-identical clusters yields
    OVERLAPPING phases (measured 乙: P1 1-12 / P2 5-50 / P3 6-6) and
    contradicts its own docstring 「等分相邻类」.  This version splits the
    min-height-ordered cluster list into 3 adjacent chunks (np.array_split),
    producing contiguous phases.  pa.three_phases is NOT modified.
    """
    labels = np.asarray(labels)
    h = np.arange(1, len(labels) + 1)
    order = sorted(set(int(v) for v in labels), key=lambda c: h[labels == c].min())
    groups = [g for g in np.array_split(np.array(order, dtype=int), 3) if len(g)]
    out = []
    for i, g in enumerate(groups, 1):
        mask = np.isin(labels, g)
        hs = h[mask]
        fr = [pa.REGION_NAMES[j] for j in range(matrix.shape[1])
              if matrix[mask, j].any()]
        out.append({"phase": i, "label": pa._phase_label(fr),
                    "range_m": [int(hs.min()), int(hs.max())],
                    "heights": [int(x) for x in hs],
                    "fractured_regions": fr,
                    "clusters": [int(x) for x in g]})
    return out


# ---------------------------------------------------------------- route 甲
def route_jia() -> dict:
    """Route-1D, verbatim from paper_alignment.build() (inputs unchanged)."""
    bending = json.loads((RES / "risk_1d_bending.json").read_text(encoding="utf-8"))
    rows = {r["bone"]: r for r in bending["rows_primary_beta4"]}
    sigma = {b: np.asarray(r["risk_by_height"], float) * float(r["sigma_c_mpa"])
             for b, r in rows.items()}
    heights = np.arange(1, 51, dtype=float)
    m1d, first, logistic = {}, {}, {}
    for cn, bone, joint, subs in pa.REGIONS:
        weak = pa._sc(subs[0][1])
        m1d[cn] = sigma[bone] >= weak
        first[cn] = _first(heights, m1d[cn])
        logistic[cn] = pa.logistic_fit(heights, m1d[cn].astype(int))
    M = np.array([[int(m1d[cn][i]) for cn in pa.REGION_NAMES]
                  for i in range(50)], dtype=int)
    clus = pa.jaccard_clusters(M, 5)
    cr = pa.cluster_ranges(clus["labels"])
    ph = three_phases_seq(clus["labels"], M)
    return dict(matrix=M, first=first, logistic=logistic,
                cluster_ranges=cr, phases=ph, matrix_by_part=m1d)


# ---------------------------------------------------------------- route 乙
def route_yi(sub_stat: str) -> dict:
    fm = json.loads(BC2_JSON.read_text(encoding="utf-8"))
    heights = np.asarray(fm["heights_m"], dtype=float)
    loads = {k: np.asarray(v, float) for k, v in fm["loads_n"].items()}
    parts = {p["bone"]: p for p in fm["parts"]}
    r1 = json.loads((RES / "bc_robust_metric_route1.json").read_text(encoding="utf-8"))
    r1_rows = {r["bone"]: r for r in r1["rows"]}

    m, first, logistic, detail = {}, {}, {}, {}
    for cn, bone, joint, subs in pa.REGIONS:
        weak = pa._sc(subs[0][1])
        if bone in pa.LONG_BONES:
            h5 = BC2_DIR / bone / f"{bone}.hdf5"
            if not h5.is_file():
                raise SystemExit(f"missing bc2 hdf5: {h5}")
            dom = r1_rows[bone]["solid_dom"]
            # pa.read_h5_domain bakes FM_DIR (legacy root) -> read bc2 root here:
            import h5py
            with h5py.File(h5, "r") as h:
                raw = h["meshes/0/nodes"][:]
                coords = np.stack([raw["x"], raw["y"], raw["z"]], 1).astype(float)
                conn = h["meshes/0/domains"][dom][:]
                sids = sorted(h["states"].keys(), key=int)
                stress = h[f"states/{sids[-1]}/element_data/stress/{dom}"][:].astype(float)
            idx = conn[:, 1:] - 1
            cent = coords[idx].mean(axis=1)
            part = pa.partition_ends_shaft(coords, cent)
            vm = pa.vm_voigt(stress)
            sel = np.median if sub_stat == "median" else (
                lambda a: np.percentile(a, 95))
            stat_ends = float(sel(vm[part["ends_mask"]]))
            stat_shaft = float(sel(vm[part["shaft_mask"]]))
            fe_solve = float(parts[bone]["fe_solve_load_n"])
            scale = loads[joint] / fe_solve
            ends_h = stat_ends * scale
            shaft_h = stat_shaft * scale
            sig_e, sig_s = pa._sc(subs[0][1]), pa._sc(subs[1][1])
            m[cn] = (ends_h >= sig_e) | (shaft_h >= sig_s)
            detail[cn] = {"domain": dom, "sub_stat": sub_stat,
                          "ends_stat_mpa": stat_ends, "shaft_stat_mpa": stat_shaft,
                          "fe_solve_n": fe_solve,
                          "ends_sigma_c": sig_e, "shaft_sigma_c": sig_s,
                          "first_ends": _first(heights, ends_h >= sig_e),
                          "first_shaft": _first(heights, shaft_h >= sig_s),
                          "n_ends": int(part["n_ends"]),
                          "load_axis": part["load_axis"]}
        else:
            vm_h = np.asarray(fm["sigma_vm_mpa"][sub_stat][bone], float)
            m[cn] = vm_h >= weak
            detail[cn] = {"domain": parts[bone]["main_domain"],
                          "sub_stat": sub_stat, "whole_sigma_c": weak,
                          "sigma_h_at_ref": float(vm_h[4])}
        first[cn] = _first(heights, m[cn])
        logistic[cn] = pa.logistic_fit(heights, m[cn].astype(int))

    M = np.array([[int(m[cn][i]) for cn in pa.REGION_NAMES] for i in range(50)],
                 dtype=int)
    clus = pa.jaccard_clusters(M, 5)
    cr = pa.cluster_ranges(clus["labels"])
    ph = three_phases_seq(clus["labels"], M)
    return dict(matrix=M, first=first, logistic=logistic, cluster_ranges=cr,
                phases=ph, detail=detail, matrix_by_part=m,
                sigma_c_weak={cn: pa._sc(pa.REGIONS[i][3][0][1])
                              for i, cn in enumerate(pa.REGION_NAMES)})


# ---------------------------------------------------------------- verdicts
def a1_verdict(phases: list[dict]) -> dict:
    b1 = float(phases[0]["range_m"][1]) if len(phases) >= 1 else None
    b2 = float(phases[1]["range_m"][1]) if len(phases) >= 2 else None
    d1 = abs(b1 - A1_TARGETS[0]) if b1 is not None else None
    d2 = abs(b2 - A1_TARGETS[1]) if b2 is not None else None
    ok = (d1 is not None and d2 is not None
          and d1 < A1_TOL and d2 < A1_TOL)
    return dict(targets_m=list(A1_TARGETS), tol_m=A1_TOL,
                ours_boundaries_m=[b1, b2], abs_delta=[d1, d2],
                verdict="PASS" if ok else "FAIL", phases=phases)


def _behav(logi: dict) -> str:
    s = logi["status"]
    return {"all_zero": "常值0", "all_one": "常值1",
            "complete_separation": "阶跃", "ok": "渐进"}.get(s, s)


def a2_verdict(logistic: dict) -> dict:
    hits, per = [], {}
    for item in A2_TARGETS_STRICT:
        bs = [_behav(logistic[cn]) for cn in item]
        hit = all(b == "阶跃" for b in bs)
        loose = any(b == "阶跃" for b in bs)
        hits.append(hit)
        per["+".join(item)] = dict(parts=item, behavior=bs, hit_strict=hit,
                                   hit_loose=loose)
    n = sum(hits)
    n_loose = sum(per["+".join(it)]["hit_loose"] for it in A2_TARGETS_STRICT)
    return dict(targets=A2_TARGETS_STRICT, hits=hits, n_hit=n,
                n_target=len(A2_TARGETS_STRICT), n_hit_loose=n_loose,
                per_target=per, verdict="PASS" if n >= 4 else "FAIL")


def a3_verdict(logistic: dict) -> dict:
    per, n_ok = {}, 0
    for cn, tau in A3_TARGETS.items():
        lg = logistic[cn]
        if lg.get("status") == "ok" and lg.get("OR"):
            rel = abs(lg["OR"] / tau - 1.0)
            ok = rel <= 0.30
            per[cn] = dict(target=tau, ours=float(lg["OR"]), rel_err=rel,
                           ci95=lg.get("ci95"), within=bool(ok))
            n_ok += int(ok)
        else:
            per[cn] = dict(target=tau, ours=None, status=lg.get("status"),
                           within=False)
    return dict(targets=A3_TARGETS, per_part=per, n_within=n_ok,
                n_target=3, verdict="PASS" if n_ok >= 2 else "FAIL")


def a4_from_matrix(yi: dict, heights: np.ndarray, fm: dict) -> dict:
    from scipy.optimize import linprog
    # Build sigma(h) per part (MPa), the same stress that feeds the matrix.
    sig = {}
    for cn, bone, joint, subs in pa.REGIONS:
        d = yi["detail"][cn]
        if bone in pa.LONG_BONES:
            # ends trajectory (driving for weak threshold); need load(h):
            loads = np.asarray(fm["loads_n"][joint], float)
            scale = loads / d["fe_solve_n"]
            sig[cn] = d["ends_stat_mpa"] * scale
        else:
            sig[cn] = np.asarray(fm["sigma_vm_mpa"]["median"][bone], float)
    names = pa.REGION_NAMES
    X = np.column_stack([sig[cn] for cn in names])          # 50 x 9
    Xc = X - X.mean(axis=0, keepdims=True)
    C = np.cov(Xc, rowvar=False)                             # 9 x 9
    w, V = np.linalg.eigh(C)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    loads2 = V[:, :2]                                        # (9,2) part loadings
    pos = {cn: [float(loads2[i, 0]), float(loads2[i, 1])]
           for i, cn in enumerate(names)}
    g1 = np.array([pos[c] for c in A4_G1])
    g2 = np.array([pos[c] for c in A4_G2])
    # LP: w.x + b >= 1 (g1); w.x + b <= -1 (g2)
    n1, n2 = len(g1), len(g2)
    A_ub, b_ub = [], []
    for x in g1:
        A_ub.append([-x[0], -x[1], -1.0]); b_ub.append(-1.0)
    for x in g2:
        A_ub.append([x[0], x[1], 1.0]); b_ub.append(-1.0)
    r = linprog(c=[0, 0, 0], A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                bounds=[(None, None)] * 3, method="highs")
    return dict(method="cov(per-part sigma(h)) -> eigh -> top-2 PC loadings; "
                       "linear separability via LP (scipy.linprog highs)",
                explained_var=[float(w[0]), float(w[1])],
                explained_var_ratio=[float(w[0] / w.sum()),
                                     float(w[1] / w.sum())],
                part_positions=pos, group1=A4_G1, group2=A4_G2,
                unlabeled=["足部"],
                separable=bool(r.success),
                verdict="PASS" if r.success else "FAIL",
                lp_message=r.message)


# ---------------------------------------------------------------- report
def md_report(jia: dict, yi: dict, fm: dict, vs: dict, elapsed: float) -> str:
    v1, v2, v3, v4 = vs["v1"], vs["v2"], vs["v3"], vs["v4"]
    v1_jia, v2_jia, v3_jia = vs["v1_jia"], vs["v2_jia"], vs["v3_jia"]
    v1_p95, v2_p95, v3_p95 = vs["v1_p95"], vs["v2_p95"], vs["v3_p95"]
    L: list[str] = []
    A = L.append
    A("# S4 对标（BC2 矩阵）— A1–A4 判定")
    A("")
    A(f"> B1.4 产出（{time.strftime('%Y-%m-%d %H:%M')}）。正式判定 = **Route 乙**"
      "（`fracture_matrix_bc2.json`，`align_pca=True`，sub-stat **median**）；"
      "Route 甲（1D risk）为交叉参照（B1.1 决策⑥）。判据原文 = "
      "`docs/methods/OpenSim_FE复现方案.md` §六。所有数字实测。")
    A("")
    A("## 0. 判定总表（Route 乙 = 正式）")
    A("")
    A("| 代号 | 判据 | 目标 | **乙（median）** | 乙(p95 敏感性) | 甲(1D) |")
    A("|---|---|---|---|---|---|")
    A(f"| A1 | 聚类边界 vs 7/16 m | \\|Δ\\|<2 m | **{v1['verdict']}** "
      f"（边界 {v1['ours_boundaries_m']}，Δ={v1['abs_delta']}） | "
      f"{md_b(v1_p95)} | {md_b(v1_jia)} |")
    A(f"| A2 | 阶跃部位集合 | ≥4/5 | **{v2['verdict']}**（{v2['n_hit']}/5"
      f"，宽松 {v2['n_hit_loose']}） | {md_b(v2_p95)} | {md_b(v2_jia)} |")
    A(f"| A3 | 渐进 OR ±30% | ≥2/3 | **{v3['verdict']}**（{v3['n_within']}/3） | "
      f"{md_b(v3_p95)} | {md_b(v3_jia)}（预期结构性失败） |")
    A(f"| A4 | 双路径 PCA 可分 | 可分 | **{v4['verdict']}** | — | —（甲无 FE 应力轨迹） |")
    A("| A5 | 首折高度（绝对） | 参考 5–12 m | 参考 | — | — |")
    A("| A6 | 跖面 BC 敏感性 ≥15% | — | ✅ 已过（降 24%，A 档前结论） | | |")
    A("")
    A("## 1. A1 细节（三阶段边界）")
    A("")
    for tag, dd in (("乙(median)", v1), ("甲", v1_jia)):
        A(f"- **{tag}**: phases = " + "；".join(
            f"P{p['phase']} {p['range_m'][0]}–{p['range_m'][1]} m"
            for p in dd.get("phases", [])))
    A("")
    A("## 2. A2 细节（逐目标行为）")
    A("")
    A("| 目标 | 部位 | 乙(median) 行为 | 命中(严格) |")
    A("|---|---|---|---|")
    for k, r in v2["per_target"].items():
        A(f"| {k} | {'+'.join(r['parts'])} | {'/'.join(r['behavior'])} | "
          f"{'✅' if r['hit_strict'] else '❌'} |")
    A("")
    A("## 3. A3 细节（OR vs [Y25]）")
    A("")
    A("| 部位 | [Y25] OR | 乙 OR | 相对误差 | ±30% 内 |")
    A("|---|---|---|---|---|")
    for cn, r in v3["per_part"].items():
        o = f"{r['ours']:.3f}" if r.get("ours") else str(r.get("status"))
        e = f"{r['rel_err']*100:.1f}%" if r.get("ours") else "—"
        A(f"| {cn} | {r['target']} | {o} | {e} | {'✅' if r['within'] else '❌'} |")
    A("")
    A("## 4. A4 细节（PCA 双路径）")
    A("")
    A(f"- 方法：{v4['method']}")
    A(f"- 解释方差比：PC1={v4['explained_var_ratio'][0]:.3f}、"
      f"PC2={v4['explained_var_ratio'][1]:.3f}")
    A(f"- LP 可分：**{v4['separable']}**（{v4['lp_message']}）")
    A(f"- 分组：G1={v4['group1']}；G2={v4['group2']}；足部未标注（§六 未列）")
    A("- 部位坐标（PC1, PC2 载荷）：")
    for cn, xy in v4["part_positions"].items():
        A(f"  - {cn}: ({xy[0]:+.4f}, {xy[1]:+.4f})")
    A("")
    A("## 5. 甲/乙矩阵首折对照")
    A("")
    A("| 部位 | 甲 1D | 乙 median |")
    A("|---|---|---|")
    for cn in pa.REGION_NAMES:
        A(f"| {cn} | {jia['first'].get(cn) or '>50'} | "
          f"{yi['first'].get(cn) or '>50'} |")
    A("")
    A(f"_生成：s4_alignment_bc2.py；耗时 {elapsed:.0f}s；数据 = "
      "risk_1d_bending.json / fracture_matrix_bc2.json + bc2 hdf5 实测。_")
    return "\n".join(L) + "\n"


def md_b(v: dict) -> str:
    return v.get("verdict", "—")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sub-stat", default="median", choices=["median", "p95"])
    args = ap.parse_args(argv)
    t0 = time.time()

    fm = json.loads(BC2_JSON.read_text(encoding="utf-8"))
    heights = np.asarray(fm["heights_m"], dtype=float)

    print("[s4] route 甲 (1D) ...", flush=True)
    jia = route_jia()
    print("[s4] route 乙 (BC2, median) ...", flush=True)
    yi = route_yi("median")
    print("[s4] route 乙 (BC2, p95 sensitivity) ...", flush=True)
    yi_p95 = route_yi("p95")

    v1 = a1_verdict(yi["phases"]); v1_jia = a1_verdict(jia["phases"])
    v1_p95 = a1_verdict(yi_p95["phases"])
    v2 = a2_verdict(yi["logistic"]); v2_jia = a2_verdict(jia["logistic"])
    v2_p95 = a2_verdict(yi_p95["logistic"])
    v3 = a3_verdict(yi["logistic"]); v3_jia = a3_verdict(jia["logistic"])
    v3_p95 = a3_verdict(yi_p95["logistic"])
    v4 = a4_from_matrix(yi, heights, fm)

    elapsed = time.time() - t0
    vs = {"v1": v1, "v2": v2, "v3": v3, "v4": v4,
          "v1_jia": v1_jia, "v2_jia": v2_jia, "v3_jia": v3_jia,
          "v1_p95": v1_p95, "v2_p95": v2_p95, "v3_p95": v3_p95}
    out = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "official_route": "乙 = fracture_matrix_bc2.json (align_pca=True, median)",
        "criteria_source": "docs/methods/OpenSim_FE复现方案.md §六",
        "A1": v1, "A2": v2, "A3": v3, "A4": v4,
        "A1_jia": v1_jia, "A2_jia": v2_jia, "A3_jia": v3_jia,
        "A1_p95": v1_p95, "A2_p95": v2_p95, "A3_p95": v3_p95,
        "first_fracture_jia_m": jia["first"],
        "first_fracture_yi_m": yi["first"],
        "logistic_yi": yi["logistic"],
        "logistic_jia": jia["logistic"],
        "cluster_ranges_yi": yi["cluster_ranges"],
        "phases_yi": yi["phases"],
        "detail_yi": yi["detail"],
        "elapsed_s": round(elapsed, 1),
    }
    BC2_OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                       encoding="utf-8")
    BC2_REPORT.write_text(md_report(jia, yi, fm, vs, elapsed),
                          encoding="utf-8")

    summary = (f"B1.4: S4 对标完成 — 正式(乙 median) A1={v1['verdict']} "
               f"A2={v2['verdict']}({v2['n_hit']}/5) A3={v3['verdict']}"
               f"({v3['n_within']}/3) A4={v4['verdict']}；"
               f"甲 A3={v3_jia['verdict']}（预期结构性失败）；"
               f"产物 S4_ALIGNMENT_BC2_REPORT.md + s4_alignment_bc2.json"
               f"（{elapsed:.0f}s）")
    log_line(summary)
    print("\n" + summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

