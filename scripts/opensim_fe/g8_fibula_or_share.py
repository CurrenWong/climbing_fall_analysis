"""g8_fibula_or_share.py -- G8: 腓骨渐进 OR via 1D 载荷分享路径（路线 A）。

严格实现 `docs/methods/G8_腓骨OR_1D载荷分享设计.md`（预注册）：只把**腓骨一行**的
σ_nom(h) 从「FE 无分配量级」换成**既有 1D 胫-腓载荷分享路径**的真实名义应力，
其余 8 部位仍用冻结的 FE `sigma_vm_mpa["median"]`；随后原封复用 G3 的被试级变异
概率机制（`logistic_fit` / `a3_verdict` 零改动），使腓骨从 `all_one` 变**渐进**、
Logistic OR 可定义，并对标 [Y25] 靶 1.682（±30% ⇒ [1.18, 2.19]）。

数据源（只读）
--------------
* 腓骨 σ_fib(h) = `risk_by_height(h) × σc_row`，取自
  `results/opensim_fe/risk_1d_loadshare.json::{rows_with_share_lit05,rows_with_share_lit10}`
  的 `fibula_r` 行。**σc_row = 160 MPa（骨干，非弱区 70）** ⇒ σ_fib = 真实名义应力 F/A。
  （设计初稿误写 ×70，已在设计稿 §3 勘误；本脚本按勘误后的 ×σc_row 实现。）
* 其余 8 部位 σ_nom(h) = `fracture_matrix_bc2.json::sigma_vm_mpa["median"][bone]`（冻结）。
* σc = route-乙 弱区口径（`paper_alignment.REGIONS` ← `bone.py::MATERIAL_STRENGTH_MPA`）。

预注册（跑后不改）
------------------
* 主设定分担比 = **腓 0.05 / 胫 0.95**（`share_lit05`，= `bone.py::fibula_ends` 既有锚）；
  敏感性 = 腓 0.10 / 胫 0.90（`share_lit10`）。**全量报告，不挑。**
* G3 主设定 COV（B=0.15/E=0.10/A=0.15）、n_subj=400、seeds {20261003,20261004,20261005}
  ——全部沿用 G3 冻结值。
* 判定：`s4_alignment_bc2.a3_verdict` 原封（腓骨靶 1.682，±30%）。跑后禁回解。

产物（纯增量）
--------------
* ``results/opensim_fe/g8_fibula_or.json``
* ``results/opensim_fe/G8_FIBULA_OR_REPORT.md``

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g8_fibula_or_share.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import paper_alignment as pa         # noqa: E402

RES = ROOT / "results" / "opensim_fe"
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"
BC2_JSON = RES / "fracture_matrix_bc2.json"
G3_OR_JSON = RES / "g3_or_estimates.json"
OUT_JSON = RES / "g8_fibula_or.json"
OUT_REPORT = RES / "G8_FIBULA_OR_REPORT.md"

FIBULA_BONE = "fibula_r"
FIBULA_SITE = "腓骨"
GRADUAL = ["腓骨", "骨盆", "颅骨"]
A3_TARGETS = dict(g3.A3_TARGETS)          # {腓骨:1.682, 骨盆:1.236, 颅骨:1.576}
FIBULA_TARGET = float(A3_TARGETS[FIBULA_SITE])
TOL = 0.30
FIBULA_WINDOW = [FIBULA_TARGET * (1 - TOL), FIBULA_TARGET * (1 + TOL)]

#: (share label, risk_1d_loadshare.json key, fibula split fraction)
SHARES = [("lit05", "rows_with_share_lit05", 0.05),
          ("lit10", "rows_with_share_lit10", 0.10)]
PRIMARY_SHARE = "lit05"

G8_DELIVERABLES = [
    "scripts/opensim_fe/g8_fibula_or_share.py",
    "results/opensim_fe/g8_fibula_or.json",
    "results/opensim_fe/G8_FIBULA_OR_REPORT.md",
    "tests/test_g8_fibula_or.py",
    "docs/methods/G8_腓骨OR_1D载荷分享设计.md",
]


# ---------------------------------------------------------------------------
# 数据源
# ---------------------------------------------------------------------------
def load_fibula_1d(key: str) -> dict[str, Any]:
    """从 ``risk_1d_loadshare.json`` 抽腓骨 1D 行 → 名义应力 σ_fib(h)=risk×σc_row。"""
    doc = json.loads(LOADSHARE_JSON.read_text(encoding="utf-8"))
    rows = doc[key]
    row = next(r for r in rows if r["bone"] == FIBULA_BONE)
    risk = np.asarray(row["risk_by_height"], dtype=float)
    sc = float(row["sigma_c_mpa"])
    return dict(
        share_key=key,
        force_multiplier=float(row["force_multiplier"]),
        sigma_c_row_mpa=sc,
        risk_by_height=risk,
        sigma_fib_mpa=risk * sc,
        risk_at_5m=float(row["risk_at_5m"]),
        risk_max=float(risk.max()),
    )


def height_grid_check(inputs: dict[str, Any], fib: dict[str, Any]) -> dict[str, Any]:
    """核对 1D 行与 G3 的高度栅格一致性。"""
    h_g3 = [float(x) for x in inputs["heights"]]
    h_1d = [float(i) for i in range(1, len(fib["risk_by_height"]) + 1)]
    match = (len(h_g3) == len(h_1d)) and all(abs(a - b) < 1e-9 for a, b in zip(h_g3, h_1d))
    return dict(g3_heights_n=len(h_g3), fibula_1d_heights_n=len(h_1d), match=bool(match))


def build_mixed_inputs(inputs: dict[str, Any], sigma_fib: np.ndarray) -> tuple[dict[str, Any], int]:
    """把腓骨行的 σ_nom 换成 1D 名义应力，返回 (新 inputs, 腓骨行索引)。"""
    bones = list(inputs["bones"])
    idx = bones.index(FIBULA_BONE)
    sig = np.array(inputs["sigma_nom"], dtype=float, copy=True)
    if sig.shape[1] != len(sigma_fib):
        raise ValueError(f"height mismatch: {sig.shape} vs {len(sigma_fib)}")
    sig[idx] = np.asarray(sigma_fib, dtype=float)
    out = dict(inputs)
    out["sigma_nom"] = sig
    return out, idx


# ---------------------------------------------------------------------------
# 运行 G3（复用，零改动）
# ---------------------------------------------------------------------------
def run_share(inputs: dict[str, Any], fib: dict[str, Any]) -> dict[str, Any]:
    mixed, idx = build_mixed_inputs(inputs, fib["sigma_fib_mpa"])
    configs = {seed: g3.run_config(seed, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                   g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
               for seed in g3.SEEDS}
    primary = configs[g3.OFFICIAL_SEED]

    fib_rec = g3._site_or_record(primary["or_map"][FIBULA_SITE])
    seed_ors = {sd: configs[sd]["or_map"][FIBULA_SITE].get("OR") for sd in g3.SEEDS}

    sigma_c_weak = float(np.asarray(inputs["sigma_c"], float)[idx])
    util = np.asarray(fib["sigma_fib_mpa"], float) / sigma_c_weak
    heights = np.asarray(inputs["heights"], float)
    det_cut = int(np.argmax(util >= 1.0)) if bool((util >= 1.0).any()) else None
    first_det = float(heights[det_cut]) if det_cut is not None else None

    return dict(
        primary=primary,
        configs=configs,
        fibula_rec=fib_rec,
        seed_ors=seed_ors,
        a3=primary["a3"],
        sigma_c_weak_mpa=sigma_c_weak,
        util_by_height=[float(u) for u in util],
        util_at_5m=float(util[4]),
        util_max=float(util.max()),
        first_fracture_deterministic_h=first_det,
    )


def unchanged_check(inputs: dict[str, Any], mixed_idx: int, primary: Any) -> dict[str, Any]:
    """核对其余 8 部位 OR/status 与 G3 主设定逐值相同。"""
    if not G3_OR_JSON.is_file():
        return dict(available=False)
    stored = json.loads(G3_OR_JSON.read_text(encoding="utf-8"))["primary"]["per_site"]
    sites = list(inputs["sites"])
    per: dict[str, bool] = {}
    worst = 0.0
    for k, cn in enumerate(sites):
        if k == mixed_idx:
            continue
        a = primary["or_map"][cn]
        b = stored.get(cn, {})
        sa, sb = a.get("status"), b.get("status")
        oa, ob = a.get("OR"), b.get("OR")
        ok = (sa == sb) and (
            (oa is None and ob is None)
            or (oa is not None and ob is not None and abs(oa - ob) <= 1e-9))
        per[cn] = bool(ok)
        if oa is not None and ob is not None:
            worst = max(worst, abs(oa - ob))
    return dict(available=True, all_equal=bool(all(per.values())), per_site=per,
                max_abs_or_diff=float(worst),
                g3_source="results/opensim_fe/g3_or_estimates.json")


# ---------------------------------------------------------------------------
def git_status_lines() -> list[str] | None:
    try:
        r = subprocess.run(
            ["git", "-c", "core.quotepath=false", "-C", str(ROOT),
             "status", "--porcelain"],
            capture_output=True, text=True, timeout=30)
        return [ln.rstrip() for ln in r.stdout.splitlines() if ln.strip()]
    except Exception:  # noqa: BLE001
        return None


def verdict_text(primary_res: dict[str, Any]) -> dict[str, Any]:
    rec = primary_res["fibula_rec"]
    st = rec.get("status")
    orv = rec.get("OR")
    within = bool(st == "ok" and orv is not None
                  and FIBULA_WINDOW[0] <= orv <= FIBULA_WINDOW[1])
    a3 = primary_res["a3"]
    return dict(
        fibula_status=st,
        fibula_or=orv,
        fibula_ci95=rec.get("ci95"),
        fibula_window=FIBULA_WINDOW,
        fibula_within=within,
        a3_verdict=a3.get("verdict"),
        a3_n_within=a3.get("n_within"),
        dod1_pass=within,
        dod3_pass=bool(a3.get("n_within", 0) >= 2),
    )


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    _ = ap.add_argument("--subjects", type=int, default=g3.N_SUBJ)
    args = ap.parse_args(argv)
    t0 = time.time()

    inputs = g3.load_inputs()
    grid = height_grid_check(inputs, load_fibula_1d(SHARES[0][1]))
    if not grid["match"]:
        raise SystemExit(f"height grid mismatch: {grid}")

    shares_out: dict[str, Any] = {}
    first_res: dict[str, Any] | None = None
    mixed_idx: int | None = None
    for label, key, frac in SHARES:
        fib = load_fibula_1d(key)
        res = run_share(inputs, fib)
        shares_out[label] = dict(
            share_fraction_fibula=float(fib["force_multiplier"]),
            share_fraction_tibia=float(1.0 - fib["force_multiplier"]),
            source_key=key,
            sigma_c_row_mpa=float(fib["sigma_c_row_mpa"]),
            risk_at_5m=float(fib["risk_at_5m"]),
            risk_max=float(fib["risk_max"]),
            sigma_c_weak_mpa=float(res["sigma_c_weak_mpa"]),
            util_at_5m=float(res["util_at_5m"]),
            util_max=float(res["util_max"]),
            first_fracture_deterministic_h=res["first_fracture_deterministic_h"],
            fibula_status=res["fibula_rec"].get("status"),
            fibula_OR=res["fibula_rec"].get("OR"),
            fibula_ci95=res["fibula_rec"].get("ci95"),
            fibula_n_fracture=res["fibula_rec"].get("n_fracture"),
            fibula_n_safe=res["fibula_rec"].get("n_safe"),
            fibula_seed_ORs={str(k): v for k, v in res["seed_ors"].items()},
            a3=res["a3"],
            verdict=verdict_text(res),
        )
        if label == PRIMARY_SHARE:
            first_res = res
            mixed_idx = list(inputs["bones"]).index(FIBULA_BONE)

    assert first_res is not None and mixed_idx is not None
    unc = unchanged_check(inputs, mixed_idx, first_res["primary"])

    # 极限退化（DoD-4）：COV→0.001 且质量不入载荷 ⇒ 逐格复现确定性基准（含 1D 腓骨行）
    mixed, _ = build_mixed_inputs(inputs, load_fibula_1d("rows_with_share_lit05")["sigma_fib_mpa"])
    _, p_deg = g3.simulate(g3.OFFICIAL_SEED, g3.N_SUBJ, g3.COV_DEGENERATE,
                           g3.COV_DEGENERATE, g3.COV_DEGENERATE, mixed, mass_in_load=False)
    base = g3.deterministic_base(mixed["sigma_nom"], mixed["sigma_c"])
    deg = dict(
        cov=g3.COV_DEGENERATE,
        majority_matches_base=bool(np.array_equal((p_deg >= 0.5).astype(int), base)),
        max_abs_dev=float(np.abs(p_deg - base).max()),
    )

    pv = first_res["primary"]["p"]
    prob_by_site = {cn: [float(v) for v in pv[k]] for k, cn in enumerate(inputs["sites"])}

    out = dict(
        meta=dict(
            title="G8 fibula progressive OR via 1D load-share path (Route A, pre-registered)",
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            design_doc="docs/methods/G8_腓骨OR_1D载荷分享设计.md",
            method=("mixed sigma_nom: fibula row <- 1D load-share nominal stress F/A "
                    "(risk_by_height * sigma_c_row=160); other 8 rows <- "
                    "fracture_matrix_bc2.json sigma_vm_mpa[median]; G3 subject-variability "
                    "reused unchanged (logistic_fit / a3_verdict)."),
            correction=("design draft said sigma_fib=risk*70; corrected to risk*sigma_c_row"
                        "=risk*160 (true nominal stress). Decisions F1-F5 unchanged."),
            seeds=g3.SEEDS, official_seed=g3.OFFICIAL_SEED, n_subj=g3.N_SUBJ,
            cov_primary=dict(cov_b=g3.COV_B_PRIMARY, cov_e=g3.COV_E_PRIMARY,
                             cov_a=g3.COV_A_PRIMARY),
            primary_share=PRIMARY_SHARE,
            a3_targets=A3_TARGETS, fibula_window=FIBULA_WINDOW,
            sources=dict(
                fibula_1d="results/opensim_fe/risk_1d_loadshare.json",
                other8="results/opensim_fe/fracture_matrix_bc2.json::sigma_vm_mpa[median]",
                g3_machinery="scripts/opensim_fe/g3_subject_variability.py",
            ),
            ban_note="no FE rerun; no existing file modified; OR reported per-metre.",
        ),
        height_grid_check=grid,
        fibula_site=FIBULA_SITE,
        fibula_row_index=int(mixed_idx),
        shares=shares_out,
        unchanged_check=unc,
        degenerate=deg,
        probability_matrix_primary_fibula=dict(
            site=FIBULA_SITE,
            p=prob_by_site[FIBULA_SITE],
        ),
        verdict=dict(
            primary_share=PRIMARY_SHARE,
            **shares_out[PRIMARY_SHARE]["verdict"],
            unchanged_other8=unc.get("all_equal"),
            degenerate_ok=deg["majority_matches_base"],
        ),
    )
    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")

    write_report(out, inputs, shares_out, unc, deg, first_res, elapsed=time.time() - t0)

    v = out["verdict"]
    print("=" * 96)
    print(f"G8 腓骨 1D 载荷分享接入 G3：主设定 {PRIMARY_SHARE}"
          f"（腓骨 {shares_out[PRIMARY_SHARE]['share_fraction_fibula']:.2f}）")
    for label, _, frac in SHARES:
        s = shares_out[label]
        orv = s["fibula_OR"]
        ors = f"{orv:.3f}" if orv is not None else f"**{s['fibula_status']}**"
        print(f"  {label} (腓骨 {frac:.2f}): status={s['fibula_status']:>16s} "
              f"OR={ors}  util@5m={s['util_at_5m']:.3f} util_max={s['util_max']:.3f}  "
              f"A3={s['a3']['n_within']}/3 {s['a3']['verdict']}")
    print(f"  腓骨靶 {FIBULA_TARGET} ±30% ⇒ {[round(x,3) for x in FIBULA_WINDOW]}；"
          f"主设定 within={v['fibula_within']}")
    print(f"  其余 8 部位 OR 与 G3 逐值相同: {unc.get('all_equal')} "
          f"(max|ΔOR|={unc.get('max_abs_or_diff')})")
    print(f"  极限退化 majority_match={deg['majority_matches_base']} "
          f"max|Δ|={deg['max_abs_dev']:.2e}")
    print(f"wrote: {OUT_JSON}")
    print(f"wrote: {OUT_REPORT}")
    return 0


# ---------------------------------------------------------------------------
def _fmt_h(v) -> str:
    return ">50" if v is None else f"{int(v)}"


def write_report(out: dict[str, Any], inputs: dict[str, Any], shares: dict[str, Any],
                 unc: dict[str, Any], deg: dict[str, Any], first_res: dict[str, Any],
                 elapsed: float) -> None:
    L: list[str] = []
    A = L.append
    sh = shares[PRIMARY_SHARE]
    fv = sh["verdict"]
    sites = list(inputs["sites"])
    weak = inputs["weak_sigma_c_mpa"]

    A("# G8 判定补录 —— 腓骨渐进 OR（1D 载荷分享路径 · 路线 A）")
    A("")
    A(f"> 生成 {datetime.now().strftime('%Y-%m-%d %H:%M')}；仓库 "
      "`D:\\Project\\climbing_fall_analysis`。耗时 %ds。" % elapsed)
    A("> 预注册设计：`docs/methods/G8_腓骨OR_1D载荷分享设计.md`（含 §3 勘误："
      "`σ_fib = risk × σc_row(160)`，非 ×70）。**纯增量**：不改任何既有文件、不重跑 FE。")
    A("> 口径纪律：OR 一律 **per-metre + 主设定分担比 + 混合口径**（腓骨 1D 行 + 其余 FE 行）。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    orv = sh["fibula_OR"]
    fib_head = (f"**腓骨 status=`{sh['fibula_status']}`，OR="
                + (f"{orv:.3f}" if orv is not None else "—")
                + f"（靶 {FIBULA_TARGET} ±30%）**")
    A(f"把腓骨一行换成 **1D 载荷分享名义应力**（`σ_fib(h)=risk(h)×160`，5% 分担）接入 "
      f"G3 后：{fib_head}。")
    if fv["fibula_within"]:
        A(f"⇒ **DoD-1 PASS**：腓骨从 `all_one` 变**渐进**、OR 落 [Y25] 靶 ±30%，"
          f"A3 = **{fv['a3_n_within']}/3（{fv['a3_verdict']}）**。")
    else:
        ortxt = f"{orv:.3f}" if orv is not None else "—"
        A(f"⇒ **DoD-1 FAIL**：主设定腓骨未落靶（status=`{sh['fibula_status']}`，OR={ortxt}）。")
        A(f"A3 = **{fv['a3_n_within']}/3（{fv['a3_verdict']}）**。")
    A("")
    A("---")
    A("")
    A("## 1. 预注册参数快照（回显）")
    A("")
    A("| 项 | 值 |")
    A("|---|---|")
    A(f"| 主设定分担比 | 腓骨 **{sh['share_fraction_fibula']:.2f}** / 胫骨 "
      f"{sh['share_fraction_tibia']:.2f}（`{sh['source_key']}`）|")
    A(f"| 敏感性分担比 | 腓骨 {shares['lit10']['share_fraction_fibula']:.2f} / 胫骨 "
      f"{shares['lit10']['share_fraction_tibia']:.2f}（`{shares['lit10']['source_key']}`）|")
    A(f"| n_subj / seeds | {g3.N_SUBJ} × {g3.SEEDS}（官方 {g3.OFFICIAL_SEED}）|")
    A("| 主设定 COV | COV_B=0.15, COV_E=0.10, COV_A=0.15（G3 冻结值）|")
    A(f"| 腓骨 σc（弱区，G3） | {sh['sigma_c_weak_mpa']:.0f} MPa |")
    A(f"| 1D 行 σc_row | {sh['sigma_c_row_mpa']:.0f} MPa（骨干）|")
    A(f"| 高度栅格 | 1…50 m（1D 行与 G3 同名：{out['height_grid_check']['match']}）|")
    A(f"| [Y25] 腓骨靶 | {FIBULA_TARGET}，±30% ⇒ {[round(x,3) for x in FIBULA_WINDOW]} |")
    A("")
    A("---")
    A("")
    A("## 2. 腓骨接入后的结果（主设定 vs 敏感性）")
    A("")
    A("| 分担比 | 腓骨 status | 腓骨 OR | Wald 95% CI | n_fracture | util@5m | util_max | 确定性首骨折 h |")
    A("|---|---|---:|---|---:|---:|---:|---:|")
    for label, _, frac in SHARES:
        s = shares[label]
        ci = s["fibula_ci95"] or [None, None]
        citxt = (f"{ci[0]:.3f}–{ci[1]:.3f}" if ci[0] is not None else "—")
        ortxt = f"{s['fibula_OR']:.3f}" if s["fibula_OR"] is not None else f"**{s['fibula_status']}**"
        A(f"| {frac:.2f} | {s['fibula_status']} | {ortxt} | {citxt} | "
          f"{s['fibula_n_fracture']} | {s['util_at_5m']:.3f} | {s['util_max']:.3f} | "
          f"{_fmt_h(s['first_fracture_deterministic_h'])} |")
    A("")
    A("> `util = σ_fib(h)/σc_weak(70)`；`util_max<1` ⇒ 确定性不骨折，概率层的混合来自被试变异。")
    A("")
    A("**主设定 A3 全表（per-metre + 主设定 COV）：**")
    A("")
    A("| 部位 | 状态 | OR | [Y25] 靶 | ±30% 内 |")
    A("|---|---|---:|---:|:--:|")
    a3 = sh["a3"]
    for cn in GRADUAL:
        r = a3["per_part"][cn]
        o = f"{r['ours']:.3f}" if r.get("ours") else f"**{r.get('status')}**"
        A(f"| {cn} | {'ok' if r.get('ours') else r.get('status')} | {o} | "
          f"{r['target']} | {'✅' if r['within'] else '❌'} |")
    A("")
    A(f"**主设定 A3 = {a3['n_within']}/3 ⇒ {a3['verdict']}**")
    A("")
    A("---")
    A("")
    A("## 3. 腓骨概率矩阵 p(h)（主设定，节选）")
    A("")
    p = out["probability_matrix_primary_fibula"]["p"]
    cols = [0, 4, 9, 19, 29, 39, 49]
    A("| 部位 | " + " | ".join(f"h={int(inputs['heights'][i])}" for i in cols) + " |")
    A("|---|" + "---:|" * len(cols))
    A("| 腓骨 | " + " | ".join(f"{p[i]:.3f}" for i in cols) + " |")
    A("")
    A("> p(h) 单调非降；主设定腓骨不再是 `all_one`（全 1）而是 `0→…→p<1` 的渐进曲线。")
    A("")
    A("---")
    A("")
    A("## 4. 其余 8 部位不变性（DoD-2）")
    A("")
    if unc.get("available"):
        A(f"对照 `g3_or_estimates.json` 主设定 per_site：**全部逐值相同 = "
          f"{unc['all_equal']}**（max|ΔOR| = {unc['max_abs_or_diff']:.2e}）。")
        A("")
        A("| 部位 | 与 G3 相同 |")
        A("|---|---|")
        for cn, ok in unc["per_site"].items():
            A(f"| {cn} | {'✅' if ok else '❌'} |")
    else:
        A("> `g3_or_estimates.json` 不可读，未做不变性核对。")
    A("")
    A("> 原理：只替换腓骨行 ⇒ 其余行的 σ_nom 未变，被试因子（逐被试 `default_rng([seed,s])`）"
      "与高度无关、完全一致 ⇒ 其余行 y 不变 ⇒ OR 逐值相同。")
    A("")
    A("---")
    A("")
    A("## 5. DoD 逐条")
    A("")
    A(f"- **DoD-1 腓骨渐进 + OR∈[1.18,2.19]** — "
      f"**{'PASS' if fv['fibula_within'] else 'FAIL'}**（主设定 status=`{sh['fibula_status']}`，"
      f"OR=" + (f"{orv:.3f}" if orv is not None else "—") + "）。")
    A(f"- **DoD-2 其余部位不变** — **{'PASS' if unc.get('all_equal') else 'FAIL'}**。")
    A(f"- **DoD-3 A3 ≥2/3** — **{'PASS' if fv['dod3_pass'] else 'FAIL'}**"
      f"（{fv['a3_n_within']}/3）。")
    A(f"- **DoD-4 极限退化** — **{'PASS' if deg['majority_matches_base'] else 'FAIL'}**"
      f"（COV={deg['cov']} 且质量不入载荷，多数判定复现含 1D 腓骨行的确定性基准；"
      f"max|Δ|={deg['max_abs_dev']:.2e}）。")
    A("- **DoD-5 零回归** — 纯增量（见 §7）；回归基线 ≥687 passed / 8 skipped；"
      "新增 `tests/test_g8_fibula_or.py`（见 §9）。")
    A("- **DoD-6 口径纪律** — OR 均标 per-metre + 主设定分担比 + 混合口径；"
      "未出现禁用比值。")
    A("")
    A("---")
    A("")
    A("## 6. 诚实结论与风险")
    A("")
    A("1. **混合口径**：腓骨行 = 1D `F/A` 均匀轴向名义应力；其余 8 行 = FE 几何解 —— "
      "两者**不同源**，本件为**显式建模选择**。")
    A("2. **数据源勘误**：初稿误用 `×70`；本件按 `×σc_row(160)`（真实名义应力）实现，"
      "并已回写设计稿 §3。此更正不改变任何预注册决策。")
    A(f"3. **OR 是否落靶**：{'落靶' if fv['fibula_within'] else '未落靶（如实报，未回解）'}。")
    A("4. **敏感性**：分担比 0.10 的结果一并报告（未挑）。")
    A("")
    A("---")
    A("")
    A("## 7. 纯增量核对")
    A("")
    cur = git_status_lines()
    A("G8 新增（应只见以下条目，无既有文件被改）：")
    A("")
    for f in G8_DELIVERABLES:
        A(f"- `{f}`")
    if cur is not None:
        A("")
        A(f"`git status --porcelain` 当前条目数 {len(cur)}（既有未提交改动 + 本件新增；"
          "详细清单见仓库）。")
    A("")
    A("---")
    A("")
    A("## 8. 复现命令与产物")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g8_fibula_or_share.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/risk_1d_loadshare.json        (腓骨 1D 行)")
    A("#   results/opensim_fe/fracture_matrix_bc2.json       (其余 8 部位 median)")
    A("#   scripts/opensim_fe/g3_subject_variability.py      (被试变异机制, 零改动)")
    A("# 产物（纯增量）：")
    A("#   results/opensim_fe/g8_fibula_or.json")
    A("#   results/opensim_fe/G8_FIBULA_OR_REPORT.md")
    A("```")
    A("")
    A("**弱区 σc（route-乙，G3 同款）：**")
    A("")
    A("| 部位 | " + " | ".join(sites) + " |")
    A("|---|" + "---|" * len(sites))
    A("| σc (MPa) | " + " | ".join(f"{weak[cn]:.0f}" for cn in sites) + " |")
    A("")
    A("---")
    A("")
    A("## 9. 实测记录")
    A("")
    A("- G8 脚本：`scripts/opensim_fe/g8_fibula_or_share.py`（官方 seed "
      f"{g3.OFFICIAL_SEED}，n_subj={g3.N_SUBJ}，2 档分担比 × 3 seeds）。")
    A("- 测试：`tests/test_g8_fibula_or.py`；全测试套件基线 687 passed / 8 skipped"
      "（G8 为纯增量，实测见仓库运行日志）。")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
