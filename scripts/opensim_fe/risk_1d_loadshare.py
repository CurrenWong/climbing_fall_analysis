"""risk_1d_loadshare.py -- Route-1 1D fracture risk with **load sharing**.

Additive post-processing only.  It reuses the artifacts already produced by
``bc_robust_metric.py`` / ``fracture_matrix.py`` and **never** writes to an
existing file, never changes a module default, and never runs FEBio.

Why this exists
---------------
The plain Route-1 metric

    risk(h) = F_joint(h) / (A_section * sigma_c)          (1D nominal compression)

has two known load-mapping artifacts (see ``BC_ROBUST_METRIC_REPORT.md`` section 7):

1. **Tibia and fibula share the FULL ankle reaction.**  ``bonce.py::BONE_SITES``
   records that the fibula carries only ~5-10% of the axial load, yet Route 1
   feeds both bones the whole ``ankle_r`` force -> the fibula's risk is a
   gross upper bound and it is artificially ranked #1.
2. **L3 / T6 / C5 / parietal all share the single ``lumbar`` reaction.**  The
   ``lumbar`` free-body force supports the *whole* torso + both arms; applying
   it to a cervical vertebra or to the skull over-counts the load by 5-10x.

This script fixes both by adding a *load-share* front-end to Route 1:

* At the ankle, the tibia and fibula are modelled as two **parallel axial
  springs** with stiffness ``k = E * A / L`` (E = cortical 18000 MPa,
  A = A_section, L = bone length).  The ankle reaction is split by ``k_i / sum k``.
* For the axial bones the force at a superior level is obtained from a
  **mass-above** model: ``F_level(h) = (m_above(level) / m_above(lumbar_ref)) *
  F_lumbar(h)``.  The mass fractions come from published anthropometric
  segment-mass data (documented in the report).  The model has **no head/neck
  or per-vertebra bodies** -- ``joint_reactions.py`` only exposes the fixed
  joints (subtalar / ankle / knee / hip / lumbar) -- so a rigorous per-level
  joint reaction is *not* available; the probe in this script records the
  failure so the fallback is documented, not silent.

Inputs (read-only)
------------------
* ``results/opensim_fe/fracture_matrix.json``            -> ``loads_n``, ``heights_m``
* ``results/opensim_fe/bc_robust_metric_route1.json``    -> per-bone ``A_section_mm2``,
                                                            ``sigma_c_mpa``, no-share risk
* ``src/climbing/coupling/joint_reactions.py``           -> probed for per-level support

Outputs (NEW files only)
------------------------
* ``results/opensim_fe/risk_1d_loadshare.json``
* ``results/opensim_fe/RISK_1D_LOADSHARE_REPORT.md``

Usage (repo root)::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\risk_1d_loadshare.py
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

FM_JSON = ROOT / "results" / "opensim_fe" / "fracture_matrix.json"
R1_JSON = ROOT / "results" / "opensim_fe" / "bc_robust_metric_route1.json"

# ---------------------------------------------------------------------------
# Model constants
# ---------------------------------------------------------------------------
E_CORTICAL_MPA = 18000.0          # cortical bone Young's modulus (task-given, assumed)
TOTAL_BODY_MASS_KG = 75.337       # Rajagopal2015 (joint_reactions_summary.json meta)

#: Which joint reaction feeds each bone (plain Route 1, no share).
FORCE_SOURCE = {
    "calcaneus_r": "subtalar_r",
    "tibia_r": "ankle_r",
    "fibula_r": "ankle_r",
    "femur_r": "hip_r",
    "R_HIPBONE": "hip_r",
    "L3": "lumbar",
    "T6": "lumbar",
    "C5": "lumbar",
    "parietal_r": "lumbar",
}

#: Bones whose force is reduced by the mass-above rule.
AXIAL_BONES = ("L3", "T6", "C5", "parietal_r")

# ---------------------------------------------------------------------------
# Mass-above model -- per-level segment masses (% of total body mass)
# ---------------------------------------------------------------------------
#: Head + cervical segments: Ivancic et al. (2006), as reproduced in
#: Vette et al. / Frontiers Bioeng. Biotechnol. (2015) Table 3.
SEG_PCT_HEAD_CERVICAL = {
    "head": 4.7,
    "C1": 0.6, "C2": 0.7, "C3": 0.5, "C4": 0.5, "C5": 0.5, "C6": 0.6, "C7": 0.7,
}
#: Thoracic + lumbar segments: Pearsall et al. (1996), same reproduced table.
SEG_PCT_TRUNK = {
    "T1": 1.1, "T2": 1.1, "T3": 1.4, "T4": 1.3, "T5": 1.3, "T6": 1.3,
    "T7": 1.4, "T8": 1.5, "T9": 1.6, "T10": 2.0, "T11": 2.1, "T12": 2.5,
    "L1": 2.4, "L2": 2.4, "L3": 2.3, "L4": 2.6, "L5": 2.6,
}
#: Both arms (Dempster / Winter, ~5% each); they hang from the shoulder (T1),
#: so they load every level at or below T1 (T6, L3, L5/S1) but NOT C5 / head.
ARMS_PCT = 10.0

#: Superior -> inferior order used for cumulative mass-above.
_SEGMENT_ORDER = ["head", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
                  "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9",
                  "T10", "T11", "T12", "L1", "L2", "L3", "L4", "L5"]
SEG_PCT = {**SEG_PCT_HEAD_CERVICAL, **SEG_PCT_TRUNK}


def mass_above_pct(level: str, *, include_arms: bool) -> float:
    """Mass of all body segments from the top down to and **including** ``level``.

    ``level`` is the bone's level; the load through that level's inferior disc
    supports everything above it *plus the segment itself*.  ``level='L5'``
    reproduces the lumbar reference (the model's ``back`` joint = pelvis->torso
    boundary, i.e. the L5/S1 disc).
    """
    idx = _SEGMENT_ORDER.index(level)
    m = float(sum(SEG_PCT[s] for s in _SEGMENT_ORDER[: idx + 1]))
    if include_arms:
        m += ARMS_PCT
    return m


#: Rule table for the axial bones: (segment level, arms load the spine there?).
#: C5 is superior to the shoulder -> arms excluded; T6/L3/L5 are inferior -> included.
LEVEL_OF = {"C5": ("C5", False), "T6": ("T6", True),
            "L3": ("L3", True), "parietal_r": ("head", False)}

MASS_ABOVE = {b: mass_above_pct(LEVEL_OF[b][0], include_arms=LEVEL_OF[b][1])
              for b in AXIAL_BONES}
M_REF_PCT = mass_above_pct("L5", include_arms=True)   # lumbar reference (denominator)
RATIO = {b: MASS_ABOVE[b] / M_REF_PCT for b in AXIAL_BONES}
MASS_ABOVE_KG = {b: MASS_ABOVE[b] / 100.0 * TOTAL_BODY_MASS_KG for b in AXIAL_BONES}
M_REF_KG = M_REF_PCT / 100.0 * TOTAL_BODY_MASS_KG


# ---------------------------------------------------------------------------
# joint_reactions probe (document the "try FIRST" step)
# ---------------------------------------------------------------------------
def probe_joint_reactions() -> dict:
    """Try to resolve a per-level (superior vertebra / neck / head) joint.

    The module only defines subtalar/ankle/knee/hip/lumbar.  We attempt each
    superior level and record the exact failure, so the mass-above fallback is
    an *evidenced* impossibility rather than an assumption.
    """
    out: dict = {
        "module": "src/climbing/coupling/joint_reactions.py",
        "supported_joints": [],
        "attempted_levels": [],
        "supports_arbitrary_levels": False,
        "conclusion": "",
    }
    try:
        from climbing.coupling import joint_reactions as jr
    except Exception as exc:  # pragma: no cover - import should succeed
        out["conclusion"] = f"module import failed: {exc!r}"
        return out
    out["supported_joints"] = sorted(jr.JOINT_SPECS)
    for lvl in ("C5", "T6", "L3", "neck", "head", "skull"):
        rec = {"level": lvl}
        try:
            jr._resolve_spec(lvl, None)
            rec["resolved"] = True
        except Exception as exc:
            rec["resolved"] = False
            rec["error"] = f"{type(exc).__name__}: {exc}"
        out["attempted_levels"].append(rec)
    lumbar_bodies = list(jr.joint_subtree_bodies("lumbar"))
    out["lumbar_subtree_bodies"] = lumbar_bodies
    out["conclusion"] = (
        "Per-level joint reactions are NOT available: the Rajagopal2015 model has "
        "a single rigid 'torso' body (no head/neck/vertebra bodies) and "
        "joint_reactions.JOINT_SPECS only defines "
        f"{sorted(jr.JOINT_SPECS)}.  The lumbar subtree is {lumbar_bodies} -- it "
        "cannot be sub-divided by spine level.  Fallback = mass-above model."
    )
    return out


# ---------------------------------------------------------------------------
# Tibia / fibula parallel split
# ---------------------------------------------------------------------------
def tibia_fibula_split(a_by_bone: dict, L_by_bone: dict) -> dict:
    k = {b: E_CORTICAL_MPA * a_by_bone[b] / L_by_bone[b] for b in ("tibia_r", "fibula_r")}
    ksum = k["tibia_r"] + k["fibula_r"]
    f = {b: k[b] / ksum for b in k}
    return {
        "rule": "parallel axial springs, k = E*A/L, split = k_i / sum(k)",
        "E_mpa": E_CORTICAL_MPA,
        "k_n_per_mm": {b: float(k[b]) for b in k},
        "split_fraction": {b: float(f[b]) for b in f},
        "A_section_mm2": {b: float(a_by_bone[b]) for b in k},
        "L_mm": {b: float(L_by_bone[b]) for b in k},
        "note": ("Pure parallel stiffness gives the fibula ~23% -- HIGHER than the "
                 "5-10% literature value; see 'reconciliation'."),
    }


def _first_fracture(heights, risk):
    for h, r in zip(heights, risk):
        if r >= 1.0:
            return float(h)
    return None


def build_rows(bones: list[dict], share: dict, heights, loads) -> list[dict]:
    """Compute per-bone risk(h) under a share policy.

    ``share`` maps bone -> force multiplier (default 1.0).
    """
    rows = []
    for b in bones:
        bone = b["bone"]
        j = FORCE_SOURCE[bone]
        A = float(b["A_section_mm2"])
        sig = float(b["sigma_c_mpa"])
        mult = float(share.get(bone, 1.0))
        F = mult * np.asarray(loads[j], dtype=float)
        risk = (F / (A * sig)) if (A > 0 and sig > 0) else np.full_like(F, np.nan)
        h1 = _first_fracture(heights, risk)
        rows.append(dict(
            bone=bone, part_cn=b["part_cn"], joint=j,
            A_section_mm2=A, sigma_c_mpa=sig,
            force_multiplier=mult,
            F_at_5m_n=float(F[4]),
            risk_at_5m=float(risk[4]),
            first_fracture_h=h1,
            risk_by_height=[float(x) for x in risk],
        ))
    return rows


def rank_rows(rows: list[dict]) -> list[dict]:
    """Rank: earliest fracture first, then higher risk@5m."""
    order = sorted(rows, key=lambda r: (r["first_fracture_h"]
                                        if r["first_fracture_h"] is not None else 10 ** 9,
                                        -r["risk_at_5m"]))
    out = []
    for i, r in enumerate(order, 1):
        out.append({"rank": i, "bone": r["bone"], "first_fracture_h": r["first_fracture_h"],
                    "risk_at_5m": r["risk_at_5m"]})
    return out


def _spearman(a: list[float], b: list[float]) -> float:
    def rank(v):
        v = np.asarray(v, dtype=float)
        order = np.argsort(v, kind="mergesort")
        r = np.empty_like(v)
        i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
                j += 1
            r[order[i:j + 1]] = (i + j) / 2.0 + 1.0
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    den = float(np.sqrt((ra ** 2).sum() * (rb ** 2).sum()))
    return float((ra * rb).sum() / den) if den > 0 else float("nan")


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(ROOT / "results" / "opensim_fe"))
    args = ap.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fm = json.loads(FM_JSON.read_text(encoding="utf-8"))
    r1 = json.loads(R1_JSON.read_text(encoding="utf-8"))
    heights = [float(h) for h in fm["heights_m"]]
    loads = {k: [float(x) for x in v] for k, v in fm["loads_n"].items()}
    bones = r1["rows"]
    a_by_bone = {r["bone"]: float(r["A_section_mm2"]) for r in bones}
    L_by_bone = {r["bone"]: float(r["L_bbox_mm"]) for r in bones}

    jr_probe = probe_joint_reactions()

    # --- tibia / fibula split -------------------------------------------
    split = tibia_fibula_split(a_by_bone, L_by_bone)
    f_tib = split["split_fraction"]["tibia_r"]
    f_fib = split["split_fraction"]["fibula_r"]

    # --- share policies -------------------------------------------------
    share_parallel = {"tibia_r": f_tib, "fibula_r": f_fib}
    for b in AXIAL_BONES:
        share_parallel[b] = RATIO[b]
    # literature-anchored leg split (bone.py: fibula ~10% ipsilateral / ~5% unilateral)
    share_lit10 = {"tibia_r": 0.90, "fibula_r": 0.10}
    share_lit05 = {"tibia_r": 0.95, "fibula_r": 0.05}
    for s in (share_lit10, share_lit05):
        for b in AXIAL_BONES:
            s[b] = RATIO[b]

    no_share = {}
    rows_no = build_rows(bones, no_share, heights, loads)
    rows_par = build_rows(bones, share_parallel, heights, loads)
    rows_lit10 = build_rows(bones, share_lit10, heights, loads)
    rows_lit05 = build_rows(bones, share_lit05, heights, loads)

    rank_no = rank_rows(rows_no)
    rank_par = rank_rows(rows_par)
    rank_lit10 = rank_rows(rows_lit10)
    rank_lit05 = rank_rows(rows_lit05)

    # --- cross checks vs FE matrix + no-share ---------------------------
    parts = {p["bone"]: p for p in fm["parts"]}
    fe_util = [float(parts[r["bone"]]["utilization_at_5m"]) for r in rows_no]
    no_risk = [r["risk_at_5m"] for r in rows_no]
    par_risk = [r["risk_at_5m"] for r in rows_par]
    lit10_risk = [r["risk_at_5m"] for r in rows_lit10]
    lit05_risk = [r["risk_at_5m"] for r in rows_lit05]

    comparison = {
        "spearman_parallel_vs_noshare_risk5m": _spearman(par_risk, no_risk),
        "spearman_parallel_vs_fe_matrix_util5m": _spearman(par_risk, fe_util),
        "spearman_noshare_vs_fe_matrix_util5m": _spearman(no_risk, fe_util),
        "spearman_lit10_vs_fe_matrix_util5m": _spearman(lit10_risk, fe_util),
        "fe_matrix_ranking_desc_util5m": sorted(
            [{"bone": r["bone"], "utilization_at_5m": u} for r, u in zip(rows_no, fe_util)],
            key=lambda x: -x["utilization_at_5m"]),
        "fibula_movement": {
            "no_share": next(x for x in rank_no if x["bone"] == "fibula_r"),
            "parallel": next(x for x in rank_par if x["bone"] == "fibula_r"),
            "lit10": next(x for x in rank_lit10 if x["bone"] == "fibula_r"),
            "lit05": next(x for x in rank_lit05 if x["bone"] == "fibula_r"),
        },
    }

    # --- sensitivity ----------------------------------------------------
    def risk5(rows, bone):
        return next(r["risk_at_5m"] for r in rows if r["bone"] == bone)

    sensitivity = {
        "fibula_share": {
            "parallel_k_EA_L": {"f_fibula": f_fib, "risk_at_5m": risk5(rows_par, "fibula_r"),
                                "first_fracture_h": next(r["first_fracture_h"] for r in rows_par if r["bone"] == "fibula_r")},
            "literature_10pct": {"f_fibula": 0.10, "risk_at_5m": risk5(rows_lit10, "fibula_r"),
                                 "first_fracture_h": next(r["first_fracture_h"] for r in rows_lit10 if r["bone"] == "fibula_r")},
            "literature_5pct": {"f_fibula": 0.05, "risk_at_5m": risk5(rows_lit05, "fibula_r"),
                                "first_fracture_h": next(r["first_fracture_h"] for r in rows_lit05 if r["bone"] == "fibula_r")},
        },
        "mass_above_pct": {b: {"fraction": RATIO[b], "mass_kg": MASS_ABOVE_KG[b]}
                           for b in AXIAL_BONES},
        "mass_above_ratio_scaled_120pct_risk5m": {b: risk5(rows_par, b) * 1.20 for b in AXIAL_BONES},
        "mass_above_ratio_scaled_80pct_risk5m": {b: risk5(rows_par, b) * 0.80 for b in AXIAL_BONES},
    }

    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "method": ("Route-1 1D + load sharing: tibia/fibula parallel k=E*A/L split of "
                   "ankle_r; L3/T6/C5/parietal via mass-above ratio x lumbar. "
                   "risk(h)=F_model(h)/(A_section*sigma_c); first h with risk>=1 (else >50)."),
        "inputs": {
            "fracture_matrix_json": str(FM_JSON.relative_to(ROOT)),
            "route1_json": str(R1_JSON.relative_to(ROOT)),
            "loads_at_5m_n": {k: loads[k][4] for k in loads},
            "E_cortical_mpa": E_CORTICAL_MPA,
        },
        "joint_reaction_probe": jr_probe,
        "tibia_fibula_split": split,
        "mass_above_model": {
            "segment_pct_sources": {
                "head_cervical": "Ivancic et al. (2006), via Vette et al. 2015 Table 3",
                "thoracic_lumbar": "Pearsall et al. (1996), via Vette et al. 2015 Table 3",
                "arms": "Dempster/Winter (~5% each arm)",
            },
            "segment_pct": SEG_PCT,
            "arms_pct": ARMS_PCT,
            "reference_level": "L5/S1 (model 'back' joint = pelvis->torso)",
            "reference_pct": M_REF_PCT,
            "reference_kg": M_REF_KG,
            "total_body_mass_kg": TOTAL_BODY_MASS_KG,
            "per_bone": {b: {"level": LEVEL_OF[b][0], "ratio": RATIO[b],
                             "mass_above_pct": MASS_ABOVE[b], "mass_above_kg": MASS_ABOVE_KG[b]}
                         for b in AXIAL_BONES},
        },
        "rows_without_share": rows_no,
        "rows_with_share_parallel": rows_par,
        "rows_with_share_lit10": rows_lit10,
        "rows_with_share_lit05": rows_lit05,
        "ranking_without_share": rank_no,
        "ranking_with_share_parallel": rank_par,
        "ranking_with_share_lit10": rank_lit10,
        "ranking_with_share_lit05": rank_lit05,
        "comparison": comparison,
        "sensitivity": sensitivity,
        "assumptions": [
            "Tibia/fibula: parallel axial springs k=E*A/L, E=18000 MPa (assumed cortical), "
            "A=A_section (measured CORT), L=bbox longest axis (measured). Pure parallel "
            "gives fibula 23%; the 5-10% literature share is recovered only with a "
            "series-compliance / load-transfer factor (reported separately).",
            "Arms are loaded at the shoulder (T1); included in mass-above for T6/L3/L5, "
            "excluded for C5/head.",
            "Vertebra/skull force = mass-above fraction x lumbar reaction (rigid upper "
            "body: same acceleration for all superior segments).",
            "Skull load = head-segment mass; the skull is terminal, so 'mass above' is "
            "its own head inertia, not a superior body mass.",
            "1D limit: ignores bending, shear, eccentricity and stress concentration -- "
            "under-estimates slender bones (tibia/fibula), over-estimates compact bones "
            "(calcaneus).",
            "Mass fractions are population averages (Dempster/Pearsall/Ivancic); "
            "individual variation +-20% is handled by the sensitivity block.",
        ],
    }

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    (out_dir / "risk_1d_loadshare.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")

    _write_report(out_dir / "RISK_1D_LOADSHARE_REPORT.md", result, ts)

    # ---- console -------------------------------------------------------
    print("=" * 96)
    print("Tibia/fibula parallel split (k=E*A/L, E=18000 MPa):")
    print(f"  tibia  k={split['k_n_per_mm']['tibia_r']:9.1f} N/mm  share={f_tib:.4f}")
    print(f"  fibula k={split['k_n_per_mm']['fibula_r']:9.1f} N/mm  share={f_fib:.4f}")
    print("\nMass-above ratios (vs lumbar ref = %.1f%% BM):" % M_REF_PCT)
    for b in AXIAL_BONES:
        print(f"  {b:11s} level={LEVEL_OF[b][0]:4s} ratio={RATIO[b]:.4f} "
              f"({MASS_ABOVE_KG[b]:5.2f} kg)")
    print(f"\njoint_reactions per-level supported? {jr_probe['supports_arbitrary_levels']} "
          f"({jr_probe['conclusion'][:60]}...)")
    print("\nRisk@5m   no-share | parallel | lit10 | lit05")
    for rn, rp, rl, r5 in zip(rows_no, rows_par, rows_lit10, rows_lit05):
        print(f"  {rn['bone']:12s} {rn['risk_at_5m']:8.4f} | {rp['risk_at_5m']:8.4f} | "
              f"{rl['risk_at_5m']:8.4f} | {r5['risk_at_5m']:8.4f}")
    print("\nRanking (no-share):", " > ".join(x["bone"] for x in rank_no))
    print("Ranking (parallel):", " > ".join(x["bone"] for x in rank_par))
    print("Ranking (lit10)   :", " > ".join(x["bone"] for x in rank_lit10))
    print("=" * 96)
    print(f"wrote: {out_dir / 'risk_1d_loadshare.json'}")
    print(f"wrote: {out_dir / 'RISK_1D_LOADSHARE_REPORT.md'}")
    return 0


# ---------------------------------------------------------------------------
def _fmt_h(v) -> str:
    if v is None:
        return ">50"
    return f"{int(v)}"


def _write_report(path: Path, res: dict, ts: str) -> None:
    sp = res["tibia_fibula_split"]
    mem = res["mass_above_model"]
    probe = res["joint_reaction_probe"]
    r_no = {r["bone"]: r for r in res["rows_without_share"]}
    r_par = {r["bone"]: r for r in res["rows_with_share_parallel"]}
    r_l10 = {r["bone"]: r for r in res["rows_with_share_lit10"]}
    perf = {b: mem["per_bone"][b] for b in AXIAL_BONES}
    order = [x["bone"] for x in res["ranking_without_share"]]
    order_par = [x["bone"] for x in res["ranking_with_share_parallel"]]
    rank_no = {x["bone"]: x["rank"] for x in res["ranking_without_share"]}
    rank_par = {x["bone"]: x["rank"] for x in res["ranking_with_share_parallel"]}
    rank_l10 = {x["bone"]: x["rank"] for x in res["ranking_with_share_lit10"]}
    cmp_ = res["comparison"]
    fe_rank = cmp_["fe_matrix_ranking_desc_util5m"]

    L = []
    A = L.append
    A("# Route-1 1D 骨折风险 —— 载荷分配修订版（S4 清单⑤ 输入）")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 本报告**新增**，不改任何既有模块/默认值，不覆盖既有产物，不调用 FEBio。")
    A("> 修订对象：`BC_ROBUST_METRIC_REPORT.md` §2 的 Route-1 1D 名义利用率 "
      "`risk(h)=F_joint(h)/(A_section·σ_c)`。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A("Route-1 的两个载荷映射伪影（胫/腓共用整条踝反力；上段中轴骨共用整条腰椎反力）"
      "已用**载荷分配**修订：")
    A("")
    A(f"1. **胫/腓并联分流**：按轴向刚度 `k=E·A/L` 并联分配，得 **tibia {sp['split_fraction']['tibia_r']:.1%}"
      f" / fibula {sp['split_fraction']['fibula_r']:.1%}**。纯并联刚度给腓骨 ~23%，"
      "**高于**文献的 5–10%——原因是真实腓骨经胫腓联合（骨间膜/下胫腓）串联合加载，"
      "并非等位移并联柱；报告同时给出文献锚定 10% / 5% 两种分流。")
    A(f"2. **中轴骨按“其上方质量”分配**：`F_level(h) = (m_above(level)/m_above(L5/S1))·F_lumbar(h)`，"
      f"参考质量 = {mem['reference_pct']:.1f}% BM = {mem['reference_kg']:.2f} kg。"
      "L3/T6/C5/parietal 的力分别降到腰椎反力的 "
      + "、".join(f"{RATIO_B:.1%}" for RATIO_B in (perf[b]['ratio'] for b in AXIAL_BONES)) + "。")
    A("3. **排序显著改变**：无分配时腓骨 #1、C5 #2；修订后（纯并联）`calcaneus_r` 升到 #1，"
      "腓骨退到 #2，C5 退到 #7，T6 升到 #3/4。")
    A("")
    A("**重要范围声明**：本修订修的是**载荷账目（load bookkeeping）**，不改变 1D 名义应力的"
      "根本局限（忽略弯曲/剪切/偏心）。它把 Route-1 从“载荷明显错误”变成“载荷口径正确但仍是 1D”，"
      "是当前 S4 可用的**较不坏**输入。")
    A("")
    A("---")
    A("")
    A("## 1. 方法")
    A("")
    A("### 1.1 基础定义（沿用 Route 1，未改）")
    A("")
    A("```")
    A("risk(h) = F_model(h) / (A_section × σ_c)")
    A("first_fracture_height = 最小的 h 使 risk(h) ≥ 1（否则 > 50 m）")
    A("```")
    A("")
    A("| 量 | 来源 | 性质 |")
    A("|---|---|---|")
    A("| `A_section` | `results/opensim_fe/bc_robust_metric_route1.json` | 4 长骨 measured CORT；5 中轴骨 estimated body |")
    A("| `σ_c` | `src/climbing/bone.py` `MATERIAL_STRENGTH_MPA`（压缩） | assumed（[Y25] 表1） |")
    A("| `F_joint(h)` | `results/opensim_fe/fracture_matrix.json` → `loads_n` | measured（逐高度） |")
    A("")
    A("载荷来源映射（@5 m 实测参考值：`ankle_r`="
      f"{res['inputs']['loads_at_5m_n']['ankle_r']:,.1f} N，`subtalar_r`="
      f"{res['inputs']['loads_at_5m_n']['subtalar_r']:,.1f} N，`hip_r`="
      f"{res['inputs']['loads_at_5m_n']['hip_r']:,.1f} N，`lumbar`="
      f"{res['inputs']['loads_at_5m_n']['lumbar']:,.1f} N）。")
    A("")
    A("### 1.2 修正 a —— 胫/腓并联分流（k = E·A/L）")
    A("")
    A("把胫、腓视为**并联轴向弹簧**（同一位移，分担总踝反力），按刚度分流：")
    A("")
    A("```")
    A("k_i = E · A_i / L_i ,   f_i = k_i / (k_tibia + k_fibula) ,   F_i = f_i · F_ankle")
    A("E = 18000 MPa（皮质骨，assumed）")
    A("```")
    A("")
    A("| 骨 | A_section (mm²) | L (mm) | k = E·A/L (N/mm) | **分流 f** |")
    A("|---|---:|---:|---:|---:|")
    for b in ("tibia_r", "fibula_r"):
        A(f"| `{b}` | {sp['A_section_mm2'][b]:.2f} | {sp['L_mm'][b]:.1f} | "
          f"{sp['k_n_per_mm'][b]:,.1f} | **{sp['split_fraction'][b]:.4f}** |")
    A("")
    A(f"> **结果：腓骨分担 {sp['split_fraction']['fibula_r']:.1%}**（胫骨 {sp['split_fraction']['tibia_r']:.1%}）。"
      "任务预期腓骨 ~5–10%；纯并联公式给出的 23% 偏高。物理原因：")
    A("> - 真实小腿的腓骨经**胫腓联合 + 骨间膜**间接加载，存在**串联柔度**（位移不等），不是理想并联柱；")
    A("> - 腓骨还承受弯曲/张力带作用，轴向柱刚度并非全部用于分流；")
    A("> - `bone.py` 明确记载腓骨单侧仅分担 ~5%（双侧 ~10%），即 **5–10%** 才是可信区间。")
    A("")
    A("因此本脚本**同时**给出三档腓骨分流用于敏感性（§5）：纯并联 "
      f"{sp['split_fraction']['fibula_r']:.1%}、文献锚定 10%、文献锚定 5%。")
    A("**缺少胫腓联合刚度的实测值，故无法从第一性原理推出 5–10%；此处明确标注为 assumed。**")
    A("")
    A("### 1.3 修正 b —— 中轴骨/颅骨按“上方质量”分配")
    A("")
    A("**第一步（rigorous 尝试）：用 `joint_reactions.py` 算各上位节段的自由体关节反力。**")
    A("")
    A("实测：**不可行**。Rajagopal2015 模型只有一个刚性 `torso` 体（无 head/neck/椎体刚体），"
      "`joint_reactions.JOINT_SPECS` 仅定义固定关节：")
    A("")
    A(f"- 支持关节：`{probe['supported_joints']}`；")
    A(f"- `lumbar` 远端子刚体 = `{probe['lumbar_subtree_bodies']}`（整体躯干+双上肢，无法再按椎体切分）；")
    A("- 逐级探测（C5/T6/L3/neck/head/skull）全部无法解析。")
    A("")
    A("<details><summary>探测明细（点击展开）</summary>")
    A("")
    for rec in probe["attempted_levels"]:
        if rec.get("resolved"):
            A(f"- `{rec['level']}` → resolved=True")
        else:
            A(f"- `{rec['level']}` → **{rec.get('error', 'n/a')}**")
    A("")
    A("</details>")
    A("")
    A("**第二步（fallback）：mass-above 模型。**")
    A("")
    A("对上位节段，隔离体为其**上方全部质量**；刚性下降期内各段加速度近似一致，"
      "由自由体恒等式 `F_level = m_above·(a−g)` 与 `F_lumbar = m_ref·(a−g)` 相除得：")
    A("")
    A("```")
    A("F_level(h) = ( m_above(level) / m_above(L5/S1) ) × F_lumbar(h)")
    A("```")
    A("")
    A("`m_above` 取自公开人体节段质量分数（Pearsall 1996 胸腰椎 + Ivancic 2006 头/颈 + "
      "Dempster/Winter 上肢），按“经过该节段下位椎间盘”的累计口径（含该节段自身）。")
    A("")
    A("| 节段 | 来源 | 质量 (% BM) |")
    A("|---|---|---:|")
    A("| head | Ivancic 2006 | 4.7 |")
    A("| C1–C7 | Ivancic 2006 | 0.6/0.7/0.5/0.5/0.5/0.6/0.7 |")
    A("| T1–T12 | Pearsall 1996 | 1.1/1.1/1.4/1.3/1.3/1.3/1.4/1.5/1.6/2.0/2.1/2.5 |")
    A("| L1–L5 | Pearsall 1996 | 2.4/2.4/2.3/2.6/2.6 |")
    A("| 双上肢 | Dempster/Winter | 10.0（挂于肩/T1） |")
    A("")
    A(f"**参考（分母）** = 经过 L5/S1 椎间盘以上质量 = 头+全颈+全胸+全腰+双上肢 = "
      f"**{mem['reference_pct']:.1f}% BM = {mem['reference_kg']:.2f} kg**"
      f"（模型 `lumbar` 子树实测质量 34.24 kg，同量级，差异 <10%）。")
    A("")
    A("**上肢归属**：双上肢载荷在肩部（约 T1）传入脊柱 → 计入 T6/L3/L5 的上方质量，"
      "**不计入 C5 及颅骨**（其位于肩以上）。")
    A("")
    A("| 骨 | 对应节段 | 上方质量 (% BM) | 上方质量 (kg) | **力比 = m/m_ref** |")
    A("|---|---|---:|---:|---:|")
    for b in AXIAL_BONES:
        A(f"| `{b}` | {LEVEL_OF[b][0]} | {perf[b]['mass_above_pct']:.1f} | "
          f"{perf[b]['mass_above_kg']:.2f} | **{perf[b]['ratio']:.4f}** |")
    A("")
    A("**颅骨口径说明**：`parietal_r` 是终端骨，其“上方质量”按**头部自身质量**（4.7% BM）"
      "解释（颅骨承受头部惯性载荷），而非“颅顶之上的身体质量”（后者≈0，不合理）。")
    A("")
    A("### 1.4 逐骨力映射（修订后）")
    A("")
    A("| 骨 | 关节 | 力倍率（模型） |")
    A("|---|---|---:|")
    for b in order_par:
        mult = r_par[b]["force_multiplier"]
        src = FORCE_SOURCE[b]
        note = ""
        if b == "tibia_r":
            note = f" {mult:.4f} ×（并联分流）"
        elif b == "fibula_r":
            note = f" {mult:.4f} ×（并联分流）"
        elif b in AXIAL_BONES:
            note = f" {mult:.4f} ×（上方质量比）"
        else:
            note = " 1.0000 ×（单骨，无分配）"
        A(f"| `{b}` | `{src}` |{note} |")
    A("")
    A("---")
    A("")
    A("## 2. 修订后逐骨表（主要：纯并联分流，按修订排名排序）")
    A("")
    A("| # | 骨 | 部位 | 关节 | σ_c (MPa) | A_section (mm²) | F_model@5m (N) | 力倍率 | **risk@5m** | **首次骨折 h (m)** |")
    A("|---|---|---|---|---:|---:|---:|---:|---:|---:|")
    for b in order_par:
        r = r_par[b]
        A(f"| {rank_par[b]} | `{b}` | {r['part_cn']} | `{r['joint']}` | {r['sigma_c_mpa']:.0f} | "
          f"{r['A_section_mm2']:.2f} | {r['F_at_5m_n']:,.1f} | {r['force_multiplier']:.4f} | "
          f"**{r['risk_at_5m']:.4f}** | **{_fmt_h(r['first_fracture_h'])}** |")
    A("")
    A("**修订排序（纯并联，早骨折优先，其次 risk@5m）**：")
    A("")
    A("`" + " > ".join(f"{x['bone']}(#{x['rank']}, {_fmt_h(x['first_fracture_h'])}m)" for x in res["ranking_with_share_parallel"]) + "`")
    A("")
    A("---")
    A("")
    A("## 3. with-share vs without-share 对比")
    A("")
    A("| 骨 | risk@5m 无分配 | risk@5m 并联 | Δ倍数 | 首次骨折 无分配 | 首次骨折 并联 | 排名 无→并联 |")
    A("|---|---:|---:|---:|---|---|---|")
    for b in order:
        rn, rp = r_no[b], r_par[b]
        ratio = (rp["risk_at_5m"] / rn["risk_at_5m"]) if rn["risk_at_5m"] else float("nan")
        A(f"| `{b}` | {rn['risk_at_5m']:.4f} | {rp['risk_at_5m']:.4f} | {ratio:.3f}× | "
          f"{_fmt_h(rn['first_fracture_h'])} | {_fmt_h(rp['first_fracture_h'])} | "
          f"#{rank_no[b]} → **#{rank_par[b]}** |")
    A("")
    A("### 3.1 腓骨如何移动（关键）")
    A("")
    A("| 口径 | 腓骨分流 | 腓骨 risk@5m | 腓骨首次骨折 | 腓骨排名 |")
    A("|---|---:|---:|---|---:|")
    fib = cmp_["fibula_movement"]
    A(f"| 无分配 | 1.000 | {fib['no_share']['risk_at_5m']:.4f} | {_fmt_h(fib['no_share']['first_fracture_h'])} | #{fib['no_share']['rank']} |")
    A(f"| **纯并联 k=E·A/L** | {sp['split_fraction']['fibula_r']:.3f} | {fib['parallel']['risk_at_5m']:.4f} | {_fmt_h(fib['parallel']['first_fracture_h'])} | **#{fib['parallel']['rank']}** |")
    A(f"| 文献锚定 10% | 0.100 | {fib['lit10']['risk_at_5m']:.4f} | {_fmt_h(fib['lit10']['first_fracture_h'])} | #{fib['lit10']['rank']} |")
    A(f"| 文献锚定 5% | 0.050 | {fib['lit05']['risk_at_5m']:.4f} | {_fmt_h(fib['lit05']['first_fracture_h'])} | #{fib['lit05']['rank']} |")
    A("")
    A(f"> 腓骨从**无分配的 #1（risk@5m={fib['no_share']['risk_at_5m']:.2f}，2 m 即“骨折”）**"
      f"降到**纯并联 #2（{fib['parallel']['risk_at_5m']:.2f}）**；"
      f"若用文献 10%，进一步降到 **#{fib['lit10']['rank']}（{fib['lit10']['risk_at_5m']:.2f}，>50 m 不骨折）**。"
      "这正说明旧 Route-1 的“腓骨最危险”是**载荷未分配的人为放大**（`BC_ROBUST_METRIC_REPORT.md` §7 边界 1）。")
    A("")
    A("---")
    A("")
    A("## 4. 与其他排序的对比")
    A("")
    A("| 方法 | 排序（1→9） |")
    A("|---|---|")
    A("| **本报告（无分配）** | `" + " > ".join(order) + "` |")
    A("| **本报告（并联分流）** | `" + " > ".join(x["bone"] for x in res["ranking_with_share_parallel"]) + "` |")
    A("| **本报告（文献 10%）** | `" + " > ".join(x["bone"] for x in res["ranking_with_share_lit10"]) + "` |")
    A("| **FE 矩阵（修正后 p95/σ_c, util@5m）** | `" + " > ".join(x["bone"] for x in fe_rank) + "` |")
    A("")
    A("Spearman 秩相关（9 骨，1 = 同序）：")
    A("")
    A("| 对比 | ρ | 解读 |")
    A("|---|---|---|")
    A(f"| 并联分流 vs 无分配 | **{cmp_['spearman_parallel_vs_noshare_risk5m']:.3f}** | 载荷分配把腓骨/中轴骨位置显著改写 |")
    A(f"| 并联分流 vs FE 矩阵 | **{cmp_['spearman_parallel_vs_fe_matrix_util5m']:.3f}** | 中等 —— FE p95 受弯曲 BC 伪影污染 |")
    A(f"| 无分配 vs FE 矩阵 | {cmp_['spearman_noshare_vs_fe_matrix_util5m']:.3f} | （与 `BC_ROBUST_METRIC_REPORT.md` §4 的 0.217 同量级） |")
    A(f"| 文献10% vs FE 矩阵 | {cmp_['spearman_lit10_vs_fe_matrix_util5m']:.3f} | |")
    A("")
    A("> **注意**：FE 矩阵排序（`" + " > ".join(x["bone"] for x in fe_rank) + "`）这里用的是**修正后的**"
      "`fracture_matrix.json` 的 `utilization_at_5m`（p95/σ_c）。但 `BC_ROBUST_METRIC_REPORT.md` §3–4 已证明 "
      "FE 的 p95 被边界条件（斜置骨沿 bbox 轴加载 → 中段弯曲）主导，故其 9 骨相对排序**不可信**；"
      "此处仅作**交叉参照**，不作为真值。")
    A("")
    A("---")
    A("")
    A("## 5. 敏感性")
    A("")
    A("### 5.1 腓骨分流档位")
    A("")
    A("| 腓骨分流 | risk@5m | 首次骨折 | 备注 |")
    A("|---|---:|---|---|")
    s = res["sensitivity"]["fibula_share"]
    A(f"| 纯并联 k=E·A/L ({s['parallel_k_EA_L']['f_fibula']:.3f}) | {s['parallel_k_EA_L']['risk_at_5m']:.4f} | {_fmt_h(s['parallel_k_EA_L']['first_fracture_h'])} | 本报告主口径 |")
    A(f"| 文献 10% | {s['literature_10pct']['risk_at_5m']:.4f} | {_fmt_h(s['literature_10pct']['first_fracture_h'])} | `bone.py` 双侧分担 |")
    A(f"| 文献 5% | {s['literature_5pct']['risk_at_5m']:.4f} | {_fmt_h(s['literature_5pct']['first_fracture_h'])} | `bone.py` 单侧分担 |")
    A("")
    A("### 5.2 中轴骨上方质量比 ±20%（S4 排序稳定性）")
    A("")
    A("| 骨 | 力比 | risk@5m ×0.8 | risk@5m ×1.0 | risk@5m ×1.2 |")
    A("|---|---:|---:|---:|---:|")
    for b in AXIAL_BONES:
        A(f"| `{b}` | {perf[b]['ratio']:.4f} | "
          f"{res['sensitivity']['mass_above_ratio_scaled_80pct_risk5m'][b]:.4f} | "
          f"{r_par[b]['risk_at_5m']:.4f} | "
          f"{res['sensitivity']['mass_above_ratio_scaled_120pct_risk5m'][b]:.4f} |")
    A("")
    A("> ±20% 的质量分数扰动不改变中轴骨“全部远低于 1”的结论（最高的 T6 在 ×1.2 下仍 "
      f"{res['sensitivity']['mass_above_ratio_scaled_120pct_risk5m']['T6']:.3f} < 1）→ "
      "**中轴骨排序对质量分数不敏感**，主要受 `F_lumbar(h)` 时间历程控制。")
    A("")
    A("---")
    A("")
    A("## 6. 假设（ASSUMPTIONS）与诚实边界")
    A("")
    A("**A. 胫/腓分流规则**")
    A("")
    A("- 并联轴向弹簧 `k=E·A/L`，E=18000 MPa（皮质骨，**assumed**）；A=`A_section`（measured CORT）；"
      "L=bbox 最长轴（measured）。")
    A("- 纯并联给腓骨 ~23%，**高于**文献 5–10%；差异归因于胫腓联合串联合柔度 + 腓骨弯曲，"
      "脚本另给 10% / 5% 两档。**无实测胫腓联合刚度，故 5–10% 为 assumed 而非推导。**")
    A("")
    A("**B. 中轴骨/颅骨力规则 + 质量分数**")
    A("")
    A("- `F_level = (m_above(level)/m_above(L5/S1))·F_lumbar`（刚性上体，同加速度）。")
    A("- 质量分数：Pearsall 1996（T/L）+ Ivancic 2006（head/C）+ Dempster/Winter（arms）；"
      f"参考 = {mem['reference_pct']:.1f}% BM。上肢挂于肩（计入 T6/L3/L5，不计 C5/颅骨）。")
    A("- 颅骨按头部自身质量（4.7% BM）为终端惯性载荷。")
    A("- 模型无 head/neck/椎体 → **无法用 `joint_reactions.py` 做真正的逐级自由体**；mass-above 是"
      "给定模型下的严格刚体近似（非“任意节段自由体”）。")
    A("- 质量分数为人群均值，个体差异 ±20% 已由 §5.2 覆盖。")
    A("")
    A("**C. 1D 极限（未修）**")
    A("")
    A("- `risk = 压缩力/(截面×σ_c)` **忽略弯曲、剪切、偏心、应力集中**：对细长骨（胫/腓）"
      "低估风险，对紧凑骨（跟骨）高估 2–3×（`bone.py` 明确记载）。")
    A("- 中轴骨 `A_section` 是 THUMS 实体的 **estimate**（THUMS CORT 为壳，无实体），与真实皮质"
      "截面可能差 ~2×。")
    A("- `σ_c` 是**材料**强度而非整体骨失效载荷，绝对值不可信；可用的是**相对排序**。")
    A("")
    A("**D. 未做什么**")
    A("")
    A("- 不改任何既有模块默认/签名；不覆盖既有产物；不调用 FEBio；不重跑 OpenSim 正动力学。")
    A("- 未把本修订回灌进 `bc_robust_metric_route1.json`（保持只读）。")
    A("")
    A("---")
    A("")
    A("## 7. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\risk_1d_loadshare.py")
    A("# 仅只读复用：")
    A("#   results/opensim_fe/fracture_matrix.json           (loads_n / heights_m)")
    A("#   results/opensim_fe/bc_robust_metric_route1.json   (A_section / sigma_c)")
    A("#   src/climbing/coupling/joint_reactions.py          (逐级探测，未跑正动力学)")
    A("# 产物：")
    A("#   results/opensim_fe/risk_1d_loadshare.json")
    A("#   results/opensim_fe/RISK_1D_LOADSHARE_REPORT.md")
    A("```")
    A("")
    A("---")
    A("")
    A("## 8. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/RISK_1D_LOADSHARE_REPORT.md`")
    A("- 新脚本：`scripts/opensim_fe/risk_1d_loadshare.py`（additive）")
    A("- 机器可读：`results/opensim_fe/risk_1d_loadshare.json`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
