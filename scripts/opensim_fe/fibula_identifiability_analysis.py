"""fibula_identifiability_analysis.py -- 腓骨 all_one 不可辨识的**只读根因取证**。

背景
----
G3（`scripts/opensim_fe/g3_subject_variability.py`，2026-10-10）把 S4 的确定性
50×9 骨折矩阵升级为被试级概率矩阵；A3 判定 **PASS 2/3**（骨盆/颅骨 OR ∈ 靶 ±30%），
唯一缺口是**腓骨 all_one**（logistic OR 无定义）。设计 §5 记为「G3 实现限制（数据问题）」。
本脚本不重跑任何 FE / 不重跑 G3 MC / 不改任何既有文件，只做**取证**：

诊断两个候选根因
----------------
* **σc 口径问题**：route-乙 腓骨 σc=70 MPa（[Y25] Table 1 腓骨两端**压缩**强度）
  是否太低 / 基准错？
* **FE 量级问题**：腓骨 FE 的 σ_nom 是否异常（单位错 / 载荷未分配 / 本身远距强度）？

方法（全部为既有产物上的**确定性算术**，无随机、无 FE 调用）
-----------------------------------------------------------
1. σ_nom(h) 取 `fracture_matrix_bc2.json::sigma_vm_mpa["median"][bone]`（= G3 的输入）。
2. σc 取 route-乙 弱区口径（`paper_alignment.REGIONS` ← `bone.py::MATERIAL_STRENGTH_MPA`）。
3. σc 灵敏性：腓骨在多大 σc 以下仍是 all_one？（临界 σc* = min_h σ_nom(h)）。
4. 姊妹对照：胫骨（同关节 `ankle_r`、同 σc=70、同 fe_solve/lrf）为何不 all_one？
5. 载荷分配反事实：把腓骨 FE 载荷乘以其生理分担比（5% / 10%，`bone.py` /
   `risk_1d_loadshare.py` 文献锚）后，是否出现「0/1 混合」使 logistic 有定义？
   —— **仅诊断，不写回、不调参凑靶**。

输入（只读）
------------
* `results/opensim_fe/fracture_matrix_bc2.json`      （BC2 σ_nom、parts、loads_n）
* `results/opensim_fe/fracture_matrix.json`          （legacy 骨干 σc 口径对照）
* `results/opensim_fe/s4_alignment_bc2.json`         （route-乙 ends/shaft 子区域统计）
* `results/opensim_fe/bc_robust_metric_route1.json`  （逐骨截面 A_section_mm2）
* `results/opensim_fe/g3_probability_matrix.json`    （G3 主设定 p[h][site]）
* `temp/fe_out/fracture_matrix_bc2/<b>/<b>_fe.json`（fe_solve/lrf/linearity）

产物（**新文件**）
------------------
* `results/opensim_fe/fibula_identifiability_analysis.json`（机器可读）
* stdout 摘要

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\fibula_identifiability_analysis.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import TypedDict

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.opensim_fe import paper_alignment as pa  # noqa: E402  (REGIONS / _sc)

RES = ROOT / "results" / "opensim_fe"
BC2_JSON = RES / "fracture_matrix_bc2.json"
LEGACY_JSON = RES / "fracture_matrix.json"
S4_JSON = RES / "s4_alignment_bc2.json"
R1_JSON = RES / "bc_robust_metric_route1.json"
G3_PROB = RES / "g3_probability_matrix.json"
BC2_DIR = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
OUT_JSON = RES / "fibula_identifiability_analysis.json"

#: 腓骨生理轴向载荷分担比（文献锚；`bone.py` fibula_ends.note / `risk_1d_loadshare.py`）。
FIBULA_SHARE_LIT = {"unilateral_05": 0.05, "ipsilateral_10": 0.10}
#: 「并行刚度」分担比（`risk_1d_loadshare.py::tibia_fibula_split`，非文献锚、仅对照）。
FIBULA_SHARE_PARALLEL = 0.23234546642114307


class Diff(TypedDict):
    bone: str
    part_cn: str
    joint: str
    sigma_c_weak_mpa: float
    sigma_c_shaft_mpa: float | None
    sigma_nom_h1_mpa: float
    sigma_nom_h5_mpa: float
    sigma_nom_h50_mpa: float
    ratio_h1: float
    ratio_h5: float
    ratio_h50: float
    ratio_min: float
    ratio_max: float
    first_fracture_h: int | None
    n_fractured_h: int
    all_one: bool


def _ratio_fields(sig: np.ndarray, sc: float) -> dict:
    r = np.asarray(sig, float) / float(sc)
    over = int((r >= 1.0).sum())
    first = next((i + 1 for i, v in enumerate(r) if v >= 1.0), None)
    return {
        "ratio_min": float(r.min()), "ratio_max": float(r.max()),
        "ratio_h1": float(r[0]), "ratio_h5": float(r[4]), "ratio_h50": float(r[-1]),
        "n_fractured_h": over, "first_fracture_h": first,
        "all_one": bool((r >= 1.0).all()),
    }


def _bone_diff(cn: str, bone: str, joint: str, subs: list[tuple[str, str]],
               sigma_bc2: dict[str, list[float]]) -> Diff:
    weak = pa._sc(subs[0][1])
    shaft = pa._sc(subs[1][1]) if len(subs) > 1 else None
    sig = np.asarray(sigma_bc2[bone], float)
    out: Diff = {
        "bone": bone, "part_cn": cn, "joint": joint,
        "sigma_c_weak_mpa": weak, "sigma_c_shaft_mpa": shaft,
        "sigma_nom_h1_mpa": float(sig[0]), "sigma_nom_h5_mpa": float(sig[4]),
        "sigma_nom_h50_mpa": float(sig[-1]),
    }
    out.update(_ratio_fields(sig, weak))  # type: ignore[typeddict-item]
    return out


def _sigma_c_sensitivity(sig: np.ndarray) -> dict:
    """腓骨 all_one 的 σc 临界值 + 若干候选 σc 下的判据（诊断用）。"""
    sig = np.asarray(sig, float)
    crit = float(sig.min())  # 恰在最小 σ_nom 之下即 h=1 安全 → 不再是 all_one
    cand = [70.0, 80.0, 160.0, 200.0, 220.0, crit, crit + 1.0, 500.0]
    rows = []
    for sc in cand:
        over = int((sig >= sc).sum())
        rows.append({"sigma_c_mpa": round(sc, 4), "n_fractured_h": over,
                     "all_one": bool(over == sig.size),
                     "first_safe_h": next((i + 1 for i, v in enumerate(sig)
                                           if v < sc), None)})
    return {"sigma_nom_min_mpa": crit, "all_one_critical_sigma_c_mpa": crit,
            "candidates": rows}


def _share_counterfactual(sig: np.ndarray, sc: float) -> dict:
    """把腓骨 σ_nom 乘以其生理分担比后的可辨识性（诊断；不写回）。

    * ``share_h1_safe``：低于此值，h=1 落到阈值下 → 出现 0/1 混合（logistic 有定义）。
    * ``share_all_safe``：低于此值，全高度都安全 → all_zero。
    """
    sig = np.asarray(sig, float)
    share_h1_safe = float(sc / sig.min())
    share_all_safe = float(sc / sig.max())
    out = {"sigma_c_mpa": sc, "share_h1_safe": share_h1_safe,
           "share_all_safe": share_all_safe, "cases": []}
    shares = {"no_share_1.0": 1.0, **FIBULA_SHARE_LIT,
              "parallel_stiffness": FIBULA_SHARE_PARALLEL}
    for name, s in shares.items():
        r = sig * float(s) / sc
        over = int((r >= 1.0).sum())
        mixed = bool(0 < over < r.size)
        out["cases"].append({
            "case": name, "share": float(s),
            "ratio_h1": float(r[0]), "ratio_h50": float(r[-1]),
            "n_fractured_h": over, "all_one": bool(over == r.size),
            "all_zero": bool(over == 0), "identifiable_0_1_mixed": mixed,
            "first_fracture_h": next((i + 1 for i, v in enumerate(r)
                                      if v >= 1.0), None),
        })
    return out


def _fe_provenance(bone: str) -> dict:
    p = BC2_DIR / bone / f"{bone}_fe.json"
    rec = json.loads(p.read_text(encoding="utf-8"))
    mesh = rec.get("mesh", {})
    return {
        "bone": bone, "status": rec.get("status"),
        "ref_load_n": rec.get("ref_load_n"),
        "fe_solve_load_n": rec.get("fe_solve_load_n"),
        "load_reduction_factor": rec.get("load_reduction_factor"),
        "n_nodes": mesh.get("n_nodes"), "by_solid": mesh.get("by_solid"),
        "linearity": rec.get("linearity"),
        "attempts": [(a.get("load_n"), a.get("ok")) for a in rec.get("attempts", [])],
    }


def build() -> dict:
    bc2 = json.loads(BC2_JSON.read_text(encoding="utf-8"))
    legacy = json.loads(LEGACY_JSON.read_text(encoding="utf-8"))
    s4 = json.loads(S4_JSON.read_text(encoding="utf-8"))
    r1 = json.loads(R1_JSON.read_text(encoding="utf-8"))
    g3 = json.loads(G3_PROB.read_text(encoding="utf-8"))

    heights = np.asarray(bc2["heights_m"], float)
    loads = {k: np.asarray(v, float) for k, v in bc2["loads_n"].items()}
    sigma_bc2 = bc2["sigma_vm_mpa"]["median"]           # G3 的 σ_nom 来源
    a_sec = {row["bone"]: row.get("A_section_mm2") for row in r1["rows"]}

    # ---- 逐部位 σ_nom/σc（BC2 median 口径 = G3 输入）-------------------
    by_part = [_bone_diff(cn, bone, joint, subs, sigma_bc2)
               for cn, bone, joint, subs in pa.REGIONS]
    by_bone = {d["bone"]: d for d in by_part}
    fib = by_bone["fibula_r"]
    tib = by_bone["tibia_r"]

    # ---- 姊妹比：同关节、同 σc、同 fe_solve/lrf --------------------------
    fib_sig = np.asarray(sigma_bc2["fibula_r"], float)
    tib_sig = np.asarray(sigma_bc2["tibia_r"], float)
    sibling = {
        "joint_shared": "ankle_r (tibia_r & fibula_r both)",
        "ankle_load_n": {"h1": float(loads["ankle_r"][0]),
                         "h5": float(loads["ankle_r"][4]),
                         "h50": float(loads["ankle_r"][-1])},
        "sigma_c_mpa_both": fib["sigma_c_weak_mpa"],
        "A_section_mm2": {"tibia_r": a_sec.get("tibia_r"),
                          "fibula_r": a_sec.get("fibula_r")},
        "sig_nom_fibula_over_tibia": {
            "h1": float(fib_sig[0] / tib_sig[0]),
            "h5": float(fib_sig[4] / tib_sig[4]),
            "h50": float(fib_sig[-1] / tib_sig[-1]),
            "min": float((fib_sig / tib_sig).min()),
            "max": float((fib_sig / tib_sig).max()),
        },
    }

    # ---- route-乙 ends/shaft（A3 判定实际口径）对照 ---------------------
    det = s4["detail_yi"]
    fib_det = det["腓骨"]
    tib_det = det["胫骨"]
    l5 = float(loads["ankle_r"][4])
    scale_h = loads["ankle_r"] / float(fib_det["fe_solve_n"])
    route_yi = {}
    for cn, d, sc_e, sc_s in (("腓骨", fib_det, 70.0, 160.0),
                              ("胫骨", tib_det, 70.0, 200.0)):
        ends_h = float(d["ends_stat_mpa"]) * scale_h
        shaft_h = float(d["shaft_stat_mpa"]) * scale_h
        route_yi[cn] = {
            "ends_stat_mpa_at_solve": float(d["ends_stat_mpa"]),
            "shaft_stat_mpa_at_solve": float(d["shaft_stat_mpa"]),
            "ends_ratio": _ratio_fields(ends_h, sc_e),
            "shaft_ratio": _ratio_fields(shaft_h, sc_s),
            "first_ends_h": d["first_ends"], "first_shaft_h": d["first_shaft"],
        }

    # ---- σc 灵敏性 + 载荷分配反事实（腓骨）------------------------------
    sigma_c_sens = _sigma_c_sensitivity(fib_sig)
    share_cf = _share_counterfactual(fib_sig, fib["sigma_c_weak_mpa"])
    # 骨干口径对照：BC2 σ_nom vs σc=160（报告 §DoD-2 所述「另一 σc 口径」）
    shaft_cf = {"sigma_c_mpa": 160.0,
                "sigma_nom_source": "BC2 median whole-domain",
                **_ratio_fields(fib_sig, 160.0)}
    legacy_fib_fe = next(p for p in legacy["parts"] if p["bone"] == "fibula_r")

    # ---- G3 主设定 p（腓骨/胫骨）---------------------------------------
    p = g3["primary"]["p"]
    g3_fib = [float(v) for v in p["腓骨"]]
    g3_tib = [float(v) for v in p["胫骨"]]

    return {
        "meta": {
            "title": "Fibula all_one identifiability -- read-only forensic analysis",
            "inputs": {
                "sigma_nom_source": "fracture_matrix_bc2.json::sigma_vm_mpa[\"median\"][bone]",
                "sigma_c_source": "paper_alignment.REGIONS weak key <- bone.py::MATERIAL_STRENGTH_MPA",
                "heights_m": [int(h) for h in heights],
            },
            "notes": [
                "确定性只读取证；零 FE 调用；零随机；不写回任何既有产物。",
                "σ_nom/σc 为越阈物理量（非禁用类）；不含 σ/σ_law 比值。",
            ],
        },
        "per_part_sigma_nom_over_sigma_c": by_part,
        "fibula": fib, "tibia": tib,
        "sigma_c_sensitivity_fibula": sigma_c_sens,
        "load_share_counterfactual_fibula": share_cf,
        "fibula_shaft_sigma_c_160": shaft_cf,
        "sibling_tibia_vs_fibula": sibling,
        "route_yi_ends_shaft": route_yi,
        "fe_provenance": {"fibula_r": _fe_provenance("fibula_r"),
                          "tibia_r": _fe_provenance("tibia_r")},
        "legacy_fibula_fe": {
            "fe_solve_load_n": legacy_fib_fe.get("fe_solve_load_n"),
            "load_reduction_factor": legacy_fib_fe.get("load_reduction_factor"),
            "sigma_vm_ref_median_mpa": (legacy_fib_fe.get("sigma_vm_ref_mpa") or {}).get("median"),
            "fe_status": legacy_fib_fe.get("fe_status"),
        },
        "g3_primary_p": {"腓骨": g3_fib, "胫骨": g3_tib},
    }


def main() -> int:
    res = build()
    OUT_JSON.write_text(json.dumps(res, indent=2, ensure_ascii=False),
                        encoding="utf-8")

    fib = res["fibula"]
    tib = res["tibia"]
    print("=" * 78)
    print("腓骨 all_one 根因取证（只读；σ_nom = BC2 median；σc = route-乙 弱区）")
    print("=" * 78)
    print(f"\n[腓骨] σc={fib['sigma_c_weak_mpa']:.0f} MPa  "
          f"σ_nom: h1={fib['sigma_nom_h1_mpa']:.3f} h5={fib['sigma_nom_h5_mpa']:.3f} "
          f"h50={fib['sigma_nom_h50_mpa']:.3f} MPa")
    print(f"       σ_nom/σc: h1={fib['ratio_h1']:.3f} h5={fib['ratio_h5']:.3f} "
          f"h50={fib['ratio_h50']:.3f}  [min={fib['ratio_min']:.3f} "
          f"max={fib['ratio_max']:.3f}]  all_one={fib['all_one']} "
          f"(n_fractured_h={fib['n_fractured_h']}/50)")
    print(f"[胫骨] σc={tib['sigma_c_weak_mpa']:.0f} MPa  "
          f"σ_nom: h1={tib['sigma_nom_h1_mpa']:.3f} h5={tib['sigma_nom_h5_mpa']:.3f} "
          f"h50={tib['sigma_nom_h50_mpa']:.3f} MPa")
    print(f"       σ_nom/σc: min={tib['ratio_min']:.3f} max={tib['ratio_max']:.3f}  "
          f"all_one={tib['all_one']} (n_fractured_h={tib['n_fractured_h']}/50)")

    print(f"\n[σc 灵敏性] 腓骨 all_one 的临界 σc* = {res['sigma_c_sensitivity_fibula']['all_one_critical_sigma_c_mpa']:.3f} MPa")
    for row in res["sigma_c_sensitivity_fibula"]["candidates"]:
        print(f"   σc={row['sigma_c_mpa']:>8.3f} -> n_frac={row['n_fractured_h']:>2} "
              f"all_one={row['all_one']} first_safe_h={row['first_safe_h']}")

    print(f"\n[载荷分配反事实] 使腓骨 h1 安全（出现 0/1 混合）的临界分担比 = "
          f"{res['load_share_counterfactual_fibula']['share_h1_safe']:.4f}；"
          f"全安全(all_zero)临界 = {res['load_share_counterfactual_fibula']['share_all_safe']:.4f}")
    for c in res["load_share_counterfactual_fibula"]["cases"]:
        print(f"   {c['case']:>20s} share={c['share']:.4f} -> ratio h1={c['ratio_h1']:.3f} "
              f"h50={c['ratio_h50']:.3f} n_frac={c['n_fractured_h']:>2} "
              f"identifiable(0/1 mixed)={c['identifiable_0_1_mixed']}")

    print(f"\n[姊妹比] σ_nom(腓骨/胫骨): h1={res['sibling_tibia_vs_fibula']['sig_nom_fibula_over_tibia']['h1']:.3f} "
          f"h5={res['sibling_tibia_vs_fibula']['sig_nom_fibula_over_tibia']['h5']:.3f} "
          f"h50={res['sibling_tibia_vs_fibula']['sig_nom_fibula_over_tibia']['h50']:.3f}")
    print(f"         A_section: tibia={res['sibling_tibia_vs_fibula']['A_section_mm2']['tibia_r']:.2f} "
          f"fibula={res['sibling_tibia_vs_fibula']['A_section_mm2']['fibula_r']:.2f} mm^2 "
          f"(ratio {res['sibling_tibia_vs_fibula']['A_section_mm2']['tibia_r']/res['sibling_tibia_vs_fibula']['A_section_mm2']['fibula_r']:.3f})")
    print(f"         ankle load: h1={res['sibling_tibia_vs_fibula']['ankle_load_n']['h1']:.1f} "
          f"h5={res['sibling_tibia_vs_fibula']['ankle_load_n']['h5']:.1f} "
          f"h50={res['sibling_tibia_vs_fibula']['ankle_load_n']['h50']:.1f} N")

    print(f"\n[FE 溯源] 腓骨: {res['fe_provenance']['fibula_r']['fe_solve_load_n']} N "
          f"lrf={res['fe_provenance']['fibula_r']['load_reduction_factor']:.3f} "
          f"linearity={res['fe_provenance']['fibula_r']['linearity']} "
          f"attempts={res['fe_provenance']['fibula_r']['attempts']}")
    print(f"          胫骨: {res['fe_provenance']['tibia_r']['fe_solve_load_n']} N "
          f"lrf={res['fe_provenance']['tibia_r']['load_reduction_factor']:.3f} "
          f"linearity_median_ratio={(res['fe_provenance']['tibia_r']['linearity'] or {}).get('ratio_median')}")
    print(f"[legacy 腓骨] fe_solve={res['legacy_fibula_fe']['fe_solve_load_n']} N "
          f"lrf={res['legacy_fibula_fe']['load_reduction_factor']:.3f}")

    print(f"\nwrote: {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
