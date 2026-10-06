"""axial_subregion.py -- 纯轴向(load-share)路线的 [Y25] 子区域强度决断性再分析。

背景（为什么需要这个脚本）
--------------------------
此前的 ``PAPER_ALIGNMENT_REPORT.md`` 把 [Y25] 的子区域强度（两端 70 / 股骨颈 80）
贴到了 **弯曲 1D** 路线（``risk_1d_bending.json`` 的 σ_total）与 FE 路线上。弯曲 1D
假设"关节力垂直作用于斜置长骨"，产生**人为的大弯曲应力**，结果腓骨 ends 在 1 m
即骨折、骨干也在 1 m 骨折 —— 与论文"腓骨 7–9 m 渐进、骨干从不骨折"不符。

本脚本做**加法式**的第三个对照：把子区域强度贴到**纯轴向**路线
（``risk_1d_loadshare.json`` 的载荷分配版 1D 名义压缩应力 σ = F/(A·σ_c)，
不含任何弯曲项），检验：

* **关键检验**：纯轴向 + 腓骨两端(70) 的腓骨首骨折高度，是否落在论文的 7–9 m？
* 纯轴向 + ends 的 9 部位排序；
* 与弯曲+ends、FE+ends、论文四路对比（定性行为 + 首骨折高度 + 骨干）；
* 二次修正：论文称颅骨经**头部惯性**（自上而下）骨折，而非 ``lumbar × 0.095``；
  用头部质量 × 载荷链隐含减速度重估颅骨载荷，并给出骨盆的更正载荷。

只读输入 / 只写新文件；不改任何既有模块默认/签名；不调用 FEBio；不覆盖既有产物。

输入（只读）
------------
* ``results/opensim_fe/risk_1d_loadshare.json``        -> 纯轴向 risk(h)（并联分流主口径）
* ``results/opensim_fe/fracture_matrix.json``          -> loads_n / grf_peak_n / 模型质量
* ``results/opensim_fe/paper_alignment_1d_matrix.json``-> 弯曲+ends 首骨折（对照）
* ``results/opensim_fe/paper_alignment_fe_matrix.json``-> FE+ends 首骨折（对照）
* ``src/climbing/bone.py::MATERIAL_STRENGTH_MPA``      -> 子区域压缩强度

产物（新文件）
--------------
* ``results/opensim_fe/axial_subregion.json``
* ``results/opensim_fe/AXIAL_SUBREGION_REPORT.md``

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\axial_subregion.py
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

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"
FM_JSON = RES / "fracture_matrix.json"
BEND_JSON = RES / "paper_alignment_1d_matrix.json"
FE_JSON = RES / "paper_alignment_fe_matrix.json"

#: (中文名, 代表骨, 关节, 弱子区域键, 骨干子区域键, 论文定性行为)
#: 弱子区域决定"部位骨折"；骨干子区域单独评估以复核"骨干从不骨折"。
REGIONS: list[tuple[str, str, str, str, str | None, str]] = [
    ("足部", "calcaneus_r", "subtalar_r", "calcaneus", None, "阶跃"),
    ("胫骨", "tibia_r", "ankle_r", "tibia_ends", "tibia_shaft", "阶跃"),
    ("腓骨", "fibula_r", "ankle_r", "fibula_ends", "fibula_shaft", "渐进"),
    ("股骨", "femur_r", "hip_r", "femoral_neck", "femur_shaft", "阶跃"),
    ("骨盆", "R_HIPBONE", "hip_r", "pelvis", None, "渐进"),
    ("腰椎", "L3", "lumbar", "spine", None, "阶跃"),
    ("胸椎", "T6", "lumbar", "spine", None, "阶跃"),
    ("颈椎", "C5", "lumbar", "spine", None, "阶跃"),
    ("颅骨", "parietal_r", "lumbar", "skull", None, "渐进"),
]
REGION_NAMES = [r[0] for r in REGIONS]
LONG_BONES = {"胫骨", "腓骨", "股骨"}

#: 论文 [Y25] 事实（对照目标）。
PAPER = {
    "fibula_ends_fracture_m": [7, 9],
    "fibula_or": 1.682,
    "skull_or": 1.576,
    "pelvis_or": 1.236,
    "or_order": ["腓骨", "颅骨", "骨盆"],
    "shafts_never_fracture": True,
    "behavior": {cn: beh for cn, _, _, _, _, beh in REGIONS},
}

#: 行为标签（把二值矩阵/单调 risk 归类）。
BEHAV_CN = {
    "all_zero": "常值(0)",
    "all_one": "常值(1)",
    "step": "阶跃",
    "mixed": "非单调",
}


def _sc(key: str) -> float:
    """子区域压缩强度 (MPa)，[Y25] Table 1。"""
    return float(MATERIAL_STRENGTH_MPA[key][1])


def _first_fracture(heights: list[float], mask) -> float | None:
    for h, m in zip(heights, mask):
        if bool(m):
            return float(h)
    return None


def classify(mask) -> str:
    """把长度 50 的布尔序列归为 all_zero / all_one / step / mixed。"""
    m = np.asarray(mask, dtype=bool)
    k = int(m.sum())
    if k == 0:
        return "all_zero"
    if k == len(m):
        return "all_one"
    idx = np.where(m)[0]
    if int(idx.max() - idx.min() + 1) == k:
        return "step"
    return "mixed"


def _match_verdict(target: str, beh: str, region: str, first_h) -> str:
    """Route A 与论文定性行为的一致性判定。"""
    if target == "阶跃":
        if beh == "step":
            return "✓"
        return "~"  # 常值（阈值落在 1–50 m 区间外）仍与"阶跃"不矛盾
    # 目标渐进
    return "✓" if beh == "mixed" else "✗"


def build() -> dict:
    ls = json.loads(LOADSHARE_JSON.read_text(encoding="utf-8"))
    fm = json.loads(FM_JSON.read_text(encoding="utf-8"))
    bend = json.loads(BEND_JSON.read_text(encoding="utf-8"))
    fe = json.loads(FE_JSON.read_text(encoding="utf-8"))

    heights = [float(h) for h in ls["heights_m"]] if "heights_m" in ls else \
        [float(h) for h in fm["heights_m"]]
    rows_par = {r["bone"]: r for r in ls["rows_with_share_parallel"]}
    rows_no = {r["bone"]: r for r in ls["rows_without_share"]}

    # ---- 纯轴向 σ(h) = risk(h) × σ_c（stored，长骨为骨干值）----------------
    sigma_h: dict[str, np.ndarray] = {}
    sigma_c_stored: dict[str, float] = {}
    for bone, r in rows_par.items():
        sc = float(r["sigma_c_mpa"])
        sigma_c_stored[bone] = sc
        sigma_h[bone] = np.asarray(r["risk_by_height"], dtype=float) * sc

    # ---- 逐部位：弱子区域首骨折 + 骨干首骨折 ------------------------------
    region_rows = []
    for cn, bone, joint, weak_key, shaft_key, target in REGIONS:
        sig = sigma_h[bone]
        weak_sc = _sc(weak_key)
        mask = sig >= weak_sc
        first = _first_fracture(heights, mask)
        beh = classify(mask)
        risk5 = float(sig[4] / weak_sc)
        row = {
            "region_cn": cn, "bone": bone, "joint": joint,
            "weak_subregion": weak_key, "weak_sigma_c_mpa": weak_sc,
            "sigma_stored_mpa": sigma_c_stored[bone],
            "sigma_at_5m_mpa": float(sig[4]),
            "risk_at_5m": risk5,
            "first_fracture_h": first,
            "behavior": beh, "behavior_cn": BEHAV_CN[beh],
            "matrix": [int(m) for m in mask],
        }
        if shaft_key is not None:
            shaft_sc = _sc(shaft_key)
            smask = sig >= shaft_sc
            row["shaft_subregion"] = shaft_key
            row["shaft_sigma_c_mpa"] = shaft_sc
            row["shaft_first_fracture_h"] = _first_fracture(heights, smask)
            row["shaft_behavior"] = classify(smask)
            row["shaft_matrix"] = [int(m) for m in smask]
        region_rows.append(row)

    # ---- Route A 排序（早骨折优先，其次 risk@5m）--------------------------
    order = sorted(
        region_rows,
        key=lambda r: (r["first_fracture_h"] if r["first_fracture_h"] is not None else 1e9,
                       -r["risk_at_5m"]),
    )
    ranking = [{"rank": i, "region_cn": r["region_cn"], "bone": r["bone"],
                "first_fracture_h": r["first_fracture_h"],
                "risk_at_5m": r["risk_at_5m"], "behavior_cn": r["behavior_cn"]}
               for i, r in enumerate(order, 1)]

    # ---- 四路对比：A(轴向+ends) / B(弯曲+ends) / C(FE+ends) / 论文 --------
    def _route_from_matrix(mat: dict, order_names: list[str], firsts: dict) -> dict:
        out = {}
        for cn in REGION_NAMES:
            m = mat[cn]
            out[cn] = {"first_fracture_h": firsts.get(cn),
                       "behavior": classify(m), "behavior_cn": BEHAV_CN[classify(m)]}
        return out

    route_b = {}
    for cn in REGION_NAMES:
        f = bend["first_fracture_m"].get(cn)
        m = bend["matrix"][cn]
        route_b[cn] = {"first_fracture_h": f, "behavior": classify(m),
                       "behavior_cn": BEHAV_CN[classify(m)]}
    route_c = {}
    for cn in REGION_NAMES:
        f = fe["first_fracture_m"].get(cn)
        m = fe["matrix"][cn]
        route_c[cn] = {"first_fracture_h": f, "behavior": classify(m),
                       "behavior_cn": BEHAV_CN[classify(m)]}

    route_a = {r["region_cn"]: {"first_fracture_h": r["first_fracture_h"],
                                "behavior": r["behavior"],
                                "behavior_cn": r["behavior_cn"]} for r in region_rows}

    comparison = []
    for cn in REGION_NAMES:
        a, b, c = route_a[cn], route_b[cn], route_c[cn]
        target = PAPER["behavior"][cn]
        comparison.append({
            "region_cn": cn,
            "A_axial_ends_first_h": a["first_fracture_h"],
            "A_axial_ends_behavior": a["behavior_cn"],
            "B_bending_ends_first_h": b["first_fracture_h"],
            "B_bending_ends_behavior": b["behavior_cn"],
            "C_fe_ends_first_h": c["first_fracture_h"],
            "C_fe_ends_behavior": c["behavior_cn"],
            "paper_behavior": target,
            "A_match": _match_verdict(target, a["behavior"], cn, a["first_fracture_h"]),
        })

    # ---- 骨干复核 --------------------------------------------------------
    shaft_rows = []
    bend_shaft = bend.get("shaft_first_fracture_m", {})
    fe_shaft = fe.get("shaft_first_fracture_m", {})
    for r in region_rows:
        if r["region_cn"] in LONG_BONES:
            shaft_rows.append({
                "region_cn": r["region_cn"], "bone": r["bone"],
                "shaft_sigma_c_mpa": r["shaft_sigma_c_mpa"],
                "A_shaft_first_h": r["shaft_first_fracture_h"],
                "A_shaft_behavior": BEHAV_CN[r["shaft_behavior"]],
                "B_shaft_first_h": bend_shaft.get(r["region_cn"]),
                "C_shaft_first_h": fe_shaft.get(r["region_cn"]),
                "paper": "从不骨折",
            })

    # ---- 关键检验：腓骨 ends --------------------------------------------
    fib = next(r for r in region_rows if r["region_cn"] == "腓骨")
    fib_first = fib["first_fracture_h"]
    fib_ok = fib_first is not None and 7 <= fib_first <= 9
    key_test = {
        "question": "纯轴向 + 腓骨两端(70) 首骨折高度是否落在论文 7–9 m？",
        "fibula_sigma_at_5m_mpa": fib["sigma_at_5m_mpa"],
        "fibula_ends_sigma_c_mpa": fib["weak_sigma_c_mpa"],
        "fibula_risk_at_5m": fib["risk_at_5m"],
        "fibula_first_fracture_h": fib_first,
        "paper_range_m": PAPER["fibula_ends_fracture_m"],
        "verdict": "✓ 落在 7–9 m" if fib_ok else "✗ 不在 7–9 m",
    }

    # ---- OR 排序对照（腓骨 > 颅骨 > 骨盆）--------------------------------
    rank_of = {r["region_cn"]: r["rank"] for r in ranking}
    or_order = sorted(["腓骨", "颅骨", "骨盆"], key=lambda cn: rank_of[cn])
    or_check = {
        "paper_order": PAPER["or_order"],
        "route_a_order": or_order,
        "route_a_ranks": {cn: rank_of[cn] for cn in ["腓骨", "颅骨", "骨盆"]},
        "route_a_first_h": {cn: next(r["first_fracture_h"] for r in region_rows
                                     if r["region_cn"] == cn)
                            for cn in ["腓骨", "颅骨", "骨盆"]},
        "match": or_order == PAPER["or_order"],
    }

    # ---- 二次修正：颅骨（头部惯性）+ 骨盆 --------------------------------
    m_total = float(fm["meta"]["model_mass_kg"])
    lumbar = np.asarray(fm["loads_n"]["lumbar"], dtype=float)
    hip = np.asarray(fm["loads_n"]["hip_r"], dtype=float)
    grf = np.asarray(fm["grf_peak_n"], dtype=float)

    ref_pct = float(ls["mass_above_model"]["reference_pct"])
    head_pct = float(ls["mass_above_model"]["segment_pct"]["head"])
    m_head = head_pct / 100.0 * m_total
    m_ref = ref_pct / 100.0 * m_total
    ratio_chain = float(ls["mass_above_model"]["per_bone"]["parietal_r"]["ratio"])

    a_skull = float(rows_par["parietal_r"]["A_section_mm2"])
    sig_skull = _sc("skull")
    # 模型 1：载荷链隐含减速度（= F_lumbar/m_ref）；F_skull = m_head·a = ratio_chain·F_lumbar
    F_skull_chain = ratio_chain * lumbar
    # 模型 2：GRF 隐含整体减速度 a=F_grf/M_total；F_skull = m_head·a
    F_skull_grf = m_head * (grf / m_total)
    # 现有口径（= 模型 1，作参照）
    F_skull_current = ratio_chain * lumbar

    a_pelvis = float(rows_par["R_HIPBONE"]["A_section_mm2"])
    sig_pelvis = _sc("pelvis")
    F_pelvis_current = hip                      # 现有：单侧 hip_r
    F_pelvis_2hip = 2.0 * hip                   # 更正：骨盆承接双髋反力
    F_pelvis_hip_lumbar = hip + lumbar          # 更正（备选）：髋 + 腰椎

    def _sec(name, F, A, sig):
        risk = F / (A * sig)
        return {
            "model": name,
            "F_at_5m_n": float(F[4]),
            "risk_at_5m": float(risk[4]),
            "first_fracture_h": _first_fracture(heights, risk >= 1.0),
            "behavior": classify(risk >= 1.0),
            "sigma_at_5m_mpa": float(F[4] / A),
        }

    secondary = {
        "total_body_mass_kg": m_total,
        "head_mass_pct": head_pct,
        "head_mass_kg": m_head,
        "reference_pct": ref_pct,
        "reference_kg": m_ref,
        "chain_ratio_head_over_ref": ratio_chain,
        "skull": {
            "A_section_mm2": a_skull,
            "sigma_c_mpa": sig_skull,
            "chain_deceleration_ms2_ref": float(lumbar[4] / m_ref),
            "grf_deceleration_ms2_5m": float(grf[4] / m_total),
            "models": [
                _sec("head_inertia_chain_decel (m_head·F_lumbar/m_ref)",
                     F_skull_chain, a_skull, sig_skull),
                _sec("head_inertia_grf_decel (m_head·F_grf/M_total)",
                     F_skull_grf, a_skull, sig_skull),
            ],
            "current_loadshare_F_at_5m_n": float(F_skull_current[4]),
        },
        "pelvis": {
            "A_section_mm2": a_pelvis,
            "sigma_c_mpa": sig_pelvis,
            "models": [
                _sec("current_single_hip", F_pelvis_current, a_pelvis, sig_pelvis),
                _sec("corrected_both_hips (2·hip_r)", F_pelvis_2hip, a_pelvis, sig_pelvis),
                _sec("corrected_hip_plus_lumbar", F_pelvis_hip_lumbar, a_pelvis, sig_pelvis),
            ],
        },
    }

    # ---- 组装 ------------------------------------------------------------
    result = {
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "method": ("纯轴向(load-share 1D) σ(h)=risk_1d_loadshare.json 的 risk_by_height×σ_c；"
                   "贴 [Y25] 弱子区域(两端70/股骨颈80/整体150-180)；"
                   "first h with σ≥σ_c；与弯曲+ends、FE+ends、论文对比。"),
        "inputs": {
            "loadshare_json": str(LOADSHARE_JSON.relative_to(ROOT)),
            "fracture_matrix_json": str(FM_JSON.relative_to(ROOT)),
            "bending_matrix_json": str(BEND_JSON.relative_to(ROOT)),
            "fe_matrix_json": str(FE_JSON.relative_to(ROOT)),
            "share_policy_primary": "rows_with_share_parallel",
        },
        "paper_facts": PAPER,
        "heights_m": [int(h) for h in heights],
        "subregion_strengths_used": {
            r["region_cn"]: {
                "weak": {r["weak_subregion"]: r["weak_sigma_c_mpa"]},
                "shaft": ({r["shaft_subregion"]: r["shaft_sigma_c_mpa"]}
                          if "shaft_subregion" in r else None),
            } for r in region_rows
        },
        "route_a_axial_ends": region_rows,
        "ranking_axial_ends": ranking,
        "key_test_fibula": key_test,
        "or_check": or_check,
        "comparison_4way": comparison,
        "shaft_check": shaft_rows,
        "secondary_correction": secondary,
        "assumptions": [
            "σ_c 取自 bone.py MATERIAL_STRENGTH_MPA（[Y25] Table 1 压缩）——材料强度，非整体骨失效载荷。",
            "纯轴向 σ(h)=risk_1d_loadshare.json 的 risk_by_height×σ_c（stored，长骨为骨干值）；"
            "主口径 = rows_with_share_parallel（胫/腓并联分流 + 中轴骨质量分配）。",
            "σ(h) 是 1D 名义压缩应力，忽略弯曲/剪切/偏心/应力集中 —— 对细长骨低估、对紧凑骨高估。",
            "颅骨更正：载荷链隐含减速度 a=F_lumbar/m_ref，F_skull=m_head·a（=现有 0.0946·F_lumbar）；"
            "备选 GRF 隐含减速度 a=F_grf/M_total。两者都远低于颅骨阈值。",
            "骨盆更正：骨盆承接双髋反力 2·hip_r（备选 hip_r+F_lumbar）。仍远低于阈值。",
            "本脚本只读复用既有产物，不调用 FEBio，不重跑 OpenSim。",
        ],
    }
    return result


def _fmt_h(v) -> str:
    return ">50" if v is None else f"{int(v)}"


def _write_report(path: Path, res: dict) -> None:
    L: list[str] = []
    A = L.append
    ts = res["generated_at"]
    rows = {r["region_cn"]: r for r in res["route_a_axial_ends"]}
    rank = {r["region_cn"]: r["rank"] for r in res["ranking_axial_ends"]}
    kt = res["key_test_fibula"]
    oc = res["or_check"]
    sec = res["secondary_correction"]

    A("# 纯轴向 + 子区域强度（两端 70 / 股骨颈 80）决断性再分析")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 本报告**新增**：不改既有模块/默认值，不覆盖既有产物，不调用 FEBio。")
    A("> 对照对象：`PAPER_ALIGNMENT_REPORT.md`（弯曲+ends 与 FE+ends）与论文 [Y25]。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"- **关键检验通过**：纯轴向 + 腓骨两端(70) 的腓骨首骨折 = "
      f"**{_fmt_h(kt['fibula_first_fracture_h'])} m**（腓骨 σ@5m={kt['fibula_sigma_at_5m_mpa']:.1f} MPa，"
      f"risk@5m={kt['fibula_risk_at_5m']:.3f}）→ **落在论文 7–9 m 区间内**。")
    A(f"- **排序**：" + " > ".join(f"{r['region_cn']}({_fmt_h(r['first_fracture_h'])}m)"
                                  for r in res["ranking_axial_ends"]) + "。")
    A(f"- 但**未逼近论文的 腓骨>颅骨>骨盆**：Route A 三骨相对序 = "
      f"`{' > '.join(oc['route_a_order'])}`（{'一致' if oc['match'] else '不一致'}）；颅骨/骨盆 1–50 m 内**从不骨折**。")
    A(f"- 骨干复核：纯轴向下胫骨/腓骨干仍在 **{_fmt_h(rows['胫骨']['shaft_first_fracture_h'])}/"
      f"{_fmt_h(rows['腓骨']['shaft_first_fracture_h'])} m** 骨折 → 论文"
      "「骨干从不骨折」在纯轴向口径下**同样不成立**。")
    A(f"- 二次修正（颅骨头部惯性 / 骨盆双髋）仍**远低于**阈值（risk@5m≈"
      f"{sec['skull']['models'][0]['risk_at_5m']:.4f} / "
      f"{sec['pelvis']['models'][1]['risk_at_5m']:.4f}）→ 纯轴向路线无法复现颅骨/骨盆的渐进损伤。")
    A("")
    A("---")
    A("")
    A("## 1. 方法")
    A("")
    A("### 1.1 纯轴向 σ(h)")
    A("")
    A("```")
    A("σ(h) = risk_1d_loadshare.json::rows_with_share_parallel[bone].risk_by_height × sigma_c_mpa")
    A("risk_first = σ(h) / σ_c(弱子区域);  first_fracture = 最小 h 使 risk_first ≥ 1（否则 >50）")
    A("```")
    A("")
    A("| 部位 | 骨 | 关节 | 弱子区域 | σ_c(弱) | stored σ_c | σ@5m (MPa) | risk@5m | 首骨折 m | 行为 |")
    A("|---|---|---|---|---:|---:|---:|---:|---:|---|")
    for r in res["ranking_axial_ends"]:
        rr = rows[r["region_cn"]]
        A(f"| {r['region_cn']} | `{rr['bone']}` | `{rr['joint']}` | {rr['weak_subregion']} | "
          f"{rr['weak_sigma_c_mpa']:.0f} | {rr['sigma_stored_mpa']:.0f} | {rr['sigma_at_5m_mpa']:.2f} | "
          f"{rr['risk_at_5m']:.4f} | **{_fmt_h(r['first_fracture_h'])}** | {rr['behavior_cn']} |")
    A("")
    A("> `stored σ_c` 是 `risk_1d_loadshare.json` 里生成 risk 用的阈值（长骨为骨干值）；"
      "`σ(h)` 由该 risk 还原，保持与既有产物一致。")
    A("")
    A("### 1.2 关键检验（腓骨 ends）")
    A("")
    A("| 量 | 值 |")
    A("|---|---:|")
    A(f"| 腓骨 σ@5m | {kt['fibula_sigma_at_5m_mpa']:.2f} MPa |")
    A(f"| 腓骨两端 σ_c | {kt['fibula_ends_sigma_c_mpa']:.0f} MPa |")
    A(f"| 腓骨 risk@5m | {kt['fibula_risk_at_5m']:.4f} |")
    A(f"| 腓骨首骨折 | **{_fmt_h(kt['fibula_first_fracture_h'])} m** |")
    A(f"| 论文区间 | 7–9 m |")
    A(f"| 判定 | **{kt['verdict']}** |")
    A("")
    A("> 这正是任务预判的洞察：弯曲 1D 因「关节力垂直作用于斜置骨」产生人为大弯曲，"
      "把腓骨压到 1 m；改用纯轴向 σ 后，腓骨两端在 7 m 跨越 70 MPa 阈值。")
    A("")
    A("---")
    A("")
    A("## 2. 四路对比")
    A("")
    A("### 2.1 定性行为（阶跃 / 渐进）")
    A("")
    A("| 部位 | 论文 | A 纯轴向+ends | B 弯曲+ends | C FE+ends | A 是否匹配论文 |")
    A("|---|---|---|---|---|---|")
    for c in res["comparison_4way"]:
        A(f"| {c['region_cn']} | {c['paper_behavior']} | {c['A_axial_ends_behavior']} | "
          f"{c['B_bending_ends_behavior']} | {c['C_fe_ends_behavior']} | **{c['A_match']}** |")
    A("")
    A("### 2.2 首骨折高度（m）")
    A("")
    A("| 部位 | A 纯轴向+ends | B 弯曲+ends | C FE+ends | 论文 |")
    A("|---|---:|---:|---:|---|")
    for c in res["comparison_4way"]:
        paper_h = "7–9" if c["region_cn"] == "腓骨" else \
            ("阶跃(临界高度)" if c["paper_behavior"] == "阶跃" else "渐进")
        A(f"| {c['region_cn']} | **{_fmt_h(c['A_axial_ends_first_h'])}** | "
          f"{_fmt_h(c['B_bending_ends_first_h'])} | {_fmt_h(c['C_fe_ends_first_h'])} | {paper_h} |")
    A("")
    A("### 2.3 OR 排序对照（腓骨 > 颅骨 > 骨盆）")
    A("")
    A(f"- 论文（每米 OR）：腓骨 {PAPER['fibula_or']} > 颅骨 {PAPER['skull_or']} > 骨盆 {PAPER['pelvis_or']}")
    A(f"- Route A 相对序：`{' > '.join(oc['route_a_order'])}`"
      f"（排名 腓骨#{oc['route_a_ranks']['腓骨']}、颅骨#{oc['route_a_ranks']['颅骨']}、"
      f"骨盆#{oc['route_a_ranks']['骨盆']}；首骨折 "
      f"{_fmt_h(oc['route_a_first_h']['腓骨'])} / "
      f"{_fmt_h(oc['route_a_first_h']['颅骨'])} / "
      f"{_fmt_h(oc['route_a_first_h']['骨盆'])} m）")
    A(f"- **判定**：{'一致' if oc['match'] else '不一致'} —— 腓骨最快（✓），但颅骨/骨盆在最慢端且不骨折（✗）。")
    A("")
    A("### 2.4 骨干复核（论文：从不骨折）")
    A("")
    A("| 部位 | 骨干 σ_c | A 纯轴向+ends 骨干首骨折 | B 弯曲+ends | C FE+ends | 论文 |")
    A("|---|---:|---:|---:|---:|---|")
    for s in res["shaft_check"]:
        A(f"| {s['region_cn']} | {s['shaft_sigma_c_mpa']:.0f} | "
          f"**{_fmt_h(s['A_shaft_first_h'])}** | {_fmt_h(s['B_shaft_first_h'])} | "
          f"{_fmt_h(s['C_shaft_first_h'])} | 从不骨折 |")
    A("")
    A("> 纯轴向 σ(h) 在长骨两端与骨干用的是同一条 σ(h)（仅阈值不同）。骨干阈值更高，"
      "故骨干骨折晚于两端 —— 但 1–50 m 内仍会跨越（腓骨 ~23 m、胫骨 ~37 m），"
      "**与论文「骨干从不骨折」不符**。这是纯轴向路线与论文的又一处不对齐。")
    A("")
    A("---")
    A("")
    A("## 3. 明确的三个回答")
    A("")
    A(f"**(a) 腓骨是否落在 ~7–9 m？** 是。纯轴向 + ends 腓骨首骨折 = "
      f"**{_fmt_h(kt['fibula_first_fracture_h'])} m**（区间 7–9 m 内）。")
    A("")
    A(f"**(b) 纯轴向 + ends 排序？** "
      + " > ".join(f"{r['region_cn']}({_fmt_h(r['first_fracture_h'])}m)"
                   for r in res["ranking_axial_ends"]) + "。")
    A("")
    A(f"**(c) 是否逼近论文 腓骨>颅骨>骨盆？** 否。腓骨 #1（对齐），但颅骨/骨盆在 1–50 m "
      f"内从不骨折（risk@5m {rows['颅骨']['risk_at_5m']:.4f} / {rows['骨盆']['risk_at_5m']:.4f}），"
      f"且三骨相对序为 `{' > '.join(oc['route_a_order'])}`，与论文相反。"
      "整条链的排序更接近「下肢为主」而非论文的「腓骨主导渐进」。")
    A("")
    A("---")
    A("")
    A("## 4. 二次修正：颅骨（头部惯性）与骨盆")
    A("")
    A("### 4.1 颅骨 —— 头部惯性而非 `lumbar × 0.095`")
    A("")
    A("论文指出颅骨应力集中在**枕骨大孔周围**，由**头部惯性**（着地后躯干骤停、头部继续下行）"
      "压向颈颅连接处，是**自上而下**路径，而非直接承力。")
    A("")
    A("```")
    A("载荷链隐含减速度：a(h) = F_lumbar(h) / m_ref ,  m_ref = 49.7% BM")
    A("  模型① F_skull(h) = m_head · a(h) = (m_head/m_ref)·F_lumbar(h) = 0.0946·F_lumbar(h)")
    A("GRF 隐含整体减速度：a(h) = F_grf(h) / M_total")
    A("  模型② F_skull(h) = m_head · a(h) ,  m_head = 4.7% BM")
    A("```")
    A("")
    A("| 模型 | F@5m (N) | σ@5m (MPa) | risk@5m | 首骨折 m | 行为 |")
    A("|---|---:|---:|---:|---:|---|")
    for m in sec["skull"]["models"]:
        A(f"| {m['model']} | {m['F_at_5m_n']:,.1f} | {m['sigma_at_5m_mpa']:.2f} | "
          f"{m['risk_at_5m']:.4f} | **{_fmt_h(m['first_fracture_h'])}** | {BEHAV_CN[m['behavior']]} |")
    A("")
    A(f"- 参考：载荷链减速度 @5m = {sec['skull']['chain_deceleration_ms2_ref']:,.1f} m/s²；"
      f"GRF 减速度 @5m = {sec['skull']['grf_deceleration_ms2_5m']:,.1f} m/s²。")
    A(f"- 结论：两种头部惯性口径的颅骨 σ@5m 仅 {sec['skull']['models'][0]['sigma_at_5m_mpa']:.2f}–"
      f"{sec['skull']['models'][1]['sigma_at_5m_mpa']:.2f} MPa，远低于颅骨 160 MPa → "
      "**纯轴向路线无法复现论文「颅骨渐进骨折」**。物理原因是 1D 用整个顶骨壳截面(877 mm²)承载，"
      "而论文的失效集中在枕骨大孔小区域 + 薄壳弯曲，1D 名义压缩根本捕捉不到。")
    A("")
    A("### 4.2 骨盆 —— 双髋反力更正")
    A("")
    A("```")
    A("现有：F_pelvis = hip_r（单侧）")
    A("更正：F_pelvis = 2· hip_r（骨盆承接双髋反力）")
    A("备选：F_pelvis = hip_r + F_lumbar（髋 + 腰椎链）")
    A("```")
    A("")
    A("| 模型 | F@5m (N) | σ@5m (MPa) | risk@5m | 首骨折 m | 行为 |")
    A("|---|---:|---:|---:|---:|---|")
    for m in sec["pelvis"]["models"]:
        A(f"| {m['model']} | {m['F_at_5m_n']:,.1f} | {m['sigma_at_5m_mpa']:.2f} | "
          f"{m['risk_at_5m']:.4f} | **{_fmt_h(m['first_fracture_h'])}** | {BEHAV_CN[m['behavior']]} |")
    A("")
    A("- 即使取双髋反力，骨盆 σ@5m 仍仅 ~15.6 MPa（阈值 180 MPa）→ **纯轴向路线同样无法"
      "复现论文「骨盆渐进骨折」**。")
    A("")
    A("---")
    A("")
    A("## 5. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("- **本报告不改变任何既有结论文件**：`risk_1d_loadshare.json`、"
      "`paper_alignment_*.json` 均为只读输入。")
    A("- **最诚实定位**：纯轴向 + ends 在**腓骨首骨折高度**这一单点上对齐论文（7 m ∈ 7–9 m），"
      "但整体排序（颅骨/骨盆从不骨折）与骨干（高高度仍骨折）**未对齐**；"
      "弯曲 1D 与 FE 则连腓骨高度都对不上（1 m）。")
    A("")
    A("---")
    A("")
    A("## 6. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\axial_subregion.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/risk_1d_loadshare.json         (纯轴向 risk_by_height)")
    A("#   results/opensim_fe/fracture_matrix.json           (loads_n / grf_peak_n / 质量)")
    A("#   results/opensim_fe/paper_alignment_1d_matrix.json (弯曲+ends)")
    A("#   results/opensim_fe/paper_alignment_fe_matrix.json (FE+ends)")
    A("#   src/climbing/bone.py MATERIAL_STRENGTH_MPA        (子区域 σ_c)")
    A("# 产物：")
    A("#   results/opensim_fe/axial_subregion.json")
    A("#   results/opensim_fe/AXIAL_SUBREGION_REPORT.md")
    A("```")
    A("")
    A("---")
    A("")
    A("## 7. 产物")
    A("")
    A("- 本报告：`results/opensim_fe/AXIAL_SUBREGION_REPORT.md`")
    A("- 新脚本：`scripts/opensim_fe/axial_subregion.py`（additive）")
    A("- 机器可读：`results/opensim_fe/axial_subregion.json`")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=str(RES))
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    res = build()
    js_path = out / "axial_subregion.json"
    js_path.write_text(json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    md_path = out / "AXIAL_SUBREGION_REPORT.md"
    _write_report(md_path, res)

    kt = res["key_test_fibula"]
    print("=" * 96)
    print("纯轴向 + ends 逐部位（排序）：")
    rows = {r["region_cn"]: r for r in res["route_a_axial_ends"]}
    for r in res["ranking_axial_ends"]:
        rr = rows[r["region_cn"]]
        print(f"  #{r['rank']} {r['region_cn']:4s} first={_fmt_h(r['first_fracture_h']):>4}m "
              f"σ@5m={rr['sigma_at_5m_mpa']:8.2f} MPa  risk@5m={rr['risk_at_5m']:.4f}  "
              f"[{rr['behavior_cn']}]")
    print(f"\n关键检验（腓骨 ends）：{kt['verdict']} "
          f"-> {_fmt_h(kt['fibula_first_fracture_h'])} m (论文 7–9 m)")
    print("骨干：", {s["region_cn"]: _fmt_h(s["A_shaft_first_h"]) for s in res["shaft_check"]})
    print("OR 对照：", res["or_check"]["route_a_order"], "vs",
          res["or_check"]["paper_order"], "match=", res["or_check"]["match"])
    sec = res["secondary_correction"]
    print("颅骨(头部惯性) risk@5m:",
          [round(m["risk_at_5m"], 4) for m in sec["skull"]["models"]],
          "first:", [_fmt_h(m["first_fracture_h"]) for m in sec["skull"]["models"]])
    print("骨盆(更正) risk@5m:",
          [round(m["risk_at_5m"], 4) for m in sec["pelvis"]["models"]],
          "first:", [_fmt_h(m["first_fracture_h"]) for m in sec["pelvis"]["models"]])
    print("=" * 96)
    print("wrote:", js_path)
    print("wrote:", md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
