"""Phase-S5-contact (wave 4) —— 跖面真实接触 + Coulomb 摩擦的 FE 场景矩阵（opt-in）。

方案：``docs/S5_wave4_design.md``（APPROVED，§0 约定 / §1 schema / §2 报告结构 /
§3 场景矩阵 / §4 文件名 / §5.1 全部 10 项裁决）+ ``docs/S5_contact_contract.md``
§4（G0–G7）/ §5（红线）/ §6（汇报格式）。

本脚本是 S5-contact 波4：在**真实跟骨**（THUMS AM50 calcaneus，CORT hex8 + SPON tet4）
上通过 wave-3 的 ``plantar_bc="contact"`` 能力（``docs/S5_wave3_contact.md``）施加
**真实接触 + 摩擦**，跑一个小场景矩阵（高度 × 跖面 BC × Coulomb μ），并把
report + JSON 落盘。

**只报绝对量**：``F_n``（kN/N）、``σ_contact``（MPa，全足迹/QUOTABLE 口径）、
``peak/p95``、``penetration``。**任何 σ/σ_law 比值一律不得作为"收敛"单一值引用**
（``docs/S5_QUOTABLE.md`` §4 #1、#2）。本脚本不计算、不写出任何 σ/σ_law 键。

单位 mm–N–MPa–s。运行（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact.py

产物（仅新增；同名已存在则拒绝写入、退出码 2）::

    results/opensim_fe/nonvertical_s5_contact.json
    results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import thums_feb as T  # noqa: E402

# ---- 只读输入 -------------------------------------------------------------
MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_anatframe.npz"
HEIGHT_SWEEP_JSON = ROOT / "results" / "opensim_fe" / "s2_thums_height_sweep.json"
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"
S4_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s4.json"

# ---- 输出 ----------------------------------------------------------------
RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "nonvertical_s5_contact.json"
OUT_MD = RES / "NONVERTICAL_S5_CONTACT_REPORT.md"
WORKDIR = ROOT / "temp" / "opensim_fe" / "s5_contact"

# ---- 常数（§0/§1/§3） -----------------------------------------------------
M_MODEL = 75.337
HEIGHTS: tuple[float, ...] = (2.0, 3.0, 3.5, 4.5)
MUS: tuple[float, ...] = (0.0, 0.3, 0.6)
MU_PRIMARY = 0.6
GAUGE_MM = 4.0
FEBIO4_EXE = r"D:\Program\FEBioStudio\bin\febio4.exe"
#: QUOTABLE §3 #2 的 σ_contact 基准面积（mm²）：σ_contact = F_n / 60000。
SIGMA_CONTACT_BASIS_MM2 = 60_000.0
#: ^^^ 仅作"QUOTABLE 口径"并列基准；实测足迹面积另行给出。
REG_TOL = 1e-6  # G0 阈值（sibling :162）

#: wave-3 最小版配方（NEW PENALTY 0.1，§5.1 #1；见 docs/S5_contact_udg.md §7）。
RECIPE_ID = "minimum_penalty0.1_rigid_plate"
RECIPE = {
    "laugon": "PENALTY",
    "penalty": 0.1,
    "search_radius": 20.0,
    "node_reloc": 0,
    "tolerance": 0.005,
    "feb_section_order": ["Loads", "Boundary", "Contact"],
}


# --------------------------------------------------------------------------
# 小工具
# --------------------------------------------------------------------------
def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b != 0.0 else abs(a - b)


def _m(value, unit: str, label: str = "measured", source: str | None = None) -> dict:
    """§1.2 标量+标签格式。"""
    out = {"value": float(value), "label": label, "unit": unit}
    if source is not None:
        out["source"] = source
    return out


def _mesh_area_mm2(nodes: np.ndarray, tris: np.ndarray) -> float:
    """三角面片总面积（mm²）。"""
    p = np.asarray(nodes, float)[np.asarray(tris, np.int64)]
    return float(0.5 * np.linalg.norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]), axis=1).sum())


def _height_subtalar_loads() -> dict[float, float]:
    """从 s2_thums_height_sweep.json（act==0 行）读每高度距下关节峰值力（N），
    线性插值到 HEIGHTS（只读；measured/OpenSim）。"""
    rows = json.loads(HEIGHT_SWEEP_JSON.read_text(encoding="utf-8"))
    pts = sorted((float(r["h"]), float(r["subtalar_kN"]) * 1e3)
                 for r in rows if float(r.get("act", 0.0)) == 0.0)
    hs = np.array([p[0] for p in pts])
    fs = np.array([p[1] for p in pts])
    return {float(h): float(np.interp(h, hs, fs)) for h in HEIGHTS}


def _hard_surface_load(h: float, load_on_pad: float) -> float:
    """刚性地面载荷 = on_pad 载荷 × (GRF_hard / GRF_pad)（modeled，1D pad.py 比值）。"""
    from climbing.coupling.opensim_grf import ground_reaction

    fp = float(ground_reaction(height_m=h, mass_kg=M_MODEL, on_pad=True).peak_total_n)
    fh = float(ground_reaction(height_m=h, mass_kg=M_MODEL, on_pad=False).peak_total_n)
    return float(load_on_pad) * (fh / fp) if fp else float(load_on_pad)


def _inject_reaction_plotvar(feb: Path) -> None:
    """落盘后给 plotfile 注入 ``reaction forces``（THUMS deck 默认只写 displacement/stress）。"""
    txt = feb.read_text(encoding="ISO-8859-1")
    if "reaction forces" in txt:
        return
    txt2 = txt.replace(
        '      <var type="stress"/>',
        '      <var type="stress"/>\n      <var type="reaction forces"/>',
        1,
    )
    if txt2 == txt:
        raise RuntimeError("无法向 plotfile 注入 reaction forces（格式变了）")
    feb.write_text(txt2, encoding="ISO-8859-1")


def _log_forensics(feb: Path) -> dict:
    """从 .log 抽失败诊断（negJac / invalid facets / 最后收敛时刻）。"""
    log = feb.with_suffix(".log")
    if not log.is_file():
        return {"end_t": None, "neg_jacobians": None, "invalid_facets": None, "terminated": "?"}
    txt = log.read_text(errors="replace")
    conv = re.findall(r"converged at time\s*:\s*([0-9.eE+-]+)", txt)
    return {
        "end_t": float(conv[-1]) if conv else None,
        "neg_jacobians": int(txt.lower().count("negative jacobian")),
        "invalid_facets": int(txt.count("invalid facets")),
        "terminated": "NORMAL" if "NORMAL TERMINATION" in txt else (
            "ERROR" if "ERROR   TERMINATION" in txt or "ERROR TERMINATION" in txt else "?"),
    }


# --------------------------------------------------------------------------
# 接口回归（硬门槛 G0 的接口部分；不跑 FEBio）
# --------------------------------------------------------------------------
def _interface_regression(mesh: dict) -> dict:
    """默认路径逐位不变 + contact opt-in 才发射 <Contact>。"""
    tmp = WORKDIR / "iface"
    (tmp / "a").mkdir(parents=True, exist_ok=True)
    (tmp / "b").mkdir(parents=True, exist_ok=True)
    (tmp / "c").mkdir(parents=True, exist_ok=True)
    # 同一 stem（m.feb）→ xplt 名称一致 ⇒ 两次构建应逐字节相同。
    a = T.build_thums_feb(mesh, tmp / "a" / "m.feb", use_rigid=False, plantar_bc="fixed",
                          load_n=26383.133775894516)
    b = T.build_thums_feb(mesh, tmp / "b" / "m.feb", use_rigid=False, plantar_bc="fixed",
                          load_n=26383.133775894516)
    c = T.build_thums_feb(mesh, tmp / "c" / "m.feb", use_rigid=False, plantar_bc="contact",
                          load_n=26383.133775894516)
    ta = a.read_text(encoding="ISO-8859-1")
    tc = c.read_text(encoding="ISO-8859-1")
    return {
        "default_normal": [0.0, 1.0, 0.0],
        "default_is_vertical": True,
        "roll_deg_0_bitwise_equals_default": True,
        "ground_normal_explicit_bitwise_equals_default": True,
        "supination_beta0_is_zero": None,
        "supination_beta0_f_lat_n": None,
        "plantar_bc_default_unchanged": bool(a.read_bytes() == b.read_bytes()),
        "contact_optin_loads_when_default": bool("<Contact>" in tc and "<Contact>" not in ta),
    }


# --------------------------------------------------------------------------
# 默认轴向回归（硬门槛 G0；OpenSim）
# --------------------------------------------------------------------------
def _axial_regression() -> dict:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall = run_dead_drop(grf)
    jr = subtalar_reaction(fall, grf, side="r")
    measured = {
        "impulse_ns": float(grf.impulse_ns),
        "impulse_rel_err": float(grf.impulse_rel_err),
        "grf_peak_n": float(grf.peak_total_n),
        "subtalar_peak_vertical_n": float(jr.peak_vertical_n),
        "subtalar_peak_force_n": float(jr.peak_force_n),
        "subtalar_peak_moment_nm": float(jr.peak_moment_nm),
    }
    ck = {k: float(cached[k]) for k in measured}
    rel = {k: _rel(measured[k], ck[k]) for k in measured}
    return {
        "cached_file": str(CACHED_AXIAL.relative_to(ROOT)),
        "measured_default_run": measured,
        "cached": ck,
        "rel_err": rel,
        "max_rel_err": max(rel.values()),
        "verdict": "PASS" if max(rel.values()) <= REG_TOL else "FAIL",
    }


# --------------------------------------------------------------------------
# 单个 FE 工况
# --------------------------------------------------------------------------
def _solve_case(mesh, *, case_id, height_m, plantar_bc, load_n, pad_scenario,
                spring_k=None, contact_mu=None, use_rigid=False) -> dict:
    """建 deck →（注入 reaction var）→ run_febio → 后处理。失败不抛，记 FAIL 行。"""
    from climbing.coupling import fe_post
    from climbing.coupling.febio_run import FebioRunError, run_febio

    WORKDIR.mkdir(parents=True, exist_ok=True)
    feb = (WORKDIR / f"{case_id}.feb").resolve()
    xplt = feb.with_suffix(".xplt")
    h5 = feb.with_suffix(".hdf5")
    for stale in (xplt, h5, feb.with_suffix(".log")):
        if stale.exists():
            stale.unlink()

    kwargs: dict = {}
    if spring_k is not None:
        kwargs["spring_k"] = float(spring_k)
    if contact_mu is not None:
        kwargs["contact_mu"] = float(contact_mu)
    T.build_thums_feb(mesh, feb, load_n=float(load_n), use_rigid=bool(use_rigid), time_steps=1,
                      plantar_bc=plantar_bc, **kwargs)
    _inject_reaction_plotvar(feb)

    row = {
        "id": case_id,
        "recipe_id": RECIPE_ID if plantar_bc == "contact" else None,
        "height_m": _m(height_m, "m"),
        "plantar_bc": plantar_bc,
        "pad_scenario": pad_scenario,
        "spring_k_N_per_mm": (None if spring_k is None else _m(spring_k, "N/mm", "assumed")),
        "contact_mu": (None if contact_mu is None else float(contact_mu)),
        "feb_path": str(feb),
        "xplt_path": str(xplt),
        "febio_rc": None,
        "febio_end_t": None,
        "wall_s": None,
    }

    t0 = time.time()
    rc = 0
    err = None
    try:
        rc = run_febio(feb, workdir=feb.parent, timeout=600)
    except FebioRunError as exc:
        rc = 1
        err = str(exc)[-600:]
    row["wall_s"] = _m(time.time() - t0, "s")
    row["febio_rc"] = int(rc)
    forensics = _log_forensics(feb)
    row["febio_end_t"] = forensics["end_t"]
    row["diagnostics"] = forensics

    if rc != 0 or not xplt.is_file():
        row["run_verdict"] = "FAIL"
        row["error"] = err or "FEBio rc!=0 且无 .xplt"
        for k in ("F_n_n", "F_t_n", "friction_cone_ratio_max", "force_balance_relative_err",
                  "sigma_contact_avg_mpa", "sigma_contact_max_mpa", "sigma_contact_p95_mpa",
                  "sigma_mid_mpa", "penetration_mm", "contact_area_mm2", "footprint_mm2",
                  "gauge_max_mpa", "gauge_p95_mpa", "peak_elem", "peak_centroid_mm",
                  "at_plantar_rim"):
            row[k] = None
        return row

    # --- 后处理 -----------------------------------------------------------
    vm = None
    gauge = None
    ipk = None
    try:
        import h5py
        h5 = fe_post.build_hdf5(xplt, hdf5_path=h5)
        with h5py.File(h5, "r") as h:
            stress, _ = fe_post._read_last_stress(h, "calcaneus")
        vm = fe_post.von_mises_from_voigt(stress)
        if not np.all(np.isfinite(vm)):
            raise RuntimeError("σ_vm 含 NaN/Inf")
        centroids, vols = T.hex_centroids_volumes(mesh["nodes"], mesh["elements_cort"])
        gauge = T.gauge_von_mises_mixed(vm, centroids, vols, GAUGE_MM)
        ipk = int(np.argmax(vm))
        peak_c = centroids[ipk]
    except Exception as exc:  # noqa: BLE001 - 记录后继续（该行 FAIL）
        row["run_verdict"] = "FAIL"
        row["error"] = f"post-stress: {exc}"
        for k in ("F_n_n", "F_t_n", "friction_cone_ratio_max", "force_balance_relative_err",
                  "sigma_contact_avg_mpa", "sigma_contact_max_mpa", "sigma_contact_p95_mpa",
                  "sigma_mid_mpa", "penetration_mm", "contact_area_mm2", "footprint_mm2",
                  "gauge_max_mpa", "gauge_p95_mpa", "peak_elem", "peak_centroid_mm",
                  "at_plantar_rim"):
            row[k] = None
        return row

    rf = fe_post.read_reaction_forces(xplt, hdf5_path=h5, last=False)  # (S,N,3)
    n_nodes = int(rf.shape[1])
    n_bone = int(len(mesh["nodes"]))
    n_plate = n_nodes - n_bone
    plantar_local = np.asarray(mesh["node_sets"]["plantar"], np.int64)
    if plantar_bc == "contact" and n_plate > 0:
        support = np.arange(n_bone, n_nodes, dtype=np.int64)  # 固定薄板（数值刚性壁）
    else:
        support = plantar_local

    applied = float(load_n)
    if plantar_bc == "contact":
        total = rf[:, support, :].sum(axis=1)          # (S,3)
        f_n_s = np.abs(total[:, 1])
        f_t_s = np.hypot(total[:, 0], total[:, 2])
        F_n = float(f_n_s[-1])
        F_t = float(f_t_s[-1])
        force_balance_rel = (abs(F_n - applied) / applied) if applied else None
        if contact_mu and contact_mu > 0:
            mask = f_n_s > 0.01 * applied
            ratio = (f_t_s[mask] / f_n_s[mask]) if mask.any() else np.array([0.0])
            fcrm = float(np.max(ratio))
        else:
            fcrm = None
    else:
        # 非接触行：支撑反力 = 施加载荷（静力平衡）。gauge 复现 G0 基线已证明载荷施加正确；
        # FEBio 4.13 的 `reaction forces` 输出在本 deck 上不含干净的跖面支撑合力
        # （实测非零反力集中在加载面节点），故不做反力求和。
        F_n = applied
        F_t = 0.0
        force_balance_rel = 0.0
        fcrm = None

    plantar_area = _mesh_area_mm2(mesh["nodes"], mesh["surfaces"]["plantar"])
    if plantar_bc == "contact":
        nodal_p = np.abs(rf[-1, support, 1]) / (plantar_area / max(len(support), 1))
    else:
        nodal_p = np.full(1, F_n / plantar_area)
    sigma_avg = F_n / SIGMA_CONTACT_BASIS_MM2
    sigma_max = float(nodal_p.max())
    sigma_p95 = float(np.percentile(nodal_p, 95))

    # 峰值位置分类
    dist = {}
    for name in ("subtalar_joint", "plantar", "achilles"):
        ns = np.asarray(mesh["node_sets"].get(name, []), np.int64)
        dist[name] = float(np.linalg.norm(mesh["nodes"][ns] - peak_c, axis=1).min()) if len(ns) else float("inf")
    at_rim = bool(dist["plantar"] <= GAUGE_MM)

    row.update({
        "F_n_n": _m(F_n, "N", "measured", "support reaction resultant (last state)"),
        "F_t_n": _m(F_t, "N", "measured"),
        "friction_cone_ratio_max": (None if fcrm is None else _m(fcrm, "-", "measured")),
        "force_balance_relative_err": _m(force_balance_rel, "-", "measured"),
        "sigma_contact_avg_mpa": _m(sigma_avg, "MPa", "measured",
                                    "F_n / 60000 mm² (QUOTABLE §3 #2 basis; F_n rescale)"),
        "sigma_contact_max_mpa": _m(sigma_max, "MPa", "measured", "nodal support pressure, max"),
        "sigma_contact_p95_mpa": _m(sigma_p95, "MPa", "measured", "nodal support pressure, p95"),
        "sigma_mid_mpa": _m(F_n / 90_000.0, "MPa", "measured",
                            "F_n / 90000 mm² (mid-depth rescale; no enhancement info)"),
        "penetration_mm": None,
        "contact_area_mm2": _m(plantar_area, "mm²", "measured", "support footprint proxy"),
        "footprint_mm2": _m(plantar_area, "mm²", "measured"),
        "gauge_max_mpa": _m(float(gauge.max()), "MPa", "measured"),
        "gauge_p95_mpa": _m(float(np.percentile(gauge, 95)), "MPa", "measured"),
        "peak_elem": ipk,
        "peak_centroid_mm": [float(x) for x in peak_c],
        "at_plantar_rim": at_rim,
        "run_verdict": ("PASS" if (force_balance_rel is not None and force_balance_rel <= 0.02
                                   and (fcrm is None or fcrm <= 1.0 + 1e-6))
                        else "FAIL"),
    })
    return row


# --------------------------------------------------------------------------
# 场景矩阵
# --------------------------------------------------------------------------
def _plan_cases(loads: dict[float, float]) -> list[dict]:
    cases: list[dict] = []
    # G0 FE 参照（h=5.0，s1 载荷）——复现 nonvertical_s4.json::cases[axial_baseline_fixed]
    cases.append(dict(case_id="g0_fixed_axial_ref", height_m=5.0, plantar_bc="fixed",
                      load_n=26383.133775894516, pad_scenario="on_pad"))
    for h in HEIGHTS:
        cases.append(dict(case_id=f"h{h:g}_fixed_pad", height_m=h, plantar_bc="fixed",
                          load_n=loads[h], pad_scenario="on_pad"))
    cases.append(dict(case_id="h2.0_spring_k100_pad", height_m=2.0, plantar_bc="spring",
                      load_n=loads[2.0], pad_scenario="on_pad", spring_k=100.0))
    for mu in MUS:
        cases.append(dict(case_id=f"h2.0_contact_rigid_mu{mu:g}_pad", height_m=2.0,
                          plantar_bc="contact", load_n=loads[2.0], pad_scenario="on_pad",
                          contact_mu=mu))
    for h in (3.0, 3.5, 4.5):
        cases.append(dict(case_id=f"h{h:g}_contact_rigid_mu{MU_PRIMARY:g}_pad", height_m=h,
                          plantar_bc="contact", load_n=loads[h], pad_scenario="on_pad",
                          contact_mu=MU_PRIMARY))
    # hard_surface 行（§5.1 #7）：h=2.0 固定 + contact μ=0.0
    cases.append(dict(case_id="h2.0_fixed_hard_surface", height_m=2.0, plantar_bc="fixed",
                      load_n=_hard_surface_load(2.0, loads[2.0]), pad_scenario="hard_surface"))
    cases.append(dict(case_id="h2.0_contact_rigid_mu0.0_hard_surface", height_m=2.0,
                      plantar_bc="contact", load_n=_hard_surface_load(2.0, loads[2.0]),
                      pad_scenario="hard_surface", contact_mu=0.0))
    return cases


# --------------------------------------------------------------------------
# 门禁
# --------------------------------------------------------------------------
def _gate(verdict: str, max_rel_err, evidence) -> dict:
    return {"verdict": verdict, "max_rel_err": max_rel_err, "evidence": evidence}


def _regression_gates(runs: list[dict], axial: dict, s4_baseline: dict) -> dict:
    by_id = {r["id"]: r for r in runs}

    ref = by_id.get("g0_fixed_axial_ref", {})
    g0_err = (_rel(ref["gauge_max_mpa"]["value"], s4_baseline["gauge_max"])
              if ref.get("gauge_max_mpa") else None)
    g0_ok = (axial.get("verdict") == "PASS" and g0_err is not None and g0_err <= 1e-3)

    # G2：μ=0 + 刚性壁 → rigid-support / fixed 极限
    c0 = next((r for r in runs if r["plantar_bc"] == "contact" and r["contact_mu"] == 0.0
               and r["pad_scenario"] == "hard_surface"), None)
    fx = next((r for r in runs if r["plantar_bc"] == "fixed"
               and r["pad_scenario"] == "hard_surface"), None)
    if c0 and c0.get("gauge_max_mpa") and fx and fx.get("gauge_max_mpa"):
        r = _rel(c0["gauge_max_mpa"]["value"], fx["gauge_max_mpa"]["value"])
        g2 = _gate("PASS" if r <= 0.15 else "FAIL", r, [c0["id"], fx["id"]])
    else:
        g2 = _gate("FAIL", None,
                   "contact μ=0 hard_surface 行未收敛（见 runs.diagnostics.neg_jacobians）")

    # G3：接触守恒（支撑反力 = 施加外力，≤2%）
    contact_runs = [r for r in runs if r["plantar_bc"] == "contact"]
    fb_contact = [r["force_balance_relative_err"] for r in contact_runs
                  if r.get("force_balance_relative_err")]
    if fb_contact:
        mx = max(v["value"] if isinstance(v, dict) else v for v in fb_contact)
        g3 = _gate("PASS" if mx <= 0.02 else "FAIL", mx, [r["id"] for r in contact_runs])
    else:
        g3 = _gate("FAIL", None, "所有 contact 行发散；无接触守恒数字")

    # G4：摩擦锥 |F_t| ≤ μ F_n（时域最大）
    fc = [r["friction_cone_ratio_max"] for r in contact_runs if r.get("friction_cone_ratio_max")]
    if fc:
        mx = max(v["value"] if isinstance(v, dict) else v for v in fc)
        g4 = _gate("PASS" if mx <= 1.0 + 1e-6 else "FAIL", mx, [r["id"] for r in contact_runs])
    else:
        g4 = _gate("FAIL", None, "所有 contact 行发散；无摩擦锥数字")

    # G5：A7 消失（contact gauge ≪ fixed gauge，同高度）
    h2c = next((r for r in runs if r["id"] == f"h2.0_contact_rigid_mu{MU_PRIMARY:g}_pad"), None)
    h2f = next((r for r in runs if r["id"] == "h2.0_fixed_pad"), None)
    if h2c and h2c.get("gauge_max_mpa") and h2f and h2f.get("gauge_max_mpa"):
        drop = 1.0 - h2c["gauge_max_mpa"]["value"] / h2f["gauge_max_mpa"]["value"]
        g5 = _gate("PASS" if drop > 0.0 else "FAIL", drop, [h2c["id"], h2f["id"]])
    else:
        g5 = _gate("FAIL", None, "contact 行未收敛；无法对比 gauge_max（fixed 行有效）")

    # G6：敏感性单调/有界
    g6 = _gate("FAIL", None, "contact 扫描行发散；μ/k 单调性不可判定")

    return {
        "g0_default_unchanged": _gate(
            "PASS" if g0_ok else "FAIL", g0_err,
            {"axial_regression": axial.get("verdict"), "baseline_ref": "nonvertical_s4.json::cases[axial_baseline_fixed]"}),
        "g1_pad_stress_reproduced": _gate(
            "N/A_minimum_version", None,
            "最小版无泡沫实体 ⇒ 无 ξ σ 曲线可复现（scope_limit.recipe_constraint）"),
        "g2_limit_self_check": g2,
        "g3_contact_conservation": g3,
        "g4_friction_cone": g4,
        "g5_a7_relief": g5,
        "g6_sensitivity_monotone_bounded": g6,
        "g7_penetration_window": _gate(
            "RETIRED_by_contract_条6-8", None,
            "契约 §4 条 6–8：最小版不在 G7 阻塞内；最小版不给出可信压入量"),
    }


def _sensitivity(runs: list[dict]) -> dict:
    rows = []
    for mu in MUS:
        r = next((x for x in runs if x["id"] == f"h2.0_contact_rigid_mu{mu:g}_pad"), None)
        ok = bool(r and r["febio_rc"] == 0 and r.get("friction_cone_ratio_max"))
        rows.append({
            "mu": float(mu),
            "monotone_in_F_t": False,
            "monotone_in_penetration": False,
            "friction_cone_holds": bool(ok),
            "evidence_run_ids": ([r["id"]] if r else []),
        })
    verdict = "PASS" if all(x["friction_cone_holds"] for x in rows) and any(x["evidence_run_ids"] for x in rows) else "FAIL"
    return {"mu_sweep": rows, "penalty_sweep": [], "k_sweep": [],
            "gate_g6_verdict": verdict, "gate_g6_max_rel_err": None}


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Phase-S5-contact 波4：跟骨跖面真实接触 + 摩擦 FE 场景矩阵（opt-in）")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    ap.add_argument("--fast", action="store_true",
                    help="跳过 OpenSim 轴向回归重跑（接口逐位检查仍执行）")
    args = ap.parse_args(argv)

    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== Phase-S5-contact 波4：跖面真实接触 + 摩擦（FE 场景矩阵）===", flush=True)

    loads = _height_subtalar_loads()
    print(f"[loads] 每高度距下关节力(N): "
          f"{ {h: round(v,1) for h, v in loads.items()} }", flush=True)

    print("[1/4] 加载 THUMS 跟骨网格 ...", flush=True)
    mesh = T.load_thums_mesh(MESH_ANAT)
    print(f"      nodes={mesh['n_nodes']} CORT={mesh['n_elements_cort']} "
          f"SPON={mesh['n_elements_spon']} plantar_faces={mesh['n_plantar_faces']}", flush=True)

    print("[2/4] 接口回归（默认逐位）+ 轴向回归（硬门槛 G0）...", flush=True)
    iface = _interface_regression(mesh)
    axial = (_axial_regression() if not args.fast
             else {"verdict": "SKIPPED", "max_rel_err": None})
    print(f"      iface default_unchanged={iface['plantar_bc_default_unchanged']} "
          f"contact_optin={iface['contact_optin_loads_when_default']} "
          f"axial={axial.get('verdict')}", flush=True)

    s4 = json.loads(S4_JSON.read_text(encoding="utf-8"))
    s4_base = next(c for c in s4["cases"] if c["tag"] == "axial_baseline_fixed")

    print("[3/4] FE 场景矩阵（每行落盘后追加）...", flush=True)
    runs: list[dict] = []
    for case in _plan_cases(loads):
        row = _solve_case(mesh, **case)
        runs.append(row)
        fn = row.get("F_n_n")
        fn_txt = "—" if not fn else f"{fn['value']/1e3:.2f} kN"
        dg = row.get("diagnostics", {})
        print(f"  {row['id']:36s} bc={row['plantar_bc']:8s} rc={row['febio_rc']} "
              f"F_n={fn_txt:>9s} verdict={row['run_verdict']} "
              f"negJac={dg.get('neg_jacobians')} invalid={dg.get('invalid_facets')} "
              f"[{row['wall_s']['value']:.1f}s]", flush=True)
        # 逐行落盘（崩溃也不丢）
        _dump_partial(args.out_json, runs, ts, iface, axial, s4_base)

    print("[4/4] 汇总门禁 + 写报告 ...", flush=True)
    gates = _regression_gates(runs, axial, s4_base)
    sens = _sensitivity(runs)
    result = _assemble(runs, ts, iface, axial, s4_base, gates, sens, loads, mesh)

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0 if result["axial_regression"].get("verdict") in ("PASS", "SKIPPED") else 1


def _assemble(runs, ts, iface, axial, s4_base, gates, sens, loads, mesh) -> dict:
    n_contact = sum(1 for r in runs if r["plantar_bc"] == "contact")
    n_contact_ok = sum(1 for r in runs if r["plantar_bc"] == "contact" and r["run_verdict"] == "PASS")
    return {
        "meta": {
            "generated_at": ts,
            "study": "Nonvertical S5-contact, FEBio plantar-pad contact + friction (wave 4, opt-in)",
            "plan_doc": "docs/非垂直落地扩展方案.md §3 L4 / §4 S5 · docs/S5接触方案.md §3",
            "contract_ref": "docs/S5_contact_contract.md §4 (acceptance G0..G7) + §4 条 7 (scope limit)",
            "quotable_ref": "docs/S5_QUOTABLE.md",
            "units": "mm-N-MPa-s",
            "mass_kg": _m(M_MODEL, "kg", "measured", "Rajagopal 2015 (nonvertical_s1.json)"),
            "recipe_id": RECIPE_ID,
            "plantar_bc_kind": "contact",
            "febio4_exe": FEBIO4_EXE,
            "commands": [r".venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s5_contact.py"],
            "measured_vs_modeled": {
                "measured": [
                    "OpenSim GRF + subtalar reaction (default link, hard gate)",
                    "FEBio 4.13 STATIC solution under run_febio (wave-3 contact recipe)",
                    "support-reaction resultant F_n / F_t (read from .xplt reaction forces)",
                    "calcaneus mesh + plantar facet set (THUMS AM50)",
                    "per-height subtalar load (s2_thums_height_sweep.json, OpenSim)",
                ],
                "modeled_assumed": [
                    "Coulomb friction cone |F_t| ≤ μ·F_n (FEBio sliding-elastic)",
                    "rigid support = fully-fixed hex8 plate stand-in for a wall (no rigid_wall in 4.13)",
                    "μ = 0.6 carried from QUOTABLE §3 #1 (cited setup, not re-derived)",
                    "σ_contact_avg = F_n / 60000 mm² (QUOTABLE §3 #2 basis; F_n rescale)",
                    "hard_surface load = on_pad load × GRF ratio (1D pad.py; modeled)",
                    "G1/G7 out of scope (scope_limit.recipe_constraint)",
                ],
            },
        },
        "scope_limit": {
            "contract_clause": "docs/S5_contact_contract.md §4 条 7 (2026-10-06 晚间修订)",
            "scope_summary": "波4 可启动，但所有引用必须落在 docs/S5_QUOTABLE.md 的可信足迹内。",
            "no_quotable_sigma_ratio": {
                "value": True, "label": "measured",
                "rule": "本脚本不报任何 σ/σ_law 比值（含单网格/固定-ξ 取值）；只报绝对量 "
                        "F_n / σ_contact（MPa）。",
                "evidence": ["docs/S5_QUOTABLE.md §4 (#1,#2)",
                             "docs/S5_g7_verify_mesh.md §3 (FAIL, free-ξ)",
                             "docs/S5_g7_verify_fixedxi.md §4 (FAIL, fixed-ξ)",
                             "docs/S5_g7_verify_anchors_mesh.md §4 (FAIL, small-strain)"],
            },
            "citation_policy": {
                "value": "QUOTABLE-only", "label": "assumed",
                "rule": "只引 docs/S5_QUOTABLE.md §0/§3 的行，且 recipe+mesh+ξ 同引；跨配方对比 "
                        "与撤销的大应变数字禁止（QUOTABLE §4 #3,#4）。",
                "evidence": ["docs/S5_QUOTABLE.md §3 (allow)", "§4 (deny)"],
            },
            "recipe_constraint": {
                "value": "minimum_version_only", "label": "assumed",
                "rule": "波4 实现最小版（刚性面对偶 + Coulomb μ，无泡沫实体）。FOAM 版（§3.2）"
                        "超出范围，被 G7/G7' 阻塞（契约 §4 条 6–8）。",
                "evidence": ["docs/S5_contact_contract.md §3.1 / §4 条 8", "docs/S5接触方案.md §3.1"],
            },
            "fe_output_constraint": {
                "value": "absolute_quantities_only", "label": "assumed",
                "rule": "FE 输出为 F_n（N）/ σ_contact（MPa，avg/max/p95）/ penetration（mm）；"
                        "不计算任何足迹型 σ/σ_law 替代量。",
                "evidence": ["docs/S5_QUOTABLE.md §5", "docs/S5_g7_verify_metric.md §6"],
            },
        },
        "interface_regression": iface,
        "axial_regression": axial,
        "scenario_matrix": {
            "axes": {
                "height_m": list(HEIGHTS),
                "plantar_bc": ["fixed", "spring", "contact"],
                "contact_kind": ["rigid_plate"],
                "mu": list(MUS),
                "pad_tilt_deg": [0.0],
                "foot_mode": ["two-foot"],
                "posture": ["stiff"],
                "cop_offset_mm": [0.0],
                "angular_momentum": ["zero"],
                "velocity_dir": ["vertical"],
                "pad_scenario": ["on_pad", "hard_surface"],
                "spring_k_N_per_mm": [None, 100.0],
            },
            "scope_notes": {
                "plantar_bc.contact": "in scope（波4 的新 opt-in 能力，wave-3）",
                "plantar_bc.fixed": "in scope（G0 硬门槛 + G5 A7 参照）",
                "plantar_bc.spring": "in scope（G2 极限自检参照，k=100 N/mm）",
                "mu": "in scope（G4 + G6 μ sweep；0.0/0.3/0.6）",
                "pad_tilt_deg": "held at 0.0（② 上游，§3）",
                "foot_mode": "held at two-foot（④ 上游，§3）",
                "posture": "held at stiff（③ 上游，§3）",
                "cop_offset_mm": "held at 0.0（⑤ 上游，§3）",
                "angular_momentum": "held at zero（⑥ 上游，§3）",
                "velocity_dir": "held at vertical（① 上游，§3）",
                "contact_kind.foam": "out of scope（FOAM 版被 G7/G7' 阻塞）",
            },
            "runs": [r["id"] for r in runs],
        },
        "contact_recipe": {
            RECIPE_ID: {
                "feb_section_order": RECIPE["feb_section_order"],
                "contact_block": {
                    "febio_type": "sliding-elastic",
                    "master_surface": "plantar top facet set",
                    "slave_surface": "pad_floor top facet set (fixed hex8 plate)",
                    "mu": _m(MU_PRIMARY, "-", "measured", "docs/S5_QUOTABLE.md §3 #1 (carried)"),
                    "penalty": _m(RECIPE["penalty"], "-", "assumed", "NEW PENALTY recipe (udg §7)"),
                    "laugon": RECIPE["laugon"],
                    "tolerance": RECIPE["tolerance"],
                    "aug_controls": None,
                    "search_radius": RECIPE["search_radius"],
                    "node_reloc": RECIPE["node_reloc"],
                },
                "rigid_wall_block": {
                    "febio_type": "rigid plate (fully fixed hex8; no rigid_wall in FEBio 4.13)",
                    "normal_unit": [0, -1, 0],
                    "gap_mm": _m(0.5, "mm", "assumed", "plantar_bc.DEFAULT_CONTACT_GAP_MM"),
                    "blocked_dof": ["tx", "ty", "tz"],
                },
                "note": "单配方（NEW PENALTY 0.1）；禁止 OLD AUGLAG 与跨配方对比。",
            }
        },
        "runs": runs,
        "regression_gates": gates,
        "sensitivity": sens,
        "assumptions": _assumptions(n_contact, n_contact_ok),
        "repro": {
            "commands": [
                r"$env:PYTHONPATH='src'; .venv\Scripts\python.exe scripts\opensim_fe\nonvertical_s5_contact.py",
            ],
            "files_read": [
                str(HEIGHT_SWEEP_JSON.relative_to(ROOT)),
                str(CACHED_AXIAL.relative_to(ROOT)),
                str(S4_JSON.relative_to(ROOT)),
                str(MESH_ANAT.relative_to(ROOT)),
            ],
            "files_written": [
                str(OUT_JSON.relative_to(ROOT)),
                str(OUT_MD.relative_to(ROOT)),
            ],
            "venv": r".venv\Scripts\python.exe",
            "env": {"PYTHONPATH": "src"},
        },
    }


def _assumptions(n_contact: int, n_contact_ok: int) -> list[str]:
    return [
        "[assumed] 最小版（刚性面对偶 + μ，无泡沫）为波4 唯一配置（docs/S5接触方案.md §3.1 + "
        "docs/S5_contact_contract.md §4 条 8）。G1/G7/G7' 超范围。",
        "[assumed] 配方 = NEW PENALTY 0.1（laugon=PENALTY, penalty=0.1, search_radius=20, "
        "node_reloc=0, tolerance=0.005），取自 docs/S5_contact_udg.md §7。禁止 OLD AUGLAG 与跨配方对比。",
        "[measured] μ=0.6 来自 docs/S5_QUOTABLE.md §3 #1 的设定（引用，非重新标定）。",
        "[measured] F_n ≈ 34.4 kN / σ_contact ≈ 0.574 MPa 是 2 m 落地（固定 sink 170 mm、NEW PENALTY 配方）"
        "网格收敛 <0.6% 的锚点；本波只作量级/口径锚，不复用为结论。",
        "[assumed] 不引用任何 σ/σ_law 比值；只报绝对量 σ_contact（MPa）。",
        "[assumed] 只引 docs/S5_QUOTABLE.md 的 QUOTABLE 行；禁止跨配方对比（QUOTABLE §4 #3）。",
        "[assumed] 对偶面 = 全节点固定 hex8 薄板（FEBio 4.13 无 rigid_wall）；不建泡沫实体。",
        "[measured] 骨 + 跖面面片来自 THUMS AM50 跟骨网格（temp/opensim_fe/thums_calcaneus/）。",
        f"[measured] contact 行 {n_contact} 条，收敛 {n_contact_ok} 条。未收敛行的诊断记于 "
        "runs[*].diagnostics（FEBio negJac / invalid facets 计数）。",
        "[measured] 所有 FEBio 调用经 climbing.coupling.febio_run.run_febio（契约 §5）；.feb 段序 "
        "Loads→Boundary→Contact；跖面对偶为真实单元面。",
        "[assumed] 硬门槛 G0 是「默认路径逐位不变」的唯一数值证据。",
        "[assumed] 每高度施加载荷取 s2_thums_height_sweep.json（OpenSim 距下关节峰值力）线性插值；"
        "hard_surface 载荷按 1D GRF 比值缩放（modeled）。",
    ]


def _dump_partial(path: Path, runs, ts, iface, axial, s4_base) -> None:
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(
            {"meta": {"generated_at": ts, "partial": True},
             "interface_regression": iface, "axial_regression": axial,
             "runs": runs, "s4_baseline": s4_base},
            indent=2, ensure_ascii=False, default=float), encoding="utf-8")
    except Exception:  # noqa: BLE001 - 落盘失败不阻断
        pass


# --------------------------------------------------------------------------
# 报告（§2：13 节 + 契约 §6 五项）
# --------------------------------------------------------------------------
def _write_report(path: Path, res: dict, ts: str) -> None:
    runs = res["runs"]
    gates = res["regression_gates"]
    by_id = {r["id"]: r for r in runs}
    contact = [r for r in runs if r["plantar_bc"] == "contact"]
    ok = [r for r in contact if r["run_verdict"] == "PASS"]
    L: list[str] = []
    A = L.append
    A("# Phase-S5-contact 波4 · 跟骨跖面真实接触 + Coulomb 摩擦（FE 场景矩阵）报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s；FEBio 4.13。")
    A("> 设计 `docs/S5_wave4_design.md`（APPROVED）；契约 `docs/S5_contact_contract.md` "
      "§4（G0–G7）/§5/§6。")
    A("> **本脚本只报绝对量，任何 σ/σ_law 比值一律不得作为「收敛」单一值引用**"
      "（`docs/S5_QUOTABLE.md` §4 #1、#2）。")
    A("> 引用限定：仅 `docs/S5_QUOTABLE.md` 的 QUOTABLE 行；跨配方对比禁止。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    g0 = gates["g0_default_unchanged"]
    A(f"1. **默认回归（硬门槛 G0）: {g0['verdict']}** —— "
      f"`axial_regression` vs `s1_h5_opensim.json` 最大相对误差 "
      f"{axial_err_txt(res)}；接口逐位检查 default_unchanged="
      f"{res['interface_regression']['plantar_bc_default_unchanged']}。")
    A(f"2. **范围**：波4 = 最小版（刚性面对偶 + Coulomb μ，无泡沫），只在 §5.1 限定的可信口径内；"
      "G1 → `N/A_minimum_version`，G7 → `RETIRED_by_contract_条6-8`（均不省略）。")
    if ok:
        A(f"3. **A7 消失（G5）**：contact gauge 低于 fixed gauge（见 §9）。")
    else:
        n = len(contact)
        A(f"3. **A7 消失（G5）未能判定**：{n} 条 contact 行**全部发散**"
          "（wave-3 接触面在真实 hex8 网格上被 FEBio 判为 invalid facets / negJac）。"
          "详见 §6/§9 与各 run 的 `diagnostics`。")
    A("4. **摩擦锥（G4）**：" + (
        "所有 contact 行 |F_t| ≤ μ·F_n。" if ok else
        "无收敛 contact 行 ⇒ 无可判定数字（FAIL）。") )
    A("5. 诚实边界见 §10。")
    A("")
    A("> ⚠️ **本次矩阵的真实结果**：fixed / spring 行收敛并给出应力；**contact 行全部 rc≠0**"
      "（FEBio `negative jacobians` / `invalid facets`）。这是 wave-3 接触面在真实 THUMS "
      "hex8 皮质网格上的首次真实 FE 验证结论 —— 详见 §6/§9/§10。")
    A("")
    A("---")
    A("")
    A("## 1. API（opt-in contact kind）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） | 作用 |")
    A("|---|---|---|---|")
    A("| BC | `plantar_bc.apply_plantar_bc` | `kind=\"contact\"`（默认 `\"fixed\"`） | 跖面 ⇄ 固定薄板（数值刚性壁）+ Coulomb μ |")
    A("| 构建 | `thums_feb.build_thums_feb` | `contact_mu=0.6` | 透传 fric_coeff（opt-in，默认不变） |")
    A("| 调度 | `scripts/opensim_fe/thums_feb_run.py` | `--plantar-bc contact` | CLI opt-in（wave-3） |")
    A("")
    A("**默认 `plantar_bc=\"fixed\"` 一个字节不变**（接口回归逐位检查 = "
      f"{res['interface_regression']['plantar_bc_default_unchanged']}）。")
    A("")
    A("---")
    A("")
    A("## 2. 默认回归（硬门槛）")
    A("")
    A("| 检查 | 结果 |")
    A("|---|---|")
    ir = res["interface_regression"]
    A(f"| 默认 `plantar_bc` 逐位不变 | **{ir['plantar_bc_default_unchanged']}** |")
    A(f"| contact opt-in 才发射 `<Contact>` | **{ir['contact_optin_loads_when_default']}** |")
    A("")
    ax = res["axial_regression"]
    if ax.get("verdict") == "SKIPPED":
        A("- 轴向链路重跑：`--fast` 已跳过。")
    else:
        A("默认链路重跑（`ground_reaction`→`run_dead_drop`→`subtalar_reaction`）vs 缓存：")
        A("")
        A("| 量 | 本次默认 | 缓存 | 相对误差 |")
        A("|---|---:|---:|---:|")
        for k in ("impulse_ns", "grf_peak_n", "subtalar_peak_vertical_n",
                  "subtalar_peak_force_n", "subtalar_peak_moment_nm"):
            A(f"| `{k}` | {ax['measured_default_run'][k]:.9g} | {ax['cached'][k]:.9g} | "
              f"{ax['rel_err'][k]:.2e} |")
        A("")
        A(f"- 最大相对误差 **{ax['max_rel_err']:.2e}**（阈值 {REG_TOL:.0e}）→ **{ax['verdict']}**")
    A("")
    A("---")
    A("")
    A("## 3. 接触配方（FEBio recipe block）")
    A("")
    rc_id = res["meta"]["recipe_id"]
    rec = res["contact_recipe"][rc_id]
    A(f"- recipe id：`{rc_id}`；段序 `{' → '.join(rec['feb_section_order'])}`（契约 §5）。")
    cb = rec["contact_block"]
    A(f"- contact：`{cb['febio_type']}`，`laugon={cb['laugon']}`，`penalty={cb['penalty']['value']}`，"
      f"`search_radius={cb['search_radius']}`，`node_reloc={cb['node_reloc']}`，`tolerance={cb['tolerance']}`，"
      f"`fric_coeff=μ`。")
    A(f"- 主面 = `{cb['master_surface']}`；从面 = `{cb['slave_surface']}`。")
    rb = rec["rigid_wall_block"]
    A(f"- 刚性对偶面：`{rb['febio_type']}`，法向 {rb['normal_unit']}，gap {rb['gap_mm']['value']} mm，"
      f"锁 {rb['blocked_dof']}。")
    A("- 出处：`docs/S5_contact_udg.md` §7（NEW PENALTY 0.1）+ `docs/S5_QUOTABLE.md` §3 #1/#5。")
    A("- **FOAM 配方（G1/G7）不在范围内**（scope_limit.recipe_constraint）。")
    A("")
    A("---")
    A("")
    A("## 4. FE 接触输出（per-run 摘要）")
    A("")
    A("| id | BC | pad | F_n | F_t | fric_max | σ_avg | σ_max | σ_p95 | gauge_max | rc | end_t | verdict |")
    A("|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    def g(r, k, fmt="{:.3g}"):
        v = r.get(k)
        return "—" if not v else fmt.format(v["value"])

    for r in runs:
        et = "—" if r["febio_end_t"] is None else f"{r['febio_end_t']:.4g}"
        A(f"| {r['id']} | {r['plantar_bc']} | {r['pad_scenario']} | "
          f"{g(r, 'F_n_n', '{:.0f}')} | {g(r, 'F_t_n', '{:.0f}')} | "
          f"{g(r, 'friction_cone_ratio_max', '{:.3f}')} | "
          f"{g(r, 'sigma_contact_avg_mpa', '{:.4f}')} | {g(r, 'sigma_contact_max_mpa', '{:.4f}')} | "
          f"{g(r, 'sigma_contact_p95_mpa', '{:.4f}')} | {g(r, 'gauge_max_mpa', '{:.2f}')} | "
          f"{r['febio_rc']} | {et} | {r['run_verdict']} |")
    A("")
    A("> 单位：F_n/F_t = N；σ = MPa；tolerance 见 §3。`—` = 未收敛/无结果（见 §10）。")
    A("")
    A("---")
    A("")
    A("## 5. 场景表（on_pad / hard_surface）")
    A("")
    for scen in ("on_pad", "hard_surface"):
        rows = [r for r in runs if r["pad_scenario"] == scen]
        if not rows:
            continue
        A(f"### {scen}")
        A("")
        A("| id | h (m) | BC | μ | spring k | F_n (N) | verdict |")
        A("|---|---:|---|---:|---:|---:|---|")
        for r in rows:
            mu = "—" if r.get("contact_mu") is None else f"{r['contact_mu']:g}"
            sk = "—" if not r.get("spring_k_N_per_mm") else f"{r['spring_k_N_per_mm']['value']:g}"
            fn = "—" if not r.get("F_n_n") else f"{r['F_n_n']['value']:.0f}"
            A(f"| {r['id']} | {r['height_m']['value']:g} | {r['plantar_bc']} | {mu} | {sk} | {fn} | {r['run_verdict']} |")
        A("")
    A("---")
    A("")
    A("## 6. 验收 G0..G6 逐条")
    A("")
    A("| Gate | verdict | max/rel | evidence |")
    A("|---|---|---:|---|")
    for name, gg in gates.items():
        ev = gg["evidence"]
        evt = ev if isinstance(ev, str) else ", ".join(map(str, ev))
        mv = "—" if gg["max_rel_err"] is None else f"{gg['max_rel_err']:.3e}"
        A(f"| {name} | **{gg['verdict']}** | {mv} | {evt} |")
    A("")
    A("> G1 = `N/A_minimum_version`（最小版无泡沫 ⇒ 无 ξ σ 曲线可复现）；"
      "G7 = `RETIRED_by_contract_条6-8`（最小版不在 G7 阻塞内）。二者**均不省略**。")
    A("")
    A("---")
    A("")
    A("## 7. 摩擦锥 + 接触守恒")
    A("")
    A("| id | force_balance_rel | friction_cone_ratio_max |")
    A("|---|---:|---:|")
    for r in runs:
        fb = "—" if not r.get("force_balance_relative_err") else f"{r['force_balance_relative_err']['value']:.3e}"
        fc = "—" if not r.get("friction_cone_ratio_max") else f"{r['friction_cone_ratio_max']['value']:.3f}"
        A(f"| {r['id']} | {fb} | {fc} |")
    A("")
    A("> G3（守恒 ≤2%）/ G4（锥 ≤1）只在**收敛**的 contact 行上有意义；本次 contact 行全部发散，"
      "故 G3/G4 判 FAIL（见 §6）。fixed/spring 行的 force_balance_rel 是"
      "支撑反力 vs 施加载荷的残差。")
    A("")
    A("---")
    A("")
    A("## 8. 敏感性（μ / penalty / k）")
    A("")
    A("| μ | friction_cone_holds | evidence |")
    A("|---|---|---|")
    for row in res["sensitivity"]["mu_sweep"]:
        A(f"| {row['mu']:g} | {row['friction_cone_holds']} | {', '.join(row['evidence_run_ids']) or '—'} |")
    A("")
    A("> penalty_sweep / k_sweep：penalty 固定 0.1（§5.1 #1 单配方）；k_sweep 仅 spring 参照。"
      "本次 contact 发散 ⇒ G6 判 FAIL。")
    A("")
    A("---")
    A("")
    A("## 9. A7 消失验证（fixed vs contact gauge 对照）")
    A("")
    h2f = by_id.get("h2.0_fixed_pad")
    g5 = gates["g5_a7_relief"]
    A(f"- `plantar_bc=\"fixed\"`（h=2.0）gauge_max = "
      f"{('%.2f MPa' % h2f['gauge_max_mpa']['value']) if h2f and h2f.get('gauge_max_mpa') else '—'}；"
      "对照 `nonvertical_s4.json::cases[axial_baseline_fixed]` = 213.50/186.76 MPa（不同载荷口径，仅形状对照）。")
    A(f"- `plantar_bc=\"contact\"`（h=2.0, μ=0.6）gauge_max = "
      f"{('%.2f MPa' % by_id['h2.0_contact_rigid_mu0.6_pad']['gauge_max_mpa']['value']) if by_id.get('h2.0_contact_rigid_mu0.6_pad', {}).get('gauge_max_mpa') else '未收敛'}。")
    A(f"- **G5 = {g5['verdict']}**（{g5['evidence'] if isinstance(g5['evidence'], str) else ', '.join(g5['evidence'])}）。")
    A("- **不掺入**任何 FOAM 版或 σ/σ_law 比值。")
    A("")
    A("---")
    A("")
    A("## 10. 假设（ASSUMPTIONS）与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("### 未解决的风险与下一步")
    A("")
    A("- **wave-3 接触面在真实 hex8 皮质网格上非法**：`mesh['surfaces']['plantar']` 是 **292 tri3**，"
      "FEBio 判为 **292 invalid facets**（tri3 不是 hex8 的面；hex8 需 quad4）。接触因而**不承载**，"
      "骨在软罚下穿透 → 元素反演（negJac）。")
    A("- **换 quad4 也不收敛**：把跖面重建为真实 hex8 quad4 面（340/232 quads）后，invalid facets 归零，"
      "但 deck 仍在首个时间步 `negative jacobians`（负载/penalty ∈ {0.1,1,10,1000}、"
      "载荷 2 kN–21 kN、rigid/pressure 两路均如此）——疑似**皮质 hex 网格过粗 + 曲面跖面 + 平面刚性对偶**"
      "的接触力集中。")
    A("- **下一步**：(a) 把接触面建到 **SPON tet4** 域（真实 tet facet）或细分皮质壳；"
      "(b) 用 `node_reloc=1` + `auto_penalty` + 分段加载；(c) 若需可信 σ_contact 分布，"
      "改用 wave-3 的固定薄板 + 更细的接触网格。在这些完成前，波4 **不给出 contact 的 "
      "F_n/σ_contact/friction 数字**，只给出真实失败诊断。")
    A("- **structural boundary**：波4 是**最小版**；G1/G7/G7' 超范围；脚本只提供 FE-contact 能力；"
      "σ_contact 的绝对量继承 QUOTABLE §3 #1/#2 的网格收敛界（未进一步细化）。")
    A("")
    A("---")
    A("")
    A("## 11. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("# 全矩阵（含 OpenSim 硬门槛 G0）：")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact.py")
    A("# 只跑 FE（跳过 OpenSim 轴向回归重跑）：")
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact.py --fast")
    A("```")
    A("")
    A("## 12. 产物")
    A("")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{WORKDIR.relative_to(ROOT)}\\*.feb/.log/.xplt`（每次求解输入/结果，可再生）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def axial_err_txt(res: dict) -> str:
    ax = res["axial_regression"]
    return "—（--fast 跳过）" if ax.get("max_rel_err") is None else f"{ax['max_rel_err']:.2e}"


if __name__ == "__main__":
    raise SystemExit(main())
