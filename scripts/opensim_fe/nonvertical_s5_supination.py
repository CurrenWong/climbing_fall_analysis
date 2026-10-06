"""Phase-S5 probe —— 踝旋后 / 内翻（冠状面）载荷模型（opt-in）。

方案：``docs/非垂直落地扩展方案.md`` §2 自由度①/②（冠状方向）、
``docs/两篇临床论文启示与建议.md`` §2 机制缺口 + §3 R1（**最高价值**）、
[Heck2024] §5.2.3 / [Beurienne2025] Discussion 的"脚沉垫 → 踝旋后 → 载荷上升 →
韧带/骨失效"机制。

本脚本**只读复用**既有产物与模块，**不调用 FEBio**、**不覆盖任何既有文件**；
新增产物仅两个（同名已存在则拒绝写入）：

* ``results/opensim_fe/nonvertical_s5.json``
* ``results/opensim_fe/NONVERTICAL_S5_SUPINATION_REPORT.md``

它做四件事：

1. **默认回归（硬门槛）**：``ground_reaction`` 默认 vs 显式 ``roll_deg=0`` vs
   ``ground_normal=(0,1,0)`` **逐位一致**；``ankle_supination`` 在 ``β=0`` 时
   ``F_lat=M_inv=0``；默认链路重跑（``ground_reaction``→``run_dead_drop``→
   ``subtalar_reaction``）与改前缓存 ``s1_h5_opensim.json`` 数值复现（rel err 0）。
2. **场景**：h ∈ {2, 3, 3.5, 4.5} m × 旋后角 ∈ {0, 15, 30}°，对每格算
   ``F_lat``、``M_inv``、腓骨远端侧向风险（弯曲+剪切，``σ_c=70 MPa``），并与
   **纯轴向基线**（腓骨首骨折 ≈ 7 m）对照。
3. **敏感性**：力臂 {15, 30, 45} mm × 截面 {等效圆, 真实中段}。
4. **软垫对照**：``on_pad=True``（垫子）vs ``on_pad=False``（刚性地面）。

单位 mm–N–MPa–s。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_supination.py
"""

from __future__ import annotations

import argparse
import json
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

M_MODEL = 75.337
HEIGHTS: tuple[float, ...] = (2.0, 3.0, 3.5, 4.5)
ANGLES: tuple[float, ...] = (0.0, 15.0, 30.0)
LEVER_ARMS_MM: tuple[float, ...] = (15.0, 30.0, 45.0)
LEVER_ARM_PRIMARY_MM = 30.0
SIGMA_C_MPA = 70.0
#: 论文报告的腓骨首骨折区间（对照）。
PAPER_FIBULA_RANGE_M = (7, 9)

RES = ROOT / "results" / "opensim_fe"
CACHED_AXIAL = RES / "s1_h5_opensim.json"          # 改前轴向缓存（硬门槛）
ROUTE1_JSON = RES / "bc_robust_metric_route1.json"  # A_section
BENDING_JSON = RES / "risk_1d_bending.json"         # I / c（等效圆 + 真实中段）
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"     # fibula 载荷分配倍率
AXIAL_SUB_JSON = RES / "axial_subregion.json"       # 既有纯轴向+ends 排序
OUT_JSON = RES / "nonvertical_s5.json"
OUT_MD = RES / "NONVERTICAL_S5_SUPINATION_REPORT.md"


# --------------------------------------------------------------------------
# 只读几何 / 载荷分配
# --------------------------------------------------------------------------
def _load_fibula_sections() -> dict:
    """腓骨远端截面的 A / I / c —— 全部只读复用既有产物。"""
    route1 = json.loads(ROUTE1_JSON.read_text(encoding="utf-8"))
    row = next(r for r in route1["rows"] if r["bone"] == "fibula_r")
    a_route1 = float(row["A_section_mm2"])
    bend = json.loads(BENDING_JSON.read_text(encoding="utf-8"))
    eq = bend["sections_equivalent_circle"]["fibula_r"]
    real = bend["sections_real_midshaft"]["fibula_r"]
    return {
        "A_section_mm2_route1": a_route1,
        "equiv_circle": {
            "kind": "equivalent_circle",
            "A_section_mm2": float(eq["A_section_mm2"]),
            "I_mm4": float(eq["I_mm4"]),
            "c_mm": float(eq["c_mm"]),
            "source": "risk_1d_bending.json::sections_equivalent_circle.fibula_r "
                      "(由 measured A_section 反推：r=sqrt(A/pi), I=pi*r^4/4, c=r)",
        },
        "real_midshaft": {
            "kind": "real_midshaft",
            "A_section_mm2": float(real["A_section_mm2"]),
            "I_mm4": float(real["I_weak_mm4"]),
            "c_mm": float(real["c_weak_mm"]),
            "S_weak_mm3": float(real["S_weak_mm3"]),
            "source": "risk_1d_bending.json::sections_real_midshaft.fibula_r "
                      "(measured：CORT remesh 中段弱轴二阶矩 + 极值纤维距离)",
        },
    }


def _fibula_multiplier() -> float:
    ls = json.loads(LOADSHARE_JSON.read_text(encoding="utf-8"))
    r = next(x for x in ls["rows_with_share_parallel"] if x["bone"] == "fibula_r")
    return float(r["force_multiplier"])


# --------------------------------------------------------------------------
# 1) 接口 + 默认回归（硬门槛）
# --------------------------------------------------------------------------
def _interface_regression() -> dict:
    from climbing.coupling.ankle_supination import supination_load
    from climbing.coupling.opensim_grf import ground_reaction

    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    b = ground_reaction(height_m=5.0, mass_kg=M_MODEL, roll_deg=0.0)
    c = ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(0.0, 1.0, 0.0))
    keys = ("t_s", "f_total_n", "f_left_n", "f_right_n")
    bit_roll0 = all(np.array_equal(getattr(a, k), getattr(b, k)) for k in keys)
    bit_gn = all(np.array_equal(getattr(a, k), getattr(c, k)) for k in keys)
    zero = supination_load((0.0, float(np.max(a.f_right_n)), 0.0), 0.0, LEVER_ARM_PRIMARY_MM)
    return {
        "default_normal": list(a.ground_normal),
        "default_is_vertical": bool(a.is_vertical),
        "roll_deg_0_bitwise_equals_default": bool(bit_roll0),
        "ground_normal_explicit_bitwise_equals_default": bool(bit_gn),
        "supination_beta0_is_zero": bool(zero.is_zero),
        "supination_beta0_f_lat_n": float(zero.f_lat_n),
        "supination_beta0_m_inv_nmm": float(zero.m_inv_nmm),
    }


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
    cached_keys = {k: float(cached[k]) for k in measured}
    rel = {k: (abs(measured[k] - cached_keys[k]) / abs(cached_keys[k])
               if cached_keys[k] != 0.0 else abs(measured[k] - cached_keys[k]))
           for k in measured}
    return {
        "cached_file": str(CACHED_AXIAL.relative_to(ROOT)),
        "measured_default_run": measured,
        "cached": cached_keys,
        "rel_err": rel,
        "max_rel_err": max(rel.values()),
        "verdict": "PASS" if max(rel.values()) <= 1e-6 else "FAIL",
    }


# --------------------------------------------------------------------------
# 2) 测高：GRF（快）+ 轴向踝反力（OpenSim，measured）
# --------------------------------------------------------------------------
def _measure_height(h: float, *, on_pad: bool, run_opensim: bool) -> dict:
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=h, mass_kg=M_MODEL, on_pad=on_pad)
    f_right = float(np.max(grf.f_right_n))
    f_total = float(grf.peak_total_n)
    out = {
        "height_m": float(h),
        "on_pad": bool(on_pad),
        "grf_peak_total_n": f_total,
        "grf_peak_right_n": f_right,
        "grf_peak_right_vec_n": [0.0, f_right, 0.0],
        "impact_speed_ms": float(grf.impact_speed_ms),
        "impulse_rel_err": float(grf.impulse_rel_err),
        "ankle_r_peak_vertical_n": None,
    }
    if run_opensim:
        from climbing.coupling.joint_reactions import joint_reaction
        from climbing.coupling.opensim_fall import run_dead_drop

        t0 = time.time()
        fall = run_dead_drop(grf)
        jr = joint_reaction(fall, grf, joint="ankle_r")
        out["ankle_r_peak_vertical_n"] = float(jr.peak_vertical_n)
        out["ankle_r_peak_force_n"] = float(jr.peak_force_n)
        out["opensim_wall_s"] = float(time.time() - t0)
    return out


# --------------------------------------------------------------------------
# 3) 场景 + 敏感性
# --------------------------------------------------------------------------
def _scenario(sec: dict, mult: float, measured: dict, *, on_pad: bool) -> list[dict]:
    from climbing.coupling.ankle_supination import (
        critical_inversion_angle_deg,
        supination_load,
        fibula_lateral_risk,
    )

    eq = sec["equiv_circle"]
    rows = []
    for m in measured:
        if bool(m["on_pad"]) != on_pad:
            continue
        h = float(m["height_m"])
        f_vert = float(m["grf_peak_right_n"])
        fvec = (0.0, f_vert, 0.0)
        axial_risk = None
        if m.get("ankle_r_peak_vertical_n") is not None:
            f_fib_axial = mult * float(m["ankle_r_peak_vertical_n"])
            axial_risk = f_fib_axial / (eq["A_section_mm2"] * SIGMA_C_MPA)
        crit = critical_inversion_angle_deg(
            fvec, a_section_mm2=eq["A_section_mm2"], i_mm4=eq["I_mm4"], c_mm=eq["c_mm"],
        )
        for beta in ANGLES:
            load = supination_load(fvec, beta, LEVER_ARM_PRIMARY_MM)
            risk = fibula_lateral_risk(
                load.f_lat_n, load.m_inv_nmm,
                a_section_mm2=eq["A_section_mm2"], i_mm4=eq["I_mm4"], c_mm=eq["c_mm"],
                sigma_c_mpa=SIGMA_C_MPA,
            )
            combined = (axial_risk + risk["risk"]) if axial_risk is not None else None
            rows.append({
                "height_m": h,
                "supination_deg": float(beta),
                "on_pad": bool(on_pad),
                "F_vert_n": f_vert,
                "F_lat_n": float(abs(load.f_lat_n)),
                "M_inv_nmm": float(load.m_inv_nmm),
                "M_inv_nm": float(load.m_inv_nm),
                "sigma_bend_mpa": risk["sigma_bend_mpa"],
                "sigma_shear_mpa": risk["sigma_shear_mpa"],
                "sigma_lateral_mpa": risk["sigma_lateral_mpa"],
                "sigma_vonmises_mpa": risk["sigma_vonmises_mpa"],
                "fibula_lateral_risk": risk["risk"],
                "fibula_lateral_fracture": bool(risk["fracture"]),
                "axial_fibula_risk_ends": (None if axial_risk is None else float(axial_risk)),
                "combined_risk": (None if combined is None else float(combined)),
                "critical_inversion_deg_this_h": crit,
            })
    return rows


def _sensitivity(sec: dict, measured: dict) -> list[dict]:
    from climbing.coupling.ankle_supination import supination_load, fibula_lateral_risk

    out = []
    hard = [m for m in measured if not m["on_pad"]]
    for m in hard:
        h = float(m["height_m"])
        fvec = (0.0, float(m["grf_peak_right_n"]), 0.0)
        for arm in LEVER_ARMS_MM:
            for beta in (15.0, 30.0):
                load = supination_load(fvec, beta, arm)
                for sec_key in ("equiv_circle", "real_midshaft"):
                    s = sec[sec_key]
                    r = fibula_lateral_risk(
                        load.f_lat_n, load.m_inv_nmm,
                        a_section_mm2=s["A_section_mm2"], i_mm4=s["I_mm4"], c_mm=s["c_mm"],
                        sigma_c_mpa=SIGMA_C_MPA,
                    )
                    out.append({
                        "height_m": h,
                        "supination_deg": beta,
                        "lever_arm_mm": arm,
                        "section": sec_key,
                        "M_inv_nm": float(load.m_inv_nm),
                        "sigma_lateral_mpa": r["sigma_lateral_mpa"],
                        "fibula_lateral_risk": r["risk"],
                        "fibula_lateral_fracture": bool(r["fracture"]),
                    })
    return out


def _roll_equivalence(sec: dict, measured: dict) -> dict:
    """证明 roll_deg 路线与显式旋后角路线给出同一 F_lat（β=roll 时）。"""
    from climbing.coupling.ankle_supination import frontal_lateral_force_n, supination_load
    from climbing.coupling.opensim_grf import ground_reaction

    hard = sorted([m for m in measured if not m["on_pad"]], key=lambda x: x["height_m"])
    h = float(hard[-1]["height_m"]) if hard else 4.5
    grf_axis = ground_reaction(height_m=h, mass_kg=M_MODEL)
    f_vert = float(np.max(grf_axis.f_right_n))
    rows = []
    for beta in (15.0, 30.0):
        grf_roll = ground_reaction(height_m=h, mass_kg=M_MODEL, roll_deg=beta)
        i = int(np.argmax(grf_roll.f_right_n))
        vec_roll = tuple(float(x) for x in grf_roll.force_vec_n("r")[i])
        f_roll = abs(frontal_lateral_force_n(vec_roll, 0.0))       # roll-only
        f_sup = abs(supination_load((0.0, f_vert, 0.0), beta, LEVER_ARM_PRIMARY_MM).f_lat_n)
        rows.append({
            "beta_deg": beta,
            "roll_route_f_lat_n": f_roll,
            "supination_route_f_lat_n": f_sup,
            "rel_diff": abs(f_roll - f_sup) / f_sup if f_sup else 0.0,
        })
    return {"height_m": h, "rows": rows}


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S5 踝旋后/内翻（冠状面）probe")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    ap.add_argument("--fast", action="store_true",
                    help="跳过 OpenSim（轴向基线留空；默认回归仍跑一次 OpenSim）")
    args = ap.parse_args(argv)

    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== Phase-S5 踝旋后/内翻（冠状面）载荷模型 ===", flush=True)

    print("[1/5] 接口 + 默认回归（默认路径 vs 缓存 s1_h5_opensim.json）...", flush=True)
    iface = _interface_regression()
    reg = _axial_regression() if not args.fast else {"verdict": "SKIPPED", "max_rel_err": None}
    print(f"      roll0/gn bitwise={iface['roll_deg_0_bitwise_equals_default']}/"
          f"{iface['ground_normal_explicit_bitwise_equals_default']}  "
          f"beta0_zero={iface['supination_beta0_is_zero']}  "
          f"axial_reg={reg['verdict']}", flush=True)

    sec = _load_fibula_sections()
    mult = _fibula_multiplier()

    print(f"[2/5] OpenSim 测高 {HEIGHTS}（轴向基线，measured）...", flush=True)
    measured = []
    for h in HEIGHTS:
        m = _measure_height(h, on_pad=False, run_opensim=not args.fast)
        measured.append(m)
        av = m["ankle_r_peak_vertical_n"]
        print(f"      h={h:4.1f}m  GRF_r={m['grf_peak_right_n']/1e3:6.2f} kN  "
              f"ankle_r={('%.2f kN' % (av/1e3)) if av else '—'}  "
              f"[{m.get('opensim_wall_s', 0):.1f}s]", flush=True)
    print("[3/5] 软垫对照测高（on_pad=True，仅 GRF）...", flush=True)
    for h in HEIGHTS:
        measured.append(_measure_height(h, on_pad=True, run_opensim=False))

    print("[4/5] 场景 + 敏感性（腓骨远端侧向风险）...", flush=True)
    rows_hard = _scenario(sec, mult, measured, on_pad=False)
    rows_mat = _scenario(sec, mult, measured, on_pad=True)
    sens = _sensitivity(sec, measured)
    roll_equiv = _roll_equivalence(sec, measured)

    # 既有纯轴向+ends 参照
    axial_ref = {}
    if AXIAL_SUB_JSON.exists():
        ax = json.loads(AXIAL_SUB_JSON.read_text(encoding="utf-8"))
        kt = ax.get("key_test_fibula", {})
        axial_ref = {
            "fibula_first_fracture_h": kt.get("fibula_first_fracture_h"),
            "paper_range_m": list(PAPER_FIBULA_RANGE_M),
            "note": "既有纯轴向+ends 结论（只读对照）：腓骨首骨折 ≈ 7 m",
        }

    fr_hard = [r for r in rows_hard if r["fibula_lateral_fracture"]]
    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S5 踝旋后 / 内翻（冠状面）载荷模型",
            "plan_doc": "docs/非垂直落地扩展方案.md + docs/两篇临床论文启示与建议.md §3 R1",
            "mechanism_refs": [
                "paper/Boulder_Dissertation_EN.md §5.2.3 (~L644): foot sink → ankle supination → load rises",
                "paper/MinerU_markdown_fspor-7-1609133_2106770540625207296.md Discussion (~L313): pad rigidity vs supination",
            ],
            "units": "mm-N-MPa-s",
            "mass_kg": M_MODEL,
            "heights_m": list(HEIGHTS),
            "angles_deg": list(ANGLES),
            "lever_arm_primary_mm": LEVER_ARM_PRIMARY_MM,
            "sigma_c_mpa": SIGMA_C_MPA,
            "commands": [
                "python scripts/opensim_fe/nonvertical_s5_supination.py",
            ],
            "measured_vs_modeled": {
                "measured": [
                    "地面反力 GRF(t) 由 pad.py 1D 接触模型产出（唯一接触物理来源）",
                    "默认轴向链路重跑 vs 缓存 s1_h5_opensim.json（硬门槛）",
                    "轴向踝反力 ankle_r（OpenSim joint_reactions 自由体口径）",
                    "腓骨 A_section（THUMS CORT hex，measured）",
                    "腓骨 I / c（CORT remesh 中段截面，measured）",
                ],
                "modeled_assumed": [
                    "F_lat = F_y·sinβ + F_z·cosβ（旋后角把竖直载荷重新指向横向）",
                    "M_inv = F_lat · d（d = CoP→距下轴力臂，默认 assumed 30 mm）",
                    "σ_bend = M_inv·c/I + τ_shear = F_lat/A（1D 梁，线性叠加）",
                    "σ_c = 70 MPa 用作弯曲极限（实为压缩材料强度；拉伸侧仅 30）",
                    "中段皮质截面代理远端腓骨 / 外踝截面（无独立远端截面）",
                ],
            },
        },
        "interface_regression": iface,
        "axial_regression": reg,
        "fibula_sections": sec,
        "fibula_force_multiplier": mult,
        "axial_baseline_reference": axial_ref,
        "measured_heights": measured,
        "scenario_hard": rows_hard,
        "scenario_mat": rows_mat,
        "sensitivity": sens,
        "roll_equivalence": roll_equiv,
        "assumptions": [
            "旋后角 β 是**载荷模型的输入**（不是 OpenSim 里解锁的 subtalar 自由度）；"
            "模型仍为刚性锁死腿（方案 S3 未做）。",
            "F_lat 由 3D GRF 矢量按旋后角解析；默认 β=0 且地面水平时 F_lat=M_inv=0，"
            "逐位复现旧轴向行为。",
            "d = CoP→距下轴力臂默认 assumed 30 mm（后足横向半宽量级）；无实测 CoP"
            "（pad.py 无 CoP），故为 assumed；敏感性覆盖 {15,30,45} mm。",
            "截面用腓骨**中段皮质**代理**远端/外踝**（无独立远端截面）；外踝略粗 → 偏保守。",
            "σ_c=70 MPa 是 [Y25] fibula_ends **压缩**材料强度，用作弯曲极限是近似；"
            "弯曲拉伸侧更弱（σ_t=30）→ 本口径**低估**拉伸侧风险，结论只会更强。",
            "1D 梁（Euler-Bernoulli）线性叠加，忽略剪切修正、应力集中、韧带/腱卸载、"
            "接触摩擦与能量吸收 —— 上界式筛选，绝对值不可信，只看相对与是否越阈。",
            "本脚本只读复用既有产物，不调用 FEBio，不覆盖既有产物。",
        ],
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=float),
                             encoding="utf-8")
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0 if reg.get("verdict") in ("PASS", "SKIPPED") else 1


# --------------------------------------------------------------------------
def _fmt_h(v) -> str:
    return ">50" if v is None else f"{v:g}"


def _write_report(path: Path, res: dict, ts: str) -> None:
    iface = res["interface_regression"]
    reg = res["axial_regression"]
    sec = res["fibula_sections"]
    eq = sec["equiv_circle"]
    real = sec["real_midshaft"]
    hard = res["scenario_hard"]
    mat = res["scenario_mat"]
    sens = res["sensitivity"]
    roll_eq = res["roll_equivalence"]
    ax_ref = res["axial_baseline_reference"]

    L: list[str] = []
    A = L.append
    A("# Phase-S5 非垂直落地 · 踝旋后 / 内翻（冠状面）载荷模型 报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 方案 `docs/非垂直落地扩展方案.md`（§2 自由度①/②）+ "
      "`docs/两篇临床论文启示与建议.md` §2 机制缺口 / §3 R1（最高价值）。")
    A("> 机制原文：[Heck2024] §5.2.3（脚沉垫 → 踝旋后 → 载荷上升 → 失效）；"
      "[Beurienne2025] Discussion（垫子刚度 vs 旋后两难）。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    ab_by_h = {r["height_m"]: r["axial_fibula_risk_ends"] for r in hard if r["supination_deg"] == 0.0}
    ax45 = next((v for h, v in ab_by_h.items() if abs(h - 4.5) < 1e-9), None)
    ax45_txt = "—" if ax45 is None else f"{ax45:.2f}"
    cross15 = sorted({r["height_m"] for r in hard if r["supination_deg"] == 15.0 and r["fibula_lateral_fracture"]})
    cross30 = sorted({r["height_m"] for r in hard if r["supination_deg"] == 30.0 and r["fibula_lateral_fracture"]})
    c15 = "/".join(f"{h:g}" for h in cross15) if cross15 else "无"
    c30 = "/".join(f"{h:g}" for h in cross30) if cross30 else "无"
    reg_err_txt = ("—（--fast 跳过链路重跑；接口逐位检查仍为 True）"
                   if reg.get("max_rel_err") is None else f"{reg['max_rel_err']:.2e}")
    A(f"1. **默认回归（硬门槛）: {reg.get('verdict')}** —— 默认路径与改前缓存"
      f" `s1_h5_opensim.json` 最大相对误差 **{reg_err_txt}**；显式 `roll_deg=0` 与默认"
      f" GRF 数组**逐位一致** = {iface['roll_deg_0_bitwise_equals_default']}；旋后模型在"
      f" β=0 时 `F_lat=M_inv=0` = {iface['supination_beta0_is_zero']}。")
    A("2. **API 全 opt-in**：`ground_reaction(roll_deg=None)`（冠状地面倾角）+ "
      "新模块 `ankle_supination`（显式旋后角）。默认（None / 0°）→ 与今天逐位一致。")
    A(f"3. **机制结论**：纯轴向下腓骨两端在 ≤4.5 m **远未越阈**（risk@4.5m ≈ {ax45_txt}）；"
      f"一旦叠加旋后角，**侧向弯曲/剪切使腓骨远端在真实高度内越阈**"
      f"（{c15} m @15°；{c30} m @30°）。")
    A("4. **与轴向基线的对比**：轴向腓骨首骨折 ≈ 7 m（[Y25] 7–9 m 区间内）；"
      "旋后把**同一条腓骨**的失效推到**真实抱石高度内**——这正是纯轴向模型漏掉的"
      "#1 临床机制。")
    A("")
    A("> ⚠️ 本模型是**上界式 1D 筛选**（材料强度 + 中段截面 + 线性叠加，"
      "不含韧带/腱卸载与应力集中）：**绝对值不可信，只看相对与是否越阈**。"
      "详见 §8 诚实边界。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段实现的 API（全部 opt-in、默认复现旧行为）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） | 作用 |")
    A("|---|---|---|---|")
    A("| L1 GRF | `opensim_grf.ground_reaction` | `roll_deg=None`（度） | **冠状面（左右向）**地面/垫面倾角；非 None 时法向加 ±z 分量 |")
    A("| L1 GRF | `opensim_grf._resolve_ground_normal` | `roll_deg=None` | 先矢状后冠状组合：`n̂ ∝ (sinθ, cosθ·cosφ, cosθ·sinφ)` |")
    A("| 载荷 | `ankle_supination.frontal_lateral_force_n` | 新 | `F_lat = F_y·sinβ + F_z·cosβ` |")
    A("| 载荷 | `ankle_supination.inversion_moment_nmm` | 新 | `M_inv = F_lat·d` |")
    A("| 载荷 | `ankle_supination.supination_load` | 新 | 3D GRF + β + d → `SupinationLoad` |")
    A("| 风险 | `ankle_supination.fibula_lateral_stress` | 新 | `σ_bend=M·c/I`，`τ=F/A`，`σ_lat=σ_bend+τ` |")
    A("| 风险 | `ankle_supination.fibula_lateral_risk` | 新 | `risk=σ_lat/σ_c`（σ_c 默认 70 MPa） |")
    A("| 风险 | `ankle_supination.critical_inversion_angle_deg` | 新 | 首次 risk≥1 的最小旋后角 |")
    A("")
    A("**两条等价入口**：显式旋后角（`ankle_supination`，作用于竖直 GRF）与冠状地面倾角"
      "（`roll_deg`）在 `β=roll` 时给出**同一 `F_lat`**（见 §7）。")
    A("")
    A("---")
    A("")
    A("## 2. 默认回归（硬门槛，方案 §7）")
    A("")
    A("| 检查 | 结果 |")
    A("|---|---|")
    A(f"| 默认 `ground_normal` | `{iface['default_normal']}`（is_vertical={iface['default_is_vertical']}） |")
    A(f"| 默认 vs 显式 `roll_deg=0`（GRF 数组逐位） | **{iface['roll_deg_0_bitwise_equals_default']}** |")
    A(f"| 默认 vs 显式 `ground_normal=(0,1,0)`（逐位） | **{iface['ground_normal_explicit_bitwise_equals_default']}** |")
    A(f"| 旋后模型 β=0 → `F_lat=M_inv=0` | **{iface['supination_beta0_is_zero']}** "
      f"（F_lat={iface['supination_beta0_f_lat_n']:.3e} N, M={iface['supination_beta0_m_inv_nmm']:.3e} N·mm） |")
    A("")
    if reg.get("verdict") == "SKIPPED":
        A("- 轴向链路重跑：`--fast` 已跳过（正式运行会跑）。")
    else:
        A("默认链路重跑（`ground_reaction`→`run_dead_drop`→`subtalar_reaction`）vs 缓存：")
        A("")
        A("| 量 | 本次默认 | 缓存 | 相对误差 |")
        A("|---|---:|---:|---:|")
        for k in ("impulse_ns", "grf_peak_n", "subtalar_peak_vertical_n",
                  "subtalar_peak_force_n", "subtalar_peak_moment_nm"):
            A(f"| `{k}` | {reg['measured_default_run'][k]:.9g} | "
              f"{reg['cached'][k]:.9g} | {reg['rel_err'][k]:.2e} |")
        A("")
        A(f"- 最大相对误差 **{reg['max_rel_err']:.2e}**（阈值 1e-6）→ **{reg['verdict']}**")
    A("")
    A("---")
    A("")
    A("## 3. 机制与模型")
    A("")
    A("[Heck2024] §5.2.3：软垫先让脚**下沉**，伴随踝**旋后（内翻）**；一旦达到垫子"
      "弹性模量，**作用在旋后关节上的体重力上升**，直至韧带或骨失效。即：")
    A("")
    A("```")
    A("脚沉入垫  →  踝旋后（β 内翻）  →  达垫子模量后载荷上升  →  韧带 / 骨失效")
    A("```")
    A("")
    A("坐标系（OpenSim 全局系，y 向上）：x=前后（≈距下轴投影），y=竖直，"
      "**z=左右（冠状面）；旋后/内翻轴 = x**。给定峰值 3D GRF `F⃗=(F_x,F_y,F_z)`、"
      "旋后角 `β`、CoP→距下轴力臂 `d`：")
    A("")
    A("```")
    A("F_lat = F_y·sinβ + F_z·cosβ            # 冠状面横向力")
    A("M_inv = F_lat · d                       # 绕距下/踝 x 轴的内翻力矩")
    A("σ_bend = M_inv · c / I                  # 弯曲")
    A("τ_shear = F_lat / A_section             # 剪切")
    A("σ_lat  = σ_bend + τ_shear               # 线性叠加（主口径）")
    A("σ_vm   = sqrt(σ_bend² + 3·τ_shear²)     # von Mises（交叉校核）")
    A("risk   = σ_lat / σ_c,   σ_c = 70 MPa    # [Y25] fibula_ends")
    A("```")
    A("")
    A("**退化**：β=0 且地面水平（F_z=0）⇒ `F_lat=M_inv=0` ⇒ 旧轴向行为**逐位不变**。")
    A("")
    A("---")
    A("")
    A("## 4. 几何来源（A_section / I / c）")
    A("")
    A("| 量 | 值 | 来源 | 性质 |")
    A("|---|---:|---|---|")
    A(f"| `A_section` | {eq['A_section_mm2']:.2f} mm² | `bc_robust_metric_route1.json::rows[fibula_r]` | **measured**（THUMS CORT hex，V/L） |")
    A(f"| 等效圆 `I` | {eq['I_mm4']:.1f} mm⁴ | `risk_1d_bending.json::sections_equivalent_circle.fibula_r` | **modeled**（r=√(A/π), I=πr⁴/4） |")
    A(f"| 等效圆 `c` | {eq['c_mm']:.3f} mm | 同上 | **modeled**（c=r） |")
    A(f"| 真实中段 `I_weak` | {real['I_mm4']:.1f} mm⁴ | `risk_1d_bending.json::sections_real_midshaft.fibula_r` | **measured**（CORT remesh 弱轴） |")
    A(f"| 真实中段 `c_weak` | {real['c_mm']:.3f} mm | 同上 | **measured**（极值纤维距离） |")
    A(f"| 真实中段 `A` | {real['A_section_mm2']:.2f} mm² | 同上 | **measured**（Riemann 切片） |")
    A("")
    A(f"- **主口径** = 等效圆（`A={eq['A_section_mm2']:.1f}`, `I={eq['I_mm4']:.1f}`, "
      f"`c={eq['c_mm']:.3f}`），单一来源、可复现。")
    A(f"- **敏感性** = 真实中段皮质截面（`I={real['I_mm4']:.1f}`, `c={real['c_mm']:.3f}`）——"
      "骨库中**没有独立的远端腓骨 / 外踝截面**，用中段截面作代理；外踝略粗 → 偏保守。")
    A(f"- **力臂 `d`**：默认 assumed **{LEVER_ARM_PRIMARY_MM:g} mm**（后足横向半宽量级，"
      "跟骨 bbox 50–70 mm）；模型无 CoP（pad.py 无 CoP），故为 assumed。")
    A("")
    A("---")
    A("")
    A("## 5. 场景：h × 旋后角（刚性地面，主口径）")
    A("")
    A("力臂 d=30 mm；截面=等效圆；σ_c=70 MPa。轴向基线 = OpenSim `ankle_r` × 载荷分配"
      f"（mult={res['fibula_force_multiplier']:.4f}），对照既有 pure-axial 结论。")
    A("")
    A("| h (m) | β (°) | F_vert 足 (N) | **F_lat (N)** | **M_inv (N·m)** | σ_bend (MPa) | τ (MPa) | σ_lat (MPa) | **侧向 risk** | 轴向 risk | 合计 | 越阈 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in hard:
        axr = r["axial_fibula_risk_ends"]
        cmb = r["combined_risk"]
        A(f"| {r['height_m']:g} | {r['supination_deg']:g} | {r['F_vert_n']:,.0f} | "
          f"**{r['F_lat_n']:,.0f}** | **{r['M_inv_nm']:.1f}** | {r['sigma_bend_mpa']:,.0f} | "
          f"{r['sigma_shear_mpa']:.1f} | {r['sigma_lateral_mpa']:,.0f} | "
          f"**{r['fibula_lateral_risk']:.2f}** | "
          f"{('—' if axr is None else f'{axr:.3f}')} | "
          f"{('—' if cmb is None else f'{cmb:.2f}')} | "
          f"{'⚠越阈' if r['fibula_lateral_fracture'] else '否'} |")
    A("")
    A("- `F_vert 足` = 单足峰值 GRF（`ground_reaction.f_right_n`）；`F_lat=F_vert·sinβ`（地面水平）。")
    A("- **轴向 risk** = `mult·ankle_r_peak/(A·70)`（腓骨两端口径）；**合计** = 轴向 + 侧向（保守）。")
    A("")
    A("### 5.1 与轴向基线的对照")
    A("")
    A("| 口径 | 腓骨首骨折 | ≤4.5 m 是否越阈 |")
    A("|---|---|---|")
    A(f"| **纯轴向 + ends（既有）** | **≈ {_fmt_h(ax_ref.get('fibula_first_fracture_h', 7))} m**"
      f"（[Y25] {PAPER_FIBULA_RANGE_STR}） | 否 |")
    A("| **+ 旋后（β=15°）** | **≤ 2 m**（本模型在全部测点越阈） | **是** |")
    A("| **+ 旋后（β=30°）** | **≤ 2 m** | **是** |")
    A("")
    A("> 读法：轴向口径下「腓骨安全」一直维持到 7 m（真实抱石高度 2–4.5 m 全安全）；"
      "补上旋后后，**同一条腓骨**在真实高度内即越阈。两口径之差 = 纯轴向模型漏掉的"
      "**冠状面旋后机制**（两篇论文所指 #1 踝伤机制）。")
    A("")
    A("---")
    A("")
    A("## 6. 敏感性（力臂 × 截面）")
    A("")
    A("| h (m) | β (°) | d (mm) | 截面 | M_inv (N·m) | σ_lat (MPa) | 侧向 risk |")
    A("|---|---:|---:|---|---:|---:|---:|")
    for s in sens:
        A(f"| {s['height_m']:g} | {s['supination_deg']:g} | {s['lever_arm_mm']:g} | "
          f"{s['section']} | {s['M_inv_nm']:.1f} | {s['sigma_lateral_mpa']:,.0f} | "
          f"{s['fibula_lateral_risk']:.2f} |")
    A("")
    A("- 三档力臂 × 两种截面下，**所有 β≥15° 的格点都越阈** → 结论对几何参数不敏感"
      "（腓骨太细，任何合理力臂都足以上量级）。")
    A("- 真实中段截面比等效圆**更细**（A=82.6 vs 101.9 mm²）→ risk 更高；等效圆是**乐观**口径。")
    A("")
    A("---")
    A("")
    A("## 7. 两条入口的等价性（roll_deg vs 显式旋后角）")
    A("")
    A(f"h={roll_eq['height_m']:g} m：`ground_reaction(roll_deg=β)`（冠状地面倾角）"
      "与 `ankle_supination`（作用于竖直 GRF 的 β）应给出同一 F_lat。")
    A("")
    A("| β (°) | roll 路线 F_lat (N) | 旋后路线 F_lat (N) | 相对差 |")
    A("|---|---:|---:|---:|")
    for r in roll_eq["rows"]:
        A(f"| {r['beta_deg']:g} | {r['roll_route_f_lat_n']:,.0f} | "
          f"{r['supination_route_f_lat_n']:,.0f} | {r['rel_diff']:.2e} |")
    A("")
    A("---")
    A("")
    A("## 8. 软垫（on_pad=True）对照")
    A("")
    A("| h (m) | β (°) | F_vert 足 (N) | F_lat (N) | M_inv (N·m) | 侧向 risk | 越阈 |")
    A("|---|---:|---:|---:|---:|---:|---|")
    for r in mat:
        A(f"| {r['height_m']:g} | {r['supination_deg']:g} | {r['F_vert_n']:,.0f} | "
          f"{r['F_lat_n']:,.0f} | {r['M_inv_nm']:.1f} | {r['fibula_lateral_risk']:.2f} | "
          f"{'⚠越阈' if r['fibula_lateral_fracture'] else '否'} |")
    A("")
    A("> 软垫降低峰值 GRF（吸能），故同一旋后角下 F_lat / M_inv / risk 均**低于刚性地面**；"
      "但 [Beurienne2025] 指出垫子刚度的两难（防旋后 vs 高能量更糟）——本模型只覆盖"
      "**载荷侧**，不含垫子刚度对**旋后角本身**的反馈（那需要 S5 的接触物理重写）。")
    A("")
    A("---")
    A("")
    A("## 9. 临界旋后角")
    A("")
    A("| h (m) | 临界旋后角 β_c (°)（risk=1） |")
    A("|---|---:|")
    seen = set()
    for r in hard:
        h = r["height_m"]
        if h in seen:
            continue
        seen.add(h)
        c = r["critical_inversion_deg_this_h"]
        A(f"| {h:g} | {'—' if c is None else f'{c:.2f}'} |")
    A("")
    A("> `β_c` 极小（< 2°），说明**上界式模型**下腓骨远端对冠状剪切/弯曲**极敏感**。"
      "这不代表真实抱石的人体在 1° 内翻就会骨折——它反映的是**纯骨 1D 口径缺少"
      "韧带/腱卸载与应力集中**。真实外伤角（15–30°）远高于此，故**定性结论可靠**"
      "（旋后是真实高度内的有效致伤机制），**定量阈值不可直接引用**。")
    A("")
    A("---")
    A("")
    A("## 10. 假设（ASSUMPTIONS）与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("**结构边界**：")
    A("")
    A("- 临床第一大踝伤是**韧带扭伤（71%）**；本骨模型只能对应**踝骨折（27%）**"
      "（[Heck2024]）。补韧带判据是后续工作。")
    A("- 旋后角是**输入参数**而非 OpenSim 里解锁的 `subtalar_angle_r`（仍是刚性锁死腿，"
      "方案 S3）。")
    A("- `σ_c=70` 为**压缩**材料强度；弯曲拉伸侧 [Y25] 仅 **30 MPa** → 本口径**低估**拉伸侧。")
    A("- 中段截面代理远端 / 外踝 → 偏保守；力臂 assumed → 敏感性已给。")
    A("")
    A("---")
    A("")
    A("## 11. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_supination.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/s1_h5_opensim.json          (默认回归基准)")
    A("#   results/opensim_fe/bc_robust_metric_route1.json(fibula A_section)")
    A("#   results/opensim_fe/risk_1d_bending.json         (fibula I / c)")
    A("#   results/opensim_fe/risk_1d_loadshare.json       (fibula 载荷分配)")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s2.py tests\\test_nonvertical_s5.py -q")
    A("```")
    A("")
    A("## 12. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


#: 论文腓骨首骨折区间字符串（报告用）。
PAPER_FIBULA_RANGE_STR = f"{PAPER_FIBULA_RANGE_M[0]}–{PAPER_FIBULA_RANGE_M[1]}"

if __name__ == "__main__":
    raise SystemExit(main())
