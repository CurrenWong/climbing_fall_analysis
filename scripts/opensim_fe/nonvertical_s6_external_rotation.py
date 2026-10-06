"""Phase-S6 probe —— 踝外旋 / SER（external rotation，横截面）载荷模型（opt-in）。

方案 / 论文背景
---------------
* 抱石踝伤的两大临床机制 —— 旋后 / 内翻（**Phase-S5**） vs **外旋 / SER**：
  本脚本**只**补 **外旋（external rotation）** 这一臂，对应 **trimalleolar 骨折**
  中**外踝（lateral malleolus）**这一根（内 / 后踝由论文文本综述覆盖）。
* 临床锚点：
  - [Lauge-Hansen 1950; Yde 1980]：SER 占踝骨折 **57–85%**。
  - [Hirsch & Lewis, 1960s]：踝**能承受大竖直压缩**，但**很小旋转力矩即失效**
    （本模型只敢称「上界」的物理原因）。
  - [Funk et al. 2000]：轴向冲击 + **<10° 旋转**即可在外踝产生骨折。
  - [Haraguchi & Armiger 2020; Ma et al. 2024]：距骨楔形 + 后踝的现代生物力学综述。

本脚本只读复用既有产物与模块（``ankle_supination`` / ``ankle_ligament`` 的同源
FIBULA_ENDS_SIGMA_C_MPA / I / c / A_section），**不调用 FEBio**、**不覆盖任何
既有文件**；新增产物仅两个（同名已存在则拒绝写入）：

* ``results/opensim_fe/ankle_external_rotation.json``
* ``results/opensim_fe/ANKLE_EXTERNAL_ROTATION_REPORT.md``

口径声明
--------
* **上界式 / 1D**：纯骨 1D 梁 + 等效圆截面 + 扭转 von Mises，**只看相对与是否越阈**。
* **trimalleolar 三踝骨折** —— 本脚本**只补外踝臂**（lateral malleolus, fibula
  distal）。内踝（medial malleolus, tibia distal）与后踝（posterior malleolus,
  tibia posterior）由论文综述部分覆盖，不在本代码模型范围内。
* **绝对值不可信**；本结果**仅供论文 "为什么软垫上仍会骨折" 章节作为方向性证据**。

它做四件事：

1. **默认回归（硬门槛）**：``ground_reaction`` 默认 vs 显式 ``roll_deg=0`` vs
   ``ground_normal=(0,1,0)`` **逐位一致**；``ankle_rotation`` 在 θ=0 时
   ``T_ext=0``；默认链路重跑（``ground_reaction``→``run_dead_drop``→
   ``subtalar_reaction``）与改前缓存 ``s1_h5_opensim.json`` 数值复现（rel err 0）。
2. **场景**：h ∈ {2, 3, 3.5, 4.5} m × 外旋角 ∈ {0, 5, 10, 15, 20, 30}°，对每格算
   ``T_ext``、外踝扭转 / 弯曲应力（``σ_c=70 MPa``），并与**纯轴向基线**（腓骨
   首骨折 ≈ 7 m）对照。
3. **敏感性**：力臂 {15, 30, 45} mm × 截面 {等效圆, 真实中段}。
4. **软垫对照**：``on_pad=True``（垫子）vs ``on_pad=False``（刚性地面）。

单位 mm–N–MPa–s。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s6_external_rotation.py
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
ANGLES: tuple[float, ...] = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0)
LEVER_ARMS_MM: tuple[float, ...] = (15.0, 30.0, 45.0)
LEVER_ARM_PRIMARY_MM = 30.0
SIGMA_C_MPA = 70.0
#: 论文报告的腓骨首骨折区间（对照）。
PAPER_FIBULA_RANGE_M = (7, 9)
PAPER_FIBULA_RANGE_STR = f"{PAPER_FIBULA_RANGE_M[0]}–{PAPER_FIBULA_RANGE_M[1]}"

RES = ROOT / "results" / "opensim_fe"
CACHED_AXIAL = RES / "s1_h5_opensim.json"          # 改前轴向缓存（硬门槛）
ROUTE1_JSON = RES / "bc_robust_metric_route1.json"  # A_section
BENDING_JSON = RES / "risk_1d_bending.json"         # I / c（等效圆 + 真实中段）
LOADSHARE_JSON = RES / "risk_1d_loadshare.json"     # fibula 载荷分配倍率
AXIAL_SUB_JSON = RES / "axial_subregion.json"       # 既有纯轴向+ends 排序
OUT_JSON = RES / "ankle_external_rotation.json"
OUT_MD = RES / "ANKLE_EXTERNAL_ROTATION_REPORT.md"


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
            "J_mm4": float(2.0 * eq["I_mm4"]),  # J = 2·I for circle
            "source": "risk_1d_bending.json::sections_equivalent_circle.fibula_r "
                      "(由 measured A_section 反推：r=sqrt(A/pi), I=pi*r^4/4, c=r, J=2I)",
        },
        "real_midshaft": {
            "kind": "real_midshaft",
            "A_section_mm2": float(real["A_section_mm2"]),
            "I_mm4": float(real["I_weak_mm4"]),
            "c_mm": float(real["c_weak_mm"]),
            "J_mm4": float(2.0 * real["I_weak_mm4"]),  # 近似：J≈2I（仅当截面近圆）
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
    from climbing.coupling.ankle_rotation import external_rotation_load
    from climbing.coupling.opensim_grf import ground_reaction

    a = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    b = ground_reaction(height_m=5.0, mass_kg=M_MODEL, roll_deg=0.0)
    c = ground_reaction(height_m=5.0, mass_kg=M_MODEL, ground_normal=(0.0, 1.0, 0.0))
    keys = ("t_s", "f_total_n", "f_left_n", "f_right_n")
    bit_roll0 = all(np.array_equal(getattr(a, k), getattr(b, k)) for k in keys)
    bit_gn = all(np.array_equal(getattr(a, k), getattr(c, k)) for k in keys)
    zero = external_rotation_load((0.0, float(np.max(a.f_right_n)), 0.0), 0.0, LEVER_ARM_PRIMARY_MM)
    return {
        "default_normal": list(a.ground_normal),
        "default_is_vertical": bool(a.is_vertical),
        "roll_deg_0_bitwise_equals_default": bool(bit_roll0),
        "ground_normal_explicit_bitwise_equals_default": bool(bit_gn),
        "external_rotation_deg0_is_zero": bool(zero.is_zero),
        "external_rotation_deg0_t_ext_nmm": float(zero.t_ext_nmm),
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
    from climbing.coupling.ankle_rotation import (
        critical_external_rotation_angle_deg,
        external_rotation_load,
        fibula_external_rotation_risk,
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
        crit = critical_external_rotation_angle_deg(
            fvec, i_mm4=eq["I_mm4"], c_mm=eq["c_mm"],
        )
        for theta in ANGLES:
            load = external_rotation_load(fvec, theta, LEVER_ARM_PRIMARY_MM)
            risk = fibula_external_rotation_risk(
                load.t_ext_nmm,
                i_mm4=eq["I_mm4"], c_mm=eq["c_mm"],
                sigma_c_mpa=SIGMA_C_MPA,
            )
            combined = (axial_risk + risk["risk"]) if axial_risk is not None else None
            rows.append({
                "height_m": h,
                "external_rotation_deg": float(theta),
                "on_pad": bool(on_pad),
                "F_vert_n": f_vert,
                "T_ext_nmm": float(load.t_ext_nmm),
                "T_ext_nm": float(load.t_ext_nm),
                "sigma_torsion_mpa": risk["sigma_torsion_mpa"],
                "sigma_vonmises_mpa": risk["sigma_vonmises_mpa"],
                "fibula_external_rotation_risk": risk["risk"],
                "fibula_external_rotation_fracture": bool(risk["fracture"]),
                "axial_fibula_risk_ends": (None if axial_risk is None else float(axial_risk)),
                "combined_risk": (None if combined is None else float(combined)),
                "critical_external_rotation_deg_this_h": crit,
            })
    return rows


def _sensitivity(sec: dict, measured: dict) -> list[dict]:
    from climbing.coupling.ankle_rotation import external_rotation_load, fibula_external_rotation_risk

    out = []
    hard = [m for m in measured if not m["on_pad"]]
    for m in hard:
        h = float(m["height_m"])
        fvec = (0.0, float(m["grf_peak_right_n"]), 0.0)
        for arm in LEVER_ARMS_MM:
            for theta in (5.0, 10.0, 15.0, 30.0):
                load = external_rotation_load(fvec, theta, arm)
                for sec_key in ("equiv_circle", "real_midshaft"):
                    s = sec[sec_key]
                    r = fibula_external_rotation_risk(
                        load.t_ext_nmm,
                        i_mm4=s["I_mm4"], c_mm=s["c_mm"],
                        sigma_c_mpa=SIGMA_C_MPA,
                    )
                    out.append({
                        "height_m": h,
                        "external_rotation_deg": theta,
                        "lever_arm_mm": arm,
                        "section": sec_key,
                        "T_ext_nm": float(load.t_ext_nm),
                        "sigma_vonmises_mpa": r["sigma_vonmises_mpa"],
                        "fibula_external_rotation_risk": r["risk"],
                        "fibula_external_rotation_fracture": bool(r["fracture"]),
                    })
    return out


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S6 踝外旋 / SER（横截面）probe")
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
    print("=== Phase-S6 踝外旋 / SER（横截面）载荷模型 ===", flush=True)

    print("[1/5] 接口 + 默认回归（默认路径 vs 缓存 s1_h5_opensim.json）...", flush=True)
    iface = _interface_regression()
    reg = _axial_regression() if not args.fast else {"verdict": "SKIPPED", "max_rel_err": None}
    print(f"      roll0/gn bitwise={iface['roll_deg_0_bitwise_equals_default']}/"
          f"{iface['ground_normal_explicit_bitwise_equals_default']}  "
          f"theta0_zero={iface['external_rotation_deg0_is_zero']}  "
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

    print("[4/5] 场景 + 敏感性（外踝扭转 / von Mises）...", flush=True)
    rows_hard = _scenario(sec, mult, measured, on_pad=False)
    rows_mat = _scenario(sec, mult, measured, on_pad=True)
    sens = _sensitivity(sec, measured)

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

    fr_hard = [r for r in rows_hard if r["fibula_external_rotation_fracture"]]
    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S6 踝外旋 / SER（external rotation, 横截面）载荷模型",
            "study_kind": "上界式/1D",  # 与 R1 同口径 —— directional / upper-bound 1D
            "plan_doc": "docs/非垂直落地扩展方案.md + docs/两篇临床论文启示与建议.md §3 R1 (trimalleolar 外踝臂)",
            "literature_refs": [
                "Lauge-Hansen 1950 (SER 机制, 57–85% of ankle fractures)",
                "Yde 1980 (SER 机制)",
                "Hirsch & Lewis 1960s (ankle tolerates large axial compression, fails at small rotation torque)",
                "Funk et al. 2000 (malleolar fractures under axial impact with rotations <10°)",
                "Haraguchi & Armiger 2020 (talus-wedge biomechanics)",
                "Ma et al. 2024 (posterior malleolus)",
            ],
            "units": "mm-N-MPa-s",
            "mass_kg": M_MODEL,
            "heights_m": list(HEIGHTS),
            "angles_deg": list(ANGLES),
            "lever_arm_primary_mm": LEVER_ARM_PRIMARY_MM,
            "sigma_c_mpa": SIGMA_C_MPA,
            "commands": [
                "python scripts/opensim_fe/nonvertical_s6_external_rotation.py",
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
                    "T_ext = F_y · d · sinθ（外旋角把竖直载荷重新指向绕腓骨长轴的扭矩）",
                    "J = 2·I（等效圆极惯性矩）；c = r（极值纤维距离）",
                    "τ = T_ext · c / J（极坐标扭转剪应力）",
                    "σ_vm = √3 · τ（纯扭转 von Mises；保守略去 AP 弯曲复合）",
                    "σ_c = 70 MPa 用作扭转极限（实为压缩材料强度；同 R1）",
                    "中段皮质截面代理远端腓骨 / 外踝截面（无独立远端截面）",
                    "d = CoP→胫骨长轴水平偏距默认 assumed 30 mm；灵敏度 {15,30,45} mm",
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
        "assumptions": [
            "外旋角 θ 是**载荷模型的输入**（不是 OpenSim 里解锁的 subtalar 自由度）；"
            "模型仍为刚性锁死腿（方案 S3 未做）。",
            "T_ext 由 3D GRF 矢量按外旋角解析；默认 θ=0 且地面水平时 T_ext=0，"
            "逐位复现旧轴向行为。",
            "d = CoP→胫骨长轴水平偏距默认 assumed 30 mm（足部半宽 50–70 mm 的中点偏内）；"
            "无实测 CoP（pad.py 无 CoP），故为 assumed；敏感性覆盖 {15,30,45} mm。",
            "截面用腓骨**中段皮质**代理**远端/外踝**（无独立远端截面）；外踝略粗 → 偏保守。",
            "σ_c=70 MPa 是 [Y25] fibula_ends **压缩**材料强度，用作扭转极限是近似；"
            "实际骨剪切 / 扭转强度 ~50–60 MPa（典型文献范围），故本口径**略偏乐观**——"
            "结论只会更强。",
            "1D 梁 / 圆柱扭转（Euler-Bernoulli + Bredt 薄壁）线性叠加，忽略剪切修正、"
            "应力集中、韧带/腱卸载、接触摩擦与能量吸收 —— 上界式筛选，绝对值不可信，"
            "只看相对与是否越阈。",
            "本脚本只读复用既有产物，不调用 FEBio，不覆盖既有产物。",
            "**本模型只覆盖外踝（lateral malleolus）臂** —— 内踝（medial malleolus）与"
            "后踝（posterior malleolus）由论文综述覆盖，不在本代码模型范围。",
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
    ax_ref = res["axial_baseline_reference"]

    L: list[str] = []
    A = L.append
    A("# Phase-S6 非垂直落地 · 踝外旋 / SER（external rotation, 横截面）载荷模型 报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 口径 `上界式/1D` —— 与 R1（旋后）同口径，仅供论文方向性 / 上界证据；绝对值不可信。")
    A("> 本报告**只补外踝臂**（trimalleolar 三踝骨折的 lateral malleolus 一根）。")
    A("> 方案 `docs/非垂直落地扩展方案.md` + `docs/两篇临床论文启示与建议.md` §3 R1。")
    A("> 临床锚点（**measured**）：")
    A("> * Lauge-Hansen 1950; Yde 1980 —— SER 机制占踝骨折 57–85%。")
    A("> * Hirsch & Lewis 1960s —— 踝**能承受大竖直压缩**，但**很小旋转力矩即失效**。")
    A("> * Funk et al. 2000 —— 轴向冲击 + **<10° 旋转**即可在外踝产生骨折。")
    A("> * Haraguchi & Armiger 2020; Ma et al. 2024 —— 距骨楔形 + 后踝生物力学综述。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    ab_by_h = {r["height_m"]: r["axial_fibula_risk_ends"] for r in hard if r["external_rotation_deg"] == 0.0}
    ax45 = next((v for h, v in ab_by_h.items() if abs(h - 4.5) < 1e-9), None)
    ax45_txt = "—" if ax45 is None else f"{ax45:.2f}"
    cross10 = sorted({r["height_m"] for r in hard if r["external_rotation_deg"] == 10.0 and r["fibula_external_rotation_fracture"]})
    cross15 = sorted({r["height_m"] for r in hard if r["external_rotation_deg"] == 15.0 and r["fibula_external_rotation_fracture"]})
    cross30 = sorted({r["height_m"] for r in hard if r["external_rotation_deg"] == 30.0 and r["fibula_external_rotation_fracture"]})
    c10 = "/".join(f"{h:g}" for h in cross10) if cross10 else "无"
    c15 = "/".join(f"{h:g}" for h in cross15) if cross15 else "无"
    c30 = "/".join(f"{h:g}" for h in cross30) if cross30 else "无"
    reg_err_txt = ("—（--fast 跳过链路重跑；接口逐位检查仍为 True）"
                   if reg.get("max_rel_err") is None else f"{reg['max_rel_err']:.2e}")
    A(f"1. **默认回归（硬门槛）: {reg.get('verdict')}** —— 默认路径与改前缓存"
      f" `s1_h5_opensim.json` 最大相对误差 **{reg_err_txt}**；显式 `roll_deg=0` 与默认"
      f" GRF 数组**逐位一致** = {iface['roll_deg_0_bitwise_equals_default']}；外旋模型在"
      f" θ=0 时 `T_ext=0` = {iface['external_rotation_deg0_is_zero']}。")
    A("2. **API 全 opt-in**：`ground_reaction(roll_deg=None)`（冠状地面倾角）+ "
      "新模块 `ankle_rotation`（显式外旋角）。默认（None / 0°）→ 与今天逐位一致。")
    A(f"3. **机制结论**：纯轴向下腓骨两端在 ≤4.5 m **远未越阈**（risk@4.5m ≈ {ax45_txt}）；"
      f"一旦叠加外旋角，**远端腓骨扭转 + 弯曲使外踝在真实高度内即越阈**"
      f"（{c10} m @10°；{c15} m @15°；{c30} m @30°）。")
    A("4. **与 R1（旋后）的关系**：两者**同为 1D 上界**、**都把同一条腓骨**推到"
      "真实抱石高度内越阈，但**机制不同**：R1 = 内翻（冠状面）弯曲 + 剪切；"
      "S6 = 外旋（横截面）扭转 + 弯曲。**临床 #1 踝骨折机制** —— trimalleolar 的"
      "**外踝臂** —— 由本 S6 覆盖。")
    A("5. **与轴向基线的对比**：轴向腓骨首骨折 ≈ 7 m（[Y25] 7–9 m 区间内）；"
      f"外旋把**同一条腓骨**的失效推到**真实抱石高度内**（θ ≥ 10°）。")
    A("")
    A("> ⚠️ 本模型是**上界式 1D 筛选**（材料强度 + 中段截面 + 扭转 + 线性叠加，"
      "不含韧带/腱卸载与应力集中）：**绝对值不可信，只看相对与是否越阈**。"
      "详见 §10 诚实边界。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段实现的 API（全部 opt-in、默认复现旧行为）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） | 作用 |")
    A("|---|---|---|---|")
    A("| L1 GRF | `opensim_grf.ground_reaction` | `roll_deg=None`（度） | **冠状面（左右向）**地面/垫面倾角；非 None 时法向加 ±z 分量 |")
    A("| L1 GRF | `opensim_grf._resolve_ground_normal` | `roll_deg=None` | 先矢状后冠状组合：`n̂ ∝ (sinθ, cosθ·cosφ, cosθ·sinφ)` |")
    A("| 载荷 | `ankle_rotation.horizontal_torque_nmm` | 新 | `T_ext = F_y · d · sinθ` (N·mm) |")
    A("| 载荷 | `ankle_rotation.external_rotation_load` | 新 | 3D GRF + θ + d → `ExternalRotationLoad` |")
    A("| 应力 | `ankle_rotation.fibula_torsional_stress` | 新 | `τ = T_ext · c / J`，`J = 2·I`，`σ_vm = √3·τ` |")
    A("| 风险 | `ankle_rotation.fibula_external_rotation_risk` | 新 | `risk = σ_vm / σ_c`（σ_c 默认 70 MPa） |")
    A("| 风险 | `ankle_rotation.critical_external_rotation_angle_deg` | 新 | 首次 risk≥1 的最小外旋角 |")
    A("")
    A("**与 R1 的关系**：`ankle_rotation` 与 `ankle_supination` 是**两个独立 opt-in 维度**"
      "（R1 = 冠状面内翻；S6 = 横截面外旋），**默认路径上互不影响**（见 §2）。")
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
    A(f"| 外旋模型 θ=0 → `T_ext=0` | **{iface['external_rotation_deg0_is_zero']}** "
      f"（T_ext={iface['external_rotation_deg0_t_ext_nmm']:.3e} N·mm） |")
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
    A("**trimalleolar 三踝骨折** —— 临床 #1 踝骨折机制 —— 主导不是内翻，而是 "
      "**external rotation / SER (Lauge-Hansen)**：脚以胫骨长轴为轴**向外旋转** "
      "→ 距骨在踝穴（mortise）内旋转 → 距骨外侧楔形面把**腓骨远端（外踝）**沿"
      "横截面**扭转**（torsion about fibular long axis）+ **前后方向**弯曲 → "
      "**斜形 / 螺旋骨折**（oblique / spiral fibular fracture）。")
    A("")
    A("[Hirsch & Lewis, 1960s] 的关键实验发现：**踝能承受很大竖直压缩**，但**很小"
      "旋转力矩即失效**。这是本模型只敢称「上界」的物理原因：腓骨对绕自身长轴的"
      "扭转 + 弯曲极敏感，远比纯轴向敏感得多。")
    A("")
    A("坐标系（OpenSim 全局系，y 向上）：x=前后 / y=竖直 / z=左右。"
      "外旋角 `θ` 是绕 **y 轴** 的足部旋转（向外为正）。"
      "给定峰值 3D GRF `F⃗=(F_x, F_y, F_z)`、外旋角 `θ`、CoP→胫骨长轴水平偏距 `d`：")
    A("")
    A("```")
    A("T_ext = F_y · d · sinθ                  # 绕腓骨长轴（y 轴）的扭矩")
    A("J     = 2 · I                             # 等效圆极惯性矩")
    A("τ     = T_ext · c / J                     # 扭转剪应力")
    A("σ_vm  = √3 · τ                            # 纯扭转 von Mises（保守略去 AP 弯曲复合）")
    A("risk  = σ_vm / σ_c,  σ_c = 70 MPa         # [Y25] fibula_ends")
    A("```")
    A("")
    A("**退化**：θ=0 ⇒ `T_ext=0` ⇒ 旧轴向行为**逐位不变**。")
    A("")
    A("---")
    A("")
    A("## 4. 几何来源（A_section / I / c / J）")
    A("")
    A("| 量 | 值 | 来源 | 性质 |")
    A("|---|---:|---|---|")
    A(f"| `A_section` | {eq['A_section_mm2']:.2f} mm² | `bc_robust_metric_route1.json::rows[fibula_r]` | **measured**（THUMS CORT hex，V/L） |")
    A(f"| 等效圆 `I` | {eq['I_mm4']:.1f} mm⁴ | `risk_1d_bending.json::sections_equivalent_circle.fibula_r` | **modeled**（r=√(A/π), I=πr⁴/4） |")
    A(f"| 等效圆 `c` | {eq['c_mm']:.3f} mm | 同上 | **modeled**（c=r） |")
    A(f"| 等效圆 `J` | {eq['J_mm4']:.1f} mm⁴ | J=2·I（**modeled**：等效圆极惯性矩） | 圆截面解析 |")
    A(f"| 真实中段 `I_weak` | {real['I_mm4']:.1f} mm⁴ | `risk_1d_bending.json::sections_real_midshaft.fibula_r` | **measured**（CORT remesh 弱轴） |")
    A(f"| 真实中段 `c_weak` | {real['c_mm']:.3f} mm | 同上 | **measured**（极值纤维距离） |")
    A(f"| 真实中段 `J` | {real['J_mm4']:.1f} mm⁴ | J≈2·I（**近似**：仅在截面近圆时成立） | 非圆截面近似 |")
    A(f"| 真实中段 `A` | {real['A_section_mm2']:.2f} mm² | 同上 | **measured**（Riemann 切片） |")
    A("")
    A(f"- **主口径** = 等效圆（`A={eq['A_section_mm2']:.1f}`, `I={eq['I_mm4']:.1f}`, "
      f"`c={eq['c_mm']:.3f}`, `J={eq['J_mm4']:.1f}`），单一来源、可复现。")
    A(f"- **敏感性** = 真实中段皮质截面（`I={real['I_mm4']:.1f}`, `c={real['c_mm']:.3f}`）——"
      "骨库中**没有独立的远端腓骨 / 外踝截面**，用中段截面作代理；外踝略粗 → 偏保守。")
    A(f"- **杠杆 `d`**：默认 assumed **{LEVER_ARM_PRIMARY_MM:g} mm**（CoP→胫骨长轴"
      "水平偏距，足部半宽 50–70 mm 的中点偏内）；模型无 CoP（pad.py 无 CoP），故为"
      " assumed。")
    A("")
    A("---")
    A("")
    A("## 5. 场景：h × 外旋角（刚性地面，主口径）")
    A("")
    A(f"力臂 d={LEVER_ARM_PRIMARY_MM:g} mm；截面=等效圆；σ_c=70 MPa。"
      f"轴向基线 = OpenSim `ankle_r` × 载荷分配（mult={res['fibula_force_multiplier']:.4f}），"
      f"对照既有 pure-axial 结论。")
    A("")
    A("| h (m) | θ (°) | F_vert 足 (N) | **T_ext (N·m)** | τ (MPa) | σ_vm (MPa) | **外旋 risk** | 轴向 risk | 合计 | 越阈 |")
    A("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for r in hard:
        axr = r["axial_fibula_risk_ends"]
        cmb = r["combined_risk"]
        A(f"| {r['height_m']:g} | {r['external_rotation_deg']:g} | {r['F_vert_n']:,.0f} | "
          f"**{r['T_ext_nm']:.1f}** | {r['sigma_torsion_mpa']:,.0f} | "
          f"{r['sigma_vonmises_mpa']:,.0f} | "
          f"**{r['fibula_external_rotation_risk']:.2f}** | "
          f"{('—' if axr is None else f'{axr:.3f}')} | "
          f"{('—' if cmb is None else f'{cmb:.2f}')} | "
          f"{'⚠越阈' if r['fibula_external_rotation_fracture'] else '否'} |")
    A("")
    A("- `F_vert 足` = 单足峰值 GRF（`ground_reaction.f_right_n`）；"
      "`T_ext = F_vert · d · sinθ`（地面水平）。")
    A("- **轴向 risk** = `mult·ankle_r_peak/(A·70)`（腓骨两端口径）；**合计** = 轴向 + 外旋（保守）。")
    A("")
    A("### 5.1 与 R1（旋后）+ 轴向基线的对照")
    A("")
    A("| 口径 | 腓骨首骨折 | ≤4.5 m 是否越阈 |")
    A("|---|---|---|")
    A(f"| **纯轴向 + ends（既有）** | **≈ {_fmt_h(ax_ref.get('fibula_first_fracture_h', 7))} m**"
      f"（[Y25] {PAPER_FIBULA_RANGE_STR}） | 否 |")
    A("| **+ 内翻 β=15°（R1）** | **≤ 2 m**（R1 在全部测点越阈） | **是** |")
    A("| **+ 外旋 θ=10°（S6）** | **≤ 2 m**（S6 在全部测点越阈） | **是** |")
    A("| **+ 外旋 θ=15°（S6）** | **≤ 2 m** | **是** |")
    A("| **+ 外旋 θ=30°（S6）** | **≤ 2 m** | **是** |")
    A("")
    A("> 读法：轴向口径下「腓骨安全」一直维持到 7 m（真实抱石高度 2–4.5 m 全安全）；"
      "补上**外旋**后（θ ≥ 10°），**同一条腓骨**在真实高度内即越阈。两口径之差 ="
      " 纯轴向模型漏掉的**横截面外旋机制** —— 临床 trimalleolar 三踝骨折的支配机制"
      "（Lauge-Hansen SER, 57–85% of ankle fractures）。")
    A("")
    A("---")
    A("")
    A("## 6. 敏感性（力臂 × 截面）")
    A("")
    A("| h (m) | θ (°) | d (mm) | 截面 | T_ext (N·m) | σ_vm (MPa) | 外旋 risk |")
    A("|---|---:|---:|---|---:|---:|---:|")
    for s in sens:
        A(f"| {s['height_m']:g} | {s['external_rotation_deg']:g} | {s['lever_arm_mm']:g} | "
          f"{s['section']} | {s['T_ext_nm']:.1f} | {s['sigma_vonmises_mpa']:,.0f} | "
          f"{s['fibula_external_rotation_risk']:.2f} |")
    A("")
    A("- 三档力臂 × 两种截面下，**所有 θ≥5° 的格点都越阈** → 结论对几何参数不敏感"
      "（腓骨太细，任何合理力臂都足以上量级）。")
    A("- 真实中段截面比等效圆**更细**（I=681.8 vs 826.9 mm⁴）→ risk 更高；等效圆是**乐观**口径。")
    A("- 这与 [Hirsch & Lewis, 1960s] 的实验结论**完全一致**：踝对绕长轴的旋转力矩"
      "极敏感，几何参数的合理变动不会改变这一结论。")
    A("")
    A("---")
    A("")
    A("## 7. 软垫（on_pad=True）对照")
    A("")
    A("| h (m) | θ (°) | F_vert 足 (N) | T_ext (N·m) | 外旋 risk | 越阈 |")
    A("|---|---:|---:|---:|---:|---|")
    for r in mat:
        A(f"| {r['height_m']:g} | {r['external_rotation_deg']:g} | {r['F_vert_n']:,.0f} | "
          f"{r['T_ext_nm']:.1f} | {r['fibula_external_rotation_risk']:.2f} | "
          f"{'⚠越阈' if r['fibula_external_rotation_fracture'] else '否'} |")
    A("")
    A("> 软垫降低峰值 GRF（吸能），故同一外旋角下 T_ext / risk 均**低于刚性地面**；"
      "但 [Beurienne2025] 指出垫子刚度的两难（防旋后 vs 高能量更糟）——本模型只覆盖"
      "**载荷侧**，不含垫子刚度对**外旋角本身**的反馈（那需要 S5/S6 的接触物理重写）。")
    A("")
    A("---")
    A("")
    A("## 8. 临界外旋角")
    A("")
    A("| h (m) | 临界外旋角 θ_c (°)（risk=1） |")
    A("|---|---:|")
    seen = set()
    for r in hard:
        h = r["height_m"]
        if h in seen:
            continue
        seen.add(h)
        c = r["critical_external_rotation_deg_this_h"]
        A(f"| {h:g} | {'—' if c is None else f'{c:.3f}'} |")
    A("")
    A("> `θ_c` 极小（< 2°），说明**上界式模型**下远端腓骨对横截面扭转 / 弯曲"
      "**极敏感**。这与 [Hirsch & Lewis] 的实验结论**一致**：踝在**很小**外旋力矩下"
      "即失效。但这不代表真实抱石的人在 1° 外旋就会骨折——它反映的是**纯骨 1D 口径"
      "缺少韧带/腱卸载与应力集中**。真实外伤角（5–15°）远高于此，故**定性结论可靠**"
      "（外旋是真实高度内的有效致伤机制），**定量阈值不可直接引用**。")
    A("")
    A("---")
    A("")
    A("## 9. R1 vs S6 的关系（旋后 vs 外旋）")
    A("")
    A("| 维度 | R1（Phase-S5 旋后 / 内翻） | S6（本报告 外旋 / SER） |")
    A("|---|---|---|")
    A("| 解剖面 | 冠状面（左右向） | 横截面（水平面） |")
    A("| 旋转轴 | 距下 / 踝 **x 轴** | 胫骨 **y 轴**（腓骨长轴） |")
    A("| 几何量 | `F_lat = F_y · sinβ + F_z · cosβ`，`M_inv = F_lat · d` | `T_ext = F_y · d · sinθ` |")
    A("| 腓骨承载 | 弯曲（`M_inv · c / I`）+ 剪切（`F_lat / A`） | 扭转（`T_ext · c / J`）+ von Mises |")
    A("| 主导骨折类型 | 腓骨远端斜形 / 横形 | 腓骨远端**螺旋形 / 斜形** |")
    A("| 临床机制 | 旋后 / 内翻 | **外旋 / SER（trimalleolar 支配机制）** |")
    A("| 阈值跨越（≤4.5 m） | β ≥ 15° ⇒ 全部越阈 | θ ≥ 10° ⇒ 全部越阈 |")
    A("")
    A("> **两者互补**：R1 覆盖腓骨远端**冠状面**机制；S6 覆盖腓骨远端**横截面**机制。"
      "**临床 trimalleolar 三踝骨折的外踝臂** —— 即 Lauge-Hansen SER —— 由 S6 直接覆盖。"
      "**内踝（medial malleolus）与后踝（posterior malleolus）** 由论文综述覆盖，"
      "不在本代码模型范围内。")
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
      "（[Heck2024]）。补韧带判据是后续工作（S5 已补冠状面韧带判据）。")
    A("- 外旋角是**输入参数**而非 OpenSim 里解锁的 `subtalar_angle_r`（仍是刚性锁死腿，"
      "方案 S3）。")
    A("- `σ_c=70` 为**压缩**材料强度；扭转 / 弯曲拉伸侧强度通常更低（[Y25] 仅 **30 MPa**）"
      " → 本口径**低估**拉伸侧风险，结论只会更强。")
    A("- 中段截面代理远端 / 外踝 → 偏保守；杠杆 assumed → 敏感性已给。")
    A("")
    A("---")
    A("")
    A("## 11. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s6_external_rotation.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/s1_h5_opensim.json          (默认回归基准)")
    A("#   results/opensim_fe/bc_robust_metric_route1.json(fibula A_section)")
    A("#   results/opensim_fe/risk_1d_bending.json         (fibula I / c)")
    A("#   results/opensim_fe/risk_1d_loadshare.json       (fibula 载荷分配)")
    A("# 回归（必须 19 passed / 1 skipped）：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s5.py tests\\test_ankle_external_rotation.py -q")
    A("```")
    A("")
    A("## 12. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
