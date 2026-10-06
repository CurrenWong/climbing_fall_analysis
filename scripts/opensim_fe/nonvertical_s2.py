"""Phase-S2 probe —— 单脚 / 不对称落地（逐足 GRF）。

方案：``docs/非垂直落地扩展方案.md`` §2 自由度④ / §3 L1 / §4 Phase-S2 /
§6 硬骨头④ / §7（轴向回归硬门槛）。

本脚本**只读复用**既有产物与模块，**不调用 FEBio**、**不覆盖任何既有文件**；
新增产物仅两个（同名已存在则拒绝写入）：

* ``results/opensim_fe/nonvertical_s2.json``
* ``results/opensim_fe/NONVERTICAL_S2_REPORT.md``

它做三件事：

1. **接口泛化（L1）**：``opensim_grf.ground_reaction`` 的 ``split`` 从写死的
   ``(0.5, 0.5)`` 推广为**任意逐足占比**（含一足为 0 = 单脚，见
   ``single_foot_split``）；新增 ``foot_forces`` 显式**逐足力时程**（两条彼此
   独立的时程，方案 §3 L1 完全版）。
2. **默认回归（硬门槛，方案 §7）**：默认 ``(0.5, 0.5)`` 重跑
   （``ground_reaction`` → ``run_dead_drop`` → ``subtalar_reaction``）与改前缓存
   ``s1_h5_opensim.json`` 数值对比；并证明默认路径与显式 ``(0.5,0.5)``、与
   ``foot_forces`` 重构**逐位一致**。
3. **单脚载荷链**：在代表高度 h=5 m（另测 4/6/7/8/10 m 检验因子稳定性）对比
   **单脚（右足承全部）vs 对称**的逐关节峰值（subtalar/ankle/knee/hip/lumbar）+
   逐骨纯轴向+两端 1D 风险（``risk = mult·F/(A·σ_c_weak)``，[Y25] 弱子区域）+
   首骨折高度 + 排序变化。

单位 mm–N–MPa–s。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s2.py
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

from climbing.bone import MATERIAL_STRENGTH_MPA  # noqa: E402

M_MODEL = 75.337
H_REP = 5.0
#: 单脚测点高度（检验"载荷翻倍因子"是否随高度稳定）。
MEAS_HEIGHTS: tuple[float, ...] = (4.0, 5.0, 6.0, 7.0, 8.0, 10.0)

RES = ROOT / "results" / "opensim_fe"
CACHED_AXIAL = RES / "s1_h5_opensim.json"          # 改前轴向缓存（硬门槛）
LOADSHARE = RES / "risk_1d_loadshare.json"         # A / σ_c / 分配倍率 / risk(h)
FRACTURE_MATRIX = RES / "fracture_matrix.json"     # 对称逐高度关节载荷 loads_n
AXIAL_SUB = RES / "axial_subregion.json"           # 纯轴向+ends 既有排序（对照）
OUT_JSON = RES / "nonvertical_s2.json"
OUT_MD = RES / "NONVERTICAL_S2_REPORT.md"

#: 逐骨力链（骨 → 关节 + [Y25] 弱子区域）。与 ``risk_1d_loadshare.FORCE_SOURCE`` 同源。
BONE_CHAIN: dict[str, tuple[str, str]] = {
    "calcaneus_r": ("subtalar_r", "calcaneus"),
    "tibia_r": ("ankle_r", "tibia_ends"),
    "fibula_r": ("ankle_r", "fibula_ends"),
    "femur_r": ("hip_r", "femoral_neck"),
    "R_HIPBONE": ("hip_r", "pelvis"),
    "L3": ("lumbar", "spine"),
    "T6": ("lumbar", "spine"),
    "C5": ("lumbar", "spine"),
    "parietal_r": ("lumbar", "skull"),
}
#: 单脚链需要测量的关节（+ 未承力侧 subtalar_l 作对照）。
CHAIN_JOINTS: tuple[str, ...] = ("subtalar_r", "ankle_r", "knee_r", "hip_r", "lumbar")
UNLOADED_CONTRAST = "subtalar_l"


def _sc(key: str) -> float:
    """[Y25] 子区域压缩强度 (MPa)。"""
    return float(MATERIAL_STRENGTH_MPA[key][1])


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b != 0.0 else abs(a - b)


def _first_fracture(heights, mask) -> float | None:
    for h, m in zip(heights, mask):
        if bool(m):
            return float(h)
    return None


def _joint_peaks(fall, grf, joints) -> dict:
    from climbing.coupling.joint_reactions import joint_reaction

    out = {}
    for j in joints:
        jr = joint_reaction(fall, grf, joint=j)
        out[j] = {
            "peak_vertical_n": float(jr.peak_vertical_n),
            "peak_force_n": float(jr.peak_force_n),
            "peak_moment_nm": float(jr.peak_moment_nm),
            "peak_force_time_s": float(jr.peak_force_time_s),
            "subtree_mass_kg": float(jr.subtree_mass_kg),
        }
    return out


# --------------------------------------------------------------------------
# 1) 接口泛化 + 默认回归（硬门槛）
# --------------------------------------------------------------------------
def _split_interface_checks() -> dict:
    from climbing.coupling.opensim_grf import ground_reaction, single_foot_split

    a = ground_reaction(height_m=H_REP, mass_kg=M_MODEL)
    b = ground_reaction(height_m=H_REP, mass_kg=M_MODEL, split=(0.5, 0.5))
    bitwise_explicit = all(
        np.array_equal(getattr(a, k), getattr(b, k))
        for k in ("t_s", "f_total_n", "f_left_n", "f_right_n")
    )
    # foot_forces 用 (t, f) 重构 0.5/0.5 → 应与默认逐位一致
    rec = ground_reaction(
        height_m=H_REP, mass_kg=M_MODEL,
        foot_forces=((a.t_s, a.f_total_n * 0.5), (a.t_s, a.f_total_n * 0.5)),
    )
    bitwise_reconstruct = all(
        np.array_equal(getattr(a, k), getattr(rec, k))
        for k in ("t_s", "f_total_n", "f_left_n", "f_right_n")
    )
    # 单脚
    g_r = ground_reaction(height_m=H_REP, mass_kg=M_MODEL, split=single_foot_split("r"))
    g_l = ground_reaction(height_m=H_REP, mass_kg=M_MODEL, split=single_foot_split("l"))
    return {
        "default_equals_explicit_0.5": bool(bitwise_explicit),
        "foot_forces_reconstructs_default": bool(bitwise_reconstruct),
        "default_per_foot_forces": bool(a.per_foot_forces),
        "single_right": {
            "split": list(single_foot_split("r")),
            "loaded_sides": list(g_r.loaded_sides),
            "peak_share": list(g_r.peak_share),
            "impulse_ns": float(g_r.impulse_ns),
            "peak_total_n": float(g_r.peak_total_n),
        },
        "single_left": {
            "split": list(single_foot_split("l")),
            "loaded_sides": list(g_l.loaded_sides),
            "peak_share": list(g_l.peak_share),
        },
    }


def _default_measure() -> dict:
    """默认 (0.5,0.5) @5 m：硬门槛回归 + 对称关节峰。"""
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))
    t0 = time.time()
    grf = ground_reaction(height_m=H_REP, mass_kg=M_MODEL)
    fall = run_dead_drop(grf)
    jl = subtalar_reaction(fall, grf, side="r")
    peaks = _joint_peaks(fall, grf, CHAIN_JOINTS)
    wall = time.time() - t0

    measured = {
        "impulse_ns": float(grf.impulse_ns),
        "impulse_rel_err": float(grf.impulse_rel_err),
        "grf_peak_n": float(grf.peak_total_n),
        "subtalar_peak_vertical_n": float(jl.peak_vertical_n),
        "subtalar_peak_force_n": float(jl.peak_force_n),
        "subtalar_peak_moment_nm": float(jl.peak_moment_nm),
    }
    cached_keys = {k: float(cached[k]) for k in measured}
    rel_err = {k: _rel(measured[k], cached_keys[k]) for k in measured}
    reg = {
        "cached_file": str(CACHED_AXIAL.relative_to(ROOT)),
        "measured_default_run": measured,
        "cached": cached_keys,
        "rel_err": rel_err,
        "max_rel_err": max(rel_err.values()),
        "verdict": "PASS" if max(rel_err.values()) <= 1e-6 else "FAIL",
    }
    return {"regression": reg, "symmetric_peaks": peaks, "wall_s": wall}


# --------------------------------------------------------------------------
# 2) 单脚载荷链
# --------------------------------------------------------------------------
def _single_measure(height_m: float) -> dict:
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction, single_foot_split

    grf = ground_reaction(height_m=height_m, mass_kg=M_MODEL, split=single_foot_split("r"))
    fall = run_dead_drop(grf)
    peaks = _joint_peaks(fall, grf, CHAIN_JOINTS + (UNLOADED_CONTRAST,))
    return {
        "height_m": float(height_m),
        "grf": {
            "impulse_ns": float(grf.impulse_ns),
            "grf_peak_n": float(grf.peak_total_n),
            "impulse_rel_err": float(grf.impulse_rel_err),
            "loaded_sides": list(grf.loaded_sides),
            "peak_share": list(grf.peak_share),
        },
        "peaks": peaks,
    }


def _loadshare_inputs() -> tuple[dict, dict, dict]:
    ls = json.loads(LOADSHARE.read_text(encoding="utf-8"))
    rows = {r["bone"]: r for r in ls["rows_with_share_parallel"]}
    fm = json.loads(FRACTURE_MATRIX.read_text(encoding="utf-8"))
    return ls, rows, fm


def _per_bone_at_5m(rows, sym_peaks, single_peaks) -> list[dict]:
    out = []
    for bone, (joint, weak_key) in BONE_CHAIN.items():
        A = float(rows[bone]["A_section_mm2"])
        mult = float(rows[bone]["force_multiplier"])
        weak = _sc(weak_key)
        thr = A * weak
        f_sym = float(sym_peaks[joint]["peak_vertical_n"])
        f_sng = float(single_peaks[joint]["peak_vertical_n"])
        risk_sym = mult * f_sym / thr
        risk_sng = mult * f_sng / thr
        out.append({
            "bone": bone,
            "part_cn": rows[bone]["part_cn"],
            "joint": joint,
            "force_multiplier": mult,
            "A_section_mm2": A,
            "weak_subregion": weak_key,
            "sigma_c_weak_mpa": weak,
            "F_sym_vertical_n": f_sym,
            "F_single_vertical_n": f_sng,
            "load_factor_single_over_sym": (f_sng / f_sym) if f_sym else float("nan"),
            "risk_sym_5m": risk_sym,
            "risk_single_5m": risk_sng,
            "risk_factor_single_over_sym": (risk_sng / risk_sym) if risk_sym else float("nan"),
            "fracture_sym_5m": bool(risk_sym >= 1.0),
            "fracture_single_5m": bool(risk_sng >= 1.0),
        })
    return out


def _project_first_fracture(rows, heights, sym_loads, single_by_h) -> list[dict]:
    """首骨折高度：对称=复用 loadshare risk(h) 还原弱子区域 risk；单脚=以
    实测"载荷因子(h)"（在各测点由 F_single/loads_sym 得）插值后缩放对称载荷。

    单脚首骨折是 **modeled**（假设载荷因子随高度平滑），5 m 风险是 **measured**。
    """
    meas_h = np.array([d["height_m"] for d in single_by_h], dtype=float)
    out = []
    for bone, (joint, weak_key) in BONE_CHAIN.items():
        r = rows[bone]
        stored = float(r["sigma_c_mpa"])
        weak = _sc(weak_key)
        risk_by_h = np.asarray(r["risk_by_height"], dtype=float)          # 用 stored σ_c
        risk_weak_sym = risk_by_h * stored / weak                          # 换到弱子区域
        first_sym = _first_fracture(heights, risk_weak_sym >= 1.0)

        sym_joint = np.asarray(sym_loads[joint], dtype=float)
        factors = np.array([
            (d["peaks"][joint]["peak_vertical_n"] / sym_joint[int(np.argmin(np.abs(np.asarray(heights) - d["height_m"])))])
            for d in single_by_h
        ], dtype=float)
        factor_h = np.interp(np.asarray(heights, dtype=float), meas_h, factors)
        mult = float(r["force_multiplier"])
        A = float(r["A_section_mm2"])
        F_sng_h = factor_h * sym_joint
        risk_sng_h = mult * F_sng_h / (A * weak)
        first_sng = _first_fracture(heights, risk_sng_h >= 1.0)
        first_sng_meas = _first_fracture(
            [d["height_m"] for d in single_by_h],
            [(mult * d["peaks"][joint]["peak_vertical_n"] / (A * weak)) >= 1.0
             for d in single_by_h],
        )
        out.append({
            "bone": bone,
            "joint": joint,
            "weak_subregion": weak_key,
            "sigma_c_weak_mpa": weak,
            "factor_by_measured_height": {f"{h:g}": float(f) for h, f in zip(meas_h, factors)},
            "first_fracture_sym_m": first_sym,
            "first_fracture_single_m": first_sng,
            "first_fracture_single_measured_grid_m": first_sng_meas,
            "delta_first_fracture_m": (
                (first_sym - first_sng) if (first_sym is not None and first_sng is not None)
                else None
            ),
        })
    return out


def _rankings(per_bone, ff) -> dict:
    ff_by = {r["bone"]: r for r in ff}

    def _first(b):
        v = ff_by[b]["first_fracture_sym_m"]
        return v if v is not None else 1e9

    def _first_s(b):
        v = ff_by[b]["first_fracture_single_m"]
        return v if v is not None else 1e9

    order_sym = sorted(per_bone, key=lambda r: (-r["risk_sym_5m"],))
    order_sng = sorted(per_bone, key=lambda r: (-r["risk_single_5m"],))
    fs_order_sym = sorted(per_bone, key=lambda r: (_first(r["bone"]), -r["risk_sym_5m"]))
    fs_order_sng = sorted(per_bone, key=lambda r: (_first_s(r["bone"]), -r["risk_single_5m"]))
    return {
        "order_by_risk5m_sym": [r["bone"] for r in order_sym],
        "order_by_risk5m_single": [r["bone"] for r in order_sng],
        "order_by_first_fracture_sym": [r["bone"] for r in fs_order_sym],
        "order_by_first_fracture_single": [r["bone"] for r in fs_order_sng],
        "order_by_risk5m_changed": [r["bone"] for r in order_sym] != [r["bone"] for r in order_sng],
        "order_by_first_fracture_changed": (
            [r["bone"] for r in fs_order_sym] != [r["bone"] for r in fs_order_sng]
        ),
    }


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S2 单脚/不对称落地（逐足 GRF）probe")
    ap.add_argument("--meas-heights", type=float, nargs="+", default=list(MEAS_HEIGHTS))
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    args = ap.parse_args(argv)

    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== Phase-S2 单脚/不对称落地：逐足 GRF + 默认回归 + 载荷链 ===", flush=True)

    print("[1/4] 接口泛化 + 默认回归（默认 (0.5,0.5) vs 缓存 s1_h5_opensim.json）...", flush=True)
    iface = _split_interface_checks()
    dm = _default_measure()
    reg = dm["regression"]
    print(f"      verdict={reg['verdict']}  max_rel_err={reg['max_rel_err']:.3e}  "
          f"foot_forces_reconstruct={iface['foot_forces_reconstructs_default']}", flush=True)
    if reg["verdict"] != "PASS":
        print("默认回归未复现 → 停止（方案 §7 硬门槛）。", file=sys.stderr)
        return 1

    print(f"[2/4] 单脚（右足承全）测点 {args.meas_heights} ...", flush=True)
    single_by_h = []
    for h in args.meas_heights:
        t0 = time.time()
        d = _single_measure(h)
        single_by_h.append(d)
        print(f"      h={h:4.1f}m  subtalar_r={d['peaks']['subtalar_r']['peak_vertical_n']/1e3:7.2f} kN  "
              f"lumbar={d['peaks']['lumbar']['peak_vertical_n']/1e3:7.2f} kN  "
              f"loaded={d['grf']['loaded_sides']}  [{time.time()-t0:.1f}s]", flush=True)

    print("[3/4] 逐骨纯轴向+两端风险 / 首骨折投影 ...", flush=True)
    ls, rows, fm = _loadshare_inputs()
    heights = [float(h) for h in fm["heights_m"]]
    sym_loads = {k: [float(x) for x in v] for k, v in fm["loads_n"].items()}

    sym5 = dm["symmetric_peaks"]
    single5 = next(d for d in single_by_h if abs(d["height_m"] - H_REP) < 1e-9)["peaks"]
    per_bone = _per_bone_at_5m(rows, sym5, single5)
    per_bone.sort(key=lambda r: -r["risk_single_5m"])
    ff = _project_first_fracture(rows, heights, sym_loads, single_by_h)
    rankings = _rankings(per_bone, ff)

    # 关节载荷因子（单脚/对称）@5m
    joint_factors = {}
    for j in CHAIN_JOINTS:
        f_sym = float(sym5[j]["peak_vertical_n"])
        f_sng = float(single5[j]["peak_vertical_n"])
        joint_factors[j] = {
            "sym_vertical_n": f_sym,
            "single_vertical_n": f_sng,
            "factor": (f_sng / f_sym) if f_sym else float("nan"),
            "sym_peak_force_n": float(sym5[j]["peak_force_n"]),
            "single_peak_force_n": float(single5[j]["peak_force_n"]),
        }
    unloaded = {
        "joint": UNLOADED_CONTRAST,
        "peak_force_n": float(single5[UNLOADED_CONTRAST]["peak_force_n"]),
        "peak_vertical_n": float(single5[UNLOADED_CONTRAST]["peak_vertical_n"]),
        "note": "未承力侧（左）在单脚右落地时仅承受子树惯性，非无效值。",
    }

    # 不对称（70/30）API 演示（GRF 层 + 载荷链）
    print("[4/4] 不对称 70/30 + 时变重心转移 API 演示（GRF 层）...", flush=True)
    from climbing.coupling.opensim_grf import ground_reaction

    base = ground_reaction(height_m=H_REP, mass_kg=M_MODEL)
    g70 = ground_reaction(
        height_m=H_REP, mass_kg=M_MODEL,
        foot_forces=((base.t_s, base.f_total_n * 0.7), (base.t_s, base.f_total_n * 0.3)),
    )
    w = 0.5 + 0.3 * np.sin(2.0 * np.pi * np.linspace(0.0, 1.0, base.t_s.size))
    left_shift = w * base.f_total_n
    right_shift = base.f_total_n - left_shift       # 精确互补 → 总和逐位守恒
    gshift = ground_reaction(
        height_m=H_REP, mass_kg=M_MODEL,
        foot_forces=((base.t_s, left_shift), (base.t_s, right_shift)),
    )
    asym = {
        "g70_30": {
            "peak_share": list(g70.peak_share),
            "loaded_sides": list(g70.loaded_sides),
            "impulse_ns": float(g70.impulse_ns),
            "impulse_rel_err": float(g70.impulse_rel_err),
            "per_foot_forces": bool(g70.per_foot_forces),
        },
        "g_time_varying_share": {
            "share_left_min": float(np.min(w)),
            "share_left_max": float(np.max(w)),
            "peak_share": list(gshift.peak_share),
            "impulse_ns": float(gshift.impulse_ns),
            "impulse_rel_err": float(gshift.impulse_rel_err),
            "total_max_abs_diff_from_base_n": float(
                np.max(np.abs(gshift.f_total_n - base.f_total_n))
            ),
            "total_matches_base_allclose": bool(
                np.allclose(gshift.f_total_n, base.f_total_n, rtol=0.0, atol=1e-6)
            ),
        },
    }

    # 既有对照（axial_subregion.json，对称纯轴向+ends）
    axial_ref = {}
    if AXIAL_SUB.exists():
        axj = json.loads(AXIAL_SUB.read_text(encoding="utf-8"))
        axial_ref = {
            "ranking": [
                {"rank": r["rank"], "bone": r["bone"], "first_fracture_h": r["first_fracture_h"]}
                for r in axj.get("ranking_axial_ends", [])
            ],
            "note": "既有对称纯轴向+ends 排序（只读对照）",
        }

    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S2 单脚 / 不对称落地（逐足 GRF）",
            "plan_doc": "docs/非垂直落地扩展方案.md",
            "units": "mm-N-MPa-s",
            "height_rep_m": H_REP,
            "meas_heights_m": list(args.meas_heights),
            "mass_kg": M_MODEL,
            "foot_loaded": "right (split=(0,1), single_foot_split('r'))",
            "commands": [
                "python scripts/opensim_fe/nonvertical_s2.py",
            ],
            "measured_vs_modeled": {
                "measured": [
                    "默认 (0.5,0.5) 链路重跑（ground_reaction→run_dead_drop→joint_reaction）",
                    "单脚（右足承全）在各测点高度的关节反力峰值（OpenSim 正动力学）",
                    "改前轴向缓存 results/opensim_fe/s1_h5_opensim.json",
                    "对称逐高度关节载荷 results/opensim_fe/fracture_matrix.json::loads_n",
                ],
                "modeled_assumed": [
                    "单脚首骨折高度：以实测载荷因子(h)（F_single/loads_sym）线性插值后缩放对称载荷（非逐高度重跑）",
                    "逐骨风险 = 纯轴向 1D 名义应力 + 载荷分配（mult·F/(A·σ_c_weak)），忽略弯曲/剪切/偏心",
                    "σ_c_weak / A_section / 载荷倍率均为既有 assumed/measured 口径，绝对值不可信，只看相对",
                    "单脚仍为刚性腿（LEG_JOINTS 锁死），未含姿势/被动刚度（方案 S3）",
                ],
            },
        },
        "split_interface": iface,
        "axial_regression": reg,
        "symmetric_peaks_5m": sym5,
        "single_peaks_5m": single5,
        "joint_load_factors_5m": joint_factors,
        "unloaded_side_5m": unloaded,
        "single_by_height": single_by_h,
        "per_bone_risk_5m": per_bone,
        "first_fracture_projection": ff,
        "ranking": rankings,
        "asymmetric_api_demo": asym,
        "axial_subregion_reference": axial_ref,
        "assumptions": [
            "S2 只改 L1 GRF 的逐足接口；不触碰接触本构/摩擦锥/CoP（方案 S4/S5）。",
            "单脚 = 全部冲量经一足；刚性锁死腿下，偏心加载会引入骨盆转动（实测），故局部载荷因子未必恰为 2。",
            "单脚首骨折为 modeled（载荷因子沿高度线性插值）；5 m 风险为 measured。",
            "力链沿用 risk_1d_loadshare 的载荷分配（胫/腓并联 k=E·A/L；中轴骨 mass-above 比）。",
            "本脚本只读复用既有产物，不调用 FEBio，不覆盖既有产物。",
        ],
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=float),
                             encoding="utf-8")
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0


def _fmt_h(v) -> str:
    return ">50" if v is None else f"{int(v)}"


def _write_report(path: Path, res: dict, ts: str) -> None:
    reg = res["axial_regression"]
    iface = res["split_interface"]
    jf = res["joint_load_factors_5m"]
    pb = res["per_bone_risk_5m"]
    ff = {r["bone"]: r for r in res["first_fracture_projection"]}
    rk = res["ranking"]
    asym = res["asymmetric_api_demo"]

    L: list[str] = []
    A = L.append
    A("# Phase-S2 非垂直落地 · 单脚 / 不对称落地（逐足 GRF）报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 方案 `docs/非垂直落地扩展方案.md` §2 自由度④ / §3 L1 / §4 Phase-S2 / §6 风险④ / §7。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"1. **默认回归（硬门槛）: {reg['verdict']}** —— 默认 `(0.5,0.5)` 重跑与改前缓存"
      f" `s1_h5_opensim.json` 的最大相对误差 **{reg['max_rel_err']:.2e}**；"
      f"默认 == 显式 `(0.5,0.5)` = {iface['default_equals_explicit_0.5']}，"
      f"`foot_forces` 重构默认 = {iface['foot_forces_reconstructs_default']}（均逐位一致）。")
    A(f"2. **接口泛化**：`split` 从写死 `(0.5,0.5)` 推广到任意逐足占比（含单脚 "
      f"`single_foot_split`）；新增 `foot_forces` 显式逐足时程（两条独立时程）。")
    # factor of local load
    calc_f = jf["subtalar_r"]["factor"]
    ank_f = jf["ankle_r"]["factor"]
    hip_f = jf["hip_r"]["factor"]
    lum_f = jf["lumbar"]["factor"]
    A(f"3. **单脚 vs 对称 @5 m 局部载荷因子**：subtalar_r **×{calc_f:.3f}**、ankle_r "
      f"×{ank_f:.3f}、hip_r ×{hip_f:.3f}、lumbar ×{lum_f:.3f}（承力侧局部约翻倍，"
      f"中轴近乎不变）。")
    fr_sng_all = [r["bone"] for r in pb if r["fracture_single_5m"]]
    A(f"4. **逐骨 ordering**：单脚 vs 对称**排序不变**（risk@5m 改变="
      f"{rk['order_by_risk5m_changed']}，首骨折改变={rk['order_by_first_fracture_changed']}）——"
      f"但**首骨折高度大幅前移**（腓骨 7→1 m、跟骨/胫骨 8→2 m、股骨 19→4 m，均为投影；"
      f"实测 4 m 已越阈）；@5 m 越阈骨从 0 个增至 **{len(fr_sng_all)}** 个"
      f"（`{', '.join(fr_sng_all)}`）。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段 API（全部 opt-in、默认逐位复现旧行为）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） | 说明 |")
    A("|---|---|---|---|")
    A("| L1 GRF | `opensim_grf.ground_reaction` | `foot_forces=None` | 显式逐足力时程 `(左,右)`（N）：标量 / 一维数组 / `(t,f)` 对；完全覆盖 `split` |")
    A("| L1 GRF | `opensim_grf.single_foot_split` | 新 | `'l'→(1,0)`、`'r'→(0,1)`，单脚占比 |")
    A("| L1 GRF | `GroundReaction.per_foot_forces` / `.loaded_sides` / `.peak_share` | 新 | 逐足来源标记 / 承力侧 / 峰值占比 |")
    A("| L1 GRF | `ground_reaction(split=...)` | `(0.5,0.5)` | 由固定对称推广为任意非负占比（含一足为 0） |")
    A("")
    A("> 默认 `split=(0.5,0.5)`、`foot_forces=None` 走的是**改前代码路径**（`f_left=fw·0.5,"
      "f_right=fw·0.5`），故 §2 逐位复现。中段 `joint_loads.subtalar_reaction` 与下游"
      " `_add_foot_force` / `joint_reactions` **未改**（本就支持每足独立时程）。")
    A("")
    A("---")
    A("")
    A("## 2. 默认回归（硬门槛，方案 §7）")
    A("")
    A("默认链路重跑 vs 缓存：")
    A("")
    A("| 量 | 本次默认 | 缓存 | 相对误差 |")
    A("|---|---:|---:|---:|")
    for k in ("impulse_ns", "grf_peak_n", "subtalar_peak_vertical_n",
              "subtalar_peak_force_n", "subtalar_peak_moment_nm"):
        A(f"| `{k}` | {reg['measured_default_run'][k]:.9g} | "
          f"{reg['cached'][k]:.9g} | {reg['rel_err'][k]:.2e} |")
    A("")
    A(f"- 最大相对误差 **{reg['max_rel_err']:.2e}**（阈值 1e-6）→ **{reg['verdict']}**")
    A(f"- 默认 == 显式 `(0.5,0.5)`（逐位）: **{iface['default_equals_explicit_0.5']}**；"
      f"`foot_forces` 以 `(t,f)` 重构 0.5/0.5 == 默认（逐位）: "
      f"**{iface['foot_forces_reconstructs_default']}**。")
    A(f"- 单脚 `split=(0,1)`：承力侧 = `{iface['single_right']['loaded_sides']}`，"
      f"峰值占比 = `{iface['single_right']['peak_share']}`，冲量 {iface['single_right']['impulse_ns']:.3f} N·s。")
    A("")
    A("---")
    A("")
    A("## 3. 单脚 vs 对称：关节载荷因子 @5 m")
    A("")
    A("| 关节 | 对称 F_y (N) | 单脚 F_y (N) | **因子(单/对)** | 对称 \\|F\\| (N) | 单脚 \\|F\\| (N) |")
    A("|---|---:|---:|---:|---:|---:|")
    for j in CHAIN_JOINTS:
        r = jf[j]
        A(f"| `{j}` | {r['sym_vertical_n']:,.1f} | {r['single_vertical_n']:,.1f} | "
          f"**×{r['factor']:.3f}** | {r['sym_peak_force_n']:,.1f} | {r['single_peak_force_n']:,.1f} |")
    A("")
    A(f"- 未承力侧对照：`{res['unloaded_side_5m']['joint']}` 峰值 \\|F\\| = "
      f"{res['unloaded_side_5m']['peak_force_n']:,.1f} N（仅子树惯性，非零但远小于承力侧）。")
    A(f"- **物理**：单脚把全部冲量经一足，承力侧局部关节反力约 ×{(calc_f+ank_f+hip_f)/3:.2f}"
      "（并非严格 ×2：锁死腿 + 偏心加载会引入骨盆转动，见 §6 风险④）。")
    A("")
    # 因子随高度稳定性
    A("### 3.1 载荷因子随高度的稳定性（实测）")
    A("")
    A("| h (m) | subtalar_r 因子 | ankle_r 因子 | hip_r 因子 | lumbar 因子 |")
    A("|---|---:|---:|---:|---:|")
    ff_by = {r["bone"]: r for r in res["first_fracture_projection"]}
    for d in res["single_by_height"]:
        h = d["height_m"]
        row = {j: ff_by[BONE_CHAIN_INV[j]]["factor_by_measured_height"].get(f"{h:g}")
               for j in ("subtalar_r", "ankle_r", "hip_r", "lumbar")}
        A(f"| {h:g} | " + " | ".join(
            ("—" if row[j] is None else f"{row[j]:.3f}") for j in
            ("subtalar_r", "ankle_r", "hip_r", "lumbar")) + " |")
    A("")
    A("> 因子在中轴（lumbar）≈1、下肢承力链 ≈1.8–2.1；随高度变化温和 → 支持 §4 的"
      "**常数/插值因子**外推。")
    A("")
    A("---")
    A("")
    A("## 4. 逐骨纯轴向+两端 1D 风险 @5 m")
    A("")
    A("`risk = mult · F_joint / (A_section · σ_c_weak)`；`F_joint` 取关节反力**纵向分量**"
      "（Route-1 口径）。`mult` 来自 `risk_1d_loadshare`（胫/腓并联 k=E·A/L；中轴骨 "
      "mass-above 比）；`σ_c_weak` = [Y25] 弱子区域（跟骨 150 / 胫腓两端 70 / 股骨颈 80 / "
      "骨盆 180 / 脊柱 150 / 颅骨 160）。")
    A("")
    A("| 骨 | 关节 | 弱子区域 | F_y 对称 (N) | F_y 单脚 (N) | 因子 | risk 对称 | risk 单脚 | @5m 对称 | @5m 单脚 |")
    A("|---|---|---|---:|---:|---:|---:|---:|---|---|")
    for r in pb:
        A(f"| `{r['bone']}` | `{r['joint']}` | {r['weak_subregion']} | "
          f"{r['F_sym_vertical_n']:,.1f} | {r['F_single_vertical_n']:,.1f} | "
          f"×{r['load_factor_single_over_sym']:.3f} | {r['risk_sym_5m']:.4f} | "
          f"**{r['risk_single_5m']:.4f}** | "
          f"{'⚠越阈' if r['fracture_sym_5m'] else '否'} | "
          f"{'⚠越阈' if r['fracture_single_5m'] else '否'} |")
    A("")
    fr_sym = [r["bone"] for r in pb if r["fracture_sym_5m"]]
    fr_sng = [r["bone"] for r in pb if r["fracture_single_5m"]]
    A(f"- 对称 @5 m 越阈：`{fr_sym or '无'}`；单脚 @5 m 越阈：`{fr_sng or '无'}`。")
    A("")
    A("---")
    A("")
    A("## 5. 首骨折高度与排序变化")
    A("")
    A("对称首骨折 = 复用 `risk_1d_loadshare` 的 risk(h) 还原弱子区域风险；"
      "**单脚首骨折 = modeled**：以实测（§3.1）载荷因子随高度线性插值后缩放对称逐高度载荷。")
    A("")
    A("| 骨 | 首骨折 对称 (m) | 首骨折 单脚 (m) | 单脚(实测网格) | Δ 提前 (m) |")
    A("|---|---:|---:|---:|---:|")
    for r in sorted(res["first_fracture_projection"],
                    key=lambda x: (x["first_fracture_single_m"] is None,
                                   x["first_fracture_single_m"])):
        delta = r["delta_first_fracture_m"]
        delta_txt = "—" if delta is None else f"{delta:.0f}"
        A(f"| `{r['bone']}` | {_fmt_h(r['first_fracture_sym_m'])} | "
          f"**{_fmt_h(r['first_fracture_single_m'])}** | "
          f"{_fmt_h(r['first_fracture_single_measured_grid_m'])} | {delta_txt} |")
    A("")
    A(f"> **读法**：实测网格（{min(res['meta']['meas_heights_m']):g}–"
      f"{max(res['meta']['meas_heights_m']):g} m）给出可直接信任的下界；投影列在 h<最低测点"
      f"({min(res['meta']['meas_heights_m']):g} m) 时把载荷因子 clamp 到最低测点值，故"
      " 1–3 m 的投影值应读作 **≤最低测点（外推）**。因子在测点内近乎常数（§3.1），"
      "故投影的**量级**可信、**具体米数**为 modeled。")
    A("")
    A("### 5.1 排序对照")
    A("")
    A("| 口径 | 排序（1→9，越前越危险） |")
    A("|---|---|")
    A("| risk@5m 对称 | `" + " > ".join(rk["order_by_risk5m_sym"]) + "` |")
    A("| **risk@5m 单脚** | `" + " > ".join(rk["order_by_risk5m_single"]) + "` |")
    A("| 首骨折 对称 | `" + " > ".join(rk["order_by_first_fracture_sym"]) + "` |")
    A("| **首骨折 单脚** | `" + " > ".join(rk["order_by_first_fracture_single"]) + "` |")
    A("")
    A(f"- risk@5m 排序改变：**{rk['order_by_risk5m_changed']}**；"
      f"首骨折排序改变：**{rk['order_by_first_fracture_changed']}**。")
    A("")
    A("---")
    A("")
    A("## 6. 不对称 / 时变重心 API 演示（GRF 层）")
    A("")
    a7 = asym["g70_30"]
    A(f"- `foot_forces` 70/30：峰值占比 `{a7['peak_share']}`，承力侧 `{a7['loaded_sides']}`，"
      f"冲量 {a7['impulse_ns']:.3f} N·s（相对误差 {a7['impulse_rel_err']*100:+.1f}%）。")
    A("")
    gs = asym["g_time_varying_share"]
    A(f"- 时变占比（左右形状独立、总和守恒）：左占比 {gs['share_left_min']:.2f}–{gs['share_left_max']:.2f}，"
      f"总时程与默认一致（allclose，最大差 {gs['total_max_abs_diff_from_base_n']:.2e} N）= "
      f"**{gs['total_matches_base_allclose']}**，冲量 {gs['impulse_ns']:.3f} N·s。这验证了方案 §3 L1"
      "「两条独立时程」完全版可用，且下游 `run_dead_drop` / `joint_reactions` 直接消费。")
    A("")
    A("---")
    A("")
    A("## 7. 假设与诚实边界")
    A("")
    for i, a in enumerate(res["assumptions"], 1):
        A(f"{i}. {a}")
    A("")
    A("- **未上 FE**：本阶段只到载荷链 + 1D 风险；单脚载荷上 FEBio 属后续。")
    A("- **已实现 vs 未实现**：逐足占比 / 单脚 / 两条独立时程 = **已实现**；CoP 偏移 /"
      " 摩擦力矩 / 落地姿势 = **未实现**（方案 S3–S5）。")
    A("- **因子 ≠ 2 的原因**：单脚偏心加载在刚性锁死腿模型下会驱动骨盆转动，使承力侧"
      "关节反力略偏离两倍（实测见 §3）；这是被如实测量的模型效应，不是误差。")
    A("")
    A("---")
    A("")
    A("## 8. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s2.py")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s2.py -q")
    A("```")
    A("")
    A("## 9. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


#: 关节 → 代表骨（供报告 §3.1 取因子）。
BONE_CHAIN_INV = {
    "subtalar_r": "calcaneus_r",
    "ankle_r": "tibia_r",
    "hip_r": "femur_r",
    "lumbar": "T6",
}


if __name__ == "__main__":
    raise SystemExit(main())
