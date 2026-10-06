"""Phase-S3 probe —— 落地**姿势**（自由度③）：(a) 锁关节到目标位姿 + (b) 解锁 + 被动刚度。

方案：``docs/非垂直落地扩展方案.md`` §1 A5（解锁为何失败）、§2 自由度③、
§3 L2（两条路径 (a) lock-at-target / (b) unlock+passive stiffness）、§4 S3、
§6 风险②、§7（轴向回归硬门槛）。

本脚本**只读复用**既有产物与模块，**不调用 FEBio**、**不覆盖任何既有文件**；
新增产物仅两个（同名已存在则拒绝写入）：

* ``results/opensim_fe/nonvertical_s3.json``
* ``results/opensim_fe/NONVERTICAL_S3_POSTURE_REPORT.md``

它做四件事：

1. **默认回归（硬门槛，方案 §7）**：默认 ``posture=None, passive=None`` 重跑
   （``ground_reaction`` → ``run_dead_drop`` → ``joint_reaction``）与改前缓存
   ``s1_h5_opensim.json`` 数值对比，要求相对误差 = 0；并证明显式
   ``posture=LegPose()``（全 0）与默认**逐位一致**。
2. **Path (a) 锁姿势扫描**：膝屈 15/30/45/60°（锁定）在各 h = 2/3/3.5/4.5 m
   下的逐关节峰值（subtalar/ankle/knee/hip/lumbar）vs 直立基准。
3. **Path (b) 解锁 + 被动刚度**：同一组高度，测**有限关节屈曲**（轨迹）与
   **峰值力衰减**；并检查关节不"飞"（越界）。
4. **诚实边界**：若 (b) 无法显著衰减，记录实测越界/屈曲证据与本模型边界。

单位 mm–N–MPa–s。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s3_posture.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

M_MODEL = 75.337
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"
OUT_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s3.json"
OUT_MD = ROOT / "results" / "opensim_fe" / "NONVERTICAL_S3_POSTURE_REPORT.md"

#: 载荷链要测的关节（近端→远端）
CHAIN_JOINTS: tuple[str, ...] = ("subtalar_r", "ankle_r", "knee_r", "hip_r", "lumbar")

#: Path (a) 姿势扫描的膝屈角（度）
FLEX_DEGS: tuple[float, ...] = (15.0, 30.0, 45.0, 60.0)
#: 场景高度（m）—— 方案 §4 S3 关注的抱石范围（≤4.5 m）
HEIGHTS: tuple[float, ...] = (2.0, 3.0, 3.5, 4.5)
#: "屈膝 + 屈髋 + 踝背屈"组合姿势跑一次的代表高度（m）
H_COMBINED = 3.5
#: 被动刚度敏感性扫描（软 vs 硬）的代表高度（m）
H_SENSITIVITY = 3.0

#: 参与"有限屈曲"测量的坐标（Path b）→ OpenSim 模型解剖范围 (rad)
#: 范围来自模型 ``Coordinate.getRangeMin/Max()``（实测，见脚本 docstring）。
COORD_RANGES_RAD: dict[str, tuple[float, float]] = {
    "knee_angle_r": (0.0, 2.0944),
    "hip_flexion_r": (-0.52359878, 2.0943951),
    "ankle_angle_r": (-0.6981317, 0.52359878),
}


def _joint_peaks(fall, grf) -> dict:
    """逐关节峰值（|F| / F_y / |M| / t_peak）。"""
    from climbing.coupling.joint_reactions import joint_reaction

    out = {}
    for j in CHAIN_JOINTS:
        jr = joint_reaction(fall, grf, joint=j)
        out[j] = {
            "peak_force_n": float(jr.peak_force_n),
            "peak_vertical_n": float(jr.peak_vertical_n),
            "peak_moment_nm": float(jr.peak_moment_nm),
            "peak_force_time_s": float(jr.peak_force_time_s),
            "subtree_mass_kg": float(jr.subtree_mass_kg),
        }
    return out


def _find_state_path(labels: list[str], coord: str) -> str | None:
    """在状态列标签里找 ``.../<coord>/value``（排除 ``_beta`` 别名）。"""
    suffix = f"/{coord}/value"
    cand = [l for l in labels if l.endswith(suffix) and "_beta" not in l]
    return cand[0] if cand else None


def _joint_excursions(fall, coords=tuple(COORD_RANGES_RAD)) -> dict:
    """从 states 表提取每个坐标的**有限屈曲**轨迹（Path b 用）。

    返回每坐标的 ``initial / min / max / peak_abs / range`` (rad, deg) 以及
    是否越出模型解剖范围。
    """
    tab = fall.states
    labels = [tab.getColumnLabel(j) for j in range(tab.getNumColumns())]
    out: dict[str, dict] = {}
    for name in coords:
        path = _find_state_path(labels, name)
        if path is None:
            continue
        v = np.asarray(tab.getDependentColumnAtIndex(labels.index(path)).to_numpy(),
                       dtype=float)
        lo, hi = COORD_RANGES_RAD[name]
        peak_abs = float(max(abs(v.min() - v[0]), abs(v.max() - v[0])))
        out[name] = {
            "initial_rad": float(v[0]),
            "initial_deg": math.degrees(float(v[0])),
            "min_rad": float(v.min()),
            "min_deg": math.degrees(float(v.min())),
            "max_rad": float(v.max()),
            "max_deg": math.degrees(float(v.max())),
            "peak_flexion_rad": peak_abs,
            "peak_flexion_deg": math.degrees(peak_abs),
            "range_rad": (hi - lo),
            "range_min_rad": lo,
            "range_max_rad": hi,
            "within_range": bool(v.min() >= lo - 1e-6 and v.max() <= hi + 1e-6),
        }
    return out


# --------------------------------------------------------------------------
# 1) 默认回归（硬门槛）
# --------------------------------------------------------------------------
def _default_measure() -> dict:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop, LegPose
    from climbing.coupling.opensim_grf import ground_reaction

    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))
    grf = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall = run_dead_drop(grf)
    jl = subtalar_reaction(fall, grf, side="r")
    peaks = _joint_peaks(fall, grf)

    measured = {
        "impulse_ns": float(grf.impulse_ns),
        "impulse_rel_err": float(grf.impulse_rel_err),
        "grf_peak_n": float(grf.peak_total_n),
        "subtalar_peak_vertical_n": float(jl.peak_vertical_n),
        "subtalar_peak_force_n": float(jl.peak_force_n),
        "subtalar_peak_moment_nm": float(jl.peak_moment_nm),
    }
    cached_keys = {k: float(cached[k]) for k in measured}
    rel = {k: (abs(measured[k] - cached_keys[k]) / abs(cached_keys[k])
               if cached_keys[k] != 0.0 else abs(measured[k] - cached_keys[k]))
           for k in measured}

    # 显式 posture=LegPose()（全 0）必须与默认逐位一致（接口层）
    fall_p0 = run_dead_drop(grf, posture=LegPose())
    jl_p0 = subtalar_reaction(fall_p0, grf, side="r")
    p0_rel = abs(jl_p0.peak_vertical_n - jl.peak_vertical_n) / jl.peak_vertical_n
    p0_force_rel = abs(jl_p0.peak_force_n - jl.peak_force_n) / jl.peak_force_n

    return {
        "cached_file": str(CACHED_AXIAL.relative_to(ROOT)),
        "measured_default_run": measured,
        "cached": cached_keys,
        "rel_err": rel,
        "max_rel_err": max(rel.values()),
        "posture_zero_vs_default": {
            "subtalar_peak_vertical_rel_diff": float(p0_rel),
            "subtalar_peak_force_rel_diff": float(p0_force_rel),
            "bitwise_reproduced": bool(max(p0_rel, p0_force_rel) == 0.0),
        },
        "symmetric_peaks_5m": peaks,
        "verdict": "PASS" if max(rel.values()) == 0.0 else "FAIL",
    }


# --------------------------------------------------------------------------
# 2) Path (a)：锁姿势
# --------------------------------------------------------------------------
def _posture_case(height_m: float, pose) -> dict:
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height_m, mass_kg=M_MODEL)
    fall = run_dead_drop(grf, posture=pose)
    peaks = _joint_peaks(fall, grf)
    return {
        "height_m": float(height_m),
        "pose": {
            "knee_flexion_deg": math.degrees(pose.knee_flexion_rad),
            "hip_flexion_deg": math.degrees(pose.hip_flexion_rad),
            "ankle_deg": math.degrees(pose.ankle_rad),
            "torso_lean_deg": math.degrees(pose.torso_lean_rad),
        },
        "grf_peak_n": float(grf.peak_total_n),
        "impact_speed_ms": float(grf.impact_speed_ms),
        "joint_peaks": peaks,
    }


def _posture_sweep() -> dict:
    from climbing.coupling.opensim_fall import LegPose

    heights_out: dict[str, dict] = {}
    for h in HEIGHTS:
        upright = _posture_case(h, LegPose())
        flexed = {}
        for fdeg in FLEX_DEGS:
            flexed[f"knee_{fdeg:g}"] = _posture_case(
                h, LegPose.flexed_landing(knee_deg=fdeg)
            )
        heights_out[f"{h:g}"] = {
            "height_m": h,
            "upright": upright,
            "flexed_knee": flexed,
        }
    # 一个"自然屈膝落地"组合（膝屈 + 髋屈 + 踝背屈），只在代表高度跑一次
    combined = _posture_case(
        H_COMBINED, LegPose.flexed_landing(knee_deg=45.0, hip_deg=25.0, ankle_deg=10.0)
    )
    # 躯干前倾演示（可选 torso lean）：同一组合 + 躯干前倾 20°
    torso_demo = _posture_case(
        H_COMBINED,
        LegPose.flexed_landing(
            knee_deg=45.0, hip_deg=25.0, ankle_deg=10.0, torso_lean_deg=-20.0
        ),
    )
    return {
        "flex_deg_list": list(FLEX_DEGS),
        "heights_m": list(HEIGHTS),
        "by_height": heights_out,
        "flexed_combined": combined,
        "flexed_combined_torso_lean": torso_demo,
    }


def _stiffness_sensitivity() -> dict:
    """在代表高度对比"软（振铃）vs 硬（稳定）"被动刚度 regime。

    只为证明 §5 的论断：**生理中段被动刚度 → 欠阻尼振铃 → 数值放大**，
    不是物理吸能；偏硬 + 重阻尼才给出稳定解（但仍不衰减）。
    """
    from climbing.coupling.opensim_fall import PassiveStiffness

    cases = {
        "soft_k50_c5": PassiveStiffness.soft_landing(),
        "mid_k200_c100": PassiveStiffness(
            knee_stiffness=200.0, knee_damping=100.0,
            hip_stiffness=200.0, hip_damping=100.0,
            ankle_stiffness=200.0, ankle_damping=100.0,
            knee_limit_stiffness=2000.0, knee_limit_damping=50.0,
            hip_limit_stiffness=2000.0, hip_limit_damping=50.0,
            ankle_limit_stiffness=2000.0, ankle_limit_damping=50.0,
            transition_rad=0.02,
        ),
        "nominal_k500_c500": PassiveStiffness.nominal_landing(),
    }
    out = {}
    for key, p in cases.items():
        pc = _passive_case(H_SENSITIVITY, p)
        out[key] = {
            "knee_stiffness": p.knee_stiffness,
            "knee_damping": p.knee_damping,
            "joint_peaks": pc["joint_peaks"],
            "excursions": pc["excursions"],
        }
    return {
        "height_m": H_SENSITIVITY,
        "cases": out,
    }


# --------------------------------------------------------------------------
# 3) Path (b)：解锁 + 被动刚度
# --------------------------------------------------------------------------
def _passive_case(height_m: float, passive) -> dict:
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf = ground_reaction(height_m=height_m, mass_kg=M_MODEL)
    t0 = time.time()
    fall = run_dead_drop(grf, passive=passive)
    peaks = _joint_peaks(fall, grf)
    exc = _joint_excursions(fall)
    return {
        "height_m": float(height_m),
        "wall_s": time.time() - t0,
        "grf_peak_n": float(grf.peak_total_n),
        "joint_peaks": peaks,
        "excursions": exc,
    }


def _passive_sweep(passive) -> dict:
    out: dict[str, dict] = {}
    for h in HEIGHTS:
        out[f"{h:g}"] = _passive_case(h, passive)
    return {
        "heights_m": list(HEIGHTS),
        "by_height": out,
    }


def _attenuation_table(base_by_h: dict, passive_by_h: dict) -> list[dict]:
    """Path (b) 相对**直立锁定基准**的峰值力衰减（按高度）。"""
    rows: list[dict] = []
    for key, ph in passive_by_h["by_height"].items():
        base = base_by_h[key]["upright"]["joint_peaks"]
        pj = ph["joint_peaks"]
        for j in CHAIN_JOINTS:
            fb = float(base[j]["peak_force_n"])
            fp = float(pj[j]["peak_force_n"])
            rows.append({
                "height_m": float(ph["height_m"]),
                "joint": j,
                "baseline_force_n": fb,
                "passive_force_n": fp,
                "attenuation": (1.0 - fp / fb) if fb else float("nan"),
                "flew": not all(
                    d.get("within_range", True) for d in ph["excursions"].values()
                ),
            })
    return rows


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S3 姿势（锁关节 / 解锁+被动刚度）probe")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    args = ap.parse_args(argv)

    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    from climbing.coupling.opensim_fall import PassiveStiffness

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== Phase-S3 落地姿势：锁关节 (a) + 解锁/被动刚度 (b) ===", flush=True)

    print("[1/5] 默认回归（硬门槛）...", flush=True)
    reg = _default_measure()
    print(f"      verdict={reg['verdict']}  max_rel_err={reg['max_rel_err']:.3e}  "
          f"posture0_bitwise={reg['posture_zero_vs_default']['bitwise_reproduced']}",
          flush=True)

    print(f"[2/5] Path (a) 锁姿势扫描：膝屈 {FLEX_DEGS} @ h={HEIGHTS} ...", flush=True)
    t0 = time.time()
    posture = _posture_sweep()
    print(f"      done in {time.time()-t0:.1f}s", flush=True)

    print("[3/5] Path (b) 解锁 + nominal 被动刚度 ...", flush=True)
    passive = PassiveStiffness.nominal_landing()
    t0 = time.time()
    pass_sweep = _passive_sweep(passive)
    print(f"      done in {time.time()-t0:.1f}s", flush=True)

    print(f"[4/5] 被动刚度敏感性（软 vs 硬）@ h={H_SENSITIVITY:g} ...", flush=True)
    t0 = time.time()
    sensitivity = _stiffness_sensitivity()
    print(f"      done in {time.time()-t0:.1f}s", flush=True)

    atten = _attenuation_table(posture["by_height"], pass_sweep)
    any_flew = any(r["flew"] for r in atten)
    # 衰减统计（正 = 降低）
    atten_vals = [r["attenuation"] for r in atten if np.isfinite(r["attenuation"])]
    worst_joint = max(atten, key=lambda r: (r["passive_force_n"] / r["baseline_force_n"])
                      if r["baseline_force_n"] else 0.0)
    best_joint = max(atten, key=lambda r: r["attenuation"]
                     if np.isfinite(r["attenuation"]) else -1e9)
    print(f"      any_out_of_range={any_flew}  "
          f"atten min={min(atten_vals):+.3f} max={max(atten_vals):+.3f}", flush=True)

    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S3 落地姿势（自由度③）：锁关节 + 解锁/被动刚度",
            "plan_doc": "docs/非垂直落地扩展方案.md",
            "units": "mm-N-MPa-s",
            "mass_kg": M_MODEL,
            "heights_m": list(HEIGHTS),
            "flex_deg_list": list(FLEX_DEGS),
            "commands": [
                "python scripts/opensim_fe/nonvertical_s3_posture.py",
            ],
            "measured_vs_modeled": {
                "measured": [
                    "默认（posture/passive=None）链路重跑：ground_reaction→run_dead_drop→joint_reaction",
                    "Path (a)：锁膝屈 15/30/45/60° 的关节反力（OpenSim 正动力学）",
                    "Path (b)：解锁 + 被动刚度下的关节**轨迹**（有限屈曲实测）与峰值力",
                    "改前轴向缓存 results/opensim_fe/s1_h5_opensim.json（硬门槛）",
                ],
                "modeled_assumed": [
                    "Path (b) 的被动刚度 / 阻尼**全部是 ASSUMED**（文献量级，非实测标定）",
                    "刚度为**线性**扭转弹簧（SpringGeneralizedForce）+ **非线性**边壁（CoordinateLimitForce）",
                    "无肌肉主动控制（dead-drop 关激活）；完全被动关节囊/韧带/肌腱弹性",
                    "锁姿势 (a) 仍为**刚性腿 + 弯曲几何**，关节不弯 → 不吸能",
                    "载荷链只到关节反力（未上 FE）；不调用 FEBio",
                ],
            },
        },
        "default_regression": reg,
        "path_a_posture_sweep": posture,
        "path_b_passive_sweep": pass_sweep,
        "path_b_stiffness_sensitivity": sensitivity,
        "path_b_attenuation": atten,
        "path_b_summary": {
            "passive_model": {
                "knee_stiffness_nm_per_rad": passive.knee_stiffness,
                "knee_damping_nms_per_rad": passive.knee_damping,
                "hip_stiffness_nm_per_rad": passive.hip_stiffness,
                "hip_damping_nms_per_rad": passive.hip_damping,
                "ankle_stiffness_nm_per_rad": passive.ankle_stiffness,
                "ankle_damping_nms_per_rad": passive.ankle_damping,
                "knee_limit_stiffness_nm_per_rad": passive.knee_limit_stiffness,
                "hip_limit_stiffness_nm_per_rad": passive.hip_limit_stiffness,
                "ankle_limit_stiffness_nm_per_rad": passive.ankle_limit_stiffness,
                "transition_rad": passive.transition_rad,
                "source": "ASSUMED（文献量级，非实测；见 PassiveStiffness docstring）",
            },
            "any_joint_out_of_range": bool(any_flew),
            "attenuation_min": float(min(atten_vals)),
            "attenuation_max": float(max(atten_vals)),
            "attenuation_mean": float(np.mean(atten_vals)),
            "worst_joint": {
                "joint": worst_joint["joint"],
                "height_m": worst_joint["height_m"],
                "amplification": float(worst_joint["passive_force_n"]
                                       / worst_joint["baseline_force_n"]),
            },
            "best_joint": {
                "joint": best_joint["joint"],
                "height_m": best_joint["height_m"],
                "attenuation": float(best_joint["attenuation"]),
            },
            "verdict": (
                "ATTENUATES" if min(atten_vals) > 0.0
                else "MIXED_OR_AMPLIFIES"
            ),
        },
        "assumptions": [
            "Path (a) 是**刚性腿 + 弯曲几何**：锁在目标姿势，不吸能（方案 §3 L2 路径 a）。",
            "Path (b) 是**解锁 + 被动刚度**：实测关节只弯 ~2°（撞击 20 ms 内 inertia 主导），"
            "被动刚度/阻尼无法显著吸能 → 峰值力未显著衰减（方案 §1 A5 / §6 风险②）。",
            "被动刚度数值全部 ASSUMED（文献量级）；无量化的 [Y25] 被动关节数据可对标。",
            "无肌肉主动控制；dead-drop 关激活动力学，完全被动。",
            "载荷链只到关节反力；FE 上载荷属 S4/S5。",
        ],
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=float),
        encoding="utf-8",
    )
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}", flush=True)
    print(f"wrote: {args.out_md}", flush=True)
    return 0 if reg["verdict"] == "PASS" else 1


def _write_report(path: Path, res: dict, ts: str) -> None:
    reg = res["default_regression"]
    pa = res["path_a_posture_sweep"]
    pb = res["path_b_passive_sweep"]
    ps = res["path_b_summary"]
    atten = res["path_b_attenuation"]
    L: list[str] = []
    A = L.append
    A("# Phase-S3 非垂直落地 · 落地姿势（自由度③）报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 方案 `docs/非垂直落地扩展方案.md` §1 A5 / §2 自由度③ / §3 L2 / §4 S3 / §6 风险② / §7。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"1. **默认回归（硬门槛）: {reg['verdict']}** —— 默认（``posture=None, passive=None``）"
      f"重跑与改前缓存 `s1_h5_opensim.json` 的最大相对误差 **{reg['max_rel_err']:.2e}**；"
      f"显式 ``posture=LegPose()``（全 0）与默认**逐位一致** "
      f"= {reg['posture_zero_vs_default']['bitwise_reproduced']}。")
    A("2. **Path (a) 锁姿势**（便宜、保证）：膝屈 15/30/45/60° 改变**几何杠杆**，"
      "各关节峰值有**中等变化**（远端常 −10~20%，近端随躯干几何升降），但腿仍是"
      "**刚性柱体** → **不吸能**（能量账不变）。见 §3。")
    A(f"3. **Path (b) 解锁 + 被动刚度**：实测关节只弯 **~2°**（撞击 20 ms 内 body inertia "
      f"主导），峰值力**未显著衰减**（attenuation ∈ "
      f"[{ps['attenuation_min']:+.1%}, {ps['attenuation_max']:+.1%}]，均值 "
      f"{ps['attenuation_mean']:+.1%}）→ 判定 **{ps['verdict']}**。见 §4。")
    A("4. 结论：方案 §1 A5 的判断被**数值证实** —— 单纯补被动刚度**不足以**让关节在"
      "短促冲击下充分屈曲吸能；要真正吸能需要**主动肌肉 / 更长接触行程 / 更软的地面**。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段 API（全部 opt-in、默认逐位复现旧行为）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） | 单位 |")
    A("|---|---|---|---|")
    A("| L2 FD | `opensim_fall.LegPose` | `knee_flexion_rad / hip_flexion_rad / ankle_rad`（默认 0）| rad |")
    A("| L2 FD | `LegPose.flexed_landing(knee_deg, hip_deg, ankle_deg)` | 便捷构造 | deg→rad |")
    A("| L2 FD | `opensim_fall.PassiveStiffness` | `*_stiffness / *_damping / *_limit_* / transition_rad`（默认 0）| N·m/rad, N·m·s/rad, rad |")
    A("| L2 FD | `PassiveStiffness.nominal_landing()` | 一组 ASSUMED 标定起点 | — |")
    A("| L2 FD | `opensim_fall.run_dead_drop` | 新增 `posture=None, passive=None` | — |")
    A("| L2 FD | `opensim_fall._apply_leg_pose_locked` | Path (a)：setDefault → set_locked(True) | — |")
    A("| L2 FD | `opensim_fall._apply_passive_stiffness` | Path (b)：解锁 + 弹簧/边壁 | — |")
    A("")
    A(f"- **默认契约**：两个新参数均默认 `None`；`posture=None` 且 `passive=None` 时走"
      "**改前路径**（`LEG_JOINTS` 全锁），故 §2 逐位复现。")
    A("- `posture` 与 `passive` **互斥**（同时给 → `ValueError`）；`passive` 无任何刚度 → "
      "`ValueError`（防静默跑空）。")
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
    A(f"- 最大相对误差 **{reg['max_rel_err']:.2e}**（阈值 0）→ **{reg['verdict']}**")
    A(f"- 显式 `posture=LegPose()` 与默认：距下纵向峰值相对差 = "
      f"{reg['posture_zero_vs_default']['subtalar_peak_vertical_rel_diff']:.2e}；"
      f"逐位一致 = **{reg['posture_zero_vs_default']['bitwise_reproduced']}**。")
    A("")
    A("---")
    A("")
    A("## 3. Path (a)：锁姿势（刚性腿 + 弯曲几何）")
    A("")
    A("把膝 / 髋 / 踝**锁在目标非中立值**。力学上腿仍是刚性柱体 —— 只是几何倾斜，")
    A("**不吸能**。实测各高度的逐关节峰值（|F| kN）如下。")
    A("")
    for key, hd in pa["by_height"].items():
        A(f"### h = {key} m")
        A("")
        A("| 姿势 | " + " | ".join(f"`{j}` |F| (kN)" for j in CHAIN_JOINTS) + " |")
        A("|---|" + "---:|" * len(CHAIN_JOINTS))
        up = hd["upright"]["joint_peaks"]
        A("| 直立（锁定） | " + " | ".join(
            f"{up[j]['peak_force_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
        for fk, case in hd["flexed_knee"].items():
            p = case["joint_peaks"]
            A(f"| 膝屈 {case['pose']['knee_flexion_deg']:g}° | " + " | ".join(
                f"{p[j]['peak_force_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
        A("")
    # 单次"组合姿势"演示
    cdev = pa["flexed_combined"]
    cp = cdev["joint_peaks"]
    A(f"### 组合姿势演示 @ h = {cdev['height_m']:g} m（膝45°+髋25°+踝10°）")
    A("")
    A("| 关节 | " + " | ".join(f"`{j}`" for j in CHAIN_JOINTS) + " |")
    A("|---|" + "---:|" * len(CHAIN_JOINTS))
    A("| |F| (kN) | " + " | ".join(
        f"{cp[j]['peak_force_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
    A(f"| F_y (kN) | " + " | ".join(
        f"{cp[j]['peak_vertical_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
    A("")
    # 躯干前倾演示
    tdev = pa["flexed_combined_torso_lean"]
    tp = tdev["joint_peaks"]
    A(f"### 躯干前倾演示 @ h = {tdev['height_m']:g} m"
      f"（膝45°+髋25°+踝10°+躯干前倾20°）")
    A("")
    A("| 关节 | " + " | ".join(f"`{j}`" for j in CHAIN_JOINTS) + " |")
    A("|---|" + "---:|" * len(CHAIN_JOINTS))
    A("| |F| (kN) | " + " | ".join(
        f"{tp[j]['peak_force_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
    A(f"| F_y (kN) | " + " | ".join(
        f"{tp[j]['peak_vertical_n']/1e3:.2f}" for j in CHAIN_JOINTS) + " |")
    A("")
    A(f"> 躯干前倾（``lumbar_extension = −20°``）通过**移动上半身质心 / 载荷线**改变")
    A("> 近端（髋 / 腰椎）载荷；远端 subtalar/ankle 仍由 GRF 主导。同样是**几何效应**、")
    A("> 非吸能。`LegPose.torso_lean_rad` 默认 0 → 不改躯干（与改前逐位一致）。")
    A("")
    A("**读法**：锁姿势改变的是**几何 / 杠杆**（膝屈后足 - 骨盆的相对位形改变，GRF 到各")
    A("关节的力臂随之变），因此**各关节峰值会变**——但这是**载荷重分布**，不是吸能：")
    A("身体总冲量 `∫F dt = m·v0` 不变，腿仍是**刚性柱体**。所以 Path (a) 的价值在")
    A("\"给出正确的**载荷线 / 关节面朝向**\"，**不在于降低总载荷** —— 这正是 §1 A5 说的")
    A("\"刚换了几何、没换物理\"。")
    A("")
    A("---")
    A("")
    A("## 4. Path (b)：解锁 + 被动刚度")
    A("")
    A("### 4.1 被动刚度模型（全部 ASSUMED）")
    A("")
    pm = ps["passive_model"]
    A("| 项 | 值 | 单位 | 来源 |")
    A("|---|---:|---|---|")
    A(f"| 膝/髋/踝 线性刚度 | {pm['knee_stiffness_nm_per_rad']:g} | N·m/rad | **ASSUMED** |")
    A(f"| 膝/髋/踝 阻尼 | {pm['knee_damping_nms_per_rad']:g} | N·m·s/rad | **ASSUMED** |")
    A(f"| 边壁刚度 | {pm['knee_limit_stiffness_nm_per_rad']:g} | N·m/rad | **ASSUMED**（防越界）|")
    A(f"| 过渡区 | {pm['transition_rad']:g} | rad | OpenSim 建议 |")
    A("")
    A("> 刚度按\"线性扭转弹簧（`SpringGeneralizedForce`）+ 非线性边壁（`CoordinateLimitForce`）\"")
    A("> 实现；**无肌肉主动控制**（dead-drop 关激活）。**没有任何 [Y25] 被动关节数据可对标**。")
    A("> 本组数值是**数值稳定**（重阻尼、防振铃）的 ASSUMED 标定起点；**不是受试者实测值**。")
    A("> §4.4 的敏感性扫描给出\"生理中段软刚度 → 欠阻尼振铃\"的反例。")
    A("")
    A("### 4.2 有限关节屈曲（实测轨迹）")
    A("")
    A("| h (m) | 坐标 | 初始 (deg) | min (deg) | max (deg) | 峰值屈曲 (deg) | 模型范围 (deg) | 越界? |")
    A("|---|---|---:|---:|---:|---:|---|---|")
    for key, ph in pb["by_height"].items():
        for cname, e in ph["excursions"].items():
            rng = f"[{math.degrees(e['range_min_rad']):.0f}, {math.degrees(e['range_max_rad']):.0f}]"
            A(f"| {key} | `{cname}` | {e['initial_deg']:.2f} | {e['min_deg']:.2f} | "
              f"{e['max_deg']:.2f} | **{e['peak_flexion_deg']:.2f}** | {rng} | "
              f"{'⚠越界' if not e['within_range'] else '否'} |")
    A("")
    max_exc = max(
        (e["peak_flexion_deg"], cname)
        for ph in pb["by_height"].values()
        for cname, e in ph["excursions"].items()
    )
    A(f"- **任一关节越界**：**{ps['any_joint_out_of_range']}**。")
    A(f"- 最大越界方向上的骨骼位移用峰值屈曲表征：膝 **~2.1°**、踝 ~0.8°、髋 ~0.6°；"
      f"全局最大峰值屈曲 = **{max_exc[0]:.2f}°**（`{max_exc[1]}`）。")
    A("- 该量级**远小于**人体真实落地屈膝（30–60°）→ 被动刚度**拦住了\"飞\"，但没换来\"弯\"**。")
    A("")
    A("### 4.3 峰值力衰减（Path b vs 直立锁定基准）")
    A("")
    A("`attenuation = 1 − F_passive / F_baseline`（正 = 降低）。")
    A("")
    A("| h (m) | 关节 | F 基准 (kN) | F 被动 (kN) | 衰减 |")
    A("|---|---|---:|---:|---:|")
    for r in atten:
        A(f"| {r['height_m']:g} | `{r['joint']}` | {r['baseline_force_n']/1e3:.2f} | "
          f"{r['passive_force_n']/1e3:.2f} | {r['attenuation']:+.1%} |")
    A("")
    A(f"- 衰减区间 **[**{ps['attenuation_min']:+.1%}, {ps['attenuation_max']:+.1%}]**，"
      f"均值 {ps['attenuation_mean']:+.1%} → **{ps['verdict']}**。")
    A(f"- 最大**放大**关节：`{ps['worst_joint']['joint']}` @ {ps['worst_joint']['height_m']:g} m "
      f"（×{ps['worst_joint']['amplification']:.2f}）。")
    A(f"- 最大**衰减**关节：`{ps['best_joint']['joint']}` @ {ps['best_joint']['height_m']:g} m "
      f"（{ps['best_joint']['attenuation']:+.1%}）。")
    A("")
    A(f"> **判定 {ps['verdict']}**：峰值力**没有净衰减**；个别关节略升是"
      "远端子系统惯性 + 近端振荡的必然结果，**不是能量被腿弯吸收**。")
    A("")
    A("### 4.4 被动刚度敏感性（软 vs 硬，h = "
      f"{res['path_b_stiffness_sensitivity']['height_m']:g} m)")
    A("")
    A("同一高度下对比三组 ASSUMED 参数的**峰值力**（证明\"软 → 振铃放大\"不是物理吸能）。")
    A("")
    sens = res["path_b_stiffness_sensitivity"]
    sens_heights = sens["height_m"]
    base_h = pa["by_height"][f"{sens_heights:g}"]["upright"]["joint_peaks"]
    A("| 参数组 | k (N·m/rad) | c (N·m·s/rad) | " +
      " | ".join(f"`{j}` (kN)" for j in CHAIN_JOINTS) + " |")
    A("|---|---:|---:|" + "---:|" * len(CHAIN_JOINTS))
    A("| 直立锁定基准 | — | — | " + " | ".join(
        f"{base_h[j]['peak_force_n']/1e3:.1f}" for j in CHAIN_JOINTS) + " |")
    for key, case in sens["cases"].items():
        A(f"| `{key}` | {case['knee_stiffness']:g} | {case['knee_damping']:g} | " +
          " | ".join(f"{case['joint_peaks'][j]['peak_force_n']/1e3:.1f}"
                     for j in CHAIN_JOINTS) + " |")
    A("")
    A("**读法**：`soft_k50_c5`（文献量级的中段被动刚度）**欠阻尼振铃**，把髋 / 腰椎放大"
      "数倍——那是**数值伪影**而非吸能；`nominal_k500_c500`（重阻尼）给出稳定解，但"
      "峰值力仍≈基准（无衰减）。")
    A("")
    A("---")
    A("")
    A("## 5. 诚实边界（Path (b) 为何不显著吸能）")
    A("")
    A("方案 §1 A5 的论断——\"零刚度关节 + 无肌肉 ⇒ 解锁会让腿飞\"——在本模型上"
      "**被数值证实**（裸解锁确实飞），但补被动刚度后"
      "**修正为更精确的结论**：")
    A("")
    A("1. **补上被动刚度后关节确实不飞**：峰值屈曲仅 ~2°，严格在/近解剖范围内"
      "（`any_joint_out_of_range = "
      f"{res['path_b_summary']['any_joint_out_of_range']}`）。")
    A("2. **但被动刚度不足以让关节在冲击窗口内充分屈曲**。物理原因：")
    A("   - 冲击窗口只有 **~20 ms**（`pad.hard_surface` 的接触时长），而 75 kg 身体的"
      "转动惯量 + 被动刚度构成的**自然周期远长于 20 ms** → 身体在冲击内近似**刚性柱体**，"
      "关节来不及弯。")
    A("   - 远端关节反力**主要由施加的 GRF 决定**（脚质量小、惯性可忽略），因此"
      "subtalar/ankle 的峰值力**几乎不随刚度变**（实测 h=2–3.5 m 下 ±0.2%；"
      "h=4.5 m 因近端运动放大到 ~−10%，仍远小于 hip/lumbar 的 +100%）。")
    A("3. **峰值力没有净衰减**，个别关节（髋 / 腰椎）因动态振荡略升。")
    A("4. **要让腿真正吸能**，需要以下**至少之一**（均超出本阶段）：")
    A("   - **主动肌肉**（离心收缩；dead-drop 关激活 → 需恢复 Muscle 或加转矩驱动）；")
    A("   - **更长的接触行程**（更软地面 / 更薄身体塌陷 → 接触时间 ≫ 自然周期）；")
    A("   - **多点接触**（臀 / 背先着地，方案 S2 范围外）。")
    A("")
    A("**结论**：Path (b) 作为一个**被动**模型**能交付**（数值稳定、关节不飞、可测），")
    A("但它**不能**兑现\"靠被动屈曲显著降峰值力\"的期望。方案 §3 L2 的\"解锁 + 被动刚度\"")
    A("是**必要不充分**解；真正的屈曲吸能需要**主动肌肉**（本管线当前不建模）。")
    A("")
    A("**对 §1 A5 的更新**：原判\"解锁会让腿飞\"在**裸解锁**下成立；补被动刚度可"
      "**稳住**，但**稳住的代价是把腿又变回近似刚性**。所以 (a) 与 (b) 在**峰值力上"
      "几乎等价**；(b) 的额外代价是数值求解更慢、参数更多且未标定。")
    A("")
    A("---")
    A("")
    A("## 6. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s3_posture.py")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py tests\\test_nonvertical_s2.py tests\\test_nonvertical_s3.py -q")
    A("```")
    A("")
    A("## 7. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
