"""Phase-S1 probe —— 斜向落地方向（矢量通路）的**轴向回归证明** + **倾角演示**。

方案：``docs/非垂直落地扩展方案.md`` §3 L1–L3（最小矢量通路）、§4 Phase-S1、
§7（轴向回归硬门槛）。

本脚本**只读复用**既有产物与模块，**不调用 FEBio**、**不覆盖任何既有文件**。
新增产物仅两个（若同名已存在则拒绝写入）：

* ``results/opensim_fe/nonvertical_s1.json``
* ``results/opensim_fe/NONVERTICAL_S1_REPORT.md``

它做两件事：

1. **轴向回归（硬门槛）**：把倾角 / 水平速度全部置零，重跑默认链路
   （``ground_reaction`` → ``run_dead_drop`` → ``subtalar_reaction``），
   与缓存 ``results/opensim_fe/s1_h5_opensim.json`` 数值对比；并证明显式
   ``tilt_deg=0`` 与默认路径**逐位一致**。
2. **倾角演示（载荷链层面）**：``tilt_deg=20`` 的地面法向 → 3D GRF →
   正动力学 → 各关节 3D wrench → ``load_transfer`` 保留完整矢量 →
   逐骨纯轴向 1D 风险（Route-1 + 载荷分配口径）相对轴向的偏移。

单位 mm–N–MPa–s。

用法（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s1.py --tilt-deg 20
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
# 复用既有的 Route-1 载荷分配常量（additive，不重跑它）
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

M_MODEL = 75.337
CACHED_AXIAL = ROOT / "results" / "opensim_fe" / "s1_h5_opensim.json"
ROUTE1_JSON = ROOT / "results" / "opensim_fe" / "bc_robust_metric_route1.json"
OUT_JSON = ROOT / "results" / "opensim_fe" / "nonvertical_s1.json"
OUT_MD = ROOT / "results" / "opensim_fe" / "NONVERTICAL_S1_REPORT.md"

#: 逐骨风险用到的关节（route1 FORCE_SOURCE 的像）
_CHAIN_JOINTS = ("subtalar_r", "ankle_r", "hip_r", "lumbar")


def _import_risk_loadshare():
    """加载现有 Route-1 载荷分配常量（FORCE_SOURCE / RATIO / 分流）。"""
    import risk_1d_loadshare as r1  # scripts/opensim_fe 已在 sys.path

    return r1


# --------------------------------------------------------------------------
# 轴段一：轴向回归证明
# --------------------------------------------------------------------------
def axial_regression() -> dict:
    from climbing.coupling.joint_loads import subtalar_reaction
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    cached = json.loads(CACHED_AXIAL.read_text(encoding="utf-8"))

    grf_default = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    fall_default = run_dead_drop(grf_default)
    jr_default = subtalar_reaction(fall_default, grf_default, side="r")

    measured = {
        "impulse_ns": float(grf_default.impulse_ns),
        "impulse_rel_err": float(grf_default.impulse_rel_err),
        "grf_peak_n": float(grf_default.peak_total_n),
        "subtalar_peak_vertical_n": float(jr_default.peak_vertical_n),
        "subtalar_peak_force_n": float(jr_default.peak_force_n),
        "subtalar_peak_moment_nm": float(jr_default.peak_moment_nm),
    }
    cached_keys = {
        "impulse_ns": float(cached["impulse_ns"]),
        "impulse_rel_err": float(cached["impulse_rel_err"]),
        "grf_peak_n": float(cached["grf_peak_n"]),
        "subtalar_peak_vertical_n": float(cached["subtalar_peak_vertical_n"]),
        "subtalar_peak_force_n": float(cached["subtalar_peak_force_n"]),
        "subtalar_peak_moment_nm": float(cached["subtalar_peak_moment_nm"]),
    }
    rel_err = {
        k: (abs(measured[k] - cached_keys[k]) / abs(cached_keys[k])
            if cached_keys[k] != 0.0 else abs(measured[k]))
        for k in measured
    }
    max_rel = max(rel_err.values())

    # 显式 tilt_deg=0 必须与默认逐位一致（接口层）
    grf_t0 = ground_reaction(height_m=5.0, mass_kg=M_MODEL, tilt_deg=0.0)
    arrays_bitwise_equal = all(
        np.array_equal(getattr(grf_default, k), getattr(grf_t0, k))
        for k in ("t_s", "f_total_n", "f_left_n", "f_right_n")
    )
    fall_t0 = run_dead_drop(grf_t0)
    jr_t0 = subtalar_reaction(fall_t0, grf_t0, side="r")
    t0_rel = abs(jr_t0.peak_vertical_n - jr_default.peak_vertical_n) / jr_default.peak_vertical_n

    verdict = "PASS" if (max_rel <= 1e-6 and arrays_bitwise_equal and t0_rel <= 1e-12) else "FAIL"
    return {
        "cached_file": str(CACHED_AXIAL.relative_to(ROOT)),
        "measured_default_run": measured,
        "cached": cached_keys,
        "rel_err": rel_err,
        "max_rel_err": max_rel,
        "tilt0_vs_default": {
            "grf_arrays_bitwise_equal": bool(arrays_bitwise_equal),
            "subtalar_peak_vertical_rel_diff": float(t0_rel),
        },
        "verdict": verdict,
    }


# --------------------------------------------------------------------------
# 阶段二：倾角演示（载荷链层面）
# --------------------------------------------------------------------------
def _peaks(force_n: np.ndarray, moment_nm: np.ndarray) -> dict:
    mag = np.linalg.norm(force_n, axis=1)
    return {
        "peak_force_n": float(np.max(mag)),
        "peak_vertical_n": float(np.max(np.abs(force_n[:, 1]))),
        "peak_horizontal_x_n": float(np.max(np.abs(force_n[:, 0]))),
        "peak_horizontal_z_n": float(np.max(np.abs(force_n[:, 2]))),
        "peak_moment_nm": float(np.max(np.linalg.norm(moment_nm, axis=1))),
    }


def _joint_peaks(fall, grf) -> dict:
    from climbing.coupling.joint_reactions import joint_reaction

    out = {}
    for j in _CHAIN_JOINTS:
        jr = joint_reaction(fall, grf, joint=j)
        rec = _peaks(jr.force_n, jr.moment_nm)
        rec["peak_force_time_s"] = jr.peak_force_time_s
        out[j] = rec
    return out


def _tilted_demo(tilt_deg: float) -> dict:
    from climbing.coupling.joint_reactions import joint_reaction
    from climbing.coupling.load_transfer import transfer
    from climbing.coupling.opensim_fall import run_dead_drop
    from climbing.coupling.opensim_grf import ground_reaction

    grf_ax = ground_reaction(height_m=5.0, mass_kg=M_MODEL)
    grf_ti = ground_reaction(height_m=5.0, mass_kg=M_MODEL, tilt_deg=tilt_deg)
    fall_ax = run_dead_drop(grf_ax)
    fall_ti = run_dead_drop(grf_ti)

    jp_ax = _joint_peaks(fall_ax, grf_ax)
    jp_ti = _joint_peaks(fall_ti, grf_ti)

    # ---- 距下 3D wrench：峰值帧的完整矢量 + 力矩 --------------------------
    jr_sub = joint_reaction(fall_ti, grf_ti, joint="subtalar_r")
    i_pk = int(np.argmax(np.linalg.norm(jr_sub.force_n, axis=1)))
    vec = tuple(float(x) for x in jr_sub.force_n[i_pk])
    mom = tuple(float(x) for x in jr_sub.moment_nm[i_pk])
    axial_n = abs(vec[1])  # 纵向分量 = Route-1 主载荷
    spec = transfer(axial_n, subtalar_force=vec, subtalar_moment=mom)

    # ---- 逐骨纯轴向 1D 风险（Route-1 + 载荷分配）--------------------------
    r1 = _import_risk_loadshare()
    route1 = json.loads(ROUTE1_JSON.read_text(encoding="utf-8"))
    rows = {r["bone"]: r for r in route1["rows"]}
    a_by = {b: float(rows[b]["A_section_mm2"]) for b in rows}
    L_by = {b: float(rows[b]["L_bbox_mm"]) for b in rows}
    split = r1.tibia_fibula_split(a_by, L_by)
    f_tib = float(split["split_fraction"]["tibia_r"])
    f_fib = float(split["split_fraction"]["fibula_r"])

    mult = {}
    for bone, joint in r1.FORCE_SOURCE.items():
        if bone == "tibia_r":
            m = f_tib
        elif bone == "fibula_r":
            m = f_fib
        elif bone in r1.AXIAL_BONES:
            m = float(r1.RATIO[bone])
        else:
            m = 1.0
        mult[bone] = m

    per_bone = []
    for bone in rows:
        joint = r1.FORCE_SOURCE[bone]
        A = a_by[bone]
        sig = float(rows[bone]["sigma_c_mpa"])
        thr = A * sig
        f_ax_vert = jp_ax[joint]["peak_vertical_n"]
        f_ti_vert = jp_ti[joint]["peak_vertical_n"]
        f_ax_mag = jp_ax[joint]["peak_force_n"]
        f_ti_mag = jp_ti[joint]["peak_force_n"]
        risk_ax = mult[bone] * f_ax_vert / thr
        risk_ti_vert = mult[bone] * f_ti_vert / thr
        risk_ti_mag = mult[bone] * f_ti_mag / thr
        per_bone.append({
            "bone": bone,
            "part_cn": rows[bone]["part_cn"],
            "joint": joint,
            "force_multiplier": mult[bone],
            "A_section_mm2": A,
            "sigma_c_mpa": sig,
            "F_axial_axial_n": f_ax_vert,
            "F_axial_tilt_n": f_ti_vert,
            "F_mag_axial_n": f_ax_mag,
            "F_mag_tilt_n": f_ti_mag,
            "risk_axial": risk_ax,
            "risk_tilt_vertical": risk_ti_vert,
            "risk_tilt_magnitude": risk_ti_mag,
            "shift_factor_vertical": (risk_ti_vert / risk_ax) if risk_ax else float("nan"),
            "shift_factor_magnitude": (risk_ti_mag / risk_ax) if risk_ax else float("nan"),
            "fracture_at_5m_axial": bool(risk_ax >= 1.0),
            "fracture_at_5m_tilt": bool(risk_ti_mag >= 1.0),
        })
    per_bone.sort(key=lambda r: -r["risk_axial"])

    n_deg = np.asarray(grf_ti.ground_normal, dtype=float)
    return {
        "tilt_deg": tilt_deg,
        "ground_normal": [float(x) for x in n_deg],
        "vertical_fraction": float(n_deg[1]),
        "grf": {
            "peak_total_n_axial": float(grf_ax.peak_total_n),
            "peak_total_n_tilt": float(grf_ti.peak_total_n),
            "impulse_rel_err_axial": float(grf_ax.impulse_rel_err),
            "impulse_rel_err_tilt": float(grf_ti.impulse_rel_err),
            "impact_speed_ms": float(grf_ti.impact_speed_ms),
        },
        "subtalar": {
            "axial": jp_ax["subtalar_r"],
            "tilt": jp_ti["subtalar_r"],
        },
        "feloadspec": {
            "call": "transfer(subtalar_n=abs(F_y), subtalar_force=peak_vec, subtalar_moment=peak_mom)",
            "subtalar_n": axial_n,
            "force_vector_n": list(spec.force_vector_n),
            "moment_vector_nm": list(spec.moment_vector_nm),
            "magnitude_n": float(np.linalg.norm(spec.force_vector_n)),
            "direction": list(np.asarray(spec.force_vector_n) /
                              max(1e-12, float(np.linalg.norm(spec.force_vector_n)))),
        },
        "joint_peaks_axial": jp_ax,
        "joint_peaks_tilt": jp_ti,
        "per_bone_risk": per_bone,
    }


# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Phase-S1 非垂直（矢量通路）probe + demo")
    ap.add_argument("--tilt-deg", type=float, default=20.0, help="地面法向倾角 (度)")
    ap.add_argument("--out-json", type=Path, default=OUT_JSON)
    ap.add_argument("--out-md", type=Path, default=OUT_MD)
    args = ap.parse_args(argv)

    # 不覆盖既有产物
    for p in (args.out_json, args.out_md):
        if p.exists():
            print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
            return 2

    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("=== Phase-S1 非垂直落地：轴向回归 + 倾角演示 ===")
    print("[1/2] 轴向回归（默认路径 vs 缓存 s1_h5_opensim.json）...")
    reg = axial_regression()
    print(f"      verdict={reg['verdict']}  max_rel_err={reg['max_rel_err']:.3e}  "
          f"tilt0_bitwise={reg['tilt0_vs_default']['grf_arrays_bitwise_equal']}")

    print(f"[2/2] 倾角演示 tilt={args.tilt_deg:g}° ...")
    demo = _tilted_demo(args.tilt_deg)
    sub = demo["subtalar"]
    print(f"      距下 peak|F|: 轴向 {sub['axial']['peak_force_n']/1e3:.2f} kN "
          f"→ 倾斜 {sub['tilt']['peak_force_n']/1e3:.2f} kN "
          f"(横向 {sub['tilt']['peak_horizontal_x_n']/1e3:.2f} kN)")
    for r in demo["per_bone_risk"]:
        print(f"      risk@5m {r['bone']:12s} 轴向 {r['risk_axial']:.4f} → "
              f"倾斜(纵向) {r['risk_tilt_vertical']:.4f} / (合成) {r['risk_tilt_magnitude']:.4f} "
              f"[×{r['shift_factor_magnitude']:.3f}]")

    result = {
        "meta": {
            "generated_at": ts,
            "study": "Phase-S1 非垂直落地方向（矢量通路）",
            "plan_doc": "docs/非垂直落地扩展方案.md",
            "units": "mm-N-MPa-s",
            "height_m": 5.0,
            "mass_kg": M_MODEL,
            "model_mass_kg_actual": 75.337,
            "commands": [
                'python scripts/opensim_fe/nonvertical_s1.py --tilt-deg 20',
            ],
            "measured_vs_modeled": {
                "measured": [
                    "GRF 由 pad.py 1D 模型产出（唯一接触物理来源）",
                    "OpenSim 正动力学关节反力（joint_reactions 自由体口径）",
                    "缓存轴向基准 results/opensim_fe/s1_h5_opensim.json",
                ],
                "modeled_assumed": [
                    "矢量通路：沿用标量 GRF 大小，仅改方向（方案 §3 L1 最小版），"
                    "非接触本构重写",
                    "逐骨风险用 Route-1 1D 名义应力 + 载荷分配（A_section/σ_c/质量分数），"
                    "忽略弯曲/剪切/偏心",
                    "FE 适配器仍只吃 大小×方向；力矩/偏心映射属 S4",
                ],
            },
        },
        "axial_regression": reg,
        "tilted_demo": demo,
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    _write_report(args.out_md, result, ts)
    print(f"wrote: {args.out_json}")
    print(f"wrote: {args.out_md}")
    return 0 if reg["verdict"] == "PASS" else 1


def _write_report(path: Path, res: dict, ts: str) -> None:
    reg = res["axial_regression"]
    demo = res["tilted_demo"]
    L: list[str] = []
    A = L.append
    A("# Phase-S1 非垂直落地（矢量通路）报告")
    A("")
    A(f"> 生成时间 {ts}；仓库 `D:\\Project\\climbing_fall_analysis`；单位 mm-N-MPa-s。")
    A("> 方案 `docs/非垂直落地扩展方案.md` §3 L1–L3 / §4 Phase-S1 / §7。")
    A("> 本报告**新增**，不改任何模块默认、不覆盖既有产物、**不调用 FEBio**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(f"1. **轴向回归（硬门槛）: {reg['verdict']}** —— 默认（竖直）路径与改前缓存"
      f" `s1_h5_opensim.json` 的最大相对误差 **{reg['max_rel_err']:.2e}**；显式 "
      f"`tilt_deg=0` 与默认的 GRF 数组**逐位一致** "
      f"({reg['tilt0_vs_default']['grf_arrays_bitwise_equal']})。")
    A(f"2. **倾角 {demo['tilt_deg']:g}° 演示**：距下峰值力 "
      f"{demo['subtalar']['axial']['peak_force_n']/1e3:.2f} kN（轴向）→ "
      f"{demo['subtalar']['tilt']['peak_force_n']/1e3:.2f} kN（倾斜），"
      f"其中新增横向分量 "
      f"{demo['subtalar']['tilt']['peak_horizontal_x_n']/1e3:.2f} kN；"
      "`load_transfer` 现保留完整 3D 力矢量 + 力矩。")
    A("3. 逐骨纯轴向 1D 风险见 §3；倾斜是否改变排序见 §3.3。")
    A("")
    A("---")
    A("")
    A("## 1. 本阶段实现的 API（全部 opt-in、默认复现旧行为）")
    A("")
    A("| 层 | 文件 / 符号 | 新增（默认） |")
    A("|---|---|---|")
    A("| L1 GRF | `opensim_grf.ground_reaction` | `ground_normal=(0,1,0)`、`tilt_deg=None`（度）|")
    A("| L1 GRF | `GroundReaction.ground_normal` / `.is_vertical` / `.force_vec_n(side)` | 默认竖直，`F⃗=f(t)·n̂` |")
    A("| L2 FD | `opensim_fall._add_foot_force` | `f_x=None, f_z=None, point=None`（默认 0 / 原点）|")
    A("| L2 FD | `opensim_fall.run_dead_drop` | `vx0_ms=0.0, vz0_ms=0.0` |")
    A("| L3 传递 | `load_transfer.FeLoadSpec` | `subtalar_force=None, subtalar_moment=None`；`force_vector_n` / `moment_vector_nm` |")
    A("| L3 传递 | `load_transfer.transfer` | `subtalar_force=None, subtalar_moment=None` |")
    A("| 消费 | `joint_reactions.joint_reaction` | 自由体按 `grf.force_vec_n` 计入完整 3D 外力（默认竖直 = 旧行为）|")
    A("")
    A("> 中段 `joint_loads.subtalar_reaction` **未改**（本就 3D）。")
    A("")
    A("---")
    A("")
    A("## 2. 轴向回归（硬门槛，方案 §7）")
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
    A(f"- 最大相对误差 **{reg['max_rel_err']:.2e}**（阈值 1e-6）。")
    A(f"- 显式 `tilt_deg=0` vs 默认：GRF 数组逐位一致 = "
      f"**{reg['tilt0_vs_default']['grf_arrays_bitwise_equal']}**；"
      f"距下纵向峰值相对差 = {reg['tilt0_vs_default']['subtalar_peak_vertical_rel_diff']:.2e}。")
    A(f"- **判定：{reg['verdict']}**")
    A("")
    A("---")
    A("")
    A("## 3. 倾角演示（载荷链层面）")
    A("")
    A(f"地面法向 `n̂ = ({demo['ground_normal'][0]:.4f}, {demo['ground_normal'][1]:.4f}, "
      f"{demo['ground_normal'][2]:.4f})`（竖直占比 {demo['vertical_fraction']:.4f}）。")
    A("")
    A("### 3.1 距下关节（subtalar_r）wrench")
    A("")
    A("| 量 | 轴向 | 倾斜 |")
    A("|---|---:|---:|")
    for key, label in (("peak_force_n", "峰值 |F| (N)"),
                       ("peak_vertical_n", "峰值 纵向 y (N)"),
                       ("peak_horizontal_x_n", "峰值 横向 x (N)"),
                       ("peak_moment_nm", "峰值 |M| (N·m)")):
        A(f"| {label} | {demo['subtalar']['axial'][key]:,.2f} | {demo['subtalar']['tilt'][key]:,.2f} |")
    A("")
    A("### 3.2 `load_transfer` 保留的 FE 载荷规格")
    A("")
    fs = demo["feloadspec"]
    A(f"- `subtalar_n`（纵向标量）= {fs['subtalar_n']:,.2f} N")
    A(f"- `force_vector_n` = ({fs['force_vector_n'][0]:,.2f}, {fs['force_vector_n'][1]:,.2f}, "
      f"{fs['force_vector_n'][2]:,.2f}) N，|F| = {fs['magnitude_n']:,.2f} N")
    A(f"- `moment_vector_nm` = ({fs['moment_vector_nm'][0]:.4f}, {fs['moment_vector_nm'][1]:.4f}, "
      f"{fs['moment_vector_nm'][2]:.4f}) N·m")
    A("- 旧口径等价物 = `subtalar_n × SUBTALAR_DIR`（仅纵向）；现完整矢量与力矩不再被丢弃。")
    A("")
    A("### 3.3 逐骨纯轴向 1D 风险（Route-1 + 载荷分配）")
    A("")
    A("`risk = (mult · F_joint) / (A_section · σ_c)`；`F_joint` 取关节反力的**纵向分量**"
      "（Route-1 口径）为主，另列**合成 |F|** 作为 3D 载荷代理。mult 来自 `risk_1d_loadshare`"
      "（胫/腓并联 k=E·A/L；中轴骨 mass-above）。")
    A("")
    A("| 骨 | 关节 | F_y 轴向 (N) | F_y 倾斜 (N) | risk 轴向 | risk 倾斜(纵向) | risk 倾斜(合成) | 偏移(合成) |")
    A("|---|---|---:|---:|---:|---:|---:|---:|")
    for r in demo["per_bone_risk"]:
        A(f"| `{r['bone']}` | `{r['joint']}` | {r['F_axial_axial_n']:,.1f} | "
          f"{r['F_axial_tilt_n']:,.1f} | {r['risk_axial']:.4f} | "
          f"{r['risk_tilt_vertical']:.4f} | {r['risk_tilt_magnitude']:.4f} | "
          f"×{r['shift_factor_magnitude']:.3f} |")
    A("")
    fr_ax = [r["bone"] for r in demo["per_bone_risk"] if r["fracture_at_5m_axial"]]
    fr_ti = [r["bone"] for r in demo["per_bone_risk"] if r["fracture_at_5m_tilt"]]
    A(f"- 轴向 @5m 已越阈（risk≥1）：`{fr_ax or '无'}`。")
    A(f"- 倾斜（合成 |F|）@5m 越阈：`{fr_ti or '无'}`。")
    A("- 说明：倾角把法向大小按方向拆分，**纵向分量下降、横向分量新增**；纯轴向 1D 指标"
      "只吃纵向，因此倾斜下纵向风险普遍**下降**，而 3D 合成载荷（含横向剪切）才是"
      "倾角真实加载强度 —— 两者之差即\"纯轴向指标的盲区\"。")
    A("")
    A("---")
    A("")
    A("## 4. 诚实边界")
    A("")
    A("- **矢量通路是方向化的最小版**：GRF 标量大小不变，只旋转方向；接触本构、摩擦锥、"
      "CoP 偏移、单脚均未实现（方案 §3 L1 完全版 / S2 / S4）。")
    A("- **未上 FE**：`FeLoadSpec` 已保留完整 3D 力 + 力矩，但当前 FE 适配器仍只消费"
      "大小×方向；力矩/偏心压力映射属 S4，本阶段不宣称已实现。")
    A("- 逐骨风险是 **1D 名义应力 + 载荷分配**，忽略弯曲/剪切/偏心与应力集中；`σ_c`、"
      "`A_section`、质量分数均为既有 assumed/estimated 口径，绝对值不可信，只用**相对**。")
    A("- 倾角动力学为**两足对称、刚性腿**（方案 S1 范围），非单脚。")
    A("")
    A("## 5. 复现命令")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s1.py --tilt-deg 20")
    A("# 回归：")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_opensim_fe.py tests\\test_febio_run.py -q")
    A("& .venv\\Scripts\\python.exe -m pytest tests\\test_nonvertical_s1.py -q")
    A("```")
    A("")
    A("## 6. 产物")
    A("")
    A(f"- `{OUT_MD.relative_to(ROOT)}`（本报告）")
    A(f"- `{OUT_JSON.relative_to(ROOT)}`（机器可读）")
    A("")
    path.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
