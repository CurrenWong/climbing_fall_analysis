"""【刚】工况 V 靶补测探针 —— RIGID_VTARGET 第 2-4 步（2026-10-08 标定轮）。

对 V3 / V6 / V8 / V9 在**真实求解器轨迹**上补测【刚】（on_pad=False）
列，并同时复跑【垫】（on_pad=True）作同参数对照，产出
``results/rigid_vtarget_probe.json``（各 part 增量合并）。

设计约束
--------
* V3/V8 用 ``posture="feet-first-stiff"``、``mass=80``（论文 §4.5 垫工况
  37.0/33.8/31.2/28.8/26.6 ms 在 mass=80 下**逐位复现**）。
* V6/V9 用 ``posture="controlled-drop"``、``mass=80``（逐字镜像
  ``tests/test_pad2d.py`` 的 ``_ankle_run``，保证与既有【垫】✅ 结果同口径）。
* 判据阈值不自设：V3=0.0125 s、V8@20=0.005 s、V6/V9 期望表照抄测试。

2026-10-08 标定轮变更
--------------------
* **求解器口径分层**：【垫】格保持默认紧容差（论文 §4.5 逐位口径）；
  【刚】格用放松容差 ``SOLVER_RIGID``（实测与紧容差 Fpk 差 <0.004%、
  fracture/util/mode 全同，见 ``temp/time_relaxed.py``），单跑 ~4× 提速。
* **并行**：所有 part 的单元格**合并进一个进程池**（``climbing.batch.
  parallel_map``，spawn 开销只付一次），worker 内按 ``kind`` 分发。
* **L0b 峰时字段**：每行补 ``tF_ms``（GRF 首峰）、``t_torso_ms``（躯干加速度
  峰），**一律限制在冲击窗内**（镜像 ``head_branch`` / ``bone.py:301`` 的
  窗口守卫：刚性地面无阻尼 ⇒ settle 事件不触发，全程轨迹的 argmax 会被
  窗后振铃污染，实测出过 966/1325 ms 的假峰）。文献对比用局部峰时，
  **不再用 settle 窗**（settle 窗含振铃尾巴，与文献"局部组织峰时
  1.8–11.9 ms"不可比）。
* **L4 V8 观测量**：新增 ``t_moment_ms``（2D 踝内翻力矩峰，= T22 Table 3
  的胫骨/韧带局部峰观测量）；V8 判据从头峰改为该列 —— 头峰归 V3 独占。
  V8a=随高度单调提前（非增），V8b=@20 m ≤5 ms。头峰列保留（``tpeak_ms_k*``）。

用法::

    python scripts/paper/rigid_vtarget_probe.py --part head   # V3/V8
    python scripts/paper/rigid_vtarget_probe.py --part v6     # V6 网格 ×2 工况
    python scripts/paper/rigid_vtarget_probe.py --part v9     # V9 内翻/外翻 ×2 工况
    python scripts/paper/rigid_vtarget_probe.py --part all
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "src"))

from climbing.batch import parallel_map  # noqa: E402
from climbing.head_branch import (  # noqa: E402
    V3_THRESHOLD_S,
    V8_THRESHOLD_AT_20M_S,
    HeadBranchParams,
    simulate_head_branch,
)
from climbing.pad import simulate_boulder_fall  # noqa: E402
from climbing.pad2d import Posture2D, simulate_boulder_fall_2d  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "rigid_vtarget_probe.json"

# ---- 求解器口径 ----------------------------------------------------------
# 【垫】= 默认紧容差（与论文 §4.5 / 全测试套件同口径，逐位不变）。
SOLVER_PAD = dict(max_ode_step=1e-4, rtol=1e-8, atol=1e-11)
# 【刚】= 放松容差（temp/time_relaxed.py 等价验证：Fpk 差 <0.004%、
# fracture/util/mode 全同；2D 单跑 ~19s → ~5s）。
SOLVER_RIGID = dict(max_ode_step=5e-4, rtol=1e-6, atol=1e-9)

# ---- V3/V8 口径（论文 §4.5 复核：feet-first-stiff + **80 kg**）----
HEAD_POSTURE = "feet-first-stiff"
HEAD_MASS_KG = 80.0
HEAD_HEIGHTS_RIGID = [1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 20.0]
HEAD_HEIGHTS_PAD = [2.0, 3.0, 5.0, 10.0, 20.0]
K_NECKS = (200.0, 1000.0, 2500.0)   # default / stiff / 标定值（文献动态颈椎
                                    # 轴向刚度上限 2.5 kN/mm，head_branch.py:115）
# V3 主判据列：刚=标定值 2500（L1 扫描：文献上限内最优）；垫=论文 §4.5
# 发表口径 200（垫侧不标定，绝对时刻本就不可移植）。
CALIB_K = {"rigid": 2500.0, "pad": 200.0}
PAD_PAPER_MS = {2.0: 37.0, 3.0: 33.8, 5.0: 31.2, 10.0: 28.8, 20.0: 26.6}
V3_HEIGHTS = [1.0, 2.0, 3.0, 5.0, 8.0]      # 同 tests/test_head_branch.py
V8_HEIGHTS = [2.0, 3.0, 5.0, 10.0, 20.0]
# V8 观测量的 2D 口径：直立刚腿 + β0=+20°（内翻），t_moment 与 β 符号无关
V8_MOMENT_BETA0_DEG = 20.0

# ---- V6/V9 口径（逐字镜像 tests/test_pad2d.py::_ankle_run）----
ANKLE_D = dict(
    mu_slide=0.4,
    ankle_d_lateral_m=0.030,
    ankle_d_medial_m=0.045,
    ligament_r_mm=22.0,
)
ANKLE_POSTURE = "controlled-drop"
ANKLE_MASS_KG = 80.0
V6_EXPECT = {          # [T22] Table 4 定性结构（照抄测试）
    (5.0, 10.0): True, (5.0, 20.0): True, (5.0, 30.0): False,
    (10.0, 10.0): True, (10.0, 20.0): True, (10.0, 30.0): True,
    (20.0, 10.0): True, (20.0, 20.0): True, (20.0, 30.0): True,
}


def _load() -> dict:
    if OUT.exists():
        return json.loads(OUT.read_text(encoding="utf-8"))
    return {}


def _save(data: dict) -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    data["meta"] = {
        "head_posture": HEAD_POSTURE,
        "head_mass_kg": HEAD_MASS_KG,
        "ankle_posture": ANKLE_POSTURE,
        "ankle_mass_kg": ANKLE_MASS_KG,
        "k_necks": list(K_NECKS),
        "calib_k": CALIB_K,
        "solver_pad": SOLVER_PAD,
        "solver_rigid": SOLVER_RIGID,
        "v8_observable_primary": "t_torso_ms (f_flex 峰 = 胫骨轴向载荷; T22 Table3 第一列同源)",
        "v8_observable_secondary": "t_moment_ms (踝内翻力矩 = 韧带序列; T22 Table3 第二列)",
        "rigid_chain_default": "kf_mult=2.0 x_max=0.015 k_stop=None (iso2+soft, L2 标定)",
        "note": "V3/V8 口径=feet-first-stiff（论文§4.5 复核）；"
                "V6/V9 口径=controlled-drop（镜像 test_pad2d._ankle_run）；"
                "垫格紧容差逐位、刚格放松容差（等价验证 <0.004%）；"
                "峰时一律限冲击窗内（head_branch 窗口守卫）",
    }
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")


# ---- 公共测量工具 ---------------------------------------------------------
def _peak_rel_ms(r, arr) -> float:
    """冲击窗内 |arr| 首峰的 (t - t0)，ms。

    必须限制在 ``meta['impact_window_s']`` 内：刚性地面无接触阻尼，
    settle 事件不触发，轨迹跑满 t_max，窗后振铃里有 50+ g 的假峰
    （argmax 全程会取到 t0+1300 ms 的假值）。与 head_branch / bone.py
    的窗口守卫同款。
    """
    t = np.asarray(r.t_s, dtype=float)
    a = np.asarray(arr, dtype=float)
    t0, t1 = r.meta["impact_window_s"]
    m = (t >= t0) & (t <= t1)
    i = int(np.argmax(np.abs(a[m])))
    return round((float(t[m][i]) - t0) * 1e3, 2)


def _window_ms(r) -> float:
    w = r.meta["impact_window_s"]
    return round((w[1] - w[0]) * 1e3, 2)


def _run_fall(h: float, on_pad: bool, posture: str, mass: float):
    solver = SOLVER_PAD if on_pad else SOLVER_RIGID
    return simulate_boulder_fall(
        height_m=h, mass_kg=mass, posture=posture, on_pad=on_pad, **solver)


def _ankle_run(height_m: float, beta0_deg: float, *, sign: float, on_pad: bool,
               posture: str | None = None):
    p = dict(ANKLE_D)
    p["ankle_beta0_rad"] = float(np.deg2rad(beta0_deg)) * float(sign)
    solver = SOLVER_PAD if on_pad else SOLVER_RIGID
    return simulate_boulder_fall_2d(
        height_m=height_m, mass_kg=ANKLE_MASS_KG,
        posture=posture or ANKLE_POSTURE,
        on_pad=on_pad, posture2d=Posture2D(**p), **solver,
    )


# ==========================================================================
# 单元格 workers（模块级 ⇒ spawn 可 pickle；单池内按 kind 分发）
# ==========================================================================
def _head_cell(task: dict) -> dict:
    h, on_pad, tag = task["h"], task["on_pad"], task["tag"]
    r = _run_fall(h, on_pad, HEAD_POSTURE, HEAD_MASS_KG)
    row = {
        "surface": tag, "h_m": h,
        "window_ms": _window_ms(r),
        "tF_ms": _peak_rel_ms(r, r.pad_force_n),
        "t_torso_ms": _peak_rel_ms(r, r.accel_torso_g),
    }
    for kn in K_NECKS:
        v = simulate_head_branch(r, HeadBranchParams(k_neck_n_per_mm=kn))
        row[f"tpeak_ms_k{int(kn)}"] = round(v.peak_time_s * 1e3, 2)
    if h in V8_HEIGHTS:
        r2 = _ankle_run(h, V8_MOMENT_BETA0_DEG, sign=+1.0, on_pad=on_pad,
                        posture=HEAD_POSTURE)
        row["t_moment_ms"] = _peak_rel_ms(r2, r2.ankle_inversion_moment_nmm)
    if tag == "pad":
        row["paper_4_5_ms"] = PAD_PAPER_MS.get(h)
    return row


def _v6_cell(task: dict) -> dict:
    h, beta, on_pad, tag = task["h"], task["beta"], task["on_pad"], task["tag"]
    r = _ankle_run(h, beta, sign=+1.0, on_pad=on_pad)
    return {"surface": tag, "h_m": h, "beta_deg": beta,
            "fracture": bool(r.ankle_fracture),
            "expect": V6_EXPECT[(h, beta)],
            "util": round(float(r.ankle_bone_utilization), 4),
            "mode": r.ankle_mode}


def _v9_cell(task: dict) -> dict:
    h, on_pad, tag = task["h"], task["on_pad"], task["tag"]
    inv = _ankle_run(h, 20.0, sign=+1.0, on_pad=on_pad)
    ev = _ankle_run(h, 20.0, sign=-1.0, on_pad=on_pad)
    strain_ok = ev.peak_ligament_strain > inv.peak_ligament_strain
    util_ok = ev.ankle_bone_utilization > inv.ankle_bone_utilization
    return {
        "surface": tag, "h_m": h,
        "inv_strain": float(inv.peak_ligament_strain),
        "ev_strain": float(ev.peak_ligament_strain),
        "inv_util": round(float(inv.ankle_bone_utilization), 4),
        "ev_util": round(float(ev.ankle_bone_utilization), 4),
        "strain_ok": strain_ok, "util_ok": util_ok,
    }


def _run_task(task: dict):
    kind = task["kind"]
    if kind == "head":
        return _head_cell(task)
    if kind == "v6":
        return _v6_cell(task)
    return _v9_cell(task)


def _build_tasks(parts: list[str]) -> list[dict]:
    tasks: list[dict] = []
    if "head" in parts:
        tasks += [dict(kind="head", h=h, on_pad=False, tag="rigid")
                  for h in HEAD_HEIGHTS_RIGID]
        tasks += [dict(kind="head", h=h, on_pad=True, tag="pad")
                  for h in HEAD_HEIGHTS_PAD]
    if "v6" in parts:
        tasks += [dict(kind="v6", h=h, beta=a, on_pad=(tag == "pad"), tag=tag)
                  for tag in ("pad", "rigid") for (h, a) in sorted(V6_EXPECT)]
    if "v9" in parts:
        tasks += [dict(kind="v9", h=h, on_pad=(tag == "pad"), tag=tag)
                  for tag in ("pad", "rigid") for h in (2.0, 3.0, 5.0)]
    return tasks


# ==========================================================================
# 结果装配（纯 CPU，进程池外）
# ==========================================================================
def _assemble_head(rows: list[dict]) -> dict:
    rows = sorted(rows, key=lambda r: (r["surface"] != "rigid", r["h_m"]))
    for row in rows:
        extra = ""
        if "t_moment_ms" in row:
            extra = f" tMom={row['t_moment_ms']:>6.2f}"
        print(f"[head] {row['surface']:>5} h={row['h_m']:>4.0f} "
              f"win={row['window_ms']:>6.2f} tF={row['tF_ms']:>6.2f} "
              f"tTorso={row['t_torso_ms']:>6.2f}{extra} "
              + " ".join(f"t{k:4.0f}={row[f'tpeak_ms_k{int(k)}']:>6.2f}"
                         for k in K_NECKS))

    verdict = {}
    for tag in ("rigid", "pad"):
        by_h = {r["h_m"]: r for r in rows if r["surface"] == tag}
        # V3：主判据 = 该面的标定 k_neck 列；全 k 列一并给出
        calib = int(CALIB_K[tag])
        v3_by_k = {
            f"k{int(k)}": {str(int(h)): (by_h[h][f"tpeak_ms_k{int(k)}"]
                                         <= V3_THRESHOLD_S * 1e3)
                           for h in V3_HEIGHTS if h in by_h}
            for k in K_NECKS
        }
        v3 = v3_by_k[f"k{calib}"]
        # V8（L4 观测量对齐）：主判据 = t_torso_ms（f_flex 峰 = 胫骨轴向
        # 载荷，= [T22] Table 3 第一列胫骨应变的同源信号）；
        # 踝矩 t_moment = 韧带序列（第二列）副判据；头峰保留为 legacy。
        seq_t = [by_h[h]["t_torso_ms"] for h in V8_HEIGHTS if h in by_h]
        mono_t = all(seq_t[i] <= seq_t[i - 1] for i in range(1, len(seq_t)))
        seq_mom = [by_h[h]["t_moment_ms"] for h in V8_HEIGHTS if h in by_h
                   and "t_moment_ms" in by_h[h]]
        mono_mom = all(seq_mom[i] <= seq_mom[i - 1] for i in range(1, len(seq_mom)))
        mono_head = {}
        for k in K_NECKS:
            seq = [by_h[h][f"tpeak_ms_k{int(k)}"] for h in V8_HEIGHTS if h in by_h]
            mono_head[f"k{int(k)}"] = all(seq[i] <= seq[i - 1]
                                          for i in range(1, len(seq)))
        at20_t = by_h.get(20.0, {}).get("t_torso_ms")
        at20_mom = by_h.get(20.0, {}).get("t_moment_ms")
        at20_head = {f"k{int(k)}": by_h[20.0][f"tpeak_ms_k{int(k)}"]
                     for k in K_NECKS if 20.0 in by_h}
        verdict[tag] = {
            "V3_calib_k": calib,
            "V3_pass_by_h": v3,
            "V3_pass_by_h_all_k": v3_by_k,
            "V3_pass": all(v3.values()) if v3 else None,
            "V8a_monotone_tibia_proxy": mono_t,          # 主判据（t_torso）
            "V8a_monotone_moment_ligament": mono_mom,    # 副判据（t_moment）
            "V8a_monotone_head_legacy": mono_head,
            "V8b_ttorso_at_20m_ms": at20_t,
            "V8b_pass": (at20_t is not None
                         and at20_t <= V8_THRESHOLD_AT_20M_S * 1e3),
            "V8b_tmoment_at_20m_ms": at20_mom,
            "V8_head_legacy_at20_ms": at20_head,
            "V8_head_legacy_pass": {key: (v <= V8_THRESHOLD_AT_20M_S * 1e3)
                                    for key, v in at20_head.items()},
        }
        print(f"[head] {tag}: V3(k{calib})={verdict[tag]['V3_pass']} "
              f"V8a(tibia)={mono_t} V8a(moment)={mono_mom} "
              f"V8b={verdict[tag]['V8b_pass']} tT@20={at20_t} "
              f"| legacy head mono={mono_head}")
    return {"rows": rows, "verdict": verdict,
            "thresholds_s": {"V3": V3_THRESHOLD_S, "V8_at20": V8_THRESHOLD_AT_20M_S}}


def _assemble_v6(cells_raw: list[dict]) -> dict:
    out = {"expect": {f"{h}/{a}": e for (h, a), e in V6_EXPECT.items()},
           "runs": {}}
    for tag in ("pad", "rigid"):
        cells = {f"{c['h_m']}/{c['beta_deg']}": c
                 for c in cells_raw if c["surface"] == tag}
        hits = sum(c["fracture"] == c["expect"] for c in cells.values())
        for key in sorted(cells, key=lambda s: tuple(map(float, s.split("/")))):
            c = cells[key]
            print(f"[v6] {tag:>5} h={c['h_m']:>4.0f} β={c['beta_deg']:>2.0f}° "
                  f"fracture={c['fracture']} (exp {c['expect']}) "
                  f"util={c['util']:.3f} mode={c['mode']}")
        flips = 0
        ks = sorted(V6_EXPECT)
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                (h0, a0), (h1, a1) = ks[i], ks[j]
                if (h0 == h1 or a0 == a1) and \
                        cells[f"{h0}/{a0}"]["fracture"] != cells[f"{h1}/{a1}"]["fracture"]:
                    flips += 1
        out["runs"][tag] = {"cells": cells, "hits": hits,
                            "total": len(V6_EXPECT), "transition_flips": flips,
                            "structure_pass": hits == len(V6_EXPECT) and flips >= 1}
        print(f"[v6] {tag}: {hits}/{len(V6_EXPECT)} 命中, flips={flips} "
              f"-> {'PASS' if out['runs'][tag]['structure_pass'] else 'FAIL'}")
    return out


def _assemble_v9(rows: list[dict]) -> dict:
    out = {"runs": {}}
    for tag in ("pad", "rigid"):
        tag_rows = sorted((r for r in rows if r["surface"] == tag),
                          key=lambda r: r["h_m"])
        ok = all(r["strain_ok"] and r["util_ok"] for r in tag_rows)
        for r in tag_rows:
            print(f"[v9] {tag:>5} h={r['h_m']:.0f} strain(inv={r['inv_strain']:.4f}"
                  f" < ev={r['ev_strain']:.4f})={r['strain_ok']} "
                  f"util(inv={r['inv_util']:.3f} < ev={r['ev_util']:.3f})"
                  f"={r['util_ok']}")
        out["runs"][tag] = {"rows": tag_rows, "pass": ok}
        print(f"[v9] {tag}: {'PASS' if ok else 'FAIL'}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", choices=["head", "v6", "v9", "all"], default="all")
    ap.add_argument("--workers", type=int, default=None,
                    help="并行进程数（默认自动：核数-2，上限 14）")
    args = ap.parse_args()

    data = _load()
    parts = (["head", "v6", "v9"] if args.part == "all" else [args.part])

    tasks = _build_tasks(parts)
    print(f"===== {len(tasks)} cells → one pool =====", flush=True)
    results = parallel_map(_run_task, tasks, workers=args.workers,
                           progress=lambda i, n: print(f"  {i}/{n}", flush=True))

    # 按 kind 分拣回各 part
    by_kind: dict[str, list] = {}
    for t_, res in zip(tasks, results):
        by_kind.setdefault(t_["kind"], []).append(res)

    for p in parts:
        print(f"\n===== part: {p} =====")
        if p == "head":
            data["V3_V8"] = _assemble_head(by_kind.get("head", []))
        elif p == "v6":
            data["V6"] = _assemble_v6(by_kind.get("v6", []))
        else:
            data["V9"] = _assemble_v9(by_kind.get("v9", []))
        _save(data)                       # 每个 part 完成即落盘
    print(f"\nwritten: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
