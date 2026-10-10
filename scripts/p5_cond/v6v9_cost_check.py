"""P5 前置侦察 · Step 3 —— 加硬垫子会不会打破已通过的 V6 / V9（代价侧）。

Step 2 证明 α=10~30 就能让 V3 过（h=2/3/5）。但 α 是有代价的：
peak_pad_n 随 α 上升（h=5: 51.6 kN → 69.5 kN @α=30）。
V6/V9 的判据直接建立在这些载荷上 ⇒ **必须实跑确认没有副作用**。

* α=1  = 发表基线 ⇒ **正对照**（不复现则本脚本结果作废）
* α=30 = V3 刚好过的档位（h=2/3 需 α=30）⇒ 代价检验对象

判据：
  V6 = 9 格全中 [T22] Table 4 定性结构（含 (5 m,30°) 必须**不**骨折）
  V9 = 外翻 strain/util 均 > 内翻（3 高度）
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.batch import parallel_map  # noqa: E402
from climbing.pad import CrashPad  # noqa: E402
from climbing.pad2d import Posture2D, simulate_boulder_fall_2d  # noqa: E402

SOLVER_PAD = dict(max_ode_step=1e-4, rtol=1e-8, atol=1e-11)
ANKLE_D = dict(mu_slide=0.4, ankle_d_lateral_m=0.030,
               ankle_d_medial_m=0.045, ligament_r_mm=22.0)
ANKLE_POSTURE = "controlled-drop"
ANKLE_MASS_KG = 80.0
V6_EXPECT = {
    (5.0, 10.0): True, (5.0, 20.0): True, (5.0, 30.0): False,
    (10.0, 10.0): True, (10.0, 20.0): True, (10.0, 30.0): True,
    (20.0, 10.0): True, (20.0, 20.0): True, (20.0, 30.0): True,
}
ALPHAS = (1.0, 30.0)
OUT = ROOT / "results" / "p5_v6v9_cost.json"
SRC_JSON = ROOT / "results" / "rigid_vtarget_probe.json"


def scaled_pad(alpha: float) -> CrashPad:
    b = CrashPad()
    return CrashPad(
        thickness_m=b.thickness_m,
        s_plateau_pa=b.s_plateau_pa * alpha,
        s_dense_pa=b.s_dense_pa * alpha,
        xi_plateau=b.xi_plateau,
        n_dense=b.n_dense,
        xi_dense_max=b.xi_dense_max,
        bottoming_modulus_pa=b.bottoming_modulus_pa * alpha,
        linear_modulus_pa=b.linear_modulus_pa,
    )


def _ankle_run(h: float, beta_deg: float, *, sign: float, alpha: float):
    p = dict(ANKLE_D)
    p["ankle_beta0_rad"] = float(np.deg2rad(beta_deg)) * float(sign)
    return simulate_boulder_fall_2d(
        height_m=h, mass_kg=ANKLE_MASS_KG, posture=ANKLE_POSTURE,
        on_pad=True, pad=scaled_pad(alpha),
        posture2d=Posture2D(**p), **SOLVER_PAD)


# ---- module-level workers (spawn picklable) ------------------------------
def _v6(task):
    h, beta, alpha = task["h"], task["beta"], task["alpha"]
    r = _ankle_run(h, beta, sign=+1.0, alpha=alpha)
    return {"alpha": alpha, "h_m": h, "beta_deg": beta,
            "fracture": bool(r.ankle_fracture),
            "expect": V6_EXPECT[(h, beta)],
            "util": round(float(r.ankle_bone_utilization), 4),
            "mode": r.ankle_mode}


def _v9(task):
    h, alpha = task["h"], task["alpha"]
    inv = _ankle_run(h, 20.0, sign=+1.0, alpha=alpha)
    ev = _ankle_run(h, 20.0, sign=-1.0, alpha=alpha)
    return {"alpha": alpha, "h_m": h,
            "inv_strain": float(inv.peak_ligament_strain),
            "ev_strain": float(ev.peak_ligament_strain),
            "inv_util": round(float(inv.ankle_bone_utilization), 4),
            "ev_util": round(float(ev.ankle_bone_utilization), 4),
            "strain_ok": bool(ev.peak_ligament_strain > inv.peak_ligament_strain),
            "util_ok": bool(ev.ankle_bone_utilization > inv.ankle_bone_utilization)}


def main() -> int:
    v6_tasks = [dict(kind="v6", h=h, beta=b, alpha=a)
                for a in ALPHAS for (h, b) in sorted(V6_EXPECT)]
    v9_tasks = [dict(kind="v9", h=h, alpha=a)
                for a in ALPHAS for h in (2.0, 3.0, 5.0)]

    print(f"[plan] V6 tasks={len(v6_tasks)}  V9 tasks={len(v9_tasks)}", flush=True)
    v6 = parallel_map(_v6, v6_tasks)
    print("[done] V6", flush=True)
    v9 = parallel_map(_v9, v9_tasks)
    print("[done] V9", flush=True)

    # ---- 正对照：alpha=1 必须复现已发表【垫】 ---------------------------
    print()
    print("=" * 74)
    print("POSITIVE CONTROL  (alpha=1 vs published pad column)")
    print("=" * 74)
    pub = json.loads(SRC_JSON.read_text(encoding="utf-8"))
    ok = True
    for row in v6:
        if row["alpha"] != 1.0:
            continue
        cell = pub["V6"]["runs"]["pad"]["cells"][f"{row['h_m']:.1f}/{row['beta_deg']:.1f}"]
        same = (row["fracture"] == cell["fracture"]
                and row["mode"] == cell["mode"]
                and abs(row["util"] - cell["util"]) < 5e-4)
        ok &= same
        print(f"  V6 {row['h_m']:>4.0f}m/{row['beta_deg']:>2.0f}deg  "
              f"pub={cell['fracture']}/{cell['mode']}/{cell['util']:<8}  "
              f"got={row['fracture']}/{row['mode']}/{row['util']:<8}  "
              f"{'PASS' if same else 'FAIL'}")
    for row in v9:
        if row["alpha"] != 1.0:
            continue
        cell = next(c for c in pub["V9"]["runs"]["pad"]["rows"]
                    if float(c["h_m"]) == row["h_m"])
        same = (abs(row["inv_strain"] - cell["inv_strain"]) < 1e-6
                and abs(row["ev_strain"] - cell["ev_strain"]) < 1e-6
                and row["strain_ok"] == cell["strain_ok"]
                and row["util_ok"] == cell["util_ok"])
        ok &= same
        print(f"  V9 {row['h_m']:>4.0f}m  strain_ok={row['strain_ok']} "
              f"util_ok={row['util_ok']}  {'PASS' if same else 'FAIL'}")
    if not ok:
        print("!! POSITIVE CONTROL FAILED -- INVALID")
        OUT.write_text(json.dumps({"v6": v6, "v9": v9, "positive_control": False},
                                  ensure_ascii=False, indent=2), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")

    # ---- 判据 + 代价 ---------------------------------------------------
    print()
    print("=" * 74)
    print(f"V6 / V9 at alpha in {ALPHAS}   (V6=9/9 must hold, V9=3/3)")
    print("=" * 74)
    verdict = {}
    for a in ALPHAS:
        rows6 = [r for r in v6 if r["alpha"] == a]
        rows9 = [r for r in v9 if r["alpha"] == a]
        hits = sum(1 for r in rows6 if r["fracture"] == r["expect"])
        v9_ok = sum(1 for r in rows9 if r["strain_ok"] and r["util_ok"])
        verdict[str(a)] = {"v6_hits": hits, "v6_total": len(rows6),
                           "v9_hits": v9_ok, "v9_total": len(rows9),
                           "v6_pass": hits == len(rows6),
                           "v9_pass": v9_ok == len(rows9)}
        print(f"\n-- alpha = {a:g} --")
        for r in rows6:
            mark = "ok  " if r["fracture"] == r["expect"] else "MISS"
            print(f"  V6 {mark} {r['h_m']:>4.0f}m/{r['beta_deg']:>2.0f}deg  "
                  f"fracture={str(r['fracture']):<5} expect={str(r['expect']):<5} "
                  f"util={r['util']:<9} mode={r['mode']}")
        for r in rows9:
            mark = "ok  " if (r["strain_ok"] and r["util_ok"]) else "MISS"
            print(f"  V9 {mark} h={r['h_m']:>4.0f}m  inv_e={r['inv_strain']:.3f} "
                  f"ev_e={r['ev_strain']:.3f}  "
                  f"inv_u={r['inv_util']:<8} ev_u={r['ev_util']:<8} "
                  f"strain_ok={r['strain_ok']} util_ok={r['util_ok']}")
        v = verdict[str(a)]
        print(f"  => V6 {v['v6_hits']}/{v['v6_total']} "
              f"{'PASS' if v['v6_pass'] else 'FAIL'}   "
              f"V9 {v['v9_hits']}/{v['v9_total']} "
              f"{'PASS' if v['v9_pass'] else 'FAIL'}")

    print()
    print("=" * 74)
    a1, a30 = verdict[str(ALPHAS[0])], verdict[str(ALPHAS[1])]
    cost = ("no regression" if (a30["v6_pass"] and a30["v9_pass"])
            else "REGRESSION -- stiffening breaks V6/V9")
    print(f"COST VERDICT (alpha {ALPHAS[0]:g} -> {ALPHAS[1]:g}): {cost}")
    print("=" * 74)

    OUT.write_text(json.dumps({"v6": v6, "v9": v9, "positive_control": True,
                               "verdict": verdict}, ensure_ascii=False,
                              indent=2), encoding="utf-8")
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
