"""论文可信度 ②③ —— MC 统计稳健性三件套 + α_pad 踝拆分实验。

为什么不直接跑 phase_p1_montecarlo.py
-------------------------------------
官方产物 ``results/p1_montecarlo.{txt,csv}``（seed=20261003）是论文表 3 的
数据源，**不可被覆盖**。本脚本只读官方 CSV 做离线重算，或经
``run_trials()`` 自跑并把结果写到 ``results/mc_stat_robustness.json``。

三个实验
--------
1. **threshold 离线扫描**（--offline）：官方 CSV 只存了 ``ankle_fracture``
   布尔值，但 ``mixed`` 模式 = 0（实测 655 none / 228 stuck / 117 sliding）
   ⇒ 每个 trial 的 ``peak_f_lig_eff`` 可由 ``peak_ligament_strain`` 精确还原：
       stuck   ⇒ eff = strain × k·L = strain × 308 N（unload=1）
       sliding ⇒ eff = 0.005 × strain × 308 N
       none    ⇒ eff = 0（β≡0 ⇒ M≡0 ⇒ f_lig≡0）
   **正对照**：T=215 重算的骨折数必须逐位等于官方 228，否则整段作废。
   扫 T ∈ {154, 165, 215, 265, 349}（±50 + 临床锁边界）。
2. **多种子稳定性**（--runs）：seed ∈ {20261004, 20261005} @ α=1，
   与官方 20261003 一起给 34.5% / 22.8% 的种子漂移。
3. **α_pad 踝骨折拆分**（--runs）：seed=20261003 固定，α ∈ {1, 10, 30}。
   α=1 兼**管线正对照**（须复现官方 345/228）。回答论文
   ``关键局限(2)`` 自认的「踝骨折拆分对垫刚度的影响未测」。

注入方式
--------
``_run_one`` 走 ``simulate_boulder_fall_2d(**task["kw"])``，kw 在**主进程**
阶段 A 从模块全局 ``FAST_KW`` 构建后随 task pickle 进 spawn worker ⇒
patch ``phase_p1_montecarlo.FAST_KW`` 加 ``pad=CrashPad(×α)`` 即可跨进程生效；
α=1 时 ``×1.0`` IEEE 精确 ⇒ 与官方默认 ``pad=None → CrashPad()`` 逐位等价。

用法::

    .venv\\Scripts\\python.exe scripts/paper/mc_stat_robustness.py --offline
    .venv\\Scripts\\python.exe scripts/paper/mc_stat_robustness.py --runs
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

OUT = ROOT / "results" / "mc_stat_robustness.json"
OFFICIAL_CSV = ROOT / "results" / "p1_montecarlo.csv"

# pad2d 失效判据常量（k·L = 14 N/mm × 22 mm）—— 与 src/climbing/pad2d.py
# 的 ATFL_STIFFNESS_N_PER_MM / ATFL_RESTING_LENGTH_MM 一致；正对照会验证。
KL_N = 14.0 * 22.0          # 308 N
UNLOAD_SLIDE = 0.005        # pad2d.ANKLE_SLIDE_UNLOAD
OFFICIAL_SEED = 20261003
OFFICIAL_N = 1000
# 官方产物计数（论文表 3 数据源）—— 正对照基准
OFFICIAL_COUNTS = {"injured": 345, "sprain": 345, "fracture": 228}

THRESHOLDS_N = (154.0, 165.0, 215.0, 265.0, 349.0)
SEEDS = (OFFICIAL_SEED, 20261004, 20261005)
ALPHAS = (1.0, 10.0, 30.0)


def wilson(k: int, n: int, z: float = 1.96) -> list[float]:
    """Wilson 二项比例 95% CI。"""
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2.0 * n)) / d
    h = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def eff_from_row(row: dict) -> float:
    """由 CSV 行精确还原 peak_f_lig_eff（依赖 mixed=0，见模块 docstring）。"""
    mode = row["ankle_mode"]
    strain = float(row["peak_ligament_strain"])
    if mode == "stuck":
        return strain * KL_N
    if mode == "sliding":
        return UNLOAD_SLIDE * strain * KL_N
    if mode == "none":
        return 0.0
    raise AssertionError(f"mixed 模式出现（离线还原不成立）: {row}")


def load_official() -> list[dict]:
    with OFFICIAL_CSV.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_json() -> dict:
    if OUT.exists():
        return json.loads(OUT.read_text(encoding="utf-8"))
    return {}


def save_json(data: dict) -> None:
    OUT.write_text(json.dumps(data, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print(f"  -> {OUT.relative_to(ROOT)}", flush=True)


# ---------------------------------------------------------------- offline --
def part_offline() -> None:
    data = load_json()
    rows = load_official()
    assert len(rows) == OFFICIAL_N, f"官方 CSV 行数 {len(rows)} != {OFFICIAL_N}"

    # ---- 正对照：T=215 离线重算必须逐位复现官方计数 ----
    off_inj = sum(int(r["ankle_injured"]) for r in rows)
    off_sp = sum(int(r["ankle_sprain"]) for r in rows)
    off_fr = sum(int(r["ankle_fracture"]) for r in rows)
    recount215 = sum(1 for r in rows if eff_from_row(r) >= 215.0)
    pc = {
        "official_counts": {"injured": off_inj, "sprain": off_sp,
                            "fracture": off_fr},
        "recount_T215": recount215,
        "match": recount215 == off_fr == OFFICIAL_COUNTS["fracture"],
    }
    print(f"[offline] official: injured={off_inj} sprain={off_sp} "
          f"fracture={off_fr}", flush=True)
    print(f"[offline] recount T=215: {recount215}  match={pc['match']}",
          flush=True)
    if not pc["match"]:
        print("[offline] POSITIVE CONTROL FAILED -- abort, do not trust scan",
              flush=True)
        data["threshold_scan"] = {"positive_control": pc, "aborted": True}
        save_json(data)
        return

    # ---- 阈值扫描 ----
    grid = []
    for t in THRESHOLDS_N:
        fr = sum(1 for r in rows if eff_from_row(r) >= t)
        inj = sum(1 for r in rows
                  if int(r["ankle_sprain"]) or eff_from_row(r) >= t)
        grid.append({"T_N": t, "fracture_n": fr,
                     "fracture_share": round(fr / OFFICIAL_N, 4),
                     "injured_n": inj,
                     "injured_share": round(inj / OFFICIAL_N, 4),
                     "ci95_fracture": wilson(fr, OFFICIAL_N)})
        print(f"  T={t:6.1f} N -> fracture {fr:3d} ({fr/10:.1f}%)  "
              f"injured {inj:3d} ({inj/10:.1f}%)", flush=True)

    # ---- 官方点估计的 Wilson CI ----
    ci = {"injured_345": wilson(OFFICIAL_COUNTS["injured"], OFFICIAL_N),
          "fracture_228": wilson(OFFICIAL_COUNTS["fracture"], OFFICIAL_N),
          "sprain_345": wilson(OFFICIAL_COUNTS["sprain"], OFFICIAL_N)}
    print(f"[offline] CI injured 34.5%: {ci['injured_345']}", flush=True)
    print(f"[offline] CI fracture 22.8%: {ci['fracture_228']}", flush=True)

    data["threshold_scan"] = {"positive_control": pc, "grid": grid,
                              "wilson_ci_official": ci}
    save_json(data)


# ------------------------------------------------------------------- runs --
def summarize(rows: list[dict]) -> dict:
    n = len(rows)
    inj = sum(int(r["ankle_injured"]) for r in rows)
    sp = sum(int(r["ankle_sprain"]) for r in rows)
    fr = sum(int(r["ankle_fracture"]) for r in rows)
    bone = sum(int(r["n_fractured_sites"]) for r in rows)
    no_rot = [r for r in rows if int(r["rotation_free"]) == 1]
    return {
        "n": n,
        "injured": inj, "injured_share": round(inj / n, 4),
        "sprain": sp, "sprain_share": round(sp / n, 4),
        "fracture": fr, "fracture_share": round(fr / n, 4),
        "bone_events": bone,
        "worst_util_max": round(max(float(r["worst_utilization"])
                                    for r in rows), 3),
        "no_rot_n": len(no_rot),
        "no_rot_injured": sum(int(r["ankle_injured"]) for r in no_rot),
        "ci95_fracture": wilson(fr, n),
        "ci95_injured": wilson(inj, n),
    }


def scaled_pad(alpha: float):
    from climbing.pad import CrashPad
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


def run_one(alpha: float, seed: int) -> dict:
    import phase_p1_montecarlo as mc
    base = dict(mc.FAST_KW)
    mc.FAST_KW = {**base, "pad": scaled_pad(alpha)}
    try:
        rows, site_rows = mc.run_trials(OFFICIAL_N, seed, fast=True,
                                        workers=None)
    finally:
        mc.FAST_KW = base
    out = summarize(rows)
    out.update({"alpha_pad": alpha, "seed": seed})
    return out


def part_runs() -> None:
    data = load_json()
    runs: list[dict] = data.get("runs", [])
    have = {(r["alpha_pad"], r["seed"]) for r in runs}

    jobs = [(a, s) for a in ALPHAS for s in (OFFICIAL_SEED,)]
    jobs += [(1.0, s) for s in SEEDS if s != OFFICIAL_SEED]

    for alpha, seed in jobs:
        if (alpha, seed) in have:
            print(f"[skip] alpha={alpha} seed={seed} (already done)",
                  flush=True)
            continue
        print(f"[run ] alpha_pad={alpha} seed={seed} n={OFFICIAL_N} ...",
              flush=True)
        res = run_one(alpha, seed)
        runs.append(res)
        data["runs"] = runs
        save_json(data)   # 增量落盘，中断不丢
        print(f"       injured={res['injured']} ({res['injured_share']:.1%}) "
              f"fracture={res['fracture']} ({res['fracture_share']:.1%}) "
              f"bone={res['bone_events']} worst_util={res['worst_util_max']}",
              flush=True)

    # ---- 正对照：α=1 / seed=20261003 必须复现官方 345/228 ----
    a1 = next((r for r in runs if r["alpha_pad"] == 1.0
               and r["seed"] == OFFICIAL_SEED), None)
    if a1:
        pc_ok = (a1["injured"] == OFFICIAL_COUNTS["injured"]
                 and a1["fracture"] == OFFICIAL_COUNTS["fracture"])
        data["alpha_pipeline_control"] = {
            "expected": OFFICIAL_COUNTS,
            "got": {"injured": a1["injured"], "fracture": a1["fracture"]},
            "match": pc_ok,
        }
        print(f"[control] alpha=1 pipeline vs official: match={pc_ok}",
              flush=True)

    # ---- 汇总视图 ----
    seeds = sorted((r for r in runs if r["alpha_pad"] == 1.0),
                   key=lambda r: r["seed"])
    if len(seeds) >= 2:
        inj = [r["injured_share"] for r in seeds]
        fr = [r["fracture_share"] for r in seeds]
        data["seed_stability"] = {
            "seeds": [r["seed"] for r in seeds],
            "injured_share_range": [min(inj), max(inj)],
            "fracture_share_range": [min(fr), max(fr)],
        }
        print(f"[seeds] injured {min(inj):.1%}..{max(inj):.1%}  "
              f"fracture {min(fr):.1%}..{max(fr):.1%}", flush=True)

    alphas = sorted((r for r in runs if r["seed"] == OFFICIAL_SEED),
                    key=lambda r: r["alpha_pad"])
    if alphas:
        data["alpha_pad_effect"] = [
            {"alpha": r["alpha_pad"], "fracture_share": r["fracture_share"],
             "injured_share": r["injured_share"],
             "sprain_share": r["sprain_share"]} for r in alphas]
        print("[alpha ] " + " | ".join(
            f"a={r['alpha_pad']:g}: fr={r['fracture_share']:.1%} "
            f"sp={r['sprain_share']:.1%}" for r in alphas), flush=True)

    save_json(data)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    ap.add_argument("--runs", action="store_true")
    args = ap.parse_args()
    if not (args.offline or args.runs):
        args.offline = args.runs = True
    if args.offline:
        part_offline()
    if args.runs:
        part_runs()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
