"""M1 全量 S1a —— 【零 FE 成本的决定性检查】share -> σ(h) -> G3 OR 是否落进 DoD-2 窗口。

**为什么这一步排在建两体 FE 之前**：M1 全量的唯一目的是把 G9 的 DoD-2
（腓骨 OR 落不进 [1.177, 2.187]）从"share 假设"层面搬到"物理"层面。
而 DoD-2 的判据链是 **share -> σ(h) 行 -> G3 概率机制 -> OR**，
其中 **σ(h) 行 = share × load(h)/A**（严格线弹性）是**解析量**：
G9 的 FE 在该链路上只负责提供 σ_ref（即 share·F/A），一旦 share 定了，
σ(h) 行就定了，**与 FE 建模细节无关**。
⇒ 所以：**先用解析 σ(h) 行喂 G3，看 OR 落不落窗**。
若落窗 ⇒ DoD-2 已由"物理 share"解决，**两体 FE 降级为验证 share 的手段（而非出口）**；
若不落窗 ⇒ 才需要两体 FE（因为那时怀疑的是"闭式 share 本身不对"）。

**口径纪律**：面积口径 3 个全扫（同 m1_precheck），因为 S0 显示 2/3 支持、1/3 不支持
⇒ 面积是主控变量，不是 share。

产物：results/opensim_fe/m1_g3_or_precheck.json + M1_G3_OR_PRECHECK.md
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fracture_matrix as fm  # noqa: E402
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402
from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "m1_g3_or_precheck.json"
OUT_REPORT = RES / "M1_G3_OR_PRECHECK.md"
SIGMA_C = 70.0
FIBULA = "腓骨"

AREAS = {"A_cort_measured_101.9": 101.9,
         "A_eff_loadflow_76.01": 76.01,
         "A_equiv_median_33": 33.0}
SHARES = [0.05, 0.10, 0.2241, 0.282, 0.30]


def or_for(sigma_h: np.ndarray) -> dict:
    """把一条腓骨 σ(h) 行喂进 G3（其余 8 部位冻结不动），取主设定腓骨 OR。"""
    inputs = g3.load_inputs()
    mixed, idx = g8.build_mixed_inputs(inputs, np.asarray(sigma_h, float))
    configs = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                 g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
               for sd in g3.SEEDS}
    primary = configs[g3.OFFICIAL_SEED]
    fr = g3._site_or_record(primary["or_map"][FIBULA])
    unc = g8.unchanged_check(inputs, idx, primary)
    return {"status": fr.get("status"), "OR": fr.get("OR"), "ci95": fr.get("ci95"),
            "n_fracture": fr.get("n_fracture"),
            "seed_ors": {str(sd): configs[sd]["or_map"][FIBULA].get("OR")
                         for sd in g3.SEEDS},
            "a3": primary["a3"],
            "unchanged_other8": unc.get("all_equal"),
            "max_abs_or_diff": unc.get("max_abs_or_diff")}


def main() -> int:
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    F_ref = float(loads[ref_idx])
    ratio = loads / F_ref

    target = float(g3.A3_TARGETS[FIBULA])
    win = [target * 0.7, target * 1.3]
    print(f"ankle_r @REF={fm.REF_H:.0f}m = {F_ref:.1f} N; 腓骨靶 OR={target} "
          f"±30% => 窗口 [{win[0]:.3f}, {win[1]:.3f}]", flush=True)

    # 正对照：冻结基线（G3 原行，腓骨不动）应**逐部位**复现已存 g3_or_estimates.json。
    # ⚠️ v1 把正对照写成"比腓骨 OR"⇒ 必然失败：冻结基线的腓骨是
    #    status="all_one"（区间内始终骨折）⇒ **结构上就没有 OR**。
    #    正确的对照对象是 (status, n_fracture, OR) 三元组的全部 9 部位。
    print("\n[C0] 正对照 —— G3 冻结原行（腓骨不动）应逐部位复现已存记录", flush=True)
    inputs = g3.load_inputs()
    base = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                              g3.COV_A_PRIMARY, inputs, g3.N_SUBJ, True)
            for sd in g3.SEEDS}
    b_map = base[g3.OFFICIAL_SEED]["or_map"]
    stored = json.loads(g8.G3_OR_JSON.read_text(encoding="utf-8"))["primary"]["per_site"]

    def _cmp(cn: str) -> dict:
        a, b = g3._site_or_record(b_map[cn]), stored.get(cn, {})
        return {"status_match": a.get("status") == b.get("status"),
                "n_fracture_match": a.get("n_fracture") == b.get("n_fracture"),
                "recomputed_OR": a.get("OR"), "stored_OR": b.get("OR"),
                "or_diff": (abs(float(a["OR"]) - float(b["OR"]))
                            if a.get("OR") is not None and b.get("OR") is not None
                            else None)}

    per = {cn: _cmp(cn) for cn in b_map}
    n_bad = [cn for cn, v in per.items() if not (v["status_match"]
                                                 and v["n_fracture_match"]
                                                 and (v["or_diff"] is None
                                                      or v["or_diff"] < 1e-9))]
    c0 = {"n_sites": len(per), "n_mismatch": len(n_bad), "mismatched": n_bad,
          "agrees": not n_bad, "per_site": per,
          "note": ("冻结基线腓骨 = all_one（无 OR）⇒ 正对照比对 "
                   "(status, n_fracture, OR) 三元组，不比腓骨 OR")}
    print(f"  {len(per)} 部位，{len(n_bad)} 个不一致 ⇒ "
          f"{'✅ 正对照通过' if c0['agrees'] else '❌ 正对照失败 ' + str(n_bad)}",
          flush=True)
    print(f"  基线腓骨: status={per['腓骨']['status_match'] and 'all_one' or per['腓骨']}"
          f"  n_fracture={b_map['腓骨'].get('n_fracture')}"
          f"  OR={b_map['腓骨'].get('OR')}（结构上不存在）", flush=True)
    if not c0["agrees"]:
        print("  ⇒ 正对照不过，**不继续**（后续 OR 全不可信）")
        OUT_JSON.write_text(json.dumps({"meta": {"generated_at": datetime.now().
                                         strftime("%Y-%m-%d %H:%M:%S")},
                                        "positive_control": c0}, indent=2,
                                       ensure_ascii=False), encoding="utf-8")
        return 2

    out: dict = {"meta": {
        "title": "M1 S1a：share->σ(h)->G3 OR 解析预检（零 FE）",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "ref_h_m": float(fm.REF_H), "ankle_ref_n": F_ref,
        "sigma_c_mpa": SIGMA_C, "fibula_target_or": target, "fibula_window": win,
        "areas_mm2": AREAS, "shares": SHARES,
        "model": ("σ(h) = share·F_ref·(load(h)/load_ref)/A  —— 严格线弹性，"
                  "与 FE 建模细节无关；FE 的唯一作用是提供 share 与 A"),
        "g3": {"n_subj": g3.N_SUBJ, "seeds": g3.SEEDS,
               "official_seed": g3.OFFICIAL_SEED,
               "cov": {"b": g3.COV_B_PRIMARY, "e": g3.COV_E_PRIMARY,
                       "a": g3.COV_A_PRIMARY}},
    }, "positive_control": c0, "areas": {}}

    for aname, area in AREAS.items():
        print("\n" + "=" * 88, flush=True)
        print(f"[{aname}] A={area:.2f} mm²", flush=True)
        rows = {}
        for sh in SHARES:
            sigma_h = sh * F_ref * ratio / area
            rec = or_for(sigma_h)
            rec.update(share=float(sh), alloc_at_ref_n=float(sh * F_ref),
                       sigma_h_ref_mpa=float(sigma_h[ref_idx]),
                       util_max=float((sigma_h / SIGMA_C).max()),
                       n_fractured_h=int((sigma_h / SIGMA_C >= 1.0).sum()),
                       first_fracture_h=next((i + 1 for i, v in
                                              enumerate(sigma_h / SIGMA_C) if v >= 1.0),
                                             None),
                       or_in_window=bool(rec["OR"] is not None
                                         and win[0] <= rec["OR"] <= win[1]))
            rows[f"{sh:.4f}"] = rec
            print(f"  share={sh:6.4f} alloc@{fm.REF_H:.0f}m={sh*F_ref:7.1f} N  "
                  f"σ@ref={sigma_h[ref_idx]:6.2f} MPa  util_max={rec['util_max']:6.2f}  "
                  f"折 {rec['n_fractured_h']:2d}/50  ⇒ status={rec['status']:<8s} "
                  f"OR={rec['OR']}  {'✅ 落窗' if rec['or_in_window'] else '❌ 未落窗'}",
                  flush=True)
        n_in = sum(r["or_in_window"] for r in rows.values())
        out["areas"][aname] = {"area_mm2": area, "rows": rows, "n_or_in_window": n_in}

    print("\n" + "=" * 88, flush=True)
    print("S1a 判读：", flush=True)
    best = None
    for aname, ad in out["areas"].items():
        print(f"  {aname}: {ad['n_or_in_window']}/{len(SHARES)} 个 share 的 OR 落窗")
        for k, r in ad["rows"].items():
            if r["or_in_window"]:
                print(f"      share={k} OR={r['OR']} ci={r['ci95']} "
                      f"A3={r['a3']['n_within']}/3 other8={r['unchanged_other8']}")
                if best is None or abs(float(r["OR"]) - target) < abs(
                        float(best[1]["OR"]) - target):
                    best = (aname, r)
    if best:
        aname, r = best
        print(f"\n  ⇒ **DoD-2 可达**：最接近靶的组合 = {aname} @ share={r['share']:.4f}"
              f" ⇒ OR={r['OR']}（靶 {target}，窗口 {[round(x,3) for x in win]}）")
        print(f"  ⇒ **两体 FE 降级为验证手段**：出口在解析层已打开；FE 只需确认"
              f"耦合下 share 能落在 {r['share']:.4f} 附近、且 A 用实测 A_eff")
    else:
        print("\n  ⇒ **零 FE 层面 DoD-2 仍不可达** ⇒ 此时才需要两体 FE"
              "（要查的是闭式 share 本身在耦合下是否成立）")
    out["verdict"] = {"doD2_reachable_in_analytic": best is not None,
                      "best": ({"area": best[0], **{k: best[1][k] for k in
                                ("share", "OR", "ci95", "status")}} if best else None)}

    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    _report(out)
    print(f"\nwrote: {OUT_JSON}\nwrote: {OUT_REPORT}", flush=True)
    return 0


def _report(out: dict) -> None:
    m, c0 = out["meta"], out["positive_control"]
    L: list[str] = []
    A = L.append
    A("# M1 全量 · S1a 解析预检 —— share→σ(h)→G3 OR")
    A("")
    A(f"> 生成 {m['generated_at']}；模型：{m['model']}")
    A(f"> **正对照 C0**：G3 冻结原行逐部位复现已存 `g3_or_estimates.json` 的 "
      f"(status, n_fracture, OR) ⇒ {c0['n_sites']} 部位中 "
      f"{c0['n_mismatch']} 个不一致 ⇒ "
      f"{'✅ 通过' if c0['agrees'] else '❌ 失败'}"
      f"（不过则不继续）。{c0.get('note','')}")
    A(f"> 腓骨靶 OR={m['fibula_target_or']}，±30% ⇒ 窗口 "
      f"{[round(x,3) for x in m['fibula_window']]}；σ_c={m['sigma_c_mpa']} MPa。")
    A("")
    A("## 0. 结论")
    A("")
    v = out["verdict"]
    if v["doD2_reachable_in_analytic"]:
        b = v["best"]
        A(f"- **DoD-2 在解析层可达**：最佳组合 **{b['area']} @ share={b['share']:.4f}** ⇒ "
          f"OR=**{b['OR']}**（status={b['status']}）。")
        A(f"- ⇒ **两体 FE 由「出口」降级为「验证手段」**：出口（OR 落窗）已由解析 share 打开；"
          f"FE 的职责收窄为确认耦合下 share ≈ {b['share']:.2f} 且用实测 A_eff 而非 median。")
    else:
        A("- **零 FE 层面 DoD-2 仍不可达** ⇒ 两体 FE 仍必需（要查闭式 share 在耦合下是否成立）。")
    A("")
    A("## 1. share × 面积口径 × OR")
    A("")
    for aname, ad in out["areas"].items():
        A(f"### {aname}（A={ad['area_mm2']:.2f} mm²）")
        A("")
        A("| share | alloc@5m (N) | σ@ref (MPa) | util_max | 折 h | status | OR | CI95 | 落窗 | A3 |")
        A("|---:|---:|---:|---:|---:|---|---:|---|:--:|---|")
        for k, r in ad["rows"].items():
            ci = r["ci95"] or [None, None]
            A(f"| {r['share']:.4f} | {r['alloc_at_ref_n']:.0f} | "
              f"{r['sigma_h_ref_mpa']:.2f} | {r['util_max']:.2f} | "
              f"{r['n_fractured_h']}/50 | {r['status']} | {r['OR']} | "
              f"{('%.3f–%.3f' % (ci[0], ci[1])) if ci[0] else '—'} | "
              f"{'✅' if r['or_in_window'] else '❌'} | "
              f"{r['a3']['n_within']}/3 |")
        A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())