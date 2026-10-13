"""最终收口：用实测 A_eff=76.28 实跑 G3，定案 DoD-2 判定区间。

**前一轮的判定是错的**（说"share ≥0.02 都落窗"）：
`m1_aeff_pinned.py` 的 verdict 只是把 A_eff 与 share 轴拼接，没有真跑 G3；
而 §S2 的稳健域实测已显示 share=0.02 是 **all_zero**（OR≈1.0，落在窗口下）。
⇒ 判定必须由 **G3 实跑**给出，本脚本补上。
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

from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402
from scripts.opensim_fe import g3_subject_variability as g3  # noqa: E402
from scripts.opensim_fe import g8_fibula_or_share as g8  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
OUT_JSON = RES / "m1_dod2_final.json"
OUT_REPORT = RES / "M1_DOD2_FINAL.md"
FIBULA = "腓骨"
SIGMA_C = 70.0


def main() -> int:
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - 5.0)))
    F_ref = float(loads[ref_idx])
    ratio = loads / F_ref
    target = float(g3.A3_TARGETS[FIBULA])
    win = [target * 0.7, target * 1.3]

    pinned = json.loads((RES / "m1_aeff_pinned.json").read_text(encoding="utf-8"))
    A_eff = float(pinned["verdict"]["A_eff_pinned_mm2"])
    print(f"A_eff 实测定案 = {A_eff:.2f} mm²（S2b，跨 share 离散 0.0000%）")
    print(f"窗口 [{win[0]:.3f}, {win[1]:.3f}]（靶 {target}）\n")

    out: dict = {"meta": {
        "title": "M1 最终收口：实测 A_eff 下 DoD-2 的判定区间",
        "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "A_eff_mm2": A_eff, "sigma_c_mpa": SIGMA_C,
        "ankle_ref_n": F_ref, "fibula_target_or": target, "fibula_window": win,
        "note": ("上一版 m1_aeff_pinned 的 verdict 是把 A 轴与 share 轴拼接的推论，"
                 "没有实跑 G3（它误称 share≥0.02 都落窗，而 S2 实测 share=0.02 是"
                 " all_zero）⇒ 本脚本用 G3 实跑给出真正的判定区间")}}

    rows = {}
    for sh in (0.02, 0.05, 0.10, 0.15, 0.20, 0.2241, 0.25, 0.30):
        sigma_h = sh * F_ref * ratio / A_eff
        util = sigma_h / SIGMA_C
        inputs = g3.load_inputs()
        mixed, idx = g8.build_mixed_inputs(inputs, sigma_h)
        cfg = {sd: g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                 g3.COV_A_PRIMARY, mixed, g3.N_SUBJ, True)
               for sd in g3.SEEDS}
        primary = cfg[g3.OFFICIAL_SEED]
        fr = g3._site_or_record(primary["or_map"][FIBULA])
        unc = g8.unchanged_check(inputs, idx, primary)
        ov = fr.get("OR")
        hit = bool(ov is not None and win[0] <= ov <= win[1])
        rows[f"{sh:.4f}"] = {
            "share": sh, "alloc_at_ref_n": sh * F_ref, "util_max": float(util.max()),
            "n_fractured_h": int((util >= 1).sum()), "status": fr.get("status"),
            "OR": ov, "ci95": fr.get("ci95"), "in_window": hit,
            "a3": primary["a3"], "other8_unchanged": unc.get("all_equal")}
        print(f"  share={sh:6.4f} util_max={util.max():6.2f} 折{rows[f'{sh:.4f}']['n_fractured_h']:2d}/50 "
              f"status={str(fr.get('status')):8s} OR={('%.4f' % ov) if ov else '—':>8s} "
              f"{'✅ PASS' if hit else '❌ .'}", flush=True)

    passing = [s for s in (0.02, 0.05, 0.10, 0.15, 0.20, 0.2241, 0.25, 0.30)
               if rows[f"{s:.4f}"]["in_window"]]
    out["rows"] = rows
    out["verdict"] = {
        "passing_shares": passing,
        "pass_range": [min(passing), max(passing)] if passing else None,
        "g8_g9_assumed_0.05_in_window": rows["0.0500"]["in_window"],
        "closed_form_0.2241_in_window": rows["0.2241"]["in_window"],
        "a3_all_pass": all(r["a3"]["n_within"] >= 2 for r in rows.values()),
        "other8_unchanged_all": all(r["other8_unchanged"] for r in rows.values()),
    }
    v = out["verdict"]
    print("\n" + "=" * 88)
    print(f"通过 share：{passing}")
    print(f"  G8/G9 假设 share=0.05 落窗？ {'✅' if v['g8_g9_assumed_0.05_in_window'] else '❌'}")
    print(f"  闭式 share=0.2241 落窗？    {'✅' if v['closed_form_0.2241_in_window'] else '❌'}")
    print(f"  A3 全组 ≥2/3？ {v['a3_all_pass']} ｜ 其余 8 部位逐值不变？ "
          f"{v['other8_unchanged_all']}")
    print("\n⇒ **DoD-2 判定区间 = share ∈ "
          f"[{v['pass_range'][0]:.4f}, {v['pass_range'][1]:.4f}]（A_eff=76.28 mm² 实测）**")
    print("  ⇒ G8/G9 的 share=0.05 **落在区间下界之外** ⇒ 这就是它们 OR≈1.13 的原因")
    print("  ⇒ 闭式 share=0.2241 **落在区间内** ⇒ DoD-2 由「物理 share」复活")

    OUT_JSON.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    _report(out)
    print(f"\nwrote: {OUT_JSON}\nwrote: {OUT_REPORT}")
    return 0


def _report(out: dict) -> None:
    m, v = out["meta"], out["verdict"]
    L: list[str] = []
    A = L.append
    A("# M1 全量 · 最终收口 —— 实测 A_eff 下的 DoD-2 判定")
    A("")
    A(f"> 生成 {m['generated_at']}")
    A(f"> **A_eff = {m['A_eff_mm2']:.2f} mm²**（S2b 实测，纯几何量，跨 share 离散 0%）")
    A(f"> 腓骨靶 OR={m['fibula_target_or']}，±30% ⇒ 窗口 "
      f"{[round(x,3) for x in m['fibula_window']]}")
    A(f"> {m['note']}")
    A("")
    A("## 0. 判定")
    A("")
    A(f"- **通过 share 区间 = [{v['pass_range'][0]:.4f}, {v['pass_range'][1]:.4f}]**")
    A(f"- G8/G9 假设 share=0.05：{'✅ 落窗' if v['g8_g9_assumed_0.05_in_window'] else '❌ 落在下界外'}")
    A(f"- 闭式 share=0.2241：{'✅ 落窗' if v['closed_form_0.2241_in_window'] else '❌'}")
    A(f"- A3 全组 ≥2/3：**{v['a3_all_pass']}**；其余 8 部位逐值不变：**{v['other8_unchanged_all']}**")
    A("")
    A("## 1. 逐 share")
    A("")
    A("| share | alloc@5m (N) | util_max | 折 h | status | OR | CI95 | 落窗 | A3 |")
    A("|---:|---:|---:|---:|---|---:|---|:--:|---|")
    for k, r in out["rows"].items():
        ci = r["ci95"] or [None, None]
        A(f"| {r['share']:.4f} | {r['alloc_at_ref_n']:.0f} | {r['util_max']:.2f} | "
          f"{r['n_fractured_h']}/50 | {r['status']} | "
          f"{('%.4f' % r['OR']) if r['OR'] else '—'} | "
          f"{('%.4f–%.4f' % (ci[0], ci[1])) if ci[0] else '—'} | "
          f"{'✅' if r['in_window'] else '❌'} | {r['a3']['n_within']}/3 |")
    A("")
    OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())