"""M1 全量 S0 前置：解析扫描 —— share≈22% 到底能不能给出"渐进"（部分 h 折、非全断）？

**为什么先做这个**：M1 全量要建两体耦合 FE（成本以小时计）。但 M1 的目的只有一个
——把 G9 的 DoD-2（OR 落不进 [1.177,2.187]）从"share 假设"层面移到"物理"层面。
若在**解析层面**（严格线弹性 σ(h)=σ_ref·load(h)/load_ref + 实测 A_eff），
share≈22% 仍然给不出渐进形态，那 M1 全量建了也不会改变 DoD-2 ⇒ **止损**。
反之则 M1 值得投。⇒ 这是"先 spike 再上马"纪律在 M1 上的应用。

**关键歧义的正面处理**：G9 用 median(fib_cort) 得出的"等效承载面积"只有 ~33 mm²，
而 M1-spike 实测 **A_cort=101.9 mm²**、G9 载荷流正对照实测 **A_eff=76.01 mm²**。
三者不一致 ⇒ 本脚本**同时**对 3 个面积口径扫描，而不是挑一个顺眼的用。
（面积口径不同 ⇒ 同一 share 落在渐进窗口的不同位置；这本身就是必须入账的敏感性。）
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fracture_matrix as fm  # noqa: E402
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
SIGMA_C = 70.0

#: 三个面积口径（mm²），来源各不相同 —— 全都扫，不挑
AREAS = {
    "A_cort_measured_101.9": 101.9,   # m1_spike: ∑V/L，实测皮质域
    "A_eff_loadflow_76.01": 76.01,    # g9 §11.7 载流积分正对照
    "A_equiv_median_33": 33.0,        # g9 median(fib_cort) 反推的等效承载面积
}


def row_stats(util: np.ndarray) -> dict:
    over = int((util >= 1.0).sum())
    return {"util_min": float(util.min()), "util_max": float(util.max()),
            "util_at_ref": None, "n_fractured_h": over,
            "frac_fractured_h": float(over / util.size),
            "first_fracture_h": next((i + 1 for i, v in enumerate(util) if v >= 1.0), None),
            "all_one": bool(over == util.size), "all_zero": bool(over == 0),
            "shape": ("all_zero" if over == 0 else
                      "all_one" if over == util.size else "progressive")}


def main() -> int:
    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    F_ref = float(loads[ref_idx])
    ratio = loads / F_ref

    print(f"ankle_r @REF={fm.REF_H:.0f}m = {F_ref:.1f} N")
    print(f"load(h) 形状: min={loads.min():.1f} @h={heights[loads.argmin()]}  "
          f"max={loads.max():.1f} @h={heights[loads.argmax()]}  "
          f"ratio 范围 [{ratio.min():.4f}, {ratio.max():.4f}]")

    out: dict = {"meta": {
        "title": "M1 S0 前置：share→σ(h) 解析渐进性扫描（严格线弹性）",
        "ref_h_m": float(fm.REF_H), "ankle_ref_n": F_ref,
        "sigma_c_mpa": SIGMA_C,
        "areas_mm2": AREAS,
        "share_hypotheses": {
            "g8_g9_assumed": 0.05,
            "literature_max": 0.30,
            "m1_spike_closed_form": None,  # 下面从 json 读
        },
        "model": ("σ(h) = σ_ref·load(h)/load_ref, σ_ref = share·F_ref/(A·1.0) "
                  "⇒ util(h) = share·F_ref·load(h)/(F_ref·A·σ_c) = share·load(h)/(A·σ_c)"
                  " —— 面积口径与 share 解耦，先扫面积"),
    }, "areas": {}}

    spike = json.loads((RES / "m1_spike_share_estimate.json").read_text(encoding="utf-8"))
    share_cf = float(spike["share_parallel"])
    out["meta"]["share_hypotheses"]["m1_spike_closed_form"] = share_cf
    print(f"m1_spike 闭式 share = {share_cf:.4f}")

    # 渐进窗口：对每个面积口径解 share 的解析区间
    shares = sorted({0.01, 0.02, 0.05, 0.10, share_cf, 0.25, 0.282, 0.30, 0.40, 0.50})
    for aname, area in AREAS.items():
        print("\n" + "=" * 84)
        print(f"[{aname}] A={area:.2f} mm²  ⇒ 压溃门槛 share = A·σ_c/F_ref = "
              f"{area * SIGMA_C / F_ref:.4f}")
        rows = {}
        for sh in shares:
            util = sh * loads / (area * SIGMA_C)
            st = row_stats(util)
            st["util_at_ref"] = float(util[ref_idx])
            st["share"] = float(sh)
            st["alloc_at_ref_n"] = float(sh * F_ref)
            rows[f"{sh:.4f}"] = st
            print(f"  share={sh:6.4f} alloc@{fm.REF_H:.0f}m={sh*F_ref:7.1f} N  "
                  f"util@ref={util[ref_idx]:6.3f}  util_max={util.max():6.3f}  "
                  f"折 {st['n_fractured_h']:2d}/{util.size} h  "
                  f"首折={str(st['first_fracture_h'] or '>' + str(util.size)):>3s}  "
                  f"⇒ {st['shape']}")
        # 解析渐进窗口（方向不可翻转，逐字推导）：
        #   util(h) = share*loads(h)/(A*sigma_c)  => util 单调正比于 loads(h)
        #   (i)  至少有一折: util_max >= 1  <=>  share*loads.max()/(A*s_c) >= 1
        #                          <=> share >= A*s_c/loads.max()   <- 下界来自分母【大】的
        #   (ii) 不是全断  : util_min <  1  <=>  share*loads.min()/(A*s_c) <  1
        #                          <=> share <  A*s_c/loads.min()   <- 上界来自分母【小】的
        #   => 窗口 = [A*s_c/loads.max(), A*s_c/loads.min())
        # v1 bug（负面对照）：把这两个界装进名字相反的变量、再反序比较，
        #   于是 2/3 的渐进窗口被判成"落在窗口外" => 给出错误的止损建议。
        share_lo = area * SIGMA_C / float(loads.max())
        share_hi = area * SIGMA_C / float(loads.min())
        prog = bool(share_lo < share_cf < share_hi)

        def _predict(sh: float) -> str:
            return ("progressive" if share_lo < sh < share_hi
                    else ("all_zero" if sh < share_lo else "all_one"))

        print(f"  => 渐进窗口 share in [{share_lo:.4f}, {share_hi:.4f})"
              f" (lo=A*s_c/loads.max={share_lo:.4f}, hi=A*s_c/loads.min={share_hi:.4f}); "
              f"closed-form share={share_cf:.4f} "
              f"{'IN-WINDOW (OK)' if prog else 'OUT-OF-WINDOW'}")
        bad = [(sh, st["shape"], _predict(sh)) for sh, st in
               ((float(k), v) for k, v in rows.items())
               if _predict(sh) != st["shape"]]
        print(f"  => 正对照（判据 vs 表内 n_fractured 形状）: "
              f"{f'OK {len(rows)}/{len(rows)} agree' if not bad else 'MISMATCH ' + str(bad)}")
        out["areas"][aname] = {"area_mm2": area,
                               "share_progressive_window": [float(share_lo), float(share_hi)],
                               "window_definition": "[A*sigma_c/loads.max, A*sigma_c/loads.min)",
                               "positive_control_rows_agree": not bad,
                               "closed_form_share_in_window": prog,
                               "rows": rows}

    # 汇总判读
    print("\n" + "=" * 84)
    print("S0 判读（决定 M1 全量是否值得建两体 FE）：")
    for aname, d in out["areas"].items():
        print(f"  {aname}: 闭式 share 在渐进窗口内？ "
              f"{'✅ 是' if d['closed_form_share_in_window'] else '❌ 否'}")
    n_in = sum(d["closed_form_share_in_window"] for d in out["areas"].values())
    print(f"\n  ⇒ {n_in}/3 个面积口径下 share≈22% 给出渐进形态")
    if n_in == 0:
        print("  ⇒ **止损建议：M1 全量（两体 FE）不会改变 DoD-2 ⇒ 不建**"
              "（σ(h)∝load(h) 是硬约束，与 FE 建模细节无关）")
    elif n_in < 3:
        print("  ⇒ **有条件投 M1 全量**：面积口径是主控，须在 M1 里用**载荷流实测 A_eff**"
              " 而非 median 口径；且 M1 的价值在于给出两体耦合下的真实 share（可 > 闭式）")
    else:
        print("  ⇒ **M1 全量值得投**：三个口径一致 ⇒ 渐进性由 share≈22% 决定，"
              "M1 只需让 FE 的 share 复现闭式值即可把 DoD-2 迁到物理层面")

    (RES / "m1_precheck_progressive.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nwrote: results/opensim_fe/m1_precheck_progressive.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())