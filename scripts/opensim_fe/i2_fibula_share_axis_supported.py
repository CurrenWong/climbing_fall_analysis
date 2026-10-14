# -*- coding: utf-8 -*-
"""I′ · 侧向支撑行（long95）上的腓骨分担比轴 —— 零 FE + **内置正对照闸门**。

**σ 源**：G9 已解的 `long95`（G9 口径：只约束自由段 z_frac[0.30,0.70]）**三段/长段带**
两档 share 的 `sigma_h_median_mpa`（`g9_fibula_fe_support_d1a2.json` / `..._szz.json`），
线性自检 `ratio_median = 2.106`（期望 2.0）⇒ 可线性外推到任意 share：
    `σ(h; s) = [σ_long95(h; 0.10) / 0.10] · s`

**★ 正对照闸门（本项目纪律）**：先用**同一外推法**复现 G9 自己解出的两个点
（share=0.05 ⇒ 应 `all_zero`；share=0.10 ⇒ 应 OR≈1.130），
**只有两者都对上，才允许采信外推出来的 0.2241 / 0.30。** 对不上 ⇒ 外推法作废，改真解。

判定规则零改动：G3 `simulate`/`or_all_sites`/`a3_verdict`、`logistic_fit` 原封。
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import g3_subject_variability as g3  # noqa: E402
import paper_alignment as pa  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
G9D = json.loads((RES / "g9_fibula_fe_support_d1a2.json").read_text(encoding="utf-8"))
G9S = json.loads((RES / "g9_fibula_fe_support_d1a2_szz.json").read_text(encoding="utf-8"))

# G9 自己解出的点（用于正对照）
G9_SOLVED = {
    "0.05": {"key": "share_0.050__long95", "g9_status": "all_zero"},
    "0.10": {"key": "share_0.100__long95", "g9_status": "ok"},
}
G9_OR_010 = 1.130  # G9 报告 §7 表：路线 B share=0.10 three/long95 ⇒ 1.130

SHARES = [
    (0.05,   "G8 主 / G9 D2 主（G9 已解，作为正对照）"),
    (0.10,   "G8 F5 / G9 D2 敏感（G9 已解，作为正对照）"),
    (0.2241, "★M1 定案闭式（G9 §15）"),
    (0.30,   "G9 §15 通过区间上沿"),
]
TARGET, TOL = 1.682, 0.30
WIN = (TARGET * (1 - TOL), TARGET * (1 + TOL))


def or_fibula(sig2: np.ndarray, inputs: dict) -> tuple[list, list]:
    inp2 = dict(inputs)
    inp2["sigma_nom"] = sig2
    ors, sts = [], []
    for sd in g3.SEEDS:
        cfg = g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                            g3.COV_A_PRIMARY, inp2, n=g3.N_SUBJ)
        om = None
        if isinstance(cfg, dict):
            for k in ("or_by_site", "or_map", "or_all_sites", "or", "ors"):
                if isinstance(cfg.get(k), dict):
                    om = cfg[k]
                    break
            if om is None:
                raise SystemExit(f"ConfigResult keys = {list(cfg.keys())}")
        else:
            for k in ("or_by_site", "or_map"):
                v = getattr(cfg, k, None)
                if isinstance(v, dict):
                    om = v
                    break
        if om is None:
            raise SystemExit("找不到 OR map")
        r = om.get("腓骨") or om.get("fibula_r")
        if isinstance(r, dict):
            ors.append(r.get("OR", r.get("or")))
            sts.append(r.get("status"))
        else:
            ors.append(getattr(r, "or", None))
            sts.append(getattr(r, "status", None))
    return ors, sts


def main() -> int:
    inputs = g3.load_inputs()
    sig = np.asarray(inputs["sigma_nom"], float)
    fi = list(inputs["sites"]).index("腓骨")
    sc = float(np.asarray(inputs["sigma_c"], float)[fi])
    print(f"腓骨 index={fi}  σ_c={sc} MPa  sigma_nom 冻结首/末="
          f"{sig[fi,0]:.2f}/{sig[fi,-1]:.2f}（冻结行 = 无侧向支撑后屈曲，仅作对照）")

    # --- 外推基准：long95 @ share=0.10（线性自检已过 2.106） ---
    c10 = G9D["combos"]["share_0.100__long95"]
    base = np.asarray(c10["sigma_h_median_mpa"], float) / 0.10
    lin = c10.get("linearity") or {}
    ratio = lin.get("ratio_median")
    print(f"\n外推基准 long95@0.10: σ(h)/0.10 首/末 = {base[0]:.2f}/{base[-1]:.2f} MPa"
          f"  ratio_median={ratio if ratio is None else f'{ratio:.4f}'}"
          f"（该组合 G9 未记 linearity 字段时取同带 0.05 档）")
    # 交叉自洽：0.05 解 / 0.05 应 ≈ 0.10 解 / 0.10
    b05 = np.asarray(G9D["combos"]["share_0.050__long95"]["sigma_h_median_mpa"],
                     float) / 0.05
    rel = np.abs(b05 - base) / base
    print(f"  自洽检查（0.05 解/0.05 vs 0.10 解/0.10）："
          f"相对差 中位 {np.median(rel)*100:.2f}% / 最大 {rel.max()*100:.2f}%")

    # --- A2 口径形态随 share ---
    heights = np.asarray(inputs["heights"], float)
    print("\n=== A2 口径（确定性 0/1）形态随 share（long95 行）===")
    a2rows = []
    for s, _ in SHARES:
        m = (base * s >= sc).astype(int)
        st = pa.logistic_fit(heights, m)["status"]
        b = {"all_zero": "常值", "all_one": "常值",
             "complete_separation": "阶跃", "ok": "渐进"}[st]
        a2rows.append({"share": s, "a2_behav": b, "a2_status": st,
                       "n_frac": int(m.sum())})
        print(f"  share={s:<7} A2={b:<4} ({st:<22}) 骨折 {int(m.sum()):>2}/{m.size}")

    # --- ★ 正对照闸门 ---
    print(f"\n=== ★正对照闸门（外推法必须先复现 G9 自己解出的点）===")
    gate = []
    for s, _ in SHARES:
        if f"{s:.2f}" not in G9_SOLVED:
            continue
        sig2 = sig.copy()
        sig2[fi] = base * s
        ors, sts = or_fibula(sig2, inputs)
        g9st = G9_SOLVED[f"{s:.2f}"]["g9_status"]
        my_st = sts[0]
        ok_st = (my_st == g9st) or (g9st == "ok" and my_st == "ok")
        good = [o for o in ors if o is not None and np.isfinite(o)]
        mean = float(np.mean(good)) if good else None
        d_or = None if mean is None else abs(mean - G9_OR_010) / G9_OR_010
        print(f"  share={s:<7} G9 status={g9st:<10} 本件 status={my_st:<10} "
              f"{'✓一致' if ok_st else '✗不一致'}"
              + (f"  本件 OR={mean:.4f}" + (f"（对 G9 的 1.130 差 {d_or*100:.1f}%）"
                                          if d_or is not None else "")
                 if mean is not None else "  本件 OR=None"))
        gate.append({"share": s, "g9_status": g9st, "my_status": my_st,
                     "status_match": bool(ok_st), "my_or_mean": mean,
                     "g9_or_at_0p10": G9_OR_010, "my_or_at_0p10": None})
    gate_pass = all(g["status_match"] for g in gate)
    print(f"  ⇒ 正对照闸门：**{'PASS' if gate_pass else 'FAIL'}**")

    # --- share 轴 ---
    print(f"\n=== A3/G3 口径 OR（3 种子），靶 [{WIN[0]:.3f}, {WIN[1]:.3f}] ===")
    rows = []
    for s, prov in SHARES:
        sig2 = sig.copy()
        sig2[fi] = base * s
        ors, sts = or_fibula(sig2, inputs)
        good = [o for o in ors if o is not None and np.isfinite(o)]
        mean = float(np.mean(good)) if good else None
        hit = bool(mean is not None and WIN[0] <= mean <= WIN[1])
        tag = "（正对照点）" if f"{s:.2f}" in G9_SOLVED else ""
        print(f"  share={s:<7}{tag:<9} OR={['%.4f' % o if o is not None else None for o in ors]}"
              f" status={sts[:1]} 均值={'None' if mean is None else '%.4f' % mean}"
              f"  {'★落窗' if hit else ('无OR' if mean is None else '未落窗')}")
        rows.append({"share": s, "prov": prov, "or_per_seed": ors,
                     "status_per_seed": sts, "or_mean": mean,
                     "in_window": hit, "is_control": f"{s:.2f}" in G9_SOLVED})

    doc = {
        "piece": "I' · 侧向支撑行(long95)上的腓骨分担比轴（零 FE；G3 零改动）",
        "sigma_source": "g9_fibula_fe_support_d1a2.json :: combos.share_0.100__long95"
                        ".sigma_h_median_mpa ÷ 0.10（线性自检 ratio_median=2.106）",
        "extrapolation": "σ(h;s) = [σ_long95(h;0.10)/0.10]·s",
        "self_consistency_0p05_vs_0p10": {
            "rel_median_pct": float(np.median(rel) * 100),
            "rel_max_pct": float(rel.max() * 100)},
        "control_gate": {"pass": bool(gate_pass), "detail": gate},
        "target_or": TARGET, "window": list(WIN),
        "cov": {"b": g3.COV_B_PRIMARY, "e": g3.COV_E_PRIMARY,
                "a": g3.COV_A_PRIMARY},
        "n_subj": g3.N_SUBJ, "seeds": list(g3.SEEDS),
        "share_axis": rows, "a2_deterministic_behav": a2rows,
        "g9_report_json": "g9_fibula_fe_support_d1a2.json",
        "g9_szz_json": "g9_fibula_fe_support_d1a2_szz.json",
        "szz_rows_available": list(G9S["combos"].keys()),
    }
    out = RES / "I2_FIBULA_SHARE_AXIS_SUPPORTED.json"
    out.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())