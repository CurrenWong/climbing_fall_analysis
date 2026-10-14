# -*- coding: utf-8 -*-
"""I · 腓骨载荷分配（分担比轴）—— 判定 G9「OR 不敏感」vs M1「share=0.2241 ⇒ OR=1.7548」。

**性质**：执行既有预注册（`G8:85 F1/F2` + `G9:D2` 的分担比轴）+ **零 FE**（σ ∝ load 已由
G9 的线性自检 ratio=2.072/2.106 佐证）。判定规则 **零改动**：G3 的 `simulate`/`or_all_sites`/
`a3_verdict`/`logistic_fit` 全部原封调用。

**关键前提（须成立）**：冻结 BC2 的腓骨行是按 **share = 1.0**（整条 `ankle_r` 反力）喂的
—— 依据 `FIBULA_FE_AUDIT_REPORT.md`：FE 把整条踝反力同时喂给胫骨与腓骨、未做胫-腓载荷分配，
腓骨 σ_nom/σc ≈ 6–44×。⇒ `σ_fib(h; s) = s · σ_fib_frozen(h)`。

**红线**：分担比**全量报告、不挑**；只取 `OFFICIAL_SEED` 之外的多种子作稳健性；
**不回解 share 去凑 OR**（与 G8/G9 同一纪律）。
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

# 分担比候选：**全部有出处，标注 provenance**；不新增无源值
SHARES = [
    (0.05,   "G8 主设定 / G9 D2 主（G9 at 5% ⇒ all_zero）"),
    (0.10,   "G8 F5 敏感性 / G9 D2 敏感（OR=1.130）"),
    (0.2241, "★M1 定案闭式 share（两端栓接轴向并联，G9 §15 / M1 全量 ⇒ OR=1.7548）"),
    (0.30,   "G9 §15 通过区间上沿 share∈[0.10,0.30]"),
]
TARGET = 1.682
TOL = 0.30
WIN = (TARGET * (1 - TOL), TARGET * (1 + TOL))


def main() -> int:
    print("=== Inputs schema ===")
    print("  ", g3.Inputs.__annotations__)
    inputs = g3.load_inputs()
    print("=== loaded ===")
    for k, v in inputs.items():
        try:
            a = np.asarray(v)
            print(f"  {k:<14} shape={a.shape} dtype={a.dtype}")
        except Exception:
            print(f"  {k:<14} {type(v).__name__}: {str(v)[:70]}")

    sig = np.asarray(inputs["sigma_nom"], float)
    sites = list(pa.REGIONS[i][0] for i in range(len(pa.REGIONS)))
    # G3 内部用的 site 顺序
    order = list(getattr(g3, "GRADUAL", []))
    fib_cn = "腓骨"
    if fib_cn not in sites:
        print(f"!! 未在 REGIONS 找到 {fib_cn}；sites={sites}")
        return 2
    fi = sites.index(fib_cn)
    print(f"\n腓骨 index={fi}  sigma_fib(h) 首/末 = "
          f"{sig[fi, 0]:.3f} / {sig[fi, -1]:.3f} MPa")
    sc = np.asarray(inputs["sigma_c"], float)[fi]
    print(f"sigma_c(腓骨) = {sc} MPa  ⇒ share=1.0 时 sigma/sigma_c = "
          f"{sig[fi, 0] / sc:.2f} … {sig[fi, -1] / sc:.2f} "
          f"（对 6–44× 的独立复核）")

    # --- A2 口径（确定性 0/1）下的形态，随 share 变化 ---
    print("\n=== A2 口径（确定性 0/1 → logistic_fit）形态随 share ===")
    heights = np.asarray(inputs["heights"], float) \
        if "heights" in inputs else np.asarray(pa.__dict__.get("_H", []), float)
    if heights.size != sig.shape[1]:
        heights = np.arange(1, sig.shape[1] + 1, dtype=float)
    rows = []
    for s, prov in SHARES:
        m = (sig[fi] * s >= sc).astype(int)
        st = pa.logistic_fit(heights, m)["status"]
        b = {"all_zero": "常值", "all_one": "常值",
             "complete_separation": "阶跃", "ok": "渐进"}[st]
        rows.append({"share": s, "prov": prov, "a2_behav": b, "a2_status": st,
                     "a2_n_frac": int(m.sum())})
        print(f"  share={s:<7} A2形态={b:<4} (status={st:<22}) "
              f"骨折格数={int(m.sum()):>2}/{m.size}")

    # --- A3/G3 口径：全量 share × 全种子 ---
    print(f"\n=== A3/G3 口径：OR（per-metre），靶 {TARGET} ±{TOL:.0%} ⇒ "
          f"[{WIN[0]:.3f}, {WIN[1]:.3f}] ===")
    print(f"  COV_B={g3.COV_B_PRIMARY} COV_E={g3.COV_E_PRIMARY} "
          f"COV_A={g3.COV_A_PRIMARY} n_subj={g3.N_SUBJ} seeds={g3.SEEDS}")
    or_rows = []
    for s, prov in SHARES:
        sig2 = sig.copy()
        sig2[fi] = sig[fi] * s
        inp2 = dict(inputs)
        inp2["sigma_nom"] = sig2
        per_seed = {}
        for sd in g3.SEEDS:
            cfg = g3.run_config(sd, g3.COV_B_PRIMARY, g3.COV_E_PRIMARY,
                                g3.COV_A_PRIMARY, inp2, n=g3.N_SUBJ)
            om = None
            if isinstance(cfg, dict):
                for k in ("or_by_site", "or_map", "or_all_sites", "or"):
                    if isinstance(cfg.get(k), dict):
                        om = cfg[k]
                        break
                if om is None:
                    print("  !! ConfigResult keys =", list(cfg.keys()))
            else:
                for k in ("or_by_site", "or_map"):
                    v = getattr(cfg, k, None)
                    if isinstance(v, dict):
                        om = v
                        break
            if om is None:
                raise SystemExit("找不到 OR map；请检查 ConfigResult 结构")
            rec = om.get(fib_cn) or om.get("fibula_r")
            if rec is None:
                raise SystemExit(f"OR map 里没有腓骨；keys={list(om.keys())}")
            per_seed[sd] = rec
        ors, sts = [], []
        for sd in g3.SEEDS:
            r = per_seed[sd]
            if isinstance(r, dict):
                o = r.get("or", r.get("odds_ratio", r.get("OR")))
                stt = r.get("status")
            else:
                o = getattr(r, "or", None)
                stt = getattr(r, "status", None)
            ors.append(o)
            sts.append(stt)
        good = [o for o in ors if o is not None and np.isfinite(o)]
        mean = float(np.mean(good)) if good else None
        hit = bool(mean is not None and WIN[0] <= mean <= WIN[1])
        print(f"  share={s:<7} OR={ors}  status={sts}  "
              f"均值={'None' if mean is None else f'{mean:.4f}'}  "
              f"{'★落窗' if hit else ('无OR' if mean is None else '未落窗')}")
        or_rows.append({"share": s, "prov": prov, "or_per_seed": ors,
                        "status_per_seed": sts, "or_mean": mean, "in_window": hit,
                        "window": list(WIN)})

    doc = {
        "piece": "I · 腓骨分担比轴（零 FE；G3/logistic_fit/a3_verdict 零改动）",
        "premise": "冻结腓骨行按 share=1.0 喂（整条 ankle_r 反力，无胫-腓载荷分配）"
                   "⇒ σ_fib(h;s) = s·σ_fib_frozen(h)",
        "target_or": TARGET, "tol": TOL, "window": list(WIN),
        "cov": {"b": g3.COV_B_PRIMARY, "e": g3.COV_E_PRIMARY,
                "a": g3.COV_A_PRIMARY},
        "n_subj": g3.N_SUBJ, "seeds": list(g3.SEEDS),
        "sigma_c_fibula": float(sc),
        "share_axis": or_rows,
        "a2_deterministic_behav": rows,
    }
    out = RES / "I_FIBULA_SHARE_AXIS.json"
    out.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n写出 {out.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())