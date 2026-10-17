"""A 件探针（FP_CONTACT 预注册）—— `docs/methods/P3_足部接触柔度_预注册.md`

为什么
------
`LEG_K_ANCHOR` 判定 `J1b` 三档全 FAIL（`F_pad`/BW 超实测靶 105–225%），诊断 =
「刚性接触尖峰伪影」。本件检验该诊断：把 `pad.py` 的**有效接触刚度**
（`hard_surface`：`k_eff = k_n × A(κ)`，实测 **548 kN/m**）换成**实测跟垫**
（Bennett & Ker 1990 / Ker 1996 in-vivo **100–175 N/mm**）⇒ `k_n ∈ {3.0e6, 2.1e6, 1.6e6}`。

**实现只用既有门 `RIGID_CONTACT["k_n_per_m"]`（`pad.py` 零改动）**。

正对照（A-J0，不过则**整份作废**）
----------------------------------
`k_n = None`（门关）必须**逐位**复现已提交的 `results/fp_anchor_probe.json`
（三档 @0.42 m 的 `F_pad` / `F_flex`）。这是对**已提交产物**的逐位核对。

判据（§5，跑前冻结）
--------------------
A-J1 `F_pad`/BW ÷ 实测靶 ∈ [0.70,1.30]（判在锚带中点 `k_n=2.1e6`，实测 `k_eff`≈128）；
A-J2 `F_pad/F_flex` 随 `k_eff` 单调不增 **且** 锚带上端三档 ≤ 2.50；
A-J3 `F_flex` 如实报告（无阈值）；A-J5 ⛔ 禁把 `k_n` 调到锚带**外**求判据变绿。

用法::

    .venv\\Scripts\\python.exe scripts/1d/fp_contact_probe.py
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.pad import RIGID_CONTACT, simulate_boulder_fall  # noqa: E402

MASS_KG = 80.0
BW = MASS_KG * 9.80665
H_MAIN = 0.42
LABELS = ("soft", "self", "stiff")
FP_TARGET_BW = {"soft": 3.38, "self": 3.70, "stiff": 5.36}
BAND = (0.70, 1.30)
RATIO_MAX = 2.50
#: 锚带中点那一档 = A-J1 的判定档
JUDGE_KN = 2.1e6
#: 实测跟垫锚（in-vivo）
ANCHOR_K_EFF_N_PER_M = (1.0e5, 1.75e5)

#: 冻结的 `k_n` 网格（§3/§7-2）
KS = (3.0e6, 2.1e6, 1.6e6)

UPSTREAM = ROOT / "results" / "fp_anchor_probe.json"
OUT = ROOT / "results" / "fp_contact_probe.json"


def _run(kn: float | None, label: str) -> dict:
    prev = RIGID_CONTACT["k_n_per_m"]
    try:
        RIGID_CONTACT["k_n_per_m"] = kn
        r = simulate_boulder_fall(height_m=H_MAIN, mass_kg=MASS_KG, on_pad=False,
                                  posture=f"anchored-{label}")
    finally:
        RIGID_CONTACT["k_n_per_m"] = prev
    fp = float(r.peak_force_n)
    ff = float(r.meta["f_flex_peak_n"])
    x = float(r.max_compression_m)
    return {
        "k_n": kn,
        "label": label,
        "f_pad_n": fp, "f_flex_n": ff, "x_m": x,
        "k_eff_n_per_m": fp / max(x, 1e-12),
        "f_pad_bw": fp / BW, "f_flex_bw": ff / BW,
        "ratio_pad_over_flex": fp / ff,
        "J1b_ratio": (fp / BW) / FP_TARGET_BW[label],
        "J1a_ratio": (ff / BW) / FP_TARGET_BW[label],
    }


def positive_control() -> tuple[bool, list[dict]]:
    """A-J0：门关时必须逐位复现已提交的 `fp_anchor_probe.json`。"""
    up = json.loads(UPSTREAM.read_text(encoding="utf-8"))
    rows, ok = [], True
    for label in LABELS:
        want = up["main_rows"][label]
        got = _run(None, label)
        good = (abs(got["f_pad_n"] - want["f_pad_n"]) < 1e-6
                and abs(got["f_flex_n"] - want["f_flex_n"]) < 1e-6)
        ok &= good
        rows.append({
            "label": label, "pass": good,
            "want": {"f_pad_n": want["f_pad_n"], "f_flex_n": want["f_flex_n"]},
            "got": {"f_pad_n": got["f_pad_n"], "f_flex_n": got["f_flex_n"]},
        })
    return ok, rows


def main() -> int:
    print("=" * 84)
    print("A-J0  POSITIVE CONTROL  (gate off must equal committed fp_anchor_probe.json)")
    print("=" * 84)
    pcok, pcrows = positive_control()
    for r in pcrows:
        print(f"  {r['label']:<6} want F_pad={r['want']['f_pad_n']/1e3:.6f} "
              f"F_flex={r['want']['f_flex_n']/1e3:.6f} | got "
              f"{r['got']['f_pad_n']/1e3:.6f} / {r['got']['f_flex_n']/1e3:.6f} "
              f"{'PASS' if r['pass'] else 'FAIL'}")
    if not pcok:
        print("!! POSITIVE CONTROL FAILED -- results INVALID, do not use")
        OUT.write_text(json.dumps({"positive_control": False, "rows": pcrows},
                                  ensure_ascii=False, indent=2), encoding="utf-8")
        return 1
    print("=> POSITIVE CONTROL PASS")

    print()
    print("=" * 84)
    print(f"O1-O4  @{H_MAIN} m hard ground; contact layer k_n x 3 anchored postures")
    print("=" * 84)
    print(f"  {'k_n':>9} {'k_eff(kN/m)':>12} | " + " | ".join(f"{l:^24}" for l in LABELS))
    print(f"  {'':>9} {'':>12} | " + " | ".join(
        f"{'F_pad/BW F_flex/BW ratio':^24}" for _ in LABELS))
    print("  " + "-" * (26 + 27 * 3))
    grid = {}
    for kn in KS:
        cells = []
        for label in LABELS:
            d = _run(kn, label)
            grid.setdefault(label, {})[f"{kn:.1e}"] = d
            cells.append(f"{d['f_pad_bw']:>7.2f} {d['f_flex_bw']:>8.2f} "
                         f"{d['ratio_pad_over_flex']:>6.2f}")
        keff = next(iter(grid[LABELS[0]][f"{kn:.1e}"].values())) if False else None
        keff = grid[LABELS[0]][f"{kn:.1e}"]["k_eff_n_per_m"]
        print(f"  {kn:>9.1e} {keff/1e3:>12.1f} | " + " | ".join(f"{c:^24}" for c in cells))

    print()
    print("=" * 84)
    print("VERDICT (pre-registered, honest binary)")
    print("=" * 84)
    judge = {label: grid[label][f"{JUDGE_KN:.1e}"] for label in LABELS}
    print(f"  A-J1 @ k_n={JUDGE_KN:.1e} (k_eff={judge['soft']['k_eff_n_per_m']/1e3:.1f} kN/m, "
          f"anchor band {ANCHOR_K_EFF_N_PER_M[0]/1e3:.0f}-{ANCHOR_K_EFF_N_PER_M[1]/1e3:.0f}):")
    a_j1 = {label: (BAND[0] <= judge[label]["J1b_ratio"] <= BAND[1]) for label in LABELS}
    for label in LABELS:
        print(f"      {label:<6} F_pad/BW={judge[label]['f_pad_bw']:>6.2f}  "
              f"ratio={judge[label]['J1b_ratio']:>5.3f}  "
              f"{'PASS' if a_j1[label] else 'FAIL'}")
    a_j1_all = all(a_j1.values())
    print(f"    A-J1 (三档) = {'PASS' if a_j1_all else 'FAIL'}")

    # A-J2: 单调（k_n 递减 ⇒ k_eff 递减 ⇒ ratio 不增）+ 锚带上端 max ratio
    keys = [f"{kn:.1e}" for kn in KS]
    mono = {label: all(grid[label][keys[i]]["ratio_pad_over_flex"]
                       >= grid[label][keys[i + 1]]["ratio_pad_over_flex"] - 1e-9
                       for i in range(len(keys) - 1)) for label in LABELS}
    top = {label: grid[label][keys[0]]["ratio_pad_over_flex"] for label in LABELS}
    a_j2 = all(mono.values()) and all(v <= RATIO_MAX for v in top.values())
    print(f"  A-J2 monotone per label: " + "  ".join(
        f"{l}={'OK' if mono[l] else 'NOT'}" for l in LABELS)
        + f" | top-of-band ratio " + "  ".join(f"{l}={top[l]:.2f}" for l in LABELS)
        + f" (<= {RATIO_MAX})")
    print(f"    A-J2 = {'PASS' if a_j2 else 'FAIL'}")

    print(f"  A-J3 F_flex (report only, no threshold):")
    for label in LABELS:
        up = json.loads(UPSTREAM.read_text(encoding="utf-8"))["main_rows"][label]
        d = judge[label]
        print(f"      {label:<6} gate-off {up['f_flex_bw']:>5.2f} BW -> "
              f"gate-on {d['f_flex_bw']:>5.2f} BW  (J1a {up['J1a_ratio']:.3f} -> "
              f"{d['J1a_ratio']:.3f})")

    payload = {
        "positive_control": True, "positive_control_rows": pcrows,
        "mass_kg": MASS_KG, "h_main_m": H_MAIN, "labels": list(LABELS),
        "fp_target_bw": FP_TARGET_BW, "band": list(BAND),
        "anchor_k_eff_n_per_m": list(ANCHOR_K_EFF_N_PER_M),
        "k_n_grid": list(KS), "judge_k_n": JUDGE_KN, "ratio_max": RATIO_MAX,
        "semantics_note": ("hard_surface 的 k_n 被乘面积：k_eff = k_n * A(kappa)；"
                           "主口径是实测 k_eff = F_pad/max_compression_m。"),
        "grid": grid,
        "verdict": {"A-J1": a_j1, "A-J1_all": a_j1_all,
                    "A-J2": a_j2, "A-J2_monotone": mono, "A-J2_top_ratio": top},
    }
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"OUT = {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
