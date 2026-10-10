"""核对：论文 §4.6末 / 出处索引 §5.1 里写的每个数字是否与 JSON 一致。

任何一条不一致 => FAIL（禁止把抄错的数字留在论文里）。
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from climbing.pad import CrashPad, hard_surface  # noqa: E402

FAIL: list[str] = []


def chk(label: str, got, want, tol: float = 0.051) -> None:
    if isinstance(want, list) and isinstance(got, list) and want and isinstance(want[0], float):
        ok = (len(got) == len(want)
              and all(abs(float(a) - float(b)) <= tol for a, b in zip(got, want)))
    elif isinstance(want, float) and isinstance(got, (int, float)):
        ok = abs(float(got) - want) <= tol
    else:
        ok = got == want
    print(f"  {'OK  ' if ok else 'FAIL'}  {label:<46} got={got!r} want={want!r}")
    if not ok:
        FAIL.append(label)


def main() -> int:
    dec = json.loads((ROOT / "results/p5_cond_decompose.json").read_text("utf-8"))
    sw = json.loads((ROOT / "results/p5_pad_stiffness_sweep.json").read_text("utf-8"))
    h20 = json.loads((ROOT / "results/p5_pad_stiffness_sweep_h20.json").read_text("utf-8"))
    cost = json.loads((ROOT / "results/p5_v6v9_cost.json").read_text("utf-8"))

    print("=" * 74)
    print("(1) 两因子分解  -- paper: 22.2/22.2/21.8 | 31.4/25.6/20.4 | 10.8/10.2/9.4")
    print("=" * 74)
    v3 = {r["h_m"]: r for r in dec["v3"]}
    chk("rigid k200 @h2/3/5", [v3[h]["rigid_k200"] for h in (2.0, 3.0, 5.0)],
        [22.2, 22.2, 21.8])
    chk("pad k2500 @h2/3/5", [v3[h]["pad_k2500"] for h in (2.0, 3.0, 5.0)],
        [31.4, 25.6, 20.4])
    chk("rigid k2500 @h2/3/5", [v3[h]["rigid_k2500"] for h in (2.0, 3.0, 5.0)],
        [10.8, 10.2, 9.4])
    chk("cond_main_pct", [round(v3[h]["cond_share_pct"], 1) for h in (2.0, 3.0, 5.0)],
        [67.6, 57.2, 46.8])
    chk("k_main_pct", [round(v3[h]["k_share_pct"], 1) for h in (2.0, 3.0, 5.0)],
        [32.4, 42.8, 53.2])
    chk("interaction_ms", [v3[h]["interaction_ms"] for h in (2.0, 3.0, 5.0)],
        [-5.8, -3.8, -1.6])
    # k 在垫上饱和 -5.6 / 在刚上 -11.4
    chk("k effect on pad @h2", round(v3[2.0]["pad_k2500"] - v3[2.0]["pad_k200"], 1), -5.6)
    chk("k effect on rigid @h2", round(v3[2.0]["rigid_k2500"] - v3[2.0]["rigid_k200"], 1), -11.4)

    print()
    print("=" * 74)
    print("(2) 垫刚度扫描 -- paper: 11.40 / 10.00 / 9.60 / 5.2 / 单调 6/6")
    print("=" * 74)
    rows = {(r["h_m"], r["alpha"]): r for r in sw["rows"]}
    chk("h5 a10 tpeak2500", rows[(5.0, 10.0)]["tpeak_ms_k2500"], 11.40)
    chk("h2 a30 tpeak2500", rows[(2.0, 30.0)]["tpeak_ms_k2500"], 10.00)
    chk("h3 a30 tpeak2500", rows[(3.0, 30.0)]["tpeak_ms_k2500"], 9.60)
    mono = all(m["t_torso_non_increasing"] and m["tpeak_non_increasing"]
               for m in sw["mono"].values())
    chk("monotone all heights", mono, True)
    r100 = next(r for r in h20["rows"] if r["alpha"] == 100.0)
    r300 = next(r for r in h20["rows"] if r["alpha"] == 300.0)
    chk("h20 a100 window_ms", r100["window_ms"], 5.2)
    chk("h20 a300 t_torso", r300["t_torso_ms"], 2.8)
    chk("h20 positive control (a1)", (next(r for r in h20["rows"]
        if r["alpha"] == 1.0)["t_torso_ms"]), 29.4)

    print()
    print("=" * 74)
    print("(3) 代价 -- paper: V6 9/9 + V9 3/3 @a30, util -15%~+27%")
    print("=" * 74)
    v = cost["verdict"]
    chk("V6 hits a30", v["30.0"]["v6_hits"], 9)
    chk("V9 hits a30", v["30.0"]["v9_hits"], 3)
    chk("V6 hits a1", v["1.0"]["v6_hits"], 9)
    chk("V9 hits a1", v["1.0"]["v9_hits"], 3)
    d = []
    for r in cost["v6"]:
        if r["alpha"] != 30.0:
            continue
        base = next(b for b in cost["v6"]
                    if b["alpha"] == 1.0 and b["h_m"] == r["h_m"]
                    and b["beta_deg"] == r["beta_deg"])
        d.append(100.0 * (r["util"] - base["util"]) / base["util"])
    lo, hi = round(min(d)), round(max(d))
    print(f"       util delta range = {lo}% .. {hi}%  (paper says -25%~+27%)")
    if not (-30 <= lo <= -20 and 20 <= hi <= 35):
        FAIL.append("util delta range")

    print()
    print("=" * 74)
    print("(4) '34x' claim -- alpha=300 foam vs hard_surface(1e7) stiffness ratio")
    print("=" * 74)
    pad300 = CrashPad(
        thickness_m=0.20,
        s_plateau_pa=CrashPad().s_plateau_pa * 300,
        s_dense_pa=CrashPad().s_dense_pa * 300,
        xi_plateau=0.15, n_dense=3.56, xi_dense_max=0.85,
        bottoming_modulus_pa=CrashPad().bottoming_modulus_pa * 300,
        linear_modulus_pa=0.0,
    )
    rigid = hard_surface(k_n_per_m=1.0e7)
    print("       x(m)   sigma_pad_a300(Pa)   sigma_rigid(Pa)   ratio")
    ratios = []
    for x in (0.002, 0.005, 0.01, 0.02):
        sp = float(pad300.stress_pa(x / pad300.thickness_m))
        sr = float(rigid.stress_pa(x / rigid.thickness_m))
        r = sp / sr
        ratios.append(r)
        print(f"       {x:<6}  {sp:>16.1f}  {sr:>16.1f}  {r:>6.1f}x")
    print(f"       => ratio range {min(ratios):.1f}x ~ {max(ratios):.1f}x")
    good = 25.0 <= min(ratios) and max(ratios) <= 34.5
    print("       paper says '26-34x' -> " + ("OK" if good else "MISMATCH"))
    if not good:
        FAIL.append("34x claim")

    print()
    print("=" * 74)
    print(f"RESULT  {'PASS' if not FAIL else 'FAIL: ' + repr(FAIL)}")
    print("=" * 74)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
