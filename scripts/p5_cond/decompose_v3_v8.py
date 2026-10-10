"""P5 前置侦察 · Step 1 —— V3/V8b 的「工况 vs k_neck」两因子分解。

输入：results/rigid_vtarget_probe.json（已有数据，**不重跑模型**）
输出：print + results/p5_cond_decompose.json

为什么做
--------
文档断言「V5/V3/V8b/踝拆分 同根因 = 工况错配 ⇒ 出路 = P5」。
但该 JSON 的 meta 显示 CALIB_K = {rigid: 2500, pad: 200} —— 两列**不仅
工况不同，k_neck 也不同**，所以 26.6 vs 8.6 ms 的差是**混合效应**。

本脚本做标准两因子分解：
  * Δ工况 = rigid − pad（在**同一** k_neck 下）
  * Δk     = k2500 − k200（在**同一**工况下）
  * interaction = Δ工况(k2500) − Δ工况(k200)

判据（可 falsify）
------------------
* **恢复率** = 该因子把 t_peak 拉回阈值的比例
    recover = (baseline − after) / (baseline − threshold)
* baseline = **pad @ k200**（论文发表口径 = 当前 ❌ 状态）
* threshold: V3 = 12.5 ms, V8b = 5.0 ms

⇒ 若「改工况」本身就能让判据 PASS ⇒ 工况错配说成立、P5 值钱
⇒ 若「改工况」只恢复部分、仍 FAIL ⇒ P5 单独救不了，须另找根因
"""
from __future__ import annotations

import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / "results" / "rigid_vtarget_probe.json"
OUT = ROOT / "results" / "p5_cond_decompose.json"

V3_THR = 12.5        # ms, head_branch.V3_THRESHOLD_S * 1000
V8B_THR = 5.0        # ms, head_branch.V8_THRESHOLD_AT_20M_S * 1000
V3_HEIGHTS = [1.0, 2.0, 3.0, 5.0, 8.0]   # V3_HEIGHTS（探针口径）


def load() -> dict:
    return json.loads(SRC.read_text(encoding="utf-8"))


def rows_by(vid: dict) -> dict:
    """{(surface, h): row}"""
    return {(r["surface"], float(r["h_m"])): r for r in vid["rows"]}


def decompose_v3(by: dict) -> list[dict]:
    """V3 用 tpeak_ms_k*（head 峰），两因子 = 工况 × k_neck。"""
    out = []
    for h in V3_HEIGHTS:
        # 垫侧只有 [2,3,5,10,20]，与 V3_HEIGHTS 取交
        base = by.get(("pad", h))
        rig = by.get(("rigid", h))
        if base is None or rig is None:
            continue
        row = {"h_m": h}
        for k in (200, 1000, 2500):
            key = f"tpeak_ms_k{k}"
            row[f"pad_k{k}"] = base[key]
            row[f"rigid_k{k}"] = rig[key]

        b, f = row["pad_k200"], row["pad_k2500"]
        r200, r2500 = row["rigid_k200"], row["rigid_k2500"]
        gap = b - V3_THR

        # 主效应（两水平均值法）+ 交互
        d_cond_k200 = r200 - b            # 只改工况
        d_cond_k2500 = r2500 - f          # 只改工况
        d_k_pad = f - b                   # 只改 k（垫上）
        d_k_rigid = r2500 - r200          # 只改 k（刚上）

        row.update(
            gap_to_thr_ms=gap,
            # 只改工况能恢复多少比例
            recover_cond_pct=100.0 * (-d_cond_k200) / gap,
            # 只改 k 能恢复多少比例
            recover_k_pct=100.0 * (-d_k_pad) / gap,
            # 两因子都改
            both_ms=r2500,
            both_pass=bool(r2500 <= V3_THR),
            cond_alone_pass=bool(r200 <= V3_THR),
            k_alone_pass=bool(f <= V3_THR),
            # 标准两因子分解（对总差 b→r2500）
            total_ms=b - r2500,
            cond_main_ms=-(d_cond_k200 + d_cond_k2500) / 2.0,
            k_main_ms=-(d_k_pad + d_k_rigid) / 2.0,
            interaction_ms=d_cond_k2500 - d_cond_k200,
        )
        tot = row["total_ms"]
        row["cond_share_pct"] = 100.0 * row["cond_main_ms"] / tot
        row["k_share_pct"] = 100.0 * row["k_main_ms"] / tot
        out.append(row)
    return out


def decompose_v8b(by: dict) -> dict:
    """V8b 用 t_torso_ms @20 m（f_flex 峰，与 k_neck 无关）。"""
    h = 20.0
    b = by[("pad", h)]["t_torso_ms"]
    r = by[("rigid", h)]["t_torso_ms"]
    gap = b - V8B_THR
    # 次要观测量（t_moment）一并记账
    b2 = by[("pad", h)].get("t_moment_ms")
    r2 = by[("rigid", h)].get("t_moment_ms")
    return {
        "h_m": h,
        "threshold_ms": V8B_THR,
        "pad_ms": b,
        "rigid_ms": r,
        "gap_to_thr_ms": gap,
        "recover_cond_ms": b - r,
        "recover_cond_pct": 100.0 * (b - r) / gap,
        "cond_alone_pass": bool(r <= V8B_THR),
        "residual_after_cond_ms": r - V8B_THR,
        # secondary observable
        "secondary_tmoment_pad_ms": b2,
        "secondary_tmoment_rigid_ms": r2,
        "secondary_pass": bool((r2 if r2 is not None else 9e9) <= V8B_THR),
    }


def main() -> None:
    data = load()
    vid = data["V3_V8"]
    by = rows_by(vid)

    v3 = decompose_v3(by)
    v8b = decompose_v8b(by)

    print("=" * 78)
    print("P5 前置侦察 Step 1 -- V3 / V8b 两因子分解（工况 x k_neck）")
    print("=" * 78)
    print(f"baseline = pad @ k200 (paper), V3 threshold = {V3_THR} ms")
    print(f"input    = {SRC.name} (no model re-run)")
    print()

    print(f"{'h(m)':>5} {'pad_k200':>9} {'pad_k2500':>10} "
          f"{'rig_k200':>9} {'rig_k2500':>10} "
          f"{'gap':>6} {'cond%':>7} {'k%':>6} {'both':>7} {'PASS':>6}")
    print("-" * 78)
    for r in v3:
        print(f"{r['h_m']:>5.1f} {r['pad_k200']:>9.1f} {r['pad_k2500']:>10.1f} "
              f"{r['rigid_k200']:>9.1f} {r['rigid_k2500']:>10.1f} "
              f"{r['gap_to_thr_ms']:>6.1f} {r['recover_cond_pct']:>7.1f} "
              f"{r['recover_k_pct']:>6.1f} {r['both_ms']:>7.1f} "
              f"{'Y' if r['both_pass'] else 'N':>6}")
    print()
    print("  cond%  = 只改工况(pad->rigid, @k200) 恢复阈值缺口的比例")
    print("  k%     = 只改 k_neck(k200->k2500, @pad) 恢复阈值缺口的比例")
    print("  both   = 两因子都改后 t_peak；PASS = <= 12.5 ms")
    print()

    print("-" * 78)
    print("单因子能否 PASS：")
    for r in v3:
        print(f"  h={r['h_m']:.1f} m  工况单改={'PASS' if r['cond_alone_pass'] else 'FAIL'}"
              f"({r['rigid_k200']:.1f} ms)   "
              f"k 单改={'PASS' if r['k_alone_pass'] else 'FAIL'}"
              f"({r['pad_k2500']:.1f} ms)   "
              f"双改={'PASS' if r['both_pass'] else 'FAIL'}({r['both_ms']:.1f} ms)")
    print()

    print("-" * 78)
    print("两因子贡献占比（总差 = pad_k200 -> rigid_k2500）:")
    print(f"{'h(m)':>5} {'total':>7} {'cond_ms':>8} {'cond%':>7} "
          f"{'k_ms':>7} {'k%':>6} {'inter_ms':>9}")
    print("-" * 78)
    for r in v3:
        print(f"{r['h_m']:>5.1f} {r['total_ms']:>7.1f} "
              f"{r['cond_main_ms']:>8.1f} {r['cond_share_pct']:>7.1f} "
              f"{r['k_main_ms']:>7.1f} {r['k_share_pct']:>6.1f} "
              f"{r['interaction_ms']:>9.2f}")
    print()

    print("=" * 78)
    print("V8b (@20 m, t_torso, threshold = 5.0 ms)")
    print("=" * 78)
    for k in ("pad_ms", "rigid_ms", "gap_to_thr_ms", "recover_cond_ms",
              "recover_cond_pct", "residual_after_cond_ms",
              "secondary_tmoment_pad_ms", "secondary_tmoment_rigid_ms"):
        print(f"  {k:<28} {v8b[k]}")
    print(f"  {'cond_alone_pass':<28} {v8b['cond_alone_pass']}")
    print(f"  {'secondary_pass':<28} {v8b['secondary_pass']}")

    payload = {"v3": v3, "v8b": v8b,
               "method": "two-factor from existing rigid_vtarget_probe.json; "
                         "no model re-run",
               "baseline": "pad @ k200 (paper published)",
               "thresholds_ms": {"v3": V3_THR, "v8b": V8B_THR}}
    OUT.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print()
    print(f"OUT = {OUT}")


if __name__ == "__main__":
    main()
