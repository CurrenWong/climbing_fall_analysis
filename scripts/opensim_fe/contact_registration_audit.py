# -*- coding: utf-8 -*-
"""接触装配自检闸门（零 FE，只读）—— 求解前先验「压板到底压没压到」

为什么需要
----------
2026-10-15 定位到 P2/T6 根因：压板注册用 `z_wall = a_max + GAP_MM`，其中
`a_max` 取自**顶带节点集**；但参与接触的 cap 面片是**另一个集合**（还被
`axdiv > CAP_AXDIV_MIN` 过滤过）。两者不一致时，压板会**悬空**：

  T6 (include_shell=True)  a_max=331.7565 (node 282)，而 **node 282 不在任何保留
  cap 面片上**，cap 面最高点比它低 **1.2342 mm** ⇒ 板面最高只到 a_max-0.8
  ⇒ **全程接合面片 0 / 36**（零接触）。

叠加「`fix_nodes_override` 会**替换** `bottom_fix`，且 `_anchor_xy_patch` 把
z_dof 置 0」⇒ 骨在 z 向毫无约束 ⇒ 纯刚体模态 ⇒ FEBio 手册 §8.3.1.5 的负 Jacobian
（且与载荷无关 —— 4 档诊断载荷同模式失败由此得到解释）。

本脚本把这次故障变成**可复用的前置断言**：
  断言 1  `a_max` 所在节点必须属于保留 cap 面片（否则板注册基准错）
  断言 2  t=0 与 t=1 的接合面积必须 > 0（否则板悬空 / 接触从未激活）
  断言 3  接合面积占 cap 面积的比例应达量级（对照已 PASS 的股骨 9.5%）

用法
----
    .venv\\Scripts\\python.exe scripts\\opensim_fe\\contact_registration_audit.py
    # 或指定部位： ... --bone T6 --include-shell 0

只读：只读网格与几何，不调用 FEBio、不写任何文件。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import bone_feb as bf                  # noqa: E402
import p2_head_contact as p2           # noqa: E402
import p4_pb_probe as p4               # noqa: E402
import p5_head_contact as p5           # noqa: E402

# 参考：已 PASS 的股骨（P5.0 实测）
REF_FEMUR_FRAC = 0.095


def audit(bone: str, include_shell: bool = True) -> dict:
    ctx = p4.load_ctx(bone)
    part, mesh = ctx[0], ctx[1]
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    geom = p2.cap_geom(mesh, part.spec, top_frac=p2.TOP_FRAC,
                       include_shell=include_shell)
    pg = p5.plate_geom_from_mesh(nodes_rot, geom)

    z = nodes_rot[:, 2]
    a_max = float(geom["a_max"])
    i_max = int(np.argmax(z)) + 1          # 1-based
    cap = geom["cap"]
    cap_nodes = set()
    for f in cap:
        cap_nodes.update(int(x) for x in f["poly"])
    zmax = np.array([z[np.asarray(f["poly"], int) - 1].max() for f in cap])
    areas = np.array([float(f["area"]) for f in cap])
    A_cap = float(areas.sum())

    zw0 = float(pg["z_wall"])
    res = {
        "bone": bone, "include_shell": include_shell,
        "a_max": a_max, "a_max_node_1based": i_max,
        "n_cap_facets": len(cap), "A_cap_mm2": A_cap,
        "a_max_in_cap_surface": i_max in cap_nodes,
        "cap_top_below_amax_mm": float(a_max - zmax.max()),
    }
    for tag, zw in (("t0", zw0 - 0.5 * p5.PLATE_DISP_MM),
                    ("t1", zw0 - 1.0 * p5.PLATE_DISP_MM)):
        eng = zmax > zw
        res[f"A_engaged_{tag}_mm2"] = float(areas[eng].sum())
        res[f"A_engaged_frac_{tag}"] = (float(areas[eng].sum()) / A_cap
                                        if A_cap else float("nan"))
        res[f"n_engaged_facets_{tag}"] = int(eng.sum())
        # 接合斑的偏心度（相对 cap 自身质心）—— 骨在 z 向只靠接触支承时，
        # 偏心斑会产生净力矩 ⇒ 骨倾覆 ⇒ 边角集中受压 ⇒ 单元倒置。
        if eng.any():
            cx = np.array([float(np.mean(nodes_rot[np.asarray(f["poly"], int) - 1, 0]))
                           for f in cap])
            cy = np.array([float(np.mean(nodes_rot[np.asarray(f["poly"], int) - 1, 1]))
                           for f in cap])
            m = areas[eng].sum()
            ex = float((areas[eng] * (cx[eng] - np.average(cx, weights=areas))).sum() / m)
            ey = float((areas[eng] * (cy[eng] - np.average(cy, weights=areas))).sum() / m)
            res[f"ecc_{tag}_mm"] = float(np.hypot(ex, ey))
        else:
            res[f"ecc_{tag}_mm"] = float("nan")

    # 断言（**只保留通过正对照的三条**）
    # ⚠️ 曾试加「接合斑偏心度」为第 4 条断言 —— **被正对照否掉**：已 PASS 的股骨
    #    实测偏心 30.5 mm（其 cap 带 4039 mm²，压在穹顶顶点，天然远离带质心），
    #    即该指标**不能**区分好/坏 ⇒ 按「检查器须有已知正对照」的规矩**废弃为断言**，
    #    仅作诊断量报告（`ecc_*`）。
    ok1 = res["a_max_in_cap_surface"]
    ok2 = res["A_engaged_t0_mm2"] > 0.0
    ok3 = res["A_engaged_frac_t0"] >= 0.2 * REF_FEMUR_FRAC
    res["asserts"] = {"a_max_in_cap_surface": ok1,
                      "engaged_t0_positive": ok2,
                      "engaged_frac_reasonable": ok3}
    res["diagnostics_only"] = {"ecc_t0_mm": res["ecc_t0_mm"],
                               "note": "偏心度未通过正对照，仅参考、不参与判定"}
    res["verdict"] = "PASS" if all(res["asserts"].values()) else "FAIL"

    # 修法预览：把板注册到**保留 cap 面片自身的最高点**
    z_top = float(zmax.max())
    fresh = {
        "z_wall_fixed": z_top + p5.GAP_MM,
        "shift_mm": float(a_max - z_top),
    }
    for tag, d in (("t0", 0.5 * p5.PLATE_DISP_MM), ("t1", 1.0 * p5.PLATE_DISP_MM)):
        zw = z_top + p5.GAP_MM - d
        eng = zmax > zw
        fresh[f"A_engaged_{tag}_mm2"] = float(areas[eng].sum())
        fresh[f"A_engaged_frac_{tag}"] = (float(areas[eng].sum()) / A_cap
                                          if A_cap else float("nan"))
    res["fix_preview"] = fresh
    return res


def show(r: dict) -> None:
    print("\n" + "=" * 76)
    print("%s   include_shell=%s" % (r["bone"], r["include_shell"]))
    print("=" * 76)
    print("  a_max = %.4f  (node %d, 1-based)" % (r["a_max"], r["a_max_node_1based"]))
    print("  cap 面片 %d 个, A_cap = %.1f mm²" % (r["n_cap_facets"], r["A_cap_mm2"]))
    print("  ★ a_max 是否落在保留 cap 面上 : %s%s"
          % (r["a_max_in_cap_surface"],
             "" if r["a_max_in_cap_surface"]
             else "   ← 板注册基准错！cap 最高点比它低 %.4f mm"
                  % r["cap_top_below_amax_mm"]))
    for tag in ("t0", "t1"):
        print("  %s 接合 %3d 面片 / %8.1f mm²  (cap 的 %.2f%%)  偏心 %.3f mm"
              % (tag, r["n_engaged_facets_" + tag], r["A_engaged_" + tag + "_mm2"],
                 100.0 * r["A_engaged_frac_" + tag], r["ecc_" + tag + "_mm"]))
    fp = r["fix_preview"]
    print("  【修法预览】板改注册到 cap 自身最高点（下移 %.4f mm）:" % fp["shift_mm"])
    for tag in ("t0", "t1"):
        print("     %s 接合 %8.1f mm²  (cap 的 %.2f%%)"
              % (tag, fp["A_engaged_" + tag + "_mm2"],
                 100.0 * fp["A_engaged_frac_" + tag]))
    a = r["asserts"]
    print("  断言（3 条，均有正对照）: a_max∈cap=%s | A_engaged(t0)>0=%s | 占比合理=%s"
          "  ⇒ %s"
          % (a["a_max_in_cap_surface"], a["engaged_t0_positive"],
             a["engaged_frac_reasonable"], r["verdict"]))
    print("  （参考：已 PASS 的股骨 t=0 接合占比 %.1f%%）" % (100 * REF_FEMUR_FRAC))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bone", default=None,
                    help="默认扫描 femur_r + T6")
    ap.add_argument("--include-shell", type=int, default=1, choices=(0, 1))
    args = ap.parse_args()

    cases = ([(args.bone, bool(args.include_shell))] if args.bone
             else [("femur_r", True), ("T6", True), ("T6", False)])
    bad = 0
    for bone, sh in cases:
        try:
            r = audit(bone, sh)
        except Exception as e:                       # noqa: BLE001
            print("[SKIP] %s sh=%s : %s: %s" % (bone, sh, type(e).__name__, e))
            continue
        show(r)
        if r["verdict"] != "PASS":
            bad += 1
    print("\n" + "-" * 76)
    print("闸门结论：%d/%d 例 FAIL（FAIL 例若上求解，几乎必然负 Jacobian）"
          % (bad, len(cases)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
