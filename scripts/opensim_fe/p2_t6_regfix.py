# -*- coding: utf-8 -*-
"""P2 · V1：**修法 1 单变量**重跑 T6 能力门（压板注册基准修正）

单变量声明
----------
**只改一处**：压板的 z 注册基准 —— 由「顶带节点集的 `a_max`」改为
「**真正参与接触的 cap 面片自身的最高点**」。其余一切（GAP_MM、PLATE_DISP_MM、
TOP_FRAC、CAP_AXDIV_MIN、penalty 公式、fric、lateral_support、3 节点锚点、
time_steps、载荷档位）**逐字不动**。

背景（`results/opensim_fe/P2_T6_ROOTCAUSE.md`）
------------------------------------------------
`p5_head_contact.plate_geom_from_mesh` 用 `z_wall = geom["a_max"] + GAP_MM`；
而 `a_max` 来自**顶带节点集**，接触面片是**另一个集合**（还被 `axdiv > 0.5` 过滤）。
T6 的 `a_max` 节点（1-based #282）**不在**任何保留 cap 面片上 ⇒ 压板悬空 1.2342 mm
⇒ t=0/t=1 接合 0/36 面片 ⇒ **全程零接触**；叠加骨在 z 向无约束 ⇒ 纯刚体模态
⇒ FEBio 手册 §8.3.1.5 的负 Jacobian。

本脚本把注册基准改为 cap 面自身最高点，其余不动，直接对比实录配置。

⚠️ 关于「修法 2（恢复底端 z 约束）」**已作废**：底端一固定，底端反力即分掉载荷
   ⇒ `C = F − R ≠ F` ⇒ `bear_frac` 判据当场失效。P5 设计**必须**让骨轴向自由、
   由接触供 z 锚位 —— 故**注册修正是唯一正确修法**（详见根因件 §3 更正）。

判据（**跑前写死**）
--------------------
  V1a  接触接合面积 > 0（修法生效的几何必要条件；零 FE 已预验 32.0 mm²）
  V1b  求解收敛（无负 Jacobian）
  V1c  接触承载 `bear_frac ∈ [0.80, 1.20]`（承预注册 §5.1）
  V1d  线性度 r² ≥ R2_MIN（承预注册）
  —— 四条全过 = **V1 PASS**；任一不过 = 如实 FAIL，不调参。

用法
----
    .venv\\Scripts\\python.exe scripts\\opensim_fe\\p2_t6_regfix.py
    # --loads 可换成别的档位；默认 [2000.0]（⇒ solve load 1220.8 N，与实录同一档）
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import g4_fe_plate as plate             # noqa: E402
import p2_head_contact as p2            # noqa: E402
import p4_pb_probe as p4                # noqa: E402
import p5_head_contact as p5            # noqa: E402

# --------------------------------------------------------------------------
# 单变量修正：压板 z 注册基准 = 保留 cap 面片自身的最高点
# --------------------------------------------------------------------------
_ORIG_PLATE_GEOM = p5.plate_geom_from_mesh
SHIFT_LOG: dict = {}


def plate_geom_regfix(nodes_rot: np.ndarray, geom: dict) -> dict:
    """`p5.plate_geom_from_mesh` 的注册修正版（只动 z 基准，其余逐字照抄）。"""
    pg = _ORIG_PLATE_GEOM(nodes_rot, geom)          # 先拿原值（含原 z_wall）
    cap = geom["cap"]
    z_cap_top = float(max(nodes_rot[np.asarray(f["poly"], int) - 1, 2].max()
                          for f in cap))
    z_wall = z_cap_top + p5.GAP_MM                  # ← 唯一改动
    shift = float(pg["z_wall"] - z_wall)
    pg["z_wall"] = z_wall
    pg["z_bot"] = z_wall - plate.PLATE_THICK_MM
    pg["z_mid"] = z_wall - 0.5 * plate.PLATE_THICK_MM
    pg["regfix"] = {"z_cap_top": z_cap_top,
                    "a_max": float(geom["a_max"]),
                    "shift_mm": shift}
    SHIFT_LOG.setdefault("calls", []).append(pg["regfix"])
    return pg


# 判据：**逐字引用官方**（`p2_head_contact` §5.1 跑前冻结），不自带副本
BEAR_LO, BEAR_HI = p2.BEAR_LO, p2.BEAR_HI
R2_MIN = p2.R2_MIN
# ⚠️ 自错记录：本脚本首版把 R2_MIN 写死成 0.90（臆造）。官方为 0.20，
#    且闸门是 **OR**：bend_ok = (r2 >= R2_MIN) or (szz_max > 0)。
#    现已改为直接引用官方常量 + 调用官方 `judge_case`，避免二次发明判据。


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--loads", nargs="*", type=float, default=None,
                    help="诊断模式：关节载荷档（×share 后为 solve load），"
                         "如 2000 ⇒ 1220.8 N。省略 = 跑**官方能力门**（4 高度）。")
    ap.add_argument("--bone", default="T6")
    ap.add_argument("--workdir", default="p2_head_contact_regfix")
    args = ap.parse_args()

    p5.plate_geom_from_mesh = plate_geom_regfix    # ← 单变量注入
    p2.BASE = ROOT / "temp" / "opensim_fe" / args.workdir   # 新目录，避免复用旧 h5

    print("=" * 92)
    print("P2 · V1  修法 1 单变量：压板 z 注册基准 = 保留 cap 面片自身最高点")
    print("  只改这一处；GAP_MM/PLATE_DISP_MM/TOP_FRAC/penalty 公式/… 逐字不动")
    print("  工作目录（新）: %s" % p2.BASE.relative_to(ROOT))
    print("  档位: %s  (solve load = L × share)" % args.loads)
    print("  判据: bear∈[%.2f,%.2f] 且 r2≥%.2f 且收敛 且 接合面积>0"
          % (BEAR_LO, BEAR_HI, R2_MIN))
    print("=" * 92, flush=True)

    t0 = time.time()
    # 求解前：先零 FE 记录注册修正量（证明 patch 生效）
    ctx = p4.load_ctx(args.bone)
    part, mesh = ctx[0], ctx[1]
    import bone_feb as bf
    geom = p2.cap_geom(mesh, part.spec, top_frac=p2.TOP_FRAC,
                       include_shell=True)
    nodes_rot = bf._align_pca_rotate(mesh["nodes"])
    pg = plate_geom_regfix(nodes_rot, geom)
    print("[regfix] a_max=%.4f  z_cap_top=%.4f  **下移 %.4f mm**"
          % (pg["regfix"]["a_max"], pg["regfix"]["z_cap_top"],
             pg["regfix"]["shift_mm"]), flush=True)

    res = p2.run_bone(args.bone, p2.HEIGHTS, force_loads=args.loads,
                      with_shell=True)
    elapsed = round(time.time() - t0, 1)

    # ---- 判定：调**官方** judge_case ----
    print("\n" + "=" * 92)
    print("V1 判定（官方 judge_case；BEAR ∈ [%.2f,%.2f]，R2_MIN=%.2f，bend=OR 条件）"
          % (BEAR_LO, BEAR_HI, R2_MIN))
    print("=" * 92)
    verdicts = []
    for r in res["rows"]:
        conv = bool(r.get("converged"))
        rec = dict(r)
        rec["converged"] = conv
        j = p2.judge_case(rec)
        print("  F=%9.1fN pen=%6.2f  conv=%-5s bear=%s  r2=%s  szzmax=%s"
              % (r["load_n"], r["penalty"], conv,
                 ("%.6f" % r["bear_frac"]) if r.get("bear_frac") is not None else "None",
                 ("%.3f" % r["r2_linear"]) if r.get("r2_linear") is not None else "None",
                 r.get("szz_max_mpa")))
        print("     bear_ok=%s bend_ok=%s ⇒ %s   %s"
              % (j["bear_ok"], j["bend_ok"], "PASS" if j["pass"] else "FAIL",
                 j["reasons"] or ""))
        if not conv:
            print("     error: %s" % str(r.get("error"))[:200])
        verdicts.append({"load_n": r["load_n"], "penalty": r["penalty"],
                         "converged": conv,
                         "bear_frac": r.get("bear_frac"),
                         "r2_linear": r.get("r2_linear"),
                         "szz_max_mpa": r.get("szz_max_mpa"),
                         "median_mpa": r.get("median_mpa"),
                         "p95_mpa": r.get("p95_mpa"),
                         "max_mpa": r.get("max_mpa"),
                         "judge": j, "error": r.get("error")})
    overall = all(v["judge"]["pass"] for v in verdicts) if verdicts else False
    n_pass = sum(1 for v in verdicts if v["judge"]["pass"])

    out = {"probe": "P2 V1 修法1 单变量（压板注册基准修正）",
           "single_variable": "z_wall 基准: a_max(顶带节点集) -> cap 面片自身最高点",
           "frozen_else": {"GAP_MM": p5.GAP_MM, "PLATE_DISP_MM": p5.PLATE_DISP_MM,
                           "TOP_FRAC": p2.TOP_FRAC, "LOAD_FRAC": p2.LOAD_FRAC,
                           "CAP_AXDIV_MIN": p5.CAP_AXDIV_MIN,
                           "penalty_formula": "50 * F_ref / 57941.3",
                           "fric": 0.0, "lateral_support": p2.LATERAL,
                           "n_anchor": p2.N_ANCHOR},
           "regfix": pg["regfix"],
           "loads_joint_n": args.loads,
           "share": res["load_share"],
           "criteria": {"bear": [BEAR_LO, BEAR_HI], "r2_min": R2_MIN},
           "rows": verdicts, "verdict": "PASS" if overall else "FAIL",
           "elapsed_s": elapsed}
    dst = ROOT / "results" / "opensim_fe" / "P2_V1_REGFIX.json"
    dst.write_text(json.dumps(out, ensure_ascii=False, indent=1), encoding="utf-8")
    print("\n  ⇒ **V1 = %s**" % out["verdict"])
    print("  产物: %s" % dst.relative_to(ROOT))
    print("  耗时: %.1f s" % elapsed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
