"""G4 ② P1 试点 · **导向 BC 重跑**（A/B 对照 §9.9 的「143/150 发散」）

预注册（**先写后跑**）
----------------------
若 §9.9 的发散根因是 **BC 屈曲**（P_cr ≈ 252 N，而试点载荷 phys=1769–7905 N、
frozen≤25266 N 全部超阈），则把**载荷端**（V7 的 ``top_load`` 带）的 x/y 自由度锁住
（= 导向压头 / 试验机导向加载）后：

* **phys（load_fraction=0.07）**：应**大面积收敛**（§9.9 是大量发散）；
* **frozen（1.0）**：载荷最高 25266 N，可能仍发散（远超任何合理承载）——
  若也大面积收敛，说明此前"不可达"完全是屈曲假象。

判定口径：与 §9.9 **同一组格**（同一高度表、同一载荷链、同一 THUMS 参数、
同一 ``steps=1``），**唯一差别 = 多一个 x/y 零位移 BC 块** → 严格 A/B。

只 ADD；复用 ``g4_fe_pilot``（载荷链/网格循环口径）/ ``g4_fe_p1b``（V7 建模）/
``g4_fe_areascan``（BC 注入），不改任何既有模块。

复现::

    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g4_fe_pilot_guided.py --hs 5,10,20,35,50
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[2]
for _p in (ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import g4_fe_areascan as asc  # noqa: E402
import g4_fe_p1b as p1b  # noqa: E402
import g4_fe_pilot as pil  # noqa: E402

OUT_RESULTS = ROOT / "results" / "opensim_fe"

#: ⚠️ **必须 align_pca=True**：pilot 原 V7 用 ``align_pca=False`` + ``load_axis=-pca``，
#: 骨长轴不在 z 上（实测 ``pca·z ≈ 0.734``）⇒ 世界 x/y 零位移里含 **~68% 的骨轴向分量**，
#: 等于把骨**部分轴向夹死**（2026-10-12 实测：R/load = 0.734 = 恰好是载荷的 z 分量占比）。
#: 改用对齐帧（骨轴 → z）后，x/y 约束才是**纯横向** ⇒ 才能叫"导向"。
V7A = dict(kind="nodal", top_frac=0.50, bot_frac=0.50, align_pca=True)
AXIS_Z = (0.0, 0.0, -1.0)          # 对齐帧下的轴向压缩


def parse_guide(s):
    """导向范围：``full`` = 整个载荷带(top 50%，= 之前的做法)；``none``；``lo:hi`` = z_frac 子段。"""
    s = str(s).strip().lower()
    if s in ("full", ""):
        return "full"
    if s in ("none", "-", "0"):
        return None
    lo, hi = s.split(":")
    return (float(lo), float(hi))


def guide_label(g) -> str:
    if g is None:
        return "none"
    if g == "full":
        return "full"
    return "%.2f-%.2f" % g


def guide_ids_for(nodes: np.ndarray, band) -> list[int]:
    """按 z_frac 取导向节点（对齐帧；z_frac=0 远端, 1 apex 端）。"""
    z = nodes[:, 2]
    zf = (z - float(z.min())) / (float(z.max()) - float(z.min()))
    lo, hi = float(band[0]), float(band[1])
    return sorted(int(i) + 1 for i in np.flatnonzero((zf >= lo) & (zf <= hi)))


def build_v7_deck(setup, *, path: Path, load_n: float, Dmax, steps: int = 1,
                  guide="full", nodes: np.ndarray | None = None) -> Path:
    """V7A deck（对齐帧）；``guide`` 指定锁 x/y 的范围（= 试验机导向压头的理想化）。"""
    p1b.build_variant(setup, path=path, load_n=float(load_n), steps=int(steps),
                      damage_Dmax=Dmax, load_axis=AXIS_Z, **V7A)
    if guide is None:
        return path
    text = path.read_text(encoding="ISO-8859-1")
    if guide == "full":
        if 'node_set="top_load"' not in text:
            raise RuntimeError("deck 里找不到 top_load（V7 载荷带）")
        text = asc.inject_boundary(text, asc.lateral_bc_xml("top_load"))
    else:
        if nodes is None:
            raise RuntimeError("guide 子段需要 nodes=")
        ids = guide_ids_for(nodes, guide)
        if not ids:
            raise RuntimeError(f"guide={guide} 选到 0 个节点")
        text = asc.inject_fragment(text, asc.nodeset_xml("guide_band", ids))
        text = asc.inject_boundary(text, asc.lateral_bc_xml("guide_band"))
    path.write_text(text, encoding="ISO-8859-1")
    return path


def solve_v7(setup, *, load_n: float, Dmax, tag: str, timeout: float,
             guide="full", nodes: np.ndarray | None = None, steps: int = 1) -> dict:
    gl = guide_label(guide)
    feb = pil.OUT_TEMP / f"gd{gl.replace('.', '')}_{tag}.feb"
    build_v7_deck(setup, path=feb, load_n=load_n, Dmax=Dmax, steps=steps,
                  guide=guide, nodes=nodes)
    rec = p1b.run_and_post(setup, feb, timeout=timeout,
                           read_damage=(Dmax is not None))
    rec.update({"tag": tag, "load_n": float(load_n), "Dmax": Dmax,
                "guide": gl})
    return rec


def axial_stats(hdf5: Path, dom: str = "fib_cort") -> dict:
    """从 h5 取**轴向应力** σ_zz 与 σ_vm —— 导向 BC 下两者不相等。

    ⚠️ 2026-10-12 实测：底部全固定 + 顶部 50% 锁 x/y ⇒ 骨在**近单轴应变**状态下变形
    （横向膨胀被约束）⇒ σ_vm 被压低：σ_vm_max/|σ_zz|_load ≈ 12.36/19.75 = **0.63**
    （ν=0.3 单轴应变理论值 0.57）⇒ **σ_vm 低估轴向应力约 1.6×**。
    骨折判据（σ_c=160 MPa 是**压缩强度**）应看 σ_zz / 主应力，故两者都记。
    """
    with h5py.File(hdf5, "r") as h5:
        sid = sorted(h5["states"].keys(), key=int)[-1]
        s = np.asarray(h5[f"states/{sid}/element_data/stress/{dom}"][:],
                       dtype=np.float64)
    vm = bf._von_mises_from_voigt(s)
    zz = s[:, 2]
    return {"vm_max": float(vm.max()), "szz_max": float(zz.max()),
            "szz_min": float(zz.min()), "szz_absmax": float(np.abs(zz).max()),
            "vm_over_abszz": (float(vm.max() / np.abs(zz).max())
                              if np.abs(zz).max() > 0 else None)}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--timeout", type=float, default=240.0)
    ap.add_argument("--hs", default="", help="逗号分隔的高度（空 = 全 1..50）")
    ap.add_argument("--fracs", default="phys", help="逗号分隔: phys,frozen")
    ap.add_argument("--dm", nargs="*", default=["0.9"],
                    help="Dmax 列表（默认 0.9；§9.9 用过 0.9/1.0；none = 纯弹性对照）")
    ap.add_argument("--tag", default="guided")
    ap.add_argument("--guide", default="full",
                    help="逗号分隔：full（整载荷带）/ none / lo:hi（z_frac 子段）")
    args = ap.parse_args()

    pil.OUT_TEMP.mkdir(parents=True, exist_ok=True)
    setup = p1b.Setup()
    nodes = setup.nodes_for(True)
    guides = [parse_guide(x) for x in args.guide.split(",") if x.strip()]
    lc = pil.load_chain()
    heights_all = lc["heights"]
    ankle_all = lc["ankle_vert"]

    # --- 正对照：z_frac 子段 [0.5,1.0] 的节点集合必须 == deck 的 top_load 集合 ---
    try:
        _tmp = pil.OUT_TEMP / "_guide_ctrl.feb"
        build_v7_deck(setup, path=_tmp, load_n=1.0, Dmax=None, guide="full")
        _t = _tmp.read_text(encoding="ISO-8859-1")
        _m = re.search(r'<NodeSet name="top_load">([^<]*)</NodeSet>', _t)
        _deck = {int(x) for x in _m.group(1).split(",") if x}
        _mine = set(guide_ids_for(nodes, (0.50, 1.00)))
        print(f"[ctrl] guide_ids(0.50,1.00)={len(_mine)} vs deck top_load={len(_deck)}"
              f" 相等={_mine == _deck}"
              f"  交集={len(_mine & _deck)}  仅deck={len(_deck - _mine)}", flush=True)
        if _mine != _deck:
            print("       ⚠️ 不等 —— 子段导向的节点集口径需复核（不阻断，仅告警）")
    except Exception as exc:  # noqa: BLE001
        print(f"[ctrl] 跳过（{type(exc).__name__}: {exc}）", flush=True)

    hs = ([float(x) for x in args.hs.split(",") if x.strip()]
          if args.hs.strip() else [float(h) for h in heights_all])

    def _frac(s: str) -> tuple[str, float]:
        s = s.strip()
        named = {"phys": pil.PHYS_FRACTION, "frozen": pil.FROZEN_FRACTION}
        if s in named:
            return s, named[s]
        return f"f{float(s):g}", float(s)

    fracs = [_frac(f) for f in args.fracs.split(",") if f.strip()]
    dms = [None if str(x).strip().lower() in ("none", "-", "elastic")
           else float(x) for x in args.dm]

    print(f"[pilot-guided] align_pca=True、load_axis=-z、top/bot=0.5、steps=1；"
          f"{len(hs)} 高度 × {len(fracs)} 组 × {len(dms)} Dmax × "
          f"{len(guides)} 导向范围；P_cr(未导向)≈252 N", flush=True)
    rows: list[dict] = []
    for guide in guides:
        gl = guide_label(guide)
        for fn, frac in fracs:
            for h in hs:
                ankle = float(ankle_all[list(heights_all).index(h)])
                load = ankle * frac
                for dm in dms:
                    tag = f"{fn}_{'el' if dm is None else 'D%g' % dm}_h{int(h):02d}"
                    t1 = time.time()
                    rec = solve_v7(setup, load_n=load, Dmax=dm, tag=tag,
                                   timeout=args.timeout, guide=guide, nodes=nodes)
                    row = {"guide": gl, "frac_name": fn,
                           "load_fraction": frac, "Dmax": dm, "h": float(h),
                           "load_n": load,
                           "ok": bool(rec.get("ok")), "rc": rec.get("rc"),
                           "sigma_max": rec.get("peak_vm_mpa"),
                           "sigma_p95": rec.get("p95_vm_mpa"),
                           "R_net_n": rec.get("R_net_n"),
                           "R_ratio": ((rec.get("R_net_n") or 0.0) / load
                                       if load else None),
                           "wall_s": rec.get("wall_clock_s"),
                           "error_reason": (rec.get("error_reason") or "")[:120]}
                    if row["ok"]:
                        try:
                            row.update(axial_stats(
                                pil.OUT_TEMP / f"gd{gl.replace('.', '')}_{tag}.hdf5"))
                        except Exception as exc:  # noqa: BLE001
                            row["axial_error"] = f"{type(exc).__name__}: {exc}"
                    rows.append(row)
                    print(f"  [guide={gl:9s} {fn:6s} "
                          f"{'el' if dm is None else 'D%g' % dm}"
                          f" h={int(h):02d}] load={load:9.1f} N  ok={row['ok']}  "
                          f"vm={_f(row['sigma_max'])}  "
                          f"|szz|max={_f(row.get('szz_absmax'))}  "
                          f"vm/|szz|={_f(row.get('vm_over_abszz'), 2)}  "
                          f"R/load={_f(row['R_ratio'], 3)}  ({time.time() - t1:.0f}s)"
                          f"  {(row['error_reason'] or '')[:30]}", flush=True)

    def rate(g) -> str:
        sub = [r for r in rows if r["guide"] == guide_label(g)]
        return f"{sum(1 for r in sub if r['ok'])}/{len(sub)}"
    print("\n[sum] 收敛率 " + "  ".join(f"{guide_label(g)}={rate(g)}"
                                      for g in guides)
          + "   （§9.9 未导向：143/150 发散）", flush=True)
    out = {"align_pca": True, "load_axis": list(AXIS_Z), "v7": V7A,
           "guide_bc": "z_frac 子段 / top_load 的 x/y zero displacement",
           "p_cr_unguided_n": 252.0, "guides": [guide_label(g) for g in guides],
           "ok_rate": {guide_label(g): rate(g) for g in guides}, "rows": rows}
    p = OUT_RESULTS / f"G4_FE_PILOT_GUIDED_{args.tag}.json"
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=str),
                 encoding="utf-8")
    print(f"[out] {p}", flush=True)
    return 0


def _f(x, nd: int = 1) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return "—"


if __name__ == "__main__":
    raise SystemExit(main())
