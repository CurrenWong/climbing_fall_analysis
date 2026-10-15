# -*- coding: utf-8 -*-
"""V1 跟骨 · 跖面接触 BC 修复探针（手册驱动）。

**动机（FEBio 手册 §8.5.3 原文）**

> FEBio often struggles with problems that have an **initial separation**. This is especially
> so when the problem is **force-driven** … In general, **when the contact interface is
> necessary to define a well-constrained problem it is best to avoid initial separation** if
> possible. **Another way around this issue is to first use a displacement to bring the
> surfaces in contact and then replace the displacement with a force.**

**现状**：``plantar_bc="contact"`` 的跖面接触是**单边约束**；``DEFAULT_CONTACT_GAP_MM = 0.5``
（**初始间隙**），且 ``s2_thums_sweep`` 只传 ``plantar_bc`` ⇒ ``contact_hold_spring_k=None``
（``plantar_bc.py`` 自己注释里的「弱 hold 弹簧」**默认关闭**）⇒ gap 闭合前骨无切向支撑
⇒ 切向刚度奇异 ⇒ Newton 发散 ⇒ 6 个皮质 hex8 负 Jacobian。

**本探针**按手册给的**三条**路径各试一遍（**不改默认**，全部 opt-in）：

| 路径 | 手册依据 | 对应变体 |
|---|---|---|
| ① **缩短初始间隙** | §8.5.3「best to avoid initial separation」 | ``gap001`` |
| ② **弱 hold 弹簧**（消刚体模态） | §8.1.3 underconstrained / §8.3.1.5 刚体模态 | ``hold*`` |
| ③ **先位移压到接触、再换力** | §8.5.3「first use a displacement … then replace」 | ``pretouch``（两步） |

外加：``swap``（primary/secondary 互换，§8.5 配方）、``ts*``（多时间步，让 ramp 更细）。

用法::

    .venv\\Scripts\\python.exe scripts\\opensim_fe\\plantar_contact_fix_probe.py --build
    .venv\\Scripts\\python.exe scripts\\opensim_fe\\plantar_contact_fix_probe.py --solve
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
sys.path.insert(0, str(_HERE))          # thums_feb / s2_thums_sweep 同目录
sys.path.insert(0, str(ROOT / "src"))

import thums_feb as T                                  # noqa: E402
import s2_thums_sweep as S2                            # noqa: E402
from climbing.coupling.load_transfer import transfer    # noqa: E402

WORKDIR_BASE = ROOT / "temp" / "opensim_fe" / "plantar_fix"
WORKDIR = WORKDIR_BASE          # main() 里按 h 重设
H_FIX = 5.0
ACT_FIX = 0.0

#: 变体表。(id, kwargs, 说明)
VARIANTS_BY_H: dict[float, list[tuple[str, dict, str]]] = {}


def _variants(kind: str = "hold_sweep") -> list[tuple[str, dict, str]]:
    """按需生成变体表（hold_k 跨 4 个数量级 + gap + swap + 参照）。"""
    if kind == "hold_sweep":
        base = [("g0_baseline", dict(),
                 "复现：gap=0.5、无 hold（＝s2_thums_sweep 的现状口径）")]
        for k in (10.0, 100.0, 1000.0, 10000.0, 100000.0):
            base.append((f"hold{_fmt(k)}", dict(contact_hold_spring_k=k),
                         f"② 弱 hold 弹簧 k={k:g} N/mm"))
        base += [
            ("gap001_hold10", dict(contact_gap_mm=0.01, contact_hold_spring_k=10.0),
             "①+② 间隙 0.01 mm ＋ hold 10"),
            ("swap_hold10", dict(contact_swap_pair=True, contact_hold_spring_k=10.0),
             "swap primary/secondary ＋ hold 10"),
            ("ref_fixed", dict(_kind="fixed"),
             "**参照**：plantar_bc=fixed（现状口径）"),
            ("ref_spring100", dict(_kind="spring", spring_k=100.0),
             "**参照**：plantar_bc=spring k=100 N/mm"),
        ]
        return base
    raise ValueError(kind)


def _fmt(k: float) -> str:
    return ("%g" % k).replace(".", "p")


def _build_one(mesh, vid: str, kwargs: dict, *, time_steps: int, spec) -> dict:
    """建 .feb + 闸门 3 自证（直接读 deck）。返回 meta。

    ``spec`` 由调用方**算一次**传入（``transfer()`` 内部会跑 OpenSim，约 25 s）。
    ``kwargs`` 里的 ``_`` 前缀键是探针私有键（``_kind`` / ``_time_steps``），
    不进 builder。
    """
    kw = {k: v for k, v in kwargs.items() if not k.startswith("_")}
    kind = kwargs.get("_kind", "contact")
    ts = int(kwargs.get("_time_steps", time_steps))

    d = (WORKDIR / vid).resolve()
    d.mkdir(parents=True, exist_ok=True)
    feb = d / "m.feb"
    for stale in (feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                  feb.with_suffix(".log")):
        if stale.exists():
            stale.unlink()

    T.build_thums_feb(
        mesh, feb,
        load_n=spec.subtalar_n, load_dir=spec.subtalar_dir,
        achilles_n=spec.achilles_n, achilles_dir=spec.achilles_dir,
        use_rigid=False, time_steps=ts,
        plantar_bc=kind, **kw,
    )

    # ---- 闸门 3：受控量必须打印（直接读 deck，不信内存） ----
    # ⚠️ 大小写：pyfebio 序列化为 <Contact>/<SurfacePair>/<NodeSet>/<DiscreteSet>
    # （首字母大写），早先的小写正则让本检查器给出**假阴性** —— 故此处加
    # ``re.I`` 并配**已知正对照**（hold 变体必须 >0，非 hold 变体必须 ==0）。
    deck = feb.read_bytes().decode("iso-8859-1", errors="replace")
    n_contact = len(re.findall(r"<Contact\b", deck, re.I))
    n_pair = len(re.findall(r"<SurfacePair\b", deck, re.I))
    n_hold = len(re.findall(r"<DiscreteSet name=\"plantar_contact_hold\"", deck,
                            re.I))
    n_hold_bc = len(re.findall(r"node_set=\"plantar_contact_hold_anchor\"", deck,
                               re.I))
    steps = len(re.findall(r"<Step\b", deck, re.I))
    meta = {
        "id": vid, "kind": kind,
        "kwargs": {k: kwargs.get(k) for k in
                   ("contact_gap_mm", "contact_hold_spring_k",
                    "contact_swap_pair", "contact_node_reloc")},
        "time_steps": ts,
        "deck_contacts": n_contact, "deck_surface_pairs": n_pair,
        "deck_hold_discrete_sets": n_hold, "deck_hold_anchor_bc": n_hold_bc,
        "deck_steps": steps,
        "f_sub_N": float(spec.subtalar_n), "f_ach_N": float(spec.achilles_n),
        "feb": str(feb),
    }
    # 闸门 3 自证（按 kind 分支，正反对照）
    if kind == "contact":
        assert n_contact >= 1, f"{vid}: deck 无 <Contact> —— 通道未生效"
        want_hold = bool((kwargs.get("contact_hold_spring_k") or 0) > 0)
        got_hold = n_hold > 0 and n_hold_bc > 0
        assert got_hold == want_hold, (
            f"{vid}: hold 弹簧请求={want_hold} 但 deck 里 {'有' if got_hold else '无'}"
            " —— 通道未生效（闸门 3）"
        )
        assert n_pair >= 1, f"{vid}: deck 无 <SurfacePair> —— 接触对未注册（闸门 1）"
    else:
        assert n_contact == 0, (
            f"{vid}: plantar_bc={kind} 却出现 <Contact> —— 通道串了（闸门 3）")
    return meta


def _neg_jac_info(log: Path) -> dict:
    """从 .log 里抓负 Jacobian 数与涉及单元。"""
    if not log.is_file():
        return {"present": False}
    txt = log.read_text(encoding="utf-8", errors="replace")
    ids: list[int] = []
    for m in re.finditer(r"[Nn]egative [Jj]acobian.*?element[s]?\s*[:=]?\s*([0-9,\s]+)",
                         txt):
        ids += [int(x) for x in re.findall(r"\d+", m.group(1))]
    n_declared = re.findall(r"(\d+)\s*negative\s+jacobians?", txt, re.I)
    return {
        "present": True,
        "n_declared": int(n_declared[-1]) if n_declared else None,
        "elements": sorted(set(ids)),
        "tail": txt[-1500:],
    }


def _solve_one(meta: dict, *, debug: str = "", post: dict | None = None) -> dict:
    """跑一个变体；返回 rc / 耗时 / 负 Jacobian / 位移量级 /（可选）gauge_max。

    ``debug``: ``""`` / ``-g`` / ``-g2``（手册 §8.2）。
    ``post``: 传给后处理的 ``{"cort_centroids":…, "cort_vol":…, "gauge_mm":…}``。
    """
    from climbing.coupling.febio_run import find_febio

    feb = Path(meta["feb"])
    exe = find_febio()
    cmd = [str(exe), "-i", str(feb)]
    if debug:
        cmd.append(debug)
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              cwd=str(feb.parent), timeout=1800)
        rc, out, err = proc.returncode, proc.stdout or "", proc.stderr or ""
    except subprocess.TimeoutExpired:
        rc, out, err = -999, "", "TIMEOUT"
    el = time.time() - t0

    rec = dict(meta)
    rec.update({"rc": int(rc), "elapsed_s": round(el, 1), "debug": debug,
                "negJ": _neg_jac_info(feb.with_suffix(".log")),
                "fatal": [l.strip() for l in (out + err).split("\n")
                          if "FATAL" in l.upper()][:4]})
    xplt = feb.with_suffix(".xplt")
    rec["xplt_exists"] = xplt.is_file()
    rec["xplt_bytes"] = xplt.stat().st_size if xplt.is_file() else 0

    # ---- 后处理：皮质 gauge_max（与 s2_thums_sweep 同口径） ----
    if rc == 0 and post and xplt.is_file():
        try:
            from climbing.coupling import fe_post
            from climbing.coupling.xplt_shim import patch_unnamed_domains
            from pyfebio import xplt as xplt_mod
            # 接触模型的 .xplt 里 FEBio 会多写一个**无名域** ⇒ pyfebio 解析器抛
            # KeyError('name')（其 nodeset 分支有空名兜底、domain 分支没有）。
            patch_unnamed_domains()
            hdf5 = feb.with_suffix(".hdf5")
            if not hdf5.is_file() or hdf5.stat().st_size == 0:
                xplt_mod.to_hdf5(str(xplt), str(hdf5))
            vm = S2._read_domain_vm(xplt, "calcaneus", hdf5)
            g = T.gauge_von_mises_mixed(vm, post["cort_centroids"],
                                        post["cort_vol"], post["gauge_mm"])
            rec["gauge_max"] = float(g.max())
            rec["vm_max"] = float(vm.max())
            rec["vm_mean"] = float(vm.mean())
            # 热点单元（gauge 峰值所在 hex8 的序号 + 质心）
            ci = int(np.argmax(g))
            rec["gauge_argmax_elem"] = ci
            rec["gauge_argmax_centroid_mm"] = [
                round(float(x), 3) for x in post["cort_centroids"][ci]]
        except Exception as exc:                      # noqa: BLE001
            rec["post_error"] = f"{type(exc).__name__}: {exc}"
    return rec


def main(argv: list[str] | None = None) -> int:
    global WORKDIR, H_FIX
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", action="store_true")
    ap.add_argument("--solve", action="store_true")
    ap.add_argument("--only", type=str, nargs="+", default=None)
    ap.add_argument("--h", type=float, default=H_FIX, help="坠落高度 (m)")
    ap.add_argument("--variants", type=str, default="hold_sweep")
    ap.add_argument("--time-steps", type=int, default=1)
    ap.add_argument("--debug", type=str, default="")
    ap.add_argument("--keep", action="store_true", help="不清理既有目录")
    args = ap.parse_args(argv)

    if not (args.build or args.solve):
        args.build = True

    H_FIX = float(args.h)
    WORKDIR = WORKDIR_BASE / ("h%g" % H_FIX)
    if WORKDIR.exists() and not args.keep:
        shutil.rmtree(WORKDIR, ignore_errors=True)
    WORKDIR.mkdir(parents=True, exist_ok=True)

    print("[load] THUMS mesh …", flush=True)
    mesh = T.load_thums_mesh()
    print(f"[load] nodes={len(mesh['nodes'])} cort={len(mesh['elements_cort'])}",
          flush=True)
    f_sub, t1 = S2._load_at(H_FIX)
    print(f"[load] load_transfer (OpenSim) …", flush=True)
    spec = transfer(f_sub, achilles_activation=ACT_FIX)
    print(f"[load] h={H_FIX} act={ACT_FIX}  subtalar={spec.subtalar_n:.1f} N "
          f"achilles={spec.achilles_n:.1f} N", flush=True)

    cort_centroids, cort_vol = T.hex_centroids_volumes(mesh["nodes"],
                                                       mesh["elements_cort"])
    gauge_mm = float(getattr(S2, "REG_LEN_MM", 4.0))
    print(f"[load] gauge_mm (ℓ) = {gauge_mm}  cort hex8 = {len(cort_centroids)}",
          flush=True)
    post = {"cort_centroids": cort_centroids, "cort_vol": cort_vol,
            "gauge_mm": gauge_mm}

    metas, rows = [], []
    for vid, kw, desc in _variants(args.variants):
        if args.only and vid not in args.only:
            continue
        print(f"\n=== {vid}  {desc}", flush=True)
        m = _build_one(mesh, vid, kw, time_steps=args.time_steps, spec=spec)
        m["desc"] = desc
        print(f"    [闸门3] <Contact>×{m['deck_contacts']} "
              f"<SurfacePair>×{m['deck_surface_pairs']} "
              f"hold_disc×{m['deck_hold_discrete_sets']} "
              f"hold_bc×{m['deck_hold_anchor_bc']} steps={m['deck_steps']} "
              f"gap={m['kwargs']['contact_gap_mm']} "
              f"hold_k={m['kwargs']['contact_hold_spring_k']} "
              f"swap={m['kwargs']['contact_swap_pair']}", flush=True)
        metas.append(m)
        if args.solve:
            r = _solve_one(m, debug=args.debug, post=post)
            nj = r["negJ"]
            gtxt = ("%.4f MPa" % r["gauge_max"]) if "gauge_max" in r else "—"
            print(f"    [solve] rc={r['rc']} {r['elapsed_s']}s  "
                  f"negJ={nj.get('n_declared')} n_elem={len(nj.get('elements') or [])} "
                  f"gauge_max={gtxt} xplt={r['xplt_bytes']}B", flush=True)
            rows.append(r)

    out = WORKDIR / "probe.json"
    out.write_text(json.dumps({"metas": metas, "rows": rows,
                               "f_sub_N": float(spec.subtalar_n),
                               "gauge_mm": gauge_mm,
                               "threshold_mpa": float(S2.THRESHOLD_MPA)},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n[out] {out}")
    if rows:
        print("\n| variant | kind | gap | hold_k | swap | ts | rc | 秒 | negJ | "
              "gauge_max | σ/σ_fixed | 热点质心 |")
        print("|---|---|---|---|---|---|---|---|---|---|---|---|")
        ref = next((r for r in rows if r["id"] == "ref_fixed"), None)
        ref_g = ref.get("gauge_max") if ref else None
        for r in rows:
            k, g = r["kwargs"], r.get("gauge_max")
            ratio = ("%.3f" % (g / ref_g)) if (g and ref_g) else "—"
            cen = r.get("gauge_argmax_centroid_mm")
            print("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                r["id"], r.get("kind"), k["contact_gap_mm"],
                k["contact_hold_spring_k"], k["contact_swap_pair"],
                r["time_steps"], r["rc"], r["elapsed_s"],
                r["negJ"].get("n_declared"),
                ("%.4f" % g) if g else "—", ratio,
                ("(%s)" % ",".join("%.1f" % x for x in cen)) if cen else "—"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
