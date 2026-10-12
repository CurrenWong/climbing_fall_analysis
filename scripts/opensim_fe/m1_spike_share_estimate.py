r"""m1_spike_share_estimate.py -- M1-spike：不做两体 FE，先用闭式并行构件 + 实测 EA
判定「耦合后腓骨的自洽分担比」，以决定是否值得投全量 M1（胫骨+腓骨+IOM）。

**为什么不用 FE 先做**：FEBio 里两体点对点耦合只能走 `Tie`/surface，而本项目已踩过
surface 被静默丢弃（`docs/notes/CASEB_CONTACT_ROOTCAUSE.md`）⇒ 新 BC 必须先有正对照；
而**分担当量**可由力学闭式给出，且不依赖任何新 BC。

**力学模型（串联-并联，两骨在重叠段 [z_prox, z_dist] 上并联）**：
    近端胫腓关节 + 远端下胫腓联合把腓骨的两端"栓"在胫骨上 ⇒ 重叠段上两骨**轴向并联**
    ⇒ 载荷按轴向刚度分配：**share = E_f·A_f / (E_f·A_f + E_t·A_t)**
    ⇒ **与几何（重叠长度）无关**（经典并联弹簧结论）⇒ 这也说明「share 沿全长均匀」
    只是单骨模型的人为假设，耦合后 share 是**由刚度比一次性决定的**。
    若再加骨间膜(IOM, Unterberg 三区弹簧)，share 只会**更高** ⇒ 本估计是**下界**。

**三种失效机制（决定该用哪个判据）**：
    ① 屈曲（孤立骨，FE 实测 P_cr=252 N，被实验带 77–370 N 佐证）
    ② 材料压溃：σ_c · A_f
    ③ 实测 p(h)：share × ankle_r(h)

**自检（判据）**：share → 1（A_f→∞）、share → 0（A_f→0）、EA 比与文献「最高 30%」对照。

复现（仓库根, PYTHONPATH=src）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\m1_spike_share_estimate.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import bone_feb as bf  # noqa: E402
import fracture_matrix as fm  # noqa: E402
from scripts.opensim_fe import fibula_fe_audit as audit  # noqa: E402

RES = ROOT / "results" / "opensim_fe"
P_CR_ISOLATED = 252.0        # N，本项目孤立腓骨 FE 屈曲阈（G4 计划 §9.13）
F_FAIL_EXP = (77.0, 370.0)   # N，实测孤立腓骨压缩失效带（Cureus 454613）
SIGMA_C_THUMS = 70.0         # MPa，G9 冻结 σ_c


def cortical_geometry(spec, mesh: dict) -> dict:
    """皮质域：体积、轴向长度、E（取材料卡）、A_cort = ΣV/L。

    连接表在 `mesh["by_solid"][mesh_name]` = {hexes, tets, pids, mat_name}；
    ⚠️ hex8 体积用仓库 `_hex8_signed_volumes`（**1-based 局部索引**），tet4 用教科书公式。
    ⚠️ `spec.solid_domains[i].topology` 只是字符串 'auto'，**不能**从它取连接表。
    """
    tot_v, zmin, zmax, mats = 0.0, np.inf, -np.inf, []
    per_dom = {}
    for sd in spec.solid_domains:
        nm = sd.mesh_name
        if "cort" not in nm.lower():
            continue
        e = mesh["by_solid"][nm]
        hexes = np.asarray(e["hexes"], int)
        tets = np.asarray(e["tets"], int)
        vols, zs = [], []
        if len(hexes):
            vols.append(np.abs(bf._hex8_signed_volumes(mesh["nodes"], hexes)))
            zs.append(mesh["nodes"][hexes - 1][:, :, 2])
        if len(tets):
            P = mesh["nodes"][tets - 1]
            vols.append(np.abs(np.einsum("ij,ij->i", np.cross(P[:, 1] - P[:, 0],
                                                          P[:, 2] - P[:, 0]),
                                         P[:, 3] - P[:, 0])) / 6.0)
            zs.append(P[:, :, 2])
        v = np.concatenate(vols) if vols else np.zeros(0)
        zz = np.concatenate([z.ravel() for z in zs]) if zs else np.zeros(0)
        mat = sd.material
        per_dom[nm] = {"n": int(len(v)), "V": float(v.sum()),
                       "E_mpa": float(getattr(mat, "elastic_mpa", 0) or 0) or None}
        tot_v += float(v.sum())
        if zz.size:
            zmin = min(zmin, float(zz.min()))
            zmax = max(zmax, float(zz.max()))
        mats.append(mat)
    L = max(zmax - zmin, 1e-9)
    A = tot_v / L
    # E 的字段名在 Material 上不固定 ⇒ 逐个试；取不到则由 main 用"两骨皮质骨同材料 ⇒ E 约掉"
    E = None
    for attr in ("elastic_mpa", "young_modulus_mpa", "E_mpa", "modulus_mpa",
                 "tangent_modulus_mpa"):
        v = getattr(mats[0], attr, None) if mats else None
        try:
            v = float(v)
        except (TypeError, ValueError):
            v = None
        if v:
            E = v
            break
    return {"V_cort_mm3": tot_v, "L_cort_mm": L, "A_cort_mm2": A, "E_mpa": E,
            "EA_N": (A * E if E else None), "per_domain": per_dom}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args(argv)
    t0 = time.time()

    frozen = audit.load_frozen()
    heights = np.asarray(frozen["heights"], int)
    loads = np.asarray(frozen["loads"]["ankle_r__vert"], float)
    ref_idx = int(np.argmin(np.abs(heights - fm.REF_H)))
    F_ref = float(loads[ref_idx])

    parts = {p.bone: p for p in fm.build_parts()}
    out: dict = {"bones": {}}
    for bone in ("fibula_r", "tibia_r"):
        p = parts[bone]
        print(f"[spike] extracting {bone} ...", flush=True)
        mesh = fm.extract_mesh(p.spec, fm.MASTER)
        for sd in p.spec.solid_domains:            # 与 g9._build_run 同样绑定材料
            sd.material = fm.read_materials(fm.MASTER, sd.pids)[sd.pids[0]]
        g = cortical_geometry(p.spec, mesh)
        out["bones"][bone] = g
        print(f"[spike] {bone}: A_cort={g['A_cort_mm2']:.1f} mm²  "
              f"L={g['L_cort_mm']:.0f} mm  E={g['E_mpa']}  "
              f"EA={g['EA_N']/1e3 if g['EA_N'] else None} kN ({time.time()-t0:.0f}s)",
              flush=True)

    fib, tib = out["bones"]["fibula_r"], out["bones"]["tibia_r"]
    EA_f, EA_t = fib["EA_N"], tib["EA_N"]
    if EA_f and EA_t:                      # 两骨 E 不同才需要 EA
        share, e_note = EA_f / (EA_f + EA_t), "EA 比"
    else:                                  # 两骨同为皮质骨 ⇒ E 约掉 ⇒ 只需面积比
        share = fib["A_cort_mm2"] / (fib["A_cort_mm2"] + tib["A_cort_mm2"])
        e_note = "两骨皮质骨同材料 ⇒ E 约掉（面积比）"
    print(f"\n[spike] A_cort: fib={fib['A_cort_mm2']:.1f}  tib={tib['A_cort_mm2']:.1f} mm²"
          f"（EA_f={EA_f}, EA_t={EA_t}）")
    print(f"[spike] 并联分担 share = {share:.4f} ({share*100:.2f}%)  [{e_note}]",
          flush=True)

    alloc = share * F_ref
    crush_f = SIGMA_C_THUMS * fib["A_cort_mm2"]
    crush_t = SIGMA_C_THUMS * tib["A_cort_mm2"]
    share_crush_f = crush_f / F_ref
    share_crush_t = crush_t / F_ref
    print(f"[spike] @REF={fm.REF_H} m: ankle_r={F_ref:.0f} N ⇒ 腓骨分配 {alloc:.0f} N"
          f"（{alloc/F_ref*100:.1f}%）")
    print(f"[spike] 容量对照：屈曲(孤立 FE)={P_CR_ISOLATED:.0f} N ｜ "
          f"压溃(σ_c·A_f)={crush_f:.0f} N ｜ 实测孤立失效带 "
          f"{F_FAIL_EXP[0]:.0f}–{F_FAIL_EXP[1]:.0f} N")
    print(f"[spike] ⇒ 该 share 下腓骨占压溃容量 {alloc/crush_f*100:.0f}%、"
          f"占屈曲容量 {alloc/P_CR_ISOLATED:.1f}×")
    print(f"[spike] 压溃门槛 share（腓骨=100% 压溃）：{share_crush_f:.3f}；"
          f"胫骨=100% 压溃：{share_crush_t:.3f}")

    # 自检
    checks = {
        "share_in_0_1": bool(0.0 < share < 1.0),
        "share_vs_assumed_0.05": share / 0.05,
        "share_vs_lit_30pct": share / 0.30,
        "A_fib_vs_prev_measured_76mm2": fib["A_cort_mm2"] / 76.0,
        "limit_EA_f_inf": 1.0, "limit_A_f_0": 0.0,
        "alloc_below_crush": bool(alloc < crush_f),
        "alloc_above_p_cr": bool(alloc > P_CR_ISOLATED),
    }
    print(f"[spike] 自检：" + " ｜ ".join(f"{k}={v}" for k, v in checks.items()))

    # share 扫描 → 分配载荷 vs 三种容量
    sweep = []
    for sh in (0.005, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20, share, 0.25, 0.30):
        a = sh * F_ref
        sweep.append({"share": float(sh), "alloc_n": a,
                      "over_p_cr": a / P_CR_ISOLATED,
                      "pct_of_crush": a / crush_f})
        print(f"        share={sh:5.3f} alloc={a:8.0f} N  ={a/P_CR_ISOLATED:6.1f}×P_cr"
              f"  {a/crush_f*100:6.1f}% 压溃")

    res = {"meta": {"title": "M1-spike 并行构件分担比估计（闭式，实测 EA）",
                    "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    "model": "重叠段轴向并联 ⇒ share = EA_f/(EA_f+EA_t)（与几何无关）；"
                             "加 IOM 只会更高 ⇒ 本值是下界",
                    "ref_h_m": float(fm.REF_H), "ankle_r_ref_n": F_ref,
                    "p_cr_isolated_n": P_CR_ISOLATED,
                    "sigma_c_thums_mpa": SIGMA_C_THUMS,
                    "f_fail_exp_n": list(F_FAIL_EXP),
                    "assumed_share_g8_g9": 0.05,
                    "literature_share_max": 0.30},
           "bones": out, "share_parallel": share, "alloc_at_ref_n": alloc,
           "crush_fibula_n": crush_f, "crush_tibia_n": crush_t,
           "share_at_crush_fibula": share_crush_f,
           "checks": checks, "sweep": sweep, "elapsed_s": time.time() - t0}
    (RES / f"m1_spike_share_estimate{args.out_suffix}.json").write_text(
        json.dumps(res, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n[spike] wrote: {RES / f'm1_spike_share_estimate{args.out_suffix}.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())