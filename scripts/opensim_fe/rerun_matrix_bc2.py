"""rerun_matrix_bc2.py -- B1.3: rerun the 9x50 fracture matrix under the NEW BC.

B1.1 decision 4): isolated roots, NEVER overwrite the old BC products:
  temp root    : temp/fe_out/fracture_matrix_bc2/     (FE caches, run.log)
  results root : results/opensim_fe/bc2/                  (renamed to
                 fracture_matrix_bc2.{json,csv,png} + FRACTURE_MATRIX_BC2_REPORT.md
                 at --finalize)
B1.1 decision 3): primary statistic = median (existing --stat semantics;
the old p95 stays the frozen legacy default in fracture_matrix.py).

How it works WITHOUT touching fracture_matrix.py on disk:
  * module-constant override: fm.OUT_TEMP / fm.OUT_RESULTS -> new roots;
  * _build_run override: same body as fracture_matrix._build_run but calls
    build_bone_feb(..., align_pca=True)  (the B1.2 opt-in, patch-test gated);
  * everything else (load chain cache, descending-load fallback, linearity
    check, compute_matrix / write_outputs / write_report) is REUSED as-is.

The load-chain cache is COPIED into the new root at --init (heights/joints
verified), so run_load_chain(reuse=True) hits the copy (~21 min saved) and
never writes to the old root.

Chunked execution (channel-resilient, one foreground call per step):
  1) --init                      copy + verify load_chain.npz, provenance json
  2) --bone <key>   (x9, any order)  solve+cache one bone (reuse=True)
  3) --finalize                  compute_matrix(median) + outputs + rename
  4) --diff                      compare vs legacy fracture_matrix.json

Run (repo root, PYTHONPATH=src):
  & .venv\\Scripts\\python.exe scripts\\opensim_fe\\rerun_matrix_bc2.py --init
  & .venv\\Scripts\\python.exe scripts\\opensim_fe\\rerun_matrix_bc2.py --bone calcaneus_r
  ...
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import fracture_matrix as fm  # noqa: E402
from bone_feb import build_bone_feb, run_bone_feb  # noqa: E402

NEW_TEMP = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"
NEW_RESULTS = ROOT / "results" / "opensim_fe" / "bc2"
OLD_RESULTS = ROOT / "results" / "opensim_fe"
B1_LOG = ROOT / "results" / "opensim_fe" / "B1_LOG.md"
BC_MODE = "align_pca=True (B1.2 opt-in; patch test P1-P5 PASS)"


def log_line(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M")
    with B1_LOG.open("a", encoding="utf-8") as fh:
        fh.write(f"[{stamp}] {msg}\n")
    print(f"[B1_LOG] {msg}", flush=True)


class AppendingLog:
    """fm.Log-compatible logger that appends across chunk invocations."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.lines: list[str] = []

    def __call__(self, msg: str) -> None:
        print(msg, flush=True)
        self.lines.append(msg)

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write("\n".join(self.lines) + "\n")
        self.lines = []


def apply_patches() -> None:
    """Override module constants + _build_run (in-memory only)."""
    fm.OUT_TEMP = NEW_TEMP
    fm.OUT_RESULTS = NEW_RESULTS

    def _build_run_bc2(part, mesh, mats, feb, load, time_steps):
        for sd in part.spec.solid_domains:
            sd.material = mats[sd.pids[0]]
        for sh in part.spec.shell_domains:
            sh.material = mats[sh.pids[0]]
        build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30,
                       bot_frac=0.30, axis=None, with_shell=True,
                       time_steps=time_steps, mesh=mesh, align_pca=True)
        run_bone_feb(feb, workdir=feb.parent, timeout=1800)

    fm._build_run = _build_run_bc2


def get_loads(log: AppendingLog) -> dict:
    loads = fm.run_load_chain(log, reuse=True)
    return loads


def select_part(key: str):
    parts = fm.build_parts()
    for p in parts:
        if key in (p.bone, p.cn):
            return p, parts
    raise SystemExit(f"unknown part {key!r}; keys: "
                     f"{[p.bone for p in parts]}")


def cmd_init() -> int:
    NEW_TEMP.mkdir(parents=True, exist_ok=True)
    NEW_RESULTS.mkdir(parents=True, exist_ok=True)
    src = fm.OUT_TEMP / "load_chain.npz"          # OLD root (pre-patch value)
    dst = NEW_TEMP / "load_chain.npz"
    if not dst.is_file():
        if not src.is_file():
            raise SystemExit(f"legacy load chain missing: {src}")
        shutil.copy2(src, dst)
    with np.load(dst) as z:
        ok = (np.array_equal(z["heights"], fm.HEIGHTS)
              and all(f"{j}__vert" in z.files for j in fm.JOINTS))
        if not ok:
            raise SystemExit("copied load_chain.npz does NOT cover "
                             f"heights {fm.HEIGHTS[0]}..{fm.HEIGHTS[-1]} x {fm.JOINTS}")
    prov = {
        "bc_mode": BC_MODE,
        "primary_stat": "median",
        "temp_root": str(NEW_TEMP),
        "results_root": str(NEW_RESULTS),
        "load_chain_copied_from": str(src),
        "load_chain_verified": True,
        "fracture_matrix_py_modified": False,
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    (NEW_TEMP / "BC2_PROVENANCE.json").write_text(
        json.dumps(prov, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(prov, indent=2, ensure_ascii=False))
    log_line("B1.3: --init 完成 — 独立根 fracture_matrix_bc2 已建、"
             "load_chain.npz 拷贝并校验（heights×joints 全覆盖）、"
             f"BC={BC_MODE}、stat=median")
    return 0


def cmd_bone(key: str, time_steps: int) -> int:
    apply_patches()
    NEW_TEMP.mkdir(parents=True, exist_ok=True)
    log = AppendingLog(NEW_TEMP / "run.log")
    part, _ = select_part(key)
    t0 = time.time()
    loads = get_loads(log)
    ref_idx = int(np.argmin(np.abs(loads["heights"] - fm.REF_H)))
    rl = float(loads[f"{part.joint}__vert"][ref_idx])
    log(f"[bc2] {part.bone}: ref load={rl/1e3:.3f} kN (joint {part.joint})")
    rec = fm.solve_bone_fe(part, rl, log, reuse=True, time_steps=time_steps)
    log.save()
    status = rec.get("status")
    err = rec.get("error", "")
    dt = time.time() - t0
    msg = (f"B1.3: bone {part.bone} ({part.cn}) -> status={status} "
           f"solve_load={rec.get('fe_solve_load_n')} "
           f"lrf={rec.get('load_reduction_factor')} ({dt:.0f}s)")
    if err:
        msg += f" error={err[:200]}"
    log_line(msg)
    print(f"\n[bc2] {part.bone}: {status} in {dt:.0f}s")
    return 0 if status in ("ok", "no_domains") else 1


def cmd_finalize(time_steps: int, stat: str) -> int:
    apply_patches()
    NEW_TEMP.mkdir(parents=True, exist_ok=True)
    NEW_RESULTS.mkdir(parents=True, exist_ok=True)
    log = AppendingLog(NEW_TEMP / "run.log")
    t_start = time.time()
    parts = fm.build_parts()

    loads = get_loads(log)
    ref_idx = int(np.argmin(np.abs(loads["heights"] - fm.REF_H)))
    ref_load_by_joint = {j: float(loads[f"{j}__vert"][ref_idx]) for j in fm.JOINTS}

    fe: dict[str, dict] = {}
    for part in parts:
        rl = ref_load_by_joint[part.joint]
        rec = fm.solve_bone_fe(part, rl, log, reuse=True, time_steps=time_steps)
        fe[part.bone] = rec
        if rec.get("status") not in ("ok", "no_domains"):
            print(f"[bc2] ⚠ {part.bone} status={rec.get('status')} "
                  f"{rec.get('error', '')[:160]}", flush=True)

    mat = fm.compute_matrix(parts, fe, loads, stat, log)

    from climbing.coupling.febio_run import find_febio, probe_febio  # noqa: E402
    exe = find_febio(verify=True)
    meta = {
        "title": "1-50 m x 9-part fracture matrix -- BC2 (align_pca=True, median)",
        "bc_mode": BC_MODE,
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "load_mode": "per-height (ground_reaction -> run_dead_drop -> joint_reaction; "
                     "NO GRF-ratio shortcut); chain COPIED from legacy cache",
        "heights_m": [int(h) for h in fm.HEIGHTS],
        "ref_height_m": fm.REF_H,
        "joints": list(fm.JOINTS),
        "primary_stat": stat,
        "stats_reported": list(fm.STATS),
        "criterion": "fracture(h) = sigma_vm(h) >= sigma_c(part)",
        "fe_load_fallback": "try ref, then 10000/5000/2500/1000 N; scale back linearly",
        "model_mass_kg": fm.MODEL_MASS_KG,
        "master_k": str(fm.MASTER.relative_to(ROOT)),
        "febio_exe": str(exe),
        "febio_version": probe_febio(exe),
        "time_steps": time_steps,
        "ref_load_by_joint_n": ref_load_by_joint,
        "total_seconds": time.time() - t_start,
    }
    fm.write_outputs(mat, parts, loads, fe, meta, log)
    fm.write_report(mat, parts, loads, fe, meta, log)
    log.save()

    # rename legacy product names -> *_bc2 (old results untouched)
    renames = [("fracture_matrix.json", "fracture_matrix_bc2.json"),
               ("fracture_matrix.csv", "fracture_matrix_bc2.csv"),
               ("fracture_matrix.png", "fracture_matrix_bc2.png"),
               ("FRACTURE_MATRIX_REPORT.md", "FRACTURE_MATRIX_BC2_REPORT.md")]
    for src, dst in renames:
        s, d = NEW_RESULTS / src, OLD_RESULTS / dst
        if s.is_file():
            if d.is_file():
                d.unlink()
            s.rename(d)
            print(f"[bc2] renamed {src} -> {dst}", flush=True)
    log_line(f"B1.3: --finalize 完成 — stat={stat}、9 骨全缓存、"
             f"产物 results/opensim_fe/fracture_matrix_bc2.{{json,csv,png}} + "
             f"FRACTURE_MATRIX_BC2_REPORT.md（耗时 {time.time()-t_start:.0f}s）")
    return 0


def _spearman(a: list, b: list) -> float:
    from scipy.stats import spearmanr
    r = spearmanr(a, b).statistic
    return float(r)


def cmd_diff(stat_new: str) -> int:
    old_p = OLD_RESULTS / "fracture_matrix.json"
    new_p = OLD_RESULTS / "fracture_matrix_bc2.json"
    if not old_p.is_file() or not new_p.is_file():
        raise SystemExit(f"missing matrix json: old={old_p.is_file()} "
                         f"new={new_p.is_file()}")
    old = json.loads(old_p.read_text(encoding="utf-8"))
    new = json.loads(new_p.read_text(encoding="utf-8"))

    def rowmap(doc, stat):
        out = {}
        for r in doc["parts"]:
            sv = (r.get("sigma_vm_ref_mpa") or {}).get(stat)
            out[r["bone"]] = {
                "part_cn": r["part_cn"],
                "first": r.get(f"first_fracture_{stat}_m"),
                "sigma_ref": sv,
                "status": r.get("fe_status"),
            }
        return out

    o95, omed = rowmap(old, "p95"), rowmap(old, "median")
    nmed = rowmap(new, stat_new)

    L: list[str] = ["# B13_DIFF — 新 BC（align_pca, median） vs 旧矩阵", ""]
    L.append("## 1. 首次骨折高度对照（m；None=≤50m 未折）")
    L.append("")
    L.append("| 部位 | 旧 p95 | 旧 median | **新 median(BC2)** | Δ(新med−旧med) |")
    L.append("|---|---|---|---|---|")
    pairs_med, pairs_new = [], []
    for bone, nr in nmed.items():
        orr, omr = o95.get(bone, {}), omed.get(bone, {})
        d = (None if (nr.get("first") is None or omr.get("first") is None)
             else nr["first"] - omr["first"])
        L.append(f"| {nr['part_cn']} (`{bone}`) | {orr.get('first')} | "
                 f"{omr.get('first')} | **{nr.get('first')}** | {d} |")
        if nr.get("first") is not None and omr.get("first") is not None:
            pairs_new.append(nr["first"]); pairs_med.append(omr["first"])
    L.append("")

    # ranks of sigma_vm_ref (both medians; also new-median vs old-p95)
    bones = [b for b in nmed if b in omed and nmed[b].get("sigma_ref") is not None
             and omed[b].get("sigma_ref") is not None]
    rho_med = (_spearman([nmed[b]["sigma_ref"] for b in bones],
                         [omed[b]["sigma_ref"] for b in bones])
               if len(bones) >= 3 else float("nan"))
    bones2 = [b for b in nmed if b in o95 and nmed[b].get("sigma_ref") is not None
              and o95[b].get("sigma_ref") is not None]
    rho_p95 = (_spearman([nmed[b]["sigma_ref"] for b in bones2],
                         [o95[b]["sigma_ref"] for b in bones2])
               if len(bones2) >= 3 else float("nan"))
    L.append("## 2. σ_vm_ref 排序 Spearman ρ（≥3 对）")
    L.append("")
    L.append(f"- 新median vs **旧median**（同统计量，纯 BC 效应）：ρ = **{rho_med:.3f}** （n={len(bones)}）")
    L.append(f"- 新median vs **旧p95**（旧生产口径）：ρ = **{rho_p95:.3f}** （n={len(bones2)}）")
    L.append("")
    agree = sum(1 for b in nmed if nmed[b].get("first") == o95.get(b, {}).get("first"))
    L.append("## 3. 首折一致性")
    L.append("")
    L.append(f"- 新median vs 旧p95 完全同值：{agree}/{len(nmed)}")
    if len(pairs_new) >= 3:
        rho_h = _spearman(pairs_new, pairs_med)
        L.append(f"- 首折高度（双方均非 None，n={len(pairs_new)}）Spearman ρ = **{rho_h:.3f}**")
    L.append("")
    L.append("_生成：rerun_matrix_bc2.py --diff；数据源为两份 JSON 实测。_")

    out_p = OLD_RESULTS / "B13_MATRIX_BC2_DIFF.md"
    out_p.write_text("\n".join(L) + "\n", encoding="utf-8")
    log_line(f"B1.3: diff 完成 — ρ(新med,旧med)={rho_med:.3f}、"
             f"ρ(新med,旧p95)={rho_p95:.3f}；产物 {out_p.name}")
    print("\n".join(L))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--init", action="store_true")
    ap.add_argument("--bone", default=None,
                    help="bone key (e.g. calcaneus_r) or Chinese part name")
    ap.add_argument("--all-fe", action="store_true")
    ap.add_argument("--finalize", action="store_true")
    ap.add_argument("--diff", action="store_true")
    ap.add_argument("--stat", default="median", choices=list(fm.STATS))
    ap.add_argument("--time-steps", type=int, default=1)
    args = ap.parse_args(argv)

    if args.init:
        return cmd_init()
    if args.bone:
        return cmd_bone(args.bone, args.time_steps)
    if args.all_fe:
        rc = 0
        for p in fm.build_parts():
            rc = max(rc, cmd_bone(p.bone, args.time_steps))
        return rc
    if args.finalize:
        return cmd_finalize(args.time_steps, args.stat)
    if args.diff:
        return cmd_diff(args.stat)
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
