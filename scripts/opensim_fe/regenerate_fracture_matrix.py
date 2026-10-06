"""Regenerate the fracture matrix from CACHED records (no FEBio / no OpenSim).

Purpose
-------
After fixing the load-reduction scaling in
``fracture_matrix.compute_matrix`` (``sigma_vm_ref`` and the per-height ``sig``
now multiply by ``load_reduction_factor``), the 9x50 matrix must be recomputed.
The normal entry point (``fracture_matrix.main``) reuses the cached per-bone FE
records *but still calls* ``find_febio(verify=True)`` + ``probe_febio`` to fill
the metadata, which launches the FEBio executable.  This driver instead replays
the matrix computation purely from cache:

    temp/opensim_fe/fracture_matrix/load_chain.npz          (per-height loads)
    temp/opensim_fe/fracture_matrix/<bone>/<bone>_fe.json   (FE solves)

and reuses the solver metadata (``febio_exe`` / ``febio_version`` / ...) already
stored in ``results/opensim_fe/fracture_matrix.json``, so FEBio is never invoked.

It writes the same products as ``main``:
    results/opensim_fe/fracture_matrix.{json,csv,png}
    results/opensim_fe/FRACTURE_MATRIX_REPORT.md

Reproduce (repo root, PYTHONPATH=src):
    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\regenerate_fracture_matrix.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parent.parent
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import fracture_matrix as fm  # noqa: E402


def load_fe_cache(parts: list[fm.Part]) -> dict[str, dict]:
    """Load each bone's cached FE record; fail loudly if a cache is missing."""
    fe: dict[str, dict] = {}
    for part in parts:
        cache = fm.OUT_TEMP / part.bone / f"{part.bone}_fe.json"
        if not cache.is_file():
            raise SystemExit(f"缺少缓存 FE 记录，无法在无 FEBio 下重放: {cache}")
        fe[part.bone] = json.loads(cache.read_text(encoding="utf-8"))
    return fe


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stat", choices=list(fm.STATS), default=fm.PRIMARY_STAT,
                    help="primary σ_vm statistic for the matrix (default p95)")
    args = ap.parse_args(argv)

    fm.OUT_TEMP.mkdir(parents=True, exist_ok=True)
    fm.OUT_RESULTS.mkdir(parents=True, exist_ok=True)
    log = fm.Log(fm.OUT_TEMP / "regen_run.log")
    t_start = time.time()

    parts = fm.build_parts()
    log(f"[regen] {len(parts)} 部位；从缓存重放（无 FEBio / 无 OpenSim）")
    loads = fm.run_load_chain(log, reuse=True)
    fe = load_fe_cache(parts)
    log(f"[regen] 载入 {len(fe)} 个缓存 FE 记录")

    mat = fm.compute_matrix(parts, fe, loads, args.stat, log)

    # Reuse solver metadata from the existing JSON so we never launch FEBio.
    prev_path = fm.OUT_RESULTS / "fracture_matrix.json"
    prev_meta: dict = {}
    if prev_path.is_file():
        prev_meta = json.loads(prev_path.read_text(encoding="utf-8")).get("meta", {})
    ref_idx = mat["ref_idx"]
    ref_load_by_joint = {j: float(loads[f"{j}__vert"][ref_idx]) for j in fm.JOINTS}
    meta = dict(prev_meta)
    meta.update({
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "primary_stat": args.stat,
        "ref_load_by_joint_n": ref_load_by_joint,
        "total_seconds": time.time() - t_start,
        "regen_note": "recomputed from cached FE records + cached load chain; "
                      "no FEBio/OpenSim run",
    })

    fm.write_outputs(mat, parts, loads, fe, meta, log)
    fm.write_report(mat, parts, loads, fe, meta, log)
    log(f"[done] 从缓存重放完成，耗时 {time.time()-t_start:.1f}s")
    log.save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
