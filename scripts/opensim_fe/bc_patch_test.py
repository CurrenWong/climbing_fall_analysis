"""bc_patch_test.py -- B1.2 positive control for the new BC (align_pca=True).

30-degree tilted square prism (10x10x40 mm, hex8 structured, E=1000 MPa,
nu=0.3) under F=1000 N end load.  Analytic uniaxial solution:
sigma_nom = F/A = 1000/100 = 10.00 MPa.

Two conditions (KICKOFF results/opensim_fe/B1_KICKOFF.md section C.1):
  * OLD BC  (align_pca=False, axis=None): negative control -- load is applied
    along the global bbox axis, which is 30 deg off the bar axis -> bending.
    Recorded only; no thresholds (it proves the checker can detect bad BC).
  * NEW BC  (align_pca=True): positive control -- bar axis aligned to +Z,
    axial loading.  Must PASS P1-P5 (see CRITERIA below).

CRITERIA (NEW BC; any FAIL -> exit 1):
  P1  mid-region (outside both builder 30% end bands) mean |sigma_zz| vs
      sigma_nom: |mean| is taken because the load is compressive (sigma_zz<0);
      rel err <= 5% (target <= 2%)
  P2  mid-region CV = std/mean of sigma_zz <= 10%
  P3  mid-region median(sigma_vm)/sigma_nom in [0.8, 1.2]  (canonical
      patch-test interior region; GLOBAL median recorded as diagnostic only
      because the end bands are 60% of the volume and legitimately deviate)
  P4  global p95(sigma_vm)/sigma_nom <= 1.5
  P5  |sum(NodalForce vectors)| == F  (rel 1e-6)  -- build-time load self-check

sigma_zz is read from the hdf5 Voigt stress (zz = component index 2), the same
source post_process uses; sigma_vm is NEVER compared against F/A.

Outputs:
  temp/fe_out/b1/patch_test.json          (both conditions, verdicts)
  temp/fe_out/b1/patch_test/{old,new}/    (.feb/.xplt/.hdf5 per condition)

Run (repo root, PYTHONPATH=src):
  & .venv\\Scripts\\python.exe scripts\\opensim_fe\\bc_patch_test.py
  # tee raw output to temp/fe_out/b1/logs/patch_test_<variant>.log
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent  # scripts/opensim_fe -> repo root
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import bone_feb as bf  # noqa: E402

# ------------------------------------------------------------------ geometry
LX = LY = 10.0          # mm cross-section
LZ = 40.0               # mm length
NX = NY = 10            # cells per cross-section side
NZ = 40                 # cells along length
TILT_DEG = 30.0         # rotate about Y to tilt the bar off the global axis
F_N = 1000.0            # N applied load
E_MPA = 1000.0
NU = 0.3
SIGMA_NOM = F_N / (LX * LY)   # = 10.00 MPa
TOP_FRAC = BOT_FRAC = 0.30    # must match fracture_matrix _build_run
PID = 1
MESH_NAME = "patch_solid"

OUT_DIR = ROOT / "temp" / "opensim_fe" / "b1"
JSON_PATH = OUT_DIR / "patch_test.json"

CRITERIA = {"P1_rel_err": 0.05, "P2_cv": 0.10,
            "P3_median_lo": 0.8, "P3_median_hi": 1.2, "P4_p95": 1.5,
            "P5_rel": 1e-6}


def build_prism() -> tuple[np.ndarray, np.ndarray]:
    """Structured hex8 prism centered at origin (span z: -20..20).

    Returns (nodes (N,3) float64, hexes (M,8) int64 1-based).
    Node index: i + (NX+1)*j + (NX+1)*(NY+1)*k (0-based).
    """
    xs = np.linspace(-LX / 2, LX / 2, NX + 1)
    ys = np.linspace(-LY / 2, LY / 2, NY + 1)
    zs = np.linspace(-LZ / 2, LZ / 2, NZ + 1)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    nodes = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)

    def nid(i: int, j: int, k: int) -> int:
        # C-order strides of meshgrid(..., indexing="ij") with shape (NX+1, NY+1, NZ+1)
        return i * (NY + 1) * (NZ + 1) + j * (NZ + 1) + k

    hexes = np.empty((NX * NY * NZ, 8), dtype=np.int64)
    e = 0
    for k in range(NZ):
        for j in range(NY):
            for i in range(NX):
                hexes[e] = [
                    nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k),
                    nid(i, j + 1, k),
                    nid(i, j, k + 1), nid(i + 1, j, k + 1),
                    nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1),
                ]
                e += 1
    return nodes, hexes + 1  # 1-based for FEBio


def tilt_y(nodes: np.ndarray, deg: float) -> np.ndarray:
    a = np.deg2rad(deg)
    c, s = np.cos(a), np.sin(a)
    R = np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])
    return nodes @ R.T


def make_mesh(nodes: np.ndarray, hexes: np.ndarray) -> dict:
    return {
        "nodes": nodes,
        "by_solid": {MESH_NAME: {
            "tets": np.zeros((0, 4), dtype=np.int64),
            "hexes": hexes,
            "pids": [PID] * len(hexes),
            "n_solid": int(len(hexes)),
            "mat_name": MESH_NAME,
        }},
        "by_shell": {},
        "n_nodes": int(len(nodes)),
    }


def make_spec() -> bf.BoneSpec:
    mat = bf.Material(pid=PID, E_mpa=E_MPA, nu=NU, src_keyword="*MAT_ELASTIC",
                      src_card=[E_MPA, NU], note="patch test E=1000 nu=0.3")
    return bf.BoneSpec(name="patch_prism", solid_domains=[bf.SolidDomainSpec(
        pids=[PID], mat_name=MESH_NAME, mesh_name=MESH_NAME, material=mat)])


def read_voigt(hdf5: Path) -> np.ndarray:
    """Last-state element Voigt stress (n_elem, 6) for our single domain."""
    import h5py
    with h5py.File(hdf5, "r") as h5:
        sids = sorted(h5["states"].keys(), key=int)
        grp = h5["states"][sids[-1]]["element_data"]["stress"]
        if MESH_NAME not in grp:
            raise RuntimeError(f"stress key {MESH_NAME!r} missing; have {list(grp)}")
        return np.asarray(grp[MESH_NAME][:], dtype=np.float64)


def mid_mask(nodes: np.ndarray, hexes: np.ndarray, axis: int) -> np.ndarray:
    """Element mask strictly outside both builder end bands (30% each)."""
    cent = nodes[hexes - 1].mean(axis=1)          # (M,3)
    lo, hi = nodes[:, axis].min(), nodes[:, axis].max()
    span = hi - lo
    z = cent[:, axis]
    return (z > lo + BOT_FRAC * span) & (z < hi - TOP_FRAC * span)


def run_condition(tag: str, align_pca: bool, nodes0: np.ndarray,
                  hexes: np.ndarray) -> dict:
    work = OUT_DIR / tag
    work.mkdir(parents=True, exist_ok=True)
    feb = work / f"patch_{tag}.feb"
    t0 = time.time()

    mesh = make_mesh(nodes0.copy(), hexes)
    spec = make_spec()
    meta = bf.build_bone_feb(spec, feb, load_n=F_N, top_frac=TOP_FRAC,
                             bot_frac=BOT_FRAC, axis=None, with_shell=False,
                             time_steps=1, mesh=mesh, align_pca=align_pca)
    # P5: |sum(NodalForce)| == F from build-time meta
    lv = np.asarray(meta["load_vec_per_node"], dtype=np.float64)
    n_top = int(meta["n_top_nodes"])
    f_sum = float(np.linalg.norm(lv * n_top))

    run = bf.run_bone_feb(feb, workdir=work, timeout=900.0)
    pp = bf.post_process(feb.with_suffix(".xplt"), feb.with_suffix(".hdf5"),
                         solid_names=[MESH_NAME], shell_specs=[], label=tag)
    dom = pp.get("domains", {}).get(f"solid_{MESH_NAME}", {})
    if not dom:
        raise RuntimeError(f"[{tag}] no domain stats: {pp.get('domains')}")

    voigt = read_voigt(feb.with_suffix(".hdf5"))
    if voigt.shape[0] != len(hexes):
        raise RuntimeError(f"[{tag}] stress n_elem {voigt.shape[0]} != {len(hexes)}"
                           " (element order/count mismatch)")
    # geometry used by THIS build (align_pca rotates a copy internally):
    nodes_used = mesh["nodes"]
    if align_pca:
        nodes_used = bf._align_pca_rotate(nodes_used)
    axis = int(meta["axis_idx"])
    mask = mid_mask(nodes_used, hexes, axis)
    sig_zz = voigt[mask, 2]
    vm = bf._von_mises_from_voigt(voigt)
    vm_mid = vm[mask]

    res = {
        "tag": tag, "align_pca": bool(align_pca),
        "axis_idx": axis, "n_top_nodes": n_top,
        "load_sum_n": f_sum, "P5_rel_err": abs(f_sum - F_N) / F_N,
        "run": run, "n_mid_elem": int(mask.sum()),
        "mid_sigma_zz_mean_mpa": float(sig_zz.mean()),
        "mid_sigma_zz_std_mpa": float(sig_zz.std()),
        # loading is compressive (sigma_zz < 0); the analytic check is on the
        # MAGNITUDE |mean sigma_zz| vs sigma_nom = F/A (spec sign slip fixed).
        "mid_sigma_zz_rel_err": float(
            abs(abs(sig_zz.mean()) - SIGMA_NOM) / SIGMA_NOM),
        "mid_sigma_zz_cv": float(sig_zz.std() / abs(sig_zz.mean())),
        # canonical patch-test region = interior (outside end bands); the
        # global median is recorded as a diagnostic only (bands = 60% of the
        # volume and legitimately deviate).
        "mid_vm_median_mpa": float(np.median(vm_mid)),
        "global_vm_median_mpa": float(np.median(vm)),
        "global_vm_p95_mpa": float(np.percentile(vm, 95)),
        "domain_stats": dom,
        "elapsed_s": round(time.time() - t0, 1),
    }
    res["P3_median_ratio"] = res["mid_vm_median_mpa"] / SIGMA_NOM
    res["global_median_ratio_diag"] = res["global_vm_median_mpa"] / SIGMA_NOM
    res["P4_p95_ratio"] = res["global_vm_p95_mpa"] / SIGMA_NOM
    return res


def judge(res: dict) -> dict:
    c = CRITERIA
    checks = {
        "P1": res["mid_sigma_zz_rel_err"] <= c["P1_rel_err"],
        "P2": res["mid_sigma_zz_cv"] <= c["P2_cv"],
        "P3": c["P3_median_lo"] <= res["P3_median_ratio"] <= c["P3_median_hi"],
        "P4": res["P4_p95_ratio"] <= c["P4_p95"],
        "P5": res["P5_rel_err"] <= c["P5_rel"],
    }
    checks["ALL"] = all(checks.values())
    return checks


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nodes0, hexes = build_prism()
    nodes_tilt = tilt_y(nodes0, TILT_DEG)
    print(f"[patch] prism {len(nodes0)} nodes / {len(hexes)} hexes, "
          f"tilt {TILT_DEG} deg, sigma_nom={SIGMA_NOM:.2f} MPa", flush=True)

    out: dict = {
        "spec": {"prism_mm": [LX, LY, LZ], "cells": [NX, NY, NZ],
                 "tilt_deg": TILT_DEG, "F_N": F_N, "E_mpa": E_MPA, "nu": NU,
                 "sigma_nom_mpa": SIGMA_NOM, "top_frac": TOP_FRAC,
                 "bot_frac": BOT_FRAC, "criteria": CRITERIA},
        "conditions": {},
    }
    rc = 0
    for tag, align in (("old", False), ("new", True)):
        print(f"\n[patch] === condition {tag} (align_pca={align}) ===",
              flush=True)
        try:
            res = run_condition(tag, align, nodes_tilt, hexes)
            if tag == "new":
                res["checks"] = judge(res)
                if not res["checks"]["ALL"]:
                    rc = 1
            out["conditions"][tag] = res
            line = (f"[patch:{tag}] axis={res['axis_idx']} "
                    f"mid sigma_zz mean={res['mid_sigma_zz_mean_mpa']:.3f} "
                    f"rel_err={res['mid_sigma_zz_rel_err']*100:.2f}% "
                    f"cv={res['mid_sigma_zz_cv']*100:.2f}% | "
                    f"mid median={res['mid_vm_median_mpa']:.3f} "
                    f"({res['P3_median_ratio']:.3f}x) p95="
                    f"{res['global_vm_p95_mpa']:.3f} ({res['P4_p95_ratio']:.3f}x) "
                    f"| glob_med_diag={res['global_median_ratio_diag']:.3f}x "
                    f"P5 rel={res['P5_rel_err']:.2e} t={res['elapsed_s']}s")
            if tag == "new":
                line += f" checks={res['checks']}"
            print(line, flush=True)
        except Exception as exc:  # noqa: BLE001
            import traceback
            out["conditions"][tag] = {
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(),
            }
            print(f"[patch:{tag}] ERROR: {type(exc).__name__}: {exc}",
                  flush=True)
            if tag == "new":
                rc = 1

    verdict = "PASS" if (
        out["conditions"].get("new", {}).get("checks", {}).get("ALL")
    ) else "FAIL"
    out["verdict"] = verdict
    JSON_PATH.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                         encoding="utf-8")
    print(f"\n[patch] verdict={verdict} -> {JSON_PATH}", flush=True)
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
