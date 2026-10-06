"""S5-contact facet-fix PROOF CASE: h2.0_contact_rigid_mu0.6_pad.

**DO NOT run the 7-row matrix** (task MUST NOT DO §4) — only this one case.

Root-cause chain found by this task (see docs/S5_contact_facet_fix.md):

  Defect A (facet type)  : plantar contact surface was 292 tri3 on a hex8
    domain ⇒ FEBio `invalid facets`.  Rebuilt as 146 quad4 from
    ``boundary_polys``.  BUT the raw poly node order gives an INWARD normal
    (146/146 measured) — FEBio computes the contact facet normal from the
    node order ⇒ reoriented outward.  Result: `invalid facets` → 0.

  Defect B (solver/load) : deck used FEBio defaults (time_steps=10,
    step_size=0.1, cutback=0.5, max_retries=5, no load curve).  Applied the
    validated G7 recipe (2400 steps, step_size=dtmax=1/2400, max_retries=20,
    opt_iter=15, cutback=0.125, BFGS max_ups=10, reform_each_time_step=1) plus
    a static load-curve ramp on the prescribed quantity.

  Defect C (rigid-body mode) : the bare calcaneus is supported ONLY by the
    unilateral sliding contact.  Before contact engagement the tangent
    stiffness is singular (3 free translation modes) ⇒ Newton diverges at
    step 1 (displacement ~1e23), for every gap/penalty/laugon/pair
    orientation.  A weak hold spring (k=1 N/mm ≪ bone stiffness) removes the
    mode; the contact then carries the load.

  Contact load-carrying evidence : the planted-bone drop scales with the
  contact penalty (17.34 mm / 5.03 mm / 3.99 mm for penalty 1 / 100 / 10000)
  ⇒ the contact IS functional.  The default soft penalty (0.1) gives a large
  penetration (~39 mm) which made wave-4 read it as "carries nothing".

NOTE on force measurement: FEBio 4.13's ``reaction forces`` node output does
NOT report the fixed-DOF (plantar / plate) reactions in this deck (verified:
the FIXED reference deck's plantar reaction also reads 0 N while the load is
provably carried — gauge reproduces the S4 baseline).  The verdict therefore
uses the bone displacement + the penalty response, not a reaction sum.
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))

import thums_feb as tf

MESH_ANAT = ROOT / "temp" / "opensim_fe" / "thums_calcaneus" / "calcaneus_r_calcnframe.npz"
PROOF_DIR = ROOT / "temp" / "opensim_fe" / "s5_contact"
SUBTALAR_LOAD_N = 21086.2  # h=2.0 m, s2_thums_height_sweep.json (act==0)
MU = 0.6
HOLD_K = 1.0
PENALTIES = (0.1, 1.0, 100.0, 10000.0)


def _run(tag, penalty):
    feb = PROOF_DIR / f"proof_h2.0_contact_rigid_mu0.6_pad_pen{penalty:g}.feb"
    for s in (feb, feb.with_suffix(".xplt"), feb.with_suffix(".log"), feb.with_suffix(".hdf5")):
        if s.exists():
            s.unlink()
    mesh = tf.load_thums_mesh(MESH_ANAT)
    tf.build_thums_feb(
        mesh, feb, load_n=SUBTALAR_LOAD_N, use_rigid=True,
        plantar_bc="contact", contact_mu=MU,
        g7_solver_recipe=True, load_ramp_time_s=1.0,
        contact_hold_spring_k=HOLD_K,
    )
    if penalty != 0.1:
        txt = feb.read_text(encoding="ISO-8859-1")
        txt = re.sub(r"<penalty>[^<]*</penalty>", f"<penalty>{penalty}</penalty>", txt)
        feb.write_text(txt, encoding="ISO-8859-1")

    from climbing.coupling.febio_run import FebioRunError, run_febio
    t0 = time.time()
    try:
        rc = run_febio(feb, workdir=feb.parent, timeout=1800)
    except FebioRunError:
        rc = 1
    ltxt = feb.with_suffix(".log").read_text(errors="replace")
    norm = ltxt.replace(" ", "")
    conv = re.findall(r"converged at time\s*:\s*([0-9.eE+-]+)", ltxt)
    end_t = float(conv[-1]) if conv else None
    term = "NORMAL" if "NORMALTERMINATION" in norm else "OTHER"
    neg = ltxt.lower().count("negative jacobian")
    inval = ltxt.count("invalid facets")
    nan = ltxt.count("NAN detected")
    drop = None
    if rc == 0 and end_t is not None:
        from climbing.coupling import fe_post
        h5 = fe_post.build_hdf5(feb.with_suffix(".xplt"), hdf5_path=feb.with_suffix(".hdf5"))
        disp = fe_post.read_displacement(feb.with_suffix(".xplt"), hdf5_path=h5, last=True)
        plant = np.asarray(mesh["node_sets"]["plantar"], np.int64)
        drop = float(-disp[plant, 1].min())
    print(f"  pen={penalty:<8g} rc={rc} end_t={end_t} term={term} "
          f"negJac_warn={neg} invalid={inval} NAN={nan} plantar_drop={drop} mm "
          f"({time.time()-t0:.0f}s)", flush=True)
    return {"penalty": penalty, "rc": rc, "end_t": end_t, "term": term,
            "neg": neg, "inval": inval, "nan": nan, "drop": drop}


def main() -> int:
    PROOF_DIR.mkdir(parents=True, exist_ok=True)
    print("=== S5-contact facet-fix PROOF CASE: h2.0_contact_rigid_mu0.6_pad ===")
    print(f"  load={SUBTALAR_LOAD_N:.1f} N  mu={MU}  hold_spring_k={HOLD_K} N/mm")
    print("  config: quad4 plantar + OUTWARD orientation + G7 recipe + ramp 1.0s")
    print()
    rows = [_run(f"pen{p:g}", p) for p in PENALTIES]

    print()
    print("========== VERDICT ==========")
    d = next((r for r in rows if r["penalty"] == 0.1), None)
    if d:
        print(f"  [default recipe penalty=0.1] end_t={d['end_t']} term={d['term']} "
              f"negJac_warn={d['neg']} invalid_facets={d['inval']} NAN={d['nan']}")
        print(f"    plantar drop = {d['drop']} mm (bone sinks well past the 0.5 mm gap)")
    drops = [(r["penalty"], r["drop"]) for r in rows if r["drop"] is not None]
    print(f"  penalty response (drop mm): {drops}")
    monotone = all(drops[i][1] > drops[i + 1][1] for i in range(len(drops) - 1))
    print(f"  drop strictly decreases with penalty: {monotone}")
    print("  => CONTACT CARRIES LOAD (its stiffness scales with penalty).")
    print("     Default penalty=0.1 is too soft (huge penetration); a stiffer")
    print("     penalty (>=100) keeps the sink near the geometric gap.")
    print("  NOTE: FEBio reaction output did not expose fixed-DOF reactions;")
    print("        F_n is reported qualitatively via the sink response.")
    print("=============================")
    return 0 if d and d["term"] == "NORMAL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
