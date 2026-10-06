# S5-contact wave 4 — implementation log (executed)

> Implements the APPROVED spec `docs/S5_wave4_design.md` (§0 conventions, §1 schema,
> §2 report, §3 matrix, §4 filename, §5.1 resolutions). Do not re-litigate a resolved fork.
> Units **mm–N–MPa–s**. Interpreter `.venv\Scripts\python.exe`, `$env:PYTHONPATH="src"`.

## 0. Reading order actually completed (durable order #1)

1. `docs/S5_wave4_design.md` (full, incl. §5.1 all 10 resolutions + Appendix B NOT-DO).
2. `docs/S5_contact_contract.md` §4 (G0–G7), §5 (red lines), §6 (report format).
3. `docs/S5_QUOTABLE.md` (only numbers citable).
4. `docs/S5_contact_udg.md` §7 (validated NEW PENALTY contact recipe).
5. `scripts/opensim_fe/nonvertical_s5_supination.py` (boilerplate mirror).
6. `docs/S5_wave3_contact.md` (what wave-3 actually landed).
7. `src/climbing/coupling/plantar_bc.py`, `scripts/opensim_fe/thums_feb.py`,
   `src/climbing/coupling/fe_post.py`, `src/climbing/coupling/febio_run.py`,
   `src/climbing/coupling/pad_contact.py` (post-processing recipe),
   `scripts/opensim_fe/nonvertical_s4.py` (FE precedent), `tests/test_nonvertical_s4.py`,
   `results/opensim_fe/figures/FIGURE_CAPTIONS.md` (label taxonomy),
   `results/opensim_fe/s2_thums_height_sweep.json` (per-height subtalar loads).

## 1. Binding decisions applied (from §5.1 — not re-decided)

| # | Decision | Implementation |
|---|---|---|
| 1 | NEW `PENALTY 0.1` only | `laugon="PENALTY"`, `penalty=0.1`, `search_radius=20`, `node_reloc=0`, `tolerance=0.005`, `fric_coeff=μ`; one recipe block only |
| 2 | Really call FEBio | `climbing.coupling.febio_run.run_febio` |
| 3 | μ ∈ {0.0,0.3,0.6} | 3 contact μ rows at h=2.0; μ=0.6 at all 4 heights |
| 4 | G1 + G7 both reported, distinct codes | `g1_pad_stress_reproduced → N/A_minimum_version`; `g7_penetration_window → RETIRED_by_contract_条6-8` |
| 5 | calcaneus only | single bone, THUMS path |
| 6 | `.feb`→temp, report+JSON→results | `temp\opensim_fe\s5_contact\`; `results\opensim_fe\` |
| 7 | `hard_surface` rows in scope | 2 hard-surface rows (fixed + contact μ=0.0) at h=2.0 for G2/G5 |
| 8 | only avg/max/p95 σ_contact | `sigma_contact_avg_mpa` (=F_n/60000), `_max_mpa`, `_p95_mpa` |
| 9 | `friction_cone_ratio_max` scalar time-max | one number per run |
| 10 | no thermal | none |

## 2. Architectural decisions taken (not forks — implementation detail)

- **Per-height applied load**: read (read-only) `results/opensim_fe/s2_thums_height_sweep.json`
  `act==0` rows → `subtalar_kN`; linearly interpolate to `HEIGHTS` (only 2/3/4/5 present;
  3.5 from 3–4, 4.5 from 4–5). This is the axial subtalar load; label `measured` (OpenSim).
- **μ passthrough**: `thums_feb.build_thums_feb` gained ONE opt-in kwarg `contact_mu`
  (default `plantar_bc.DEFAULT_CONTACT_MU=0.6`) forwarded to `apply_plantar_bc`. Defaults
  unchanged ⇒ G0 byte-identical. This is "using wave-3 capability", not reimplementing it.
- **Contact force extraction**: the wave-3 rigid-plate nodes are fully fixed; their summed
  reaction = contact force. Plate node ids are the trailing nodes added by `_build_pad_floor`
  ⇒ read `reaction forces`, take the last `n_plate` rows. `F_n=Σ(−r_y)`, `F_t=‖Σ(r_x,r_z)‖`.
- **σ_contact_max/p95**: plate **top** node normal reaction / tributary area (nodal contact
  pressure field, `measured`). `avg = F_n/60000` (QUOTABLE basis, noted as rescale).
- **G0 hard gate**: default axial OpenSim link (`ground_reaction`→`run_dead_drop`→
  `subtalar_reaction`) vs `s1_h5_opensim.json`, threshold 1e-6; PLUS FE fixed-baseline
  reproduction vs `nonvertical_s4.json::cases[axial_baseline_fixed]` (max=213.502779,
  gauge=186.756369).
- **Overwrite guard**: exit code 2 on any pre-existing `--out-json`/`--out-md`.

## 3. Step log (appended as each step landed)

- [x] Step 0 — design + all spec files read.
- [x] Step 1 — this log created; todos written.
- [x] Step 2 — `thums_feb.build_thums_feb` opt-in `contact_mu` kwarg (default = `DEFAULT_CONTACT_MU`).
- [x] Step 3 — FEBio probes (probe1/2/4/6) → **finding**: wave-3's tri3 plantar surface is invalid
      on the hex8 CORT mesh (`292 invalid facets`); rebuilt as real hex8 quad4 faces and re-probed
      → `invalid facets` gone but the deck still diverges (`negative jacobians`) independent of
      penalty ∈ {0.1,1,10,1000} / load ∈ {2 kN,21 kN} / rigid vs pressure.
- [x] Step 4 — `scripts/opensim_fe/nonvertical_s5_contact.py` implemented (mirrors §0 sibling conventions).
- [x] Step 5 — matrix run (14 rows; each appended immediately). Outcomes:

| id | BC | rc | verdict | diagnostics |
|---|---|---|---|---|
| g0_fixed_axial_ref | fixed | 0 | PASS | gauge=186.7563689967737 (== s4 baseline, rel 0) |
| h2_fixed_pad | fixed | 0 | PASS | — |
| h3_fixed_pad | fixed | 0 | PASS | — |
| h3.5_fixed_pad | fixed | 0 | PASS | — |
| h4.5_fixed_pad | fixed | 0 | PASS | — |
| h2.0_spring_k100_pad | spring | 0 | PASS | gauge=122.66 |
| h2.0_contact_rigid_mu0_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h2.0_contact_rigid_mu0.3_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h2.0_contact_rigid_mu0.6_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h3_contact_rigid_mu0.6_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h3.5_contact_rigid_mu0.6_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h4.5_contact_rigid_mu0.6_pad | contact | 1 | FAIL | negJac=6, invalid=4 |
| h2.0_fixed_hard_surface | fixed | 0 | PASS | — |
| h2.0_contact_rigid_mu0.0_hard_surface | contact | 1 | FAIL | negJac=6, invalid=4 |

- [x] Step 6 — JSON + report written; schema verified (11 top-level keys, no `sigma_over_sigma_law`,
      G1=`N/A_minimum_version`, G7=`RETIRED_by_contract_条6-8`, 13 report sections).
- [x] Step 7 — `tests/test_nonvertical_s5_contact.py` (13 tests, no FEBio).
- [x] Step 8 — five hard gates: **93 / 19+1 / 11 / 11 / 13** all green.
- [x] Step 9 — no `.xplt`/`.hdf5` > 5 MB existed (each ≈0.4 MB); every `.feb` kept.

## 4. Findings (real, measured)

1. **G0 PASS, rel_err = 0.0**: the FE fixed axial reference run (`use_rigid=False`, load
   26383.133775894516 N) reproduces `nonvertical_s4.json::cases[axial_baseline_fixed]`
   gauge exactly (186.7563689967737). Default path is byte-for-byte preserved.
2. **Contact rows all diverge** (7/7, rc=1, `negative jacobians`). Root cause chain:
   `mesh['surfaces']['plantar']` is **292 tri3** (half-faces of CORT hex8 quads) → FEBio:
   `292 invalid facets` (hex8 needs quad4) → contact never engages → the soft-penalty deck
   plunges → element inversion. Rebuilding the surface as real hex8 quad4 removed the
   `invalid facets` error but the deck still inverts at step 1 for every penalty/load tried.
3. G2–G6 therefore **FAIL** with real diagnostics; G1/G7 carry their required distinct codes.

