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

---

## 5. CORRECTION (2026-10-06 晚，wave-4 → wave-5 facet/solver fix)

> Appended by the follow-up task. **Wave-4 §4 findings #2/#3 above are amended below;
> the measurements in §3 (matrix outcomes) and §4 #1 (G0 PASS) stand.**
> Full plan + evidence: `docs\S5_contact_facet_fix.md`.

### 5.1 What wave-4 got right vs what it missed

* **Right**: the plantar surface was `292 tri3` on a **hex8** domain → FEBio `invalid facets`
  → contact never engaged. Rebuilding it as `quad4` removed the error.
* **Missed #1 (orientation)**: the raw `boundary_polys` **node order** yields an
  **INWARD** facet normal (measured 146/146 with `normal.y > 0`). FEBio computes a contact
  facet's normal from the **node order** (right-hand rule), *not* from the stored normal
  vector that `_polygon_centroids_normals` flips. An inward master surface makes the gap
  sign wrong ⇒ contact never closes. Fixed by reordering each emitted quad so its
  right-hand-rule normal matches the stored **outward** normal.
* **Missed #2 (solver/load)**: the deck used FEBio defaults (`time_steps=10`,
  `step_size=0.1`, `cutback=0.5`, `max_retries=5`, no `<LoadData>`). The validated G7
  recipe (`docs\S5_contact_udg.md` §7: 2400 steps, `dtmax=1/2400`, `cutback=0.125`,
  `max_retries=20`, `opt_iter=15`, BFGS `max_ups=10`, `reform_each_time_step=1`, plus a
  static load-curve ramp) is now applied (opt-in `g7_solver_recipe=True`,
  `load_ramp_time_s`). The **fixed** reference deck reproduces `end_t=1.0`, `negJac=0`.
* **Missed #3 (rigid-body mode)**: the bare calcaneus is supported **only** by the
  unilateral sliding contact. Before engagement the tangent stiffness is singular (3 free
  translations) ⇒ Newton diverges at step 1 (`displacement ≈ 2.3e25`), invariant to gap
  (0.5 → 1e-6), `node_reloc`, pair orientation, `penalty`, `laugon`, load,
  `search_radius`, `two_pass`, solver `symmetric_stiffness`, and ghost X/Z fixation.
  A weak hold spring (`contact_hold_spring_k`, k=1 N/mm ≪ bone stiffness) removes the
  mode; the contact then runs.

### 5.2 Corrected verdict on wave-4 §4 #2 ("contact never engages")

**The contact DOES carry load — wave-4 read it as "carries nothing" because the default
soft penalty (`0.1`, calibrated for a FOAM pad) gives a huge penetration on a 15000 MPa
cortical bone vs a rigid plate.** Measured on the single proof case
`h2.0_contact_rigid_mu0.6_pad` (2 m, μ=0.6), 2400 G7 steps, hold spring k=1:

| penalty | rc | end_t | negJac warnings | invalid facets | NAN | plantar sink (mm) |
|---|---|---|---|---|---|---|
| 0.1 (recipe default) | 0 | 1.0 | 75 | 0 | 0 | 39.16 |
| 1 | 0 | 1.0 | 6 | 0 | 0 | 17.34 |
| 100 | 0 | 1.0 | **0** | 0 | 0 | 5.03 |
| 10000 | 0 | 1.0 | **0** | 0 | 0 | 3.99 |

The sink **decreases monotonically with penalty** ⇒ the contact stiffness is active ⇒ the
contact carries load. `penalty ≥ 100` also drives the negative-jacobian warnings to **0**.

### 5.3 Measurement caveat (a hidden defect in the wave-4 `F_n` extraction)

FEBio 4.13's `reaction forces` **node output does not report the fixed-DOF (plantar/plate)
reactions in this deck**: the **fixed** reference deck's plantar reaction also reads
`0.00 N` while its gauge reproduces the S4 baseline bit-for-bit (load provably carried).
Therefore `nonvertical_s5_contact.py`'s plate-reaction `F_n` extraction (§2 "Contact force
extraction") is **unreliable as written** and must move to a contact-pressure output or a
prescribed-displacement reaction path before any 2 m `F_n` is quoted.

### 5.4 Files changed by this correction

| Action | File |
|---|---|
| changed | `src\climbing\coupling\plantar_bc.py` (`quad4` + outward orientation + `swap_pair` + `hold_spring_k`; `_apply_spring` anchor-id fix) |
| changed | `scripts\opensim_fe\thums_feb.py` (`_build_thums_plantar_quads`; G7 recipe; static ramp; `model.control_ = None`; contact opt-ins) |
| added | `scripts\opensim_fe\proof_s5_facet_fix.py` |
| added | `tests\test_s5_contact_facet_fix.py` (15 tests) |
| added | `docs\S5_contact_facet_fix.md` |

**Defaults unchanged**: `plantar_bc="fixed"`, `pad_domain_type=""`, `pad_hg=None`,
`DEFAULT_PENALTY=1.0`. Every new kwarg is opt-in. Hard gates: **93 / 19+1 / 11 / 11 / 13 +
15 new**.


---

## 9. 波5 补注（2026-10-07，APPEND-ONLY）—— `docs/S5_wave5_matrix.md`

> 本节由波5追加，不改动上文任何结论。上文「contact 行全部发散」是**波4 默认 `<Control>` +
> `use_rigid=False`** 下的实测，仍然成立。

波5 把 commit `315c373` 的修复（quad4 + G7 配方 + opt-in hold 弹簧）接进真实矩阵，并发现
一个上文未识别的**剩余阻塞**：

1. 真实矩阵（`calcaneus_r_anatframe.npz`，跖面 **137 quad4 / 0 tri3**）在
   **`use_rigid=True`（幽灵刚体加载）** + **G7 配方** + **hold k=1 N/mm** 下，
   **7 条先前失败 contact 行全部 CONVERGED+CARRYING**（rc=0、end_t=1.0）。
2. **hold 弹簧必需**：同一 case 去 hold ⇒ 2–3 s 内发散（接触闭合前 3 刚体模态）。
3. **`use_rigid=False` 的 PressureLoad 路径仍发散**（接触闭合 t≈0.42 处 FEBio 收敛判据退化
   → `Max nr of reformations reached`）。波4 矩阵 contact 行正是该路径 ⇒ 这是矩阵与 proof
   （`use_rigid=True`）的**真实差异之一**（另一差异是 mesh 帧：proof calcnframe / 矩阵 anatframe）。

产物：`results/opensim_fe/nonvertical_s5_contact_wave5.json` + `NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md`。
门禁：**93 / 19+1 / 11 / 11 / 13 / 15** 全绿（未动任何既有 `results/opensim_fe/*`）。
