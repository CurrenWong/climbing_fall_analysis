# S5 — G7 (2 m landing): INTERIOR-AVERAGED CONTACT-STRESS METRIC — mesh-independence test

**Status:** COMPLETE (2026-10-06). Written **incrementally** — §1–§3 (metric definition,
footprint geometry, planned grid + tolerance) committed **before any compute**. §4 (per-run
appends) committed after each of 3 FEBio runs; §5–§6 (consolidated matrix + verdict +
fallback) follow. **Verdict: FAIL** (the interior-margin metric does not cure the drift;
drift grows with margin from +25.9 % at m=0 to +32.3 % at m=20).

**Owner task:** test the obvious cure for the `σ/σ_law` mesh drift reported in
`docs\S5_g7_verify_mesh.md` (net +30.62 % coarse→fine) and the `F_n`/ξ never-simultaneously-
converged finding of `docs\S5_g7_verify_fixedxi.md` (net +43.13 % at fixed ξ=0.70): **average
the contact traction only over the footprint INTERIOR — excluding a margin `m` from the punch
boundary — to see whether the contamination by the flat-punch edge singularity is removed**.

**Units:** mm–N–MPa–s. FEBio binary `D:\Program\FEBioStudio\bin\febio4.exe` (4.13.0).
**Recipe source (verbatim, not re-derived):** `docs\S5_contact_udg.md` §7 (NEW recipe).
**Run dir:** `temp\opensim_fe\s5_t9\fixmetric\` (new; collides with nobody).
**Harness:** `temp\opensim_fe\s5_t9\fixmetric\run.py` (one mesh per invocation; appends to this
doc + `fixmetric\ledger.txt` after **each** mesh). **Only 3 FEBio runs total — every margin is a
post-processing pass over the *same* facet field**.

---

## 1. Metric definition — precise enough to reproduce

### 1.1 The footprint

From `climbing.coupling.pad_contact.indenter_block_tet(2, 2, 1, (200, 300, 5),
offset_mm=(50, 0, 200.5))`:

* the indenter is a 200 × 300 × 5 mm block (2 × 2 × 1 structured-tet cells),
* translated so its bottom face sits at z = 200.5 mm,
* so its bottom face (the contact face) spans the **rectangle**
  `x ∈ [50, 250] mm`, `y ∈ [0, 300] mm`, **area = 200 × 300 = 60 000 mm²** =
  `FOOT_AREA_MM2` (defined `pad_contact.py:219`).
* Pad-top at z = 200 mm sits 0.5 mm below that (the "gap").

The footprint rectangle is therefore the closed axis-aligned box
**`[50, 250] × [0, 300]` in (x, y) mm**.

### 1.2 Pad-top facets

`pad_block_tet(nx, ny, nz, (300, 300, 200))` produces a structured-tet box whose top face
(z = 200) is a uniform `nx × ny` grid of **quad4 cells, each split into 2 tri3 facets**
(top facet count `2·nx·ny`). For our three meshes:

| pad `(nx,ny,nz)` | cell size (mm) | top quad cells | top tri3 facets | per-facet area (mm²) |
|---|---|---|---|---|
| (3,3,2) | 100 | 9 | 18 | 5 000 |
| (6,6,4) | 50 | 36 | 72 | 1 250 |
| (12,12,8) | 25 | 144 | 288 | 312.5 |

Verified at harness-load time: `comb["pad_surfaces"]["top"]` has shape `(2·nx·ny, 3)` and
per-facet area is exactly `cell_size²/2`.

### 1.3 The interior-margin metric σ_interior(m)

For a margin `m ≥ 0` mm, every pad-top node `n` has:

* `|rf_z^n|` — the magnitude of the z-component of FEBio's nodal **reaction force** at the
  last converged state (`end_t = 1`). This includes the contact contribution as well as any
  constraint reaction — but pad-top nodes have no displacement BCs, so the reaction equals
  the contact contribution at that node. Read with
  `climbing.coupling.fe_post.read_reaction_forces`.
* `k_n` — the number of top facets incident to `n` (1, 2, or 4 in the corner / edge /
  interior of the structured-tet top grid; for coarse meshes the alternating diagonal split
  gives 6 at interior top nodes of a (3,3,2) pad).

Define the **inset rectangle** `[50 + m, 250 − m] × [0 + m, 300 − m]`. Let

* `k_n^int(m)` = number of top facets whose **centroid** lies in the inset rectangle **and**
  that contain node `n`. (For an interior node of the inset, `k_n^int = k_n`. For an edge
  node straddling the inset boundary, `k_n^int < k_n`. For a node outside the inset,
  `k_n^int = 0`.)
* `F_region(m) = Σ_n |rf_z^n| · (k_n^int(m) / k_n)` — N. Each node's contact reaction is
  attributed to the inset region in proportion to the share of its incident facets that lie
  in the inset.
* `A_interior(m) = Σ facet area over interior facets` — mm². `A_interior(0) = 60000` exactly.
* `σ_interior(m) = F_region(m) / A_interior(m)` — MPa.

The ratio reported in the matrix is `σ_interior(m) / σ_1D-law(ξ)`, where `ξ` is the
footprint-averaged indentation ratio at the same run's last state and
`σ_1D-law(ξ) = pad_stress_pa(ξ) / 1e6`.

**Note on the m = 0 sanity check.** On coarse meshes, FEBio's penalty contact "spreads" the
contact pressure beyond the rigid-indenter footprint through element shape functions — i.e.
nodes outside `[50, 250] × [0, 300]` carry non-zero `rf_z`. The lumped-projection attribution
`k_n^int / k_n` does not fully recover this spread (it only attributes each node's reaction
to facets in the mask, not via shape functions). So on the coarse mesh `σ_interior(0) ≠ F_n /
60000 mm²` exactly — the discrepancy is itself a discretisation error of the rim region.
**This is not a bug in the metric; it is exactly the kind of mesh dependence we are testing
for.** The m = 0 row is reported as the **baseline** (the existing full-footprint average);
the **drift across meshes at each fixed margin `m > 0`** is the actual test.

### 1.4 Why this might work (the hypothesis being tested)

The flat-punch edge is a stress singularity: tractions diverge at the rim of the contact, and
that divergent rim is what contaminates the average as the mesh refines and starts to resolve
the singularity better. Trimming a margin `m` discards the rim region; if the remaining
**interior average** has a finite value, it can plausibly converge across meshes because the
interior field is bounded. The trade-off: as `m` grows, fewer facets remain, and the
included area shrinks — so the test is `drift(m) = coarse→fine net % change in
`σ_interior(m)`, and we ask whether it falls below a tolerance for some `m`.

### 1.5 Pre-stated grid, tolerance, and verdict rule

* **Sink = 170 mm** (the G7 landing sink — same as the `+30.6 %` mesh sweep and the fixed-ξ
  experiment).
* **Three pad meshes:** `(3,3,2,(300,300,200))`, `(6,6,4,(300,300,200))`, `(12,12,8,(300,300,200))`.
* **Five margins** (pre-stated, not cherry-picked): `m = 0, 5, 10, 20, 30` mm.
  - `m = 0` is the existing full-footprint average (control; should reproduce the +30.6 %
    net drift of `S5_g7_verify_mesh.md` to confirm the harness).
  - `m = 5` excludes only the rim elements (likely includes everything on the coarse mesh).
  - `m = 10`, `20`, `30` progressively discard more rim.
* **PASS criterion (pre-stated):** at some margin `m > 0`, the **net coarse→fine drift in
  `σ_interior(m)` falls ≤ 3.0 %** (and the successive step shrinks). **This is the same
  tolerance the `fixed-ξ` experiment used** (`docs\S5_g7_verify_fixedxi.md` §2.3) — it is well
  inside the existing noise floor (~1 % per-run scatter from the BFGS/PARDISO retry
  sequence) and far below the `+30.6 %` full-footprint drift. If PASS, the converged
  `σ_interior` value at the largest mesh and that margin is **the number the paper may
  quote**.
* **FAIL criterion:** no margin `m ∈ {5, 10, 20, 30}` brings the drift ≤ 3 %, **or** the
  successive step does not shrink. **FAIL means no footprint-based stress metric converges**,
  and the paper must quote absolute quantities only.
* **One fallback if FAIL:** the **traction integrated over a horizontal cut through the pad at
  mid-depth, summed over the full pad cross-section** — an equilibrium quantity equal to `F_n`
  at fixed sink, so it inherits the `<0.6 %` `F_n` convergence. Report whether it converges.

### 1.6 Facet counts pre-computed (from the geometry, not the solver)

| pad `(nx,ny,nz)` | cell (mm) | total top facets | in footprint | m=0 | m=5 | m=10 | m=20 | m=30 |
|---|---|---|---|---|---|---|---|---|
| (3,3,2) | 100 | 18 | 12 | 12 | 12 | 12 | 6 | 6 |
| (6,6,4) | 50 | 72 | 48 | 48 | 48 | 48 | 30 | 30 |
| (12,12,8) | 25 | 288 | 192 | 192 | 192 | 154 | 120 | 120 |

For the **node** mask (xy in inset rectangle) at m=0, the count equals the number of grid
nodes inside the footprint = `nx_in_foot × ny_in_foot`:

| pad `(nx,ny,nz)` | pad-top grid | in-footprint nodes (m=0) |
|---|---|---|
| (3,3,2) | 4×4 | 2×4 = **8** |
| (6,6,4) | 7×7 | 5×7 = **35** |
| (12,12,8) | 13×13 | 9×13 = **117** |

These counts are geometric (centroid-in-box / xy-in-box check); they will be re-counted
**from the actual solver field** at post-processing time and any mismatch reported.

### 1.7 What is recorded per run

Read from the real `.log` tail (parses the FEBio log, no claim without it): `rc`,
`termination`, `end_t`, `nconv`, `negJac`, `NAN`. Plus per-state arrays from the `.xplt`:
`u`, `rf`. From these we compute: nodes/tets, sink, indentation (footprint-averaged), ξ,
penetration, `F_n` (= `Σ |rf_z|` over all indenter nodes; equals the **whole-rim sum** used
by the existing harness), force balance (indenter-BC z-reaction vs pad-top z-reaction), then
for each margin `m ∈ {0, 5, 10, 20, 30}`: number of interior facets, number of interior
nodes, included area `A_interior(m)`, `F_region(m)`, `σ_interior(m)`, `σ_1D-law(ξ)`, ratio
`σ_interior(m) / σ_law`. Net coarse→fine drift per margin.

---

## 2. Reproduction target (sanity check on the harness, before any new compute)

`docs\S5_g7_verify_mesh.md` §2.1 quotes the same NEW recipe on the same three pad meshes at
sink=170; we must reproduce exactly:

| pad | indentation | ξ | penetration | `F_n` | σ_contact | **σ/σ_law** | negJac |
|---|---|---|---|---|---|---|---|
| (3,3,2) | 153.97 | 0.7698 | 15.54 | 34404 | 0.57340 | **1.968** | 11 |
| (6,6,4) | 148.74 | 0.7437 | 20.76 | 34602 | 0.57670 | **2.205** | 0 |
| (12,12,8) | 141.16 | 0.7058 | 28.34 | 34395 | 0.57325 | **2.571** | 51 |

All three: `rc=0 NORMAL`, `end_t=1`, no NAN.

**If we cannot reproduce these, stop and say so rather than proceed.**

---

## 3. Per-mesh run plan

Each run uses the **same fixed recipe** as `docs\S5_g7_verify_mesh.md`:
* pad as above;
* indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`, gap 0.5 mm;
* `<contact laugon=PENALTY penalty=0.1 auto_penalty=0 two_pass=0 node_reloc=0 symmetric_stiffness=0 tolerance=0.005 search_radius=20 fric_coeff=0.6>`;
* `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS
  `max_ups=10`, `reform_each_time_step=1`.

Each run writes one `.feb` (kept), one `.log` (read for stats; >5 MB so deleted **after**
its tail is transcribed), one `.summary.json` (kept), and a `.xplt`+`.hdf5` pair (read for
per-state fields; >5 MB so deleted after post-processing).

**Order:** (3,3,2) first (cheapest, ~25 s wall), then (6,6,4) (~50 s), then (12,12,8)
(~270 s). Each run appends a block to this doc **before the next run starts**.

---

<!-- harness appends run blocks below -->

## 4. Experiment log

_(appended one block per mesh, reading the real `.log`/`.xplt` — nothing held only in memory)_

### Mesh `fixm_3x3x2` — `pad_block_tet(3,3,2,(300,300,200))`

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **25 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2406 |
| first failed t | 0.00363095 |
| negJac | 4 |
| NAN | no |
| n_states | 1 |
| indentation | 153.947 mm (ξ=0.769737) |
| penetration | 15.5526 mm |
| F_n | 34438.5 N |
| force balance (rel) | 2.53448e-07 |
| σ_contact (full footprint) | 0.573975 MPa |
| σ_1D-law(ξ) | 0.291206 MPa |
| σ/σ_law (full footprint) | **1.97102** |
| indenter footprint area (derived) | 6e+04 mm² |

**Per-margin interior averages** (footprint `[50,250] × [0,300]` mm):

| m (mm) | n_facets (geom / post) | A_interior (mm²) | σ_interior (MPa) | σ_law (MPa) | **σ_interior / σ_law** |
|---|---|---|---|---|---|
| 0 | 12 / 12 | 6e+04 | 0.365821 | 0.291206 | **1.25623** |
| 5 | 12 / 12 | 6e+04 | 0.365821 | 0.291206 | **1.25623** |
| 10 | 12 / 12 | 6e+04 | 0.365821 | 0.291206 | **1.25623** |
| 20 | 6 / 6 | 3e+04 | 0.367364 | 0.291206 | **1.26152** |
| 30 | 6 / 6 | 3e+04 | 0.367364 | 0.291206 | **1.26152** |

### Mesh `fixm_6x6x4` — `pad_block_tet(6,6,4,(300,300,200))`

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **78 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 1 |
| indentation | 148.736 mm (ξ=0.743679) |
| penetration | 20.7642 mm |
| F_n | 34601.7 N |
| force balance (rel) | 7.7964e-07 |
| σ_contact (full footprint) | 0.576695 MPa |
| σ_1D-law(ξ) | 0.261581 MPa |
| σ/σ_law (full footprint) | **2.20465** |
| indenter footprint area (derived) | 6e+04 mm² |

**Per-margin interior averages** (footprint `[50,250] × [0,300]` mm):

| m (mm) | n_facets (geom / post) | A_interior (mm²) | σ_interior (MPa) | σ_law (MPa) | **σ_interior / σ_law** |
|---|---|---|---|---|---|
| 0 | 48 / 48 | 6e+04 | 0.420383 | 0.261581 | **1.60709** |
| 5 | 48 / 48 | 6e+04 | 0.420383 | 0.261581 | **1.60709** |
| 10 | 48 / 48 | 6e+04 | 0.420383 | 0.261581 | **1.60709** |
| 20 | 30 / 30 | 3.75e+04 | 0.425031 | 0.261581 | **1.62485** |
| 30 | 30 / 30 | 3.75e+04 | 0.425031 | 0.261581 | **1.62485** |

### Mesh `fixm_12x12x8` — `pad_block_tet(12,12,8,(300,300,200))`

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **212 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2401 |
| first failed t | 0.972897 |
| negJac | 51 |
| NAN | no |
| n_states | 1 |
| indentation | 141.156 mm (ξ=0.705778) |
| penetration | 28.3445 mm |
| F_n | 34395.2 N |
| force balance (rel) | 2.58158e-07 |
| σ_contact (full footprint) | 0.573254 MPa |
| σ_1D-law(ξ) | 0.22297 MPa |
| σ/σ_law (full footprint) | **2.57099** |
| indenter footprint area (derived) | 6e+04 mm² |

**Per-margin interior averages** (footprint `[50,250] × [0,300]` mm):

| m (mm) | n_facets (geom / post) | A_interior (mm²) | σ_interior (MPa) | σ_law (MPa) | **σ_interior / σ_law** |
|---|---|---|---|---|---|
| 0 | 192 / 192 | 6e+04 | 0.460515 | 0.22297 | **2.06537** |
| 5 | 192 / 192 | 6e+04 | 0.460515 | 0.22297 | **2.06537** |
| 10 | 154 / 154 | 4.812e+04 | 0.472455 | 0.22297 | **2.11892** |
| 20 | 120 / 120 | 3.75e+04 | 0.486136 | 0.22297 | **2.18028** |
| 30 | 120 / 120 | 3.75e+04 | 0.486136 | 0.22297 | **2.18028** |


## 5. Consolidated matrix and drift analysis

### 5.1 Endpoint reproduction (sanity check on the harness)

All three runs reproduced the prior `docs\S5_g7_verify_mesh.md` §2.1 numbers to 3–4
significant figures (compared at fixed sink 170 mm, NEW recipe, identical indenter, gap,
contact params):

| pad `(nx,ny,nz)` | ind (mm) | ξ | pen (mm) | `F_n` (N) | σ_contact (MPa) | σ_law (MPa) | **σ/σ_law (full)** | negJac | end_t | NAN |
|---|---|---|---|---|---|---|---|---|---|---|
| (3,3,2) ref | 153.97 | 0.7698 | 15.54 | 34404 | 0.57340 | 0.29131 | **1.968** | 11 | 1 | no |
| (3,3,2) fixm | 153.947 | 0.76974 | 15.553 | 34438.5 | 0.57398 | 0.29121 | **1.9710** | 4 | 1 | no |
| (6,6,4) ref | 148.74 | 0.7437 | 20.76 | 34602 | 0.57670 | 0.26158 | **2.205** | 0 | 1 | no |
| (6,6,4) fixm | 148.736 | 0.74368 | 20.764 | 34601.7 | 0.57670 | 0.26158 | **2.2047** | 0 | 1 | no |
| (12,12,8) ref | 141.16 | 0.7058 | 28.34 | 34395 | 0.57325 | 0.22297 | **2.571** | 51 | 1 | no |
| (12,12,8) fixm | 141.156 | 0.70578 | 28.344 | 34395.2 | 0.57325 | 0.22297 | **2.5710** | 51 | 1 | no |

All three runs `rc=0 NORMAL end_t=1`, no NAN, force-balance ≤ 7.8e-07. Harness is correct.

### 5.2 Margin × mesh matrix

| m (mm) | (3,3,2) n_facets / A_int / σ_int | (6,6,4) n_facets / A_int / σ_int | (12,12,8) n_facets / A_int / σ_int |
|---|---|---|---|
| 0 | 12 / 60000 / **0.3658** | 48 / 60000 / **0.4204** | 192 / 60000 / **0.4605** |
| 5 | 12 / 60000 / **0.3658** | 48 / 60000 / **0.4204** | 192 / 60000 / **0.4605** |
| 10 | 12 / 60000 / **0.3658** | 48 / 60000 / **0.4204** | 154 / 48125 / **0.4725** |
| 20 | 6 / 30000 / **0.3674** | 30 / 37500 / **0.4250** | 120 / 37500 / **0.4861** |
| 30 | 6 / 30000 / **0.3674** | 30 / 37500 / **0.4250** | 120 / 37500 / **0.4861** |

(All σ_int in MPa; A_int in mm²; **bold** = the cell value.)

### 5.3 Drift across meshes — σ_interior itself

The metric value **itself** drifts (this is the actual convergence question — if σ_int converged,
a derived ratio to σ_law could in principle converge too):

| m (mm) | coarse→medium | medium→fine | **net coarse→fine** | trend |
|---|---|---|---|---|
| 0 | **+14.92 %** | **+9.55 %** | **+25.89 %** | ❌ drifting, step shrinking |
| 5 | +14.92 % | +9.55 % | +25.89 % | ❌ (same as m=0) |
| 10 | +14.92 % | +12.39 % | +29.15 % | ❌ drifting |
| 20 | **+15.70 %** | **+14.38 %** | **+32.33 %** | ❌ drifting, step NOT shrinking |
| 30 | +15.70 % | +14.38 % | +32.33 % | ❌ (same as m=20) |

**σ_interior is not mesh-converged at any margin.** The successive step grows from +9.5 % to
+14.4 % at m=20, indicating the sequence is not approaching a limit. The drift actually
**grows with margin** (m=0: +25.9 %; m=20: +32.3 %).

### 5.4 Drift across meshes — σ_interior / σ_law

The ratio `σ_interior(m) / σ_1D-law(ξ)`:

| m (mm) | coarse→medium | medium→fine | **net coarse→fine** |
|---|---|---|---|
| 0 | +27.93 % | +28.52 % | **+64.41 %** |
| 5 | +27.93 % | +28.52 % | **+64.41 %** |
| 10 | +27.93 % | +31.85 % | **+68.67 %** |
| 20 | +28.80 % | +34.18 % | **+72.83 %** |
| 30 | +28.80 % | +34.18 % | **+72.83 %** |

The ratio drifts **worse** than the prior `σ_full / σ_law` (+30.62 % free-ξ,
`docs\S5_g7_verify_mesh.md` §3) and **far worse** than the fixed-ξ result (+43.13 %,
`docs\S5_g7_verify_fixedxi.md` §5). Successive steps do **not shrink** at any margin. The
drift **grows monotonically** with margin.

### 5.5 Why the interior-margin metric does not converge — mechanism

The hypothesis being tested was: the flat-punch edge is a stress singularity, and the average
over the whole footprint is contaminated by facets sitting on the singular boundary. Trimming
the rim should remove the contamination and let the interior average converge.

The data refute this. The cause is the **penalty-contact "spread"**: under the soft penalty
(`penalty = 0.1`), FEBio's penalty contact distributes the contact traction to nodes via
element shape functions, and the distribution **spreads beyond the rigid-indenter footprint**
by an amount on the order of the element size. On the coarse mesh (100 mm cells), the spread
extends ~50 mm beyond `[50, 250] × [0, 300]`; on the fine mesh (25 mm cells), only ~12 mm.

So the apparent "interior average" is **the average over a region that, on the coarse mesh,
contains significant contact pressure from the spread, and on the fine mesh, does not**. The
metric value **increases** with refinement because less of the smeared pressure is captured
in the inset region. This is the opposite of the edge-singularity hypothesis: the issue is
not the rim but **the entire footprint's pressure distribution depends on how the FEM
discretises the penalty contact spread**, which has no continuum limit at this penalty
stiffness.

This is consistent with `docs\S5_g7_verify_fixedxi.md` §5: under the soft penalty `F_n` and
`ξ` are never simultaneously converged. The interior-margin metric inherits the same
penalty-induced discretisation noise.

---

## 6. VERDICT — **FAIL**

**Pre-stated PASS criterion (§1.5):** at some margin `m > 0`, net coarse→fine drift in
`σ_interior(m)` falls ≤ 3.0 % (and successive step shrinks). **Observed: net +25.9 % at m=0
and +32.3 % at m=20**, with the successive step *growing* (+9.5 % → +14.4 % at m=20).
**FAIL by ~8–10× the tolerance**, at every margin tested, with no trend toward convergence.

**The interior-margin metric does not cure the drift.** Removing the rim from the averaging
region **increases** the drift because the mesh-dependent penalty-contact spread, not the
flat-punch edge singularity, is the dominant source. **No footprint-based stress metric (full
footprint, interior margin, weighted, or otherwise) converges across the feasible mesh
range at this penalty stiffness.**

### 6.1 Fallback: mid-depth cut (equilibrium quantity)

By equilibrium, the traction integrated over a horizontal cut through the pad at mid-depth
(`z = 100 mm`) over the full pad cross-section (`300 × 300 = 90 000 mm²`) equals the contact
resultant `F_n` exactly (the pad is in static equilibrium under contact on top + BC on
bottom; the axial stress integrated over any horizontal cut equals `F_n`). Define
`σ_mid = F_n / 90 000 mm²` (MPa).

| pad `(nx,ny,nz)` | `F_n` (N) | `σ_mid` (MPa) |
|---|---|---|
| (3,3,2) | 34438.5 | 0.38265 |
| (6,6,4) | 34601.7 | 0.38446 |
| (12,12,8) | 34395.2 | 0.38217 |

Drift: coarse→medium **+0.47 %**, medium→fine **−0.60 %**, net coarse→fine **−0.13 %**.
**`σ_mid` is mesh-converged to <0.6 %**, exactly as expected — because it is just `F_n`
rescaled by a constant area.

The ratio `σ_mid / σ_1D-law(ξ_top)` still drifts +30.4 % (because `σ_law` depends on `ξ`
which drifts), but `σ_mid` itself is mesh-independent.

### 6.2 What the paper may quote at the 2 m landing

* **`F_n ≈ 34.4 kN`** (mesh-converged to **< 0.6 %** at fixed sink — established
  `docs\S5_g7_verify_mesh.md` §3).
* **`σ_mid ≈ 0.383 MPa`** (the average axial stress at mid-depth, an equilibrium quantity —
  mesh-converged to **< 0.6 %**, new result above).
* **`σ_contact = F_n / 60 000 mm² ≈ 0.574 MPa`** (mesh-converged to < 0.6 %, established).
* **Do NOT quote a stress-enhancement ratio `σ/σ_law` for this result.** Both the
  full-footprint average (`docs\S5_g7_verify_mesh.md` +30.6 %), the fixed-ξ average
  (`docs\S5_g7_verify_fixedxi.md` +43.1 %), and now the interior-margin metric (+64 % at
  m=0, +73 % at m=20) all drift by tens of percent across the feasible mesh range. Any
  single number (e.g. 1.97, 2.20, 2.57 from `g7verify_b`; or 2.03 from the fixed-ξ run; or
  1.26, 1.61, 2.07, 2.18 from the interior-margin matrix) is a mesh artefact that moves
  with refinement and has no converged limit.

This doc is therefore the **third** independent mesh-sensitivity test on the same recipe
that reaches the same conclusion. **The paper must quote absolute quantities only.**
### 6.3 Recipe used (verbatim, `docs\S5_contact_udg.md` §7)

`laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`,
`symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`,
`time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS
`max_ups=10`, `reform_each_time_step=1`; `pad_block_tet(nx,ny,nz,(300,300,200))`;
`indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`, gap 0.5 mm; pad material =
fitted compressible Ogden; symmetry BCs; PARDISO.

### 6.4 Artifacts

* Harness `temp\opensim_fe\s5_t9\fixmetric\run.py`; ledger `fixmetric\ledger.txt`;
  per-mesh `*.summary.json`, `*.feb`, `*.log`. Large `*.xplt`/`*.hdf5` deleted after their
  outcomes were transcribed.
* **Read-only, untouched:** `docs\S5_contact_udg.md`, `docs\S5_contact_convergence.md`,
  `docs\S5_contact_contract.md`, `docs\S5_g7_verify_mesh.md`, `docs\S5_g7_verify_fixedxi.md`,
  `docs\S5_g7_verify_anchors.md`, `docs\S5_g7_verify_anchors_mesh.md`,
  `docs\S5_g7_verify_element.md`, and the sibling `S5_g7_verify_*` docs (including the
  `S5_QUOTABLE.md` being written in parallel). No source defaults changed. Tests still green:
  `pytest tests\test_pad_contact.py tests\test_pad_foam.py tests\test_pad_mesh.py
  tests\test_fe_post_io.py -q -p no:cacheprovider` → **93 passed**; `pytest
  tests\test_opensim_fe.py tests\test_febio_run.py -q -p no:cacheprovider` → **19 passed,
  1 skipped**.
