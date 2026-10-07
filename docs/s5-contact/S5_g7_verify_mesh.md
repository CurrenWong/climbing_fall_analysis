# S5 — G7 (2 m landing) contact result: mesh-convergence verification

**Status:** COMPLETE (2026-10-06) — 3-point mesh sweep done (coarse/medium/fine), each block appended
as it finished; verdict in §3: **FAIL — the G7 `σ/σ_law` is not mesh-converged; only the resultant
`F_n ≈ 34.4 kN` is (< 0.6 %).**
**Owner task:** establish whether the currently-proposed G7 number (from a **single `tet4` pad
mesh**) is **mesh-converged** at sink = 170 mm, so it can be quoted.
**Units:** mm–N–MPa–s. FEBio binary `D:\Program\FEBioStudio\bin\febio4.exe` (4.13.0).
**Recipe source (do not re-derive):** `docs\S5_contact_udg.md` §7 (NEW recipe).
**Run dir:** `temp\opensim_fe\s5_t9\g7verify_b\` (new; collides with nothing).
**Harness:** `temp\opensim_fe\s5_t9\g7verify_b\g7verify.py` (one mesh per invocation; appends
to this doc + `ledger.txt` after **each** mesh).

---

## 1. Current G7 result — setup transcription (no compute)

### 1.1 Deck

Produced by `temp\opensim_fe\s5_t9\udg_contact\udg_contact_sweep2.py` (group `tet`), deck family
`udg_contact\n7_tet4_s170.feb` — byte-identical twins `n7b_tet4_s170.feb` (reproducibility run),
`n7c_tet4_s170_fb.feb`, `n7d_tet4_s170_fb.feb` (force-balance runs).

### 1.2 Mesh — exact args

```python
pad_block_tet(6, 6, 4, (300.0, 300.0, 200.0))          # S5 pad; default size_mm
indenter_block_tet(2, 2, 1, (200.0, 300.0, 5.0),
                   offset_mm=(50.0, 0.0, 200.5))       # fixed indenter; bottom z = 200.5
```

* gap `gap_mm = 0.5` (pad top z = 200.0 → indenter bottom z = 200.5).
* pad `<SolidDomain name="pad_tets" mat="pad"/>` — **default (tet4) domain**, *not* `udg-hex`.
* pad material: compressible **Ogden** (`k = 1`, `pressure_model = 2`), emitted at MPa scale.

### 1.3 Node / element counts (verified by counting the deck `n7_tet4_s170.feb`)

| body | mesh | nodes | elements | top facets |
|---|---|---|---|---|
| pad | `pad_block_tet(6,6,4)` | (6+1)(6+1)(4+1) = **245** | 6·6·6·4 = **864** tet4 | 2·6·6 = 72 tri3 |
| indenter | `indenter_block_tet(2,2,1)` | 3·3·2 = **18** | 6·2·2·1 = **24** tet4 | 2·2·2 = 8 tri3 (top + bottom) |
| **combined deck** | | **263** (grep `<node id=` = 263 ✓) | **888** (grep `<elem id=` = 888 ✓) | pad_top + ind_top + ind_bot tri3 = 88 ✓ |

### 1.4 Recipe (NEW — `docs\S5_contact_udg.md` §7, quoted verbatim, not re-derived)

Contact block:

```xml
<contact name="indenter_pad" surface_pair="indenter_pad_pair" type="sliding-elastic">
  <laugon>PENALTY</laugon>      <penalty>0.1</penalty>      <auto_penalty>0</auto_penalty>
  <two_pass>0</two_pass>        <node_reloc>0</node_reloc>
  <symmetric_stiffness>0</symmetric_stiffness>
  <tolerance>0.005</tolerance>  <search_radius>20</search_radius>
  <fric_coeff>0.6</fric_coeff>
</contact>
```

Control / stepper: `time_steps=2400`, `step_size = dtmax = 1/2400`, `cutback=0.125`,
`max_retries=20`, `opt_iter=15`, BFGS `max_ups=10`, `reform_each_time_step=1`, PARDISO,
`max_refs=25`. Load: `displacement`, `sinkage_mm=170`, `lateral_mm=(0,0)`. BC: pad bottom z
fixed; pad x=0 face x fixed; pad y=0 face y fixed; all indenter nodes x/y/z prescribed.

### 1.5 Current single-mesh numbers (read from `.log` tail + `.xplt` at last state)

`n7_tet4_s170` / `n7b_tet4_s170` (×2 identical) → `n7c_tet4_s170_fb` / `n7d_tet4_s170_fb`
(force balance): **`rc=0`, NORMAL TERMINATION, `end_t=1`, 2400/2400 steps, `negJac=0`, no NAN**,
wall-clock ≈ 40–46 s, `n_states=2401`.

| quantity | value |
|---|---|
| indentation | **148.7 mm** (ξ = **0.7437**) |
| penetration | **20.76 mm** |
| contact force `F_n` | **34 600 N** |
| force balance (indenter-BC vs pad-top reaction) | **7.80e-07** |
| σ_contact | **0.5767 MPa** |
| σ_1D-law(ξ) | **0.2616 MPa** |
| **σ/σ_law** | **2.205** |

> This is the number under test. The question: does it survive mesh refinement at sink = 170 mm?

---

## 2. Mesh sweep at G7 sink = 170 mm (NEW recipe)

**Design.** Indenter **fixed** across the whole sweep (`indenter_block_tet(2,2,1,(200,300,5),
offset_mm=(50,0,200.5))`, gap 0.5 mm) — **only the pad mesh changes**, so exactly one variable
moves. Chosen densities (coarse → medium → fine, element size 100 → 50 → 25 mm):

| # | `pad_block_tet(nx,ny,nz, size_mm)` | pad nodes | pad tets | note |
|---|---|---|---|---|
| 1 | `(3, 3, 2, (300,300,200))` | 48 | 108 | very coarse (100 mm cells) |
| 2 | `(6, 6, 4, (300,300,200))` | 245 | 864 | **= the current single-mesh G7** (anchor) |
| 3 | `(12, 12, 8, (300,300,200))` | 1521 | 6912 | very fine (25 mm cells) |

Each block below is appended by the harness **immediately after that mesh finishes** (durability:
nothing is held only in memory). A mesh that fails to converge is recorded honestly and is **not**
treated as data.

<!-- harness appends mesh blocks below -->
### Mesh `g7v_b_3x3x2` — `pad_block_tet(3,3,2,(300,300,200))`

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **23 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2416 |
| first failed t | 0.0036998 |
| negJac | 11 |
| NAN | no |
| n_states | 2417 |
| indentation | 153.965 mm (ξ=0.769825) |
| penetration | 15.535 mm |
| F_n | 34404.2 N |
| force balance (rel) | 3.95635e-09 |
| σ_contact | 0.573404 MPa |
| σ_1D-law | 0.291311 MPa |
| **σ/σ_law** | **1.96836** |

### Mesh `g7v_b_6x6x4` — `pad_block_tet(6,6,4,(300,300,200))`

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **50 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation | 148.736 mm (ξ=0.743679) |
| penetration | 20.7642 mm |
| F_n | 34601.7 N |
| force balance (rel) | 7.7964e-07 |
| σ_contact | 0.576695 MPa |
| σ_1D-law | 0.261581 MPa |
| **σ/σ_law** | **2.20465** |

### Mesh `g7v_b_12x12x8` — `pad_block_tet(12,12,8,(300,300,200))`

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **272 s**.

| metric | value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2401 |
| first failed t | 0.972897 |
| negJac | 51 |
| NAN | no |
| n_states | 2402 |
| indentation | 141.156 mm (ξ=0.705778) |
| penetration | 28.3445 mm |
| F_n | 34395.2 N |
| force balance (rel) | 2.58158e-07 |
| σ_contact | 0.573254 MPa |
| σ_1D-law | 0.22297 MPa |
| **σ/σ_law** | **2.57099** |

### 2.1 Consolidated sweep table

| # | pad `(nx,ny,nz)` | pad nodes/tets | wall (s) | `end_t` | nconv | negJac | NAN | indentation (mm) | ξ | penetration (mm) | `F_n` (N) | force-balance | σ_contact (MPa) | σ_law (MPa) | **σ/σ_law** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | (3,3,2) coarse | 48 / 108 | 23 | 1 | 2416 | 11 | no | 153.965 | 0.76983 | 15.535 | 34404.2 | 3.96e-09 | 0.573404 | 0.291311 | **1.96836** |
| 2 | (6,6,4) medium | 245 / 864 | 50 | 1 | 2400 | 0 | no | 148.736 | 0.74368 | 20.764 | 34601.7 | 7.80e-07 | 0.576695 | 0.261581 | **2.20465** |
| 3 | (12,12,8) fine | 1521 / 6912 | 272 | 1 | 2401 | 51 | no | 141.156 | 0.70578 | 28.344 | 34395.2 | 2.58e-07 | 0.573254 | 0.222970 | **2.57099** |

All three meshes terminated `rc=0 / NORMAL / end_t=1` with no NAN. The medium `(6,6,4)` run is the
**exact reproduction** of the current single-mesh G7 result (`ind 148.736` vs `148.7`; `F_n 34601.7`
vs `34 600`; `σ/σ_law 2.20465` vs `2.205`; force-balance `7.80e-07` identical) — so the harness
agrees with the prior result to 4–5 significant figures. **Note the negative-Jacobian events:**
coarse `negJac=11`, medium `negJac=0`, fine `negJac=51` — element distortion is *not* eliminated by
refinement.

---

## 3. Convergence judgement

**Drift between successive refinements** (relative to the coarser of each pair):

| quantity | coarse→medium | medium→fine | net coarse→fine | trend |
|---|---|---|---|---|
| `F_n` | **+0.574 %** | **−0.597 %** | −0.026 % | ✅ converged (<0.6 %) |
| σ_contact | +0.574 % | −0.597 % | −0.026 % | ✅ converged (<0.6 %) |
| force balance (rel) | 4e-9 → 8e-7 | 8e-7 → 3e-7 | all ≈ machine zero | ✅ (identity check holds) |
| indentation | **−3.40 %** | **−5.10 %** | **−8.32 %** | ❌ drifting, step growing |
| penetration | **+33.66 %** | **+36.51 %** | **+82.46 %** | ❌ drifting, step growing |
| σ_1D-law (via ξ) | −10.21 % | −14.76 % | −23.46 % | ❌ drifting |
| **σ/σ_law** | **+12.00 %** | **+16.62 %** | **+30.62 %** | ❌ **drifting, step growing** |

### VERDICT — **FAIL** (the single-mesh G7 number is NOT mesh-converged)

* **`σ/σ_law` is not converged.** It moves **+12.00 %** from coarse→medium and **+16.62 %** from
  medium→fine — i.e. the successive step does **not shrink** (it grows), so the sequence is *not*
  approaching a limit in this range. Net coarse→fine is **+30.62 %** (1.968 → 2.205 → 2.571).
  A quantity whose refinement step grows from 12 % to 17 % is not "within a few percent" of its
  continuum value.
* **The deformation field is not converged.** Indentation drifts **−3.40 % → −5.10 %**
  (net −8.32 %); penetration — its complement at fixed 170 mm travel and 0.5 mm gap — drifts
  **+33.7 % → +36.5 %** (net +82.5 %). Again the step grows with refinement.
* **Only the resultant force is mesh-insensitive.** `F_n ≈ 34.4 kN` and `σ_contact = F_n / 60000 mm²`
  vary by **< 0.6 %** across a 64× element-count change. The quoted *force* is robust; the quoted
  *ratio* is not.
* **Refinement does not stabilise the element response.** `negJac` goes `11 → 0 → 51`: the finest
  mesh already has 51 negative-Jacobian events (recovered by retries — it still reached `end_t=1`,
  but the element-distortion regime is active).

**Consequence for the paper.** The G7 `σ/σ_law = 2.205` from the single `(6,6,4)` `tet4` mesh
**must not be quoted as a mesh-converged result**. The final value keeps climbing with refinement
(net +30.6 % and still rising), so the single-mesh number understates it by an unbounded amount at
this resolution. **What *can* be quoted at this load level is the resultant contact force
`F_n ≈ 34.4 kN` (mesh-converged to < 0.6 %)**, not the flat-punch enhancement ratio.

**Why (mechanism, not over-claim).** The flat-punch edge is a stress singularity, so the averaged
pad-top indentation converges slowly (no fixed order), while the global equilibrium resultant does
not depend on resolving the singularity — hence `F_n` is stable but `indentation`/ξ (and therefore
`σ/σ_law`) are not. This is consistent with the sibling observation
(`docs\S5_contact_convergence.md` §3.3) that tet4 + this Ogden foam sits outside its robust domain
at `ξ ≳ 0.5`; here the coarse/fine meshes both trigger negative-Jacobian events.


