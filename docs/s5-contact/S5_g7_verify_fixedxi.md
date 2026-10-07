# S5 — G7 (2 m landing): FIXED-ξ mesh test — is `σ/σ_law` quotable at all?

**Status:** COMPLETE (2026-10-06). Written **incrementally** — §1/§2 (transcription +
pre-stated protocol) were committed **before any compute**; §3 was appended after **every single
run** (including the off-target coarse seed); §4/§5 consolidated last.
**Verdict (§5): FAIL — at FIXED ξ≈0.70 the ratio drifts +43.13 % net (worse than the +30.62 %
free-ξ sweep), so the ratio is intrinsically mesh-dependent; quote absolute quantities only.**
**Owner task:** run the **decisive experiment** that separates the two candidate explanations for
the `+30.6 %` ratio drift seen in `docs\S5_g7_verify_mesh.md`:

* **(H1) ξ artefact** — the drift is *entirely* because the soft penalty makes penetration
  mesh-dependent, so ξ moved (0.7698 → 0.7437 → 0.7058) and `σ_1D-law(ξ)` moved with it; if ξ is
  **held fixed**, the ratio stops drifting → the ratio **is** quotable.
* **(H2) intrinsic mesh dependence** — the ratio still drifts when ξ is held fixed → the ratio is
  genuinely mesh-dependent and the paper may quote **absolute quantities only**.

**Units:** mm–N–MPa–s. FEBio `D:\Program\FEBioStudio\bin\febio4.exe` (4.13.0).
**Recipe source (do not re-derive):** `docs\S5_contact_udg.md` §7 (NEW recipe).
**Run dir:** `temp\opensim_fe\s5_t9\fixb\` (new; collides with nobody).
**Harness:** `temp\opensim_fe\s5_t9\fixb\run.py` (one mesh+one sink per invocation; appends to
this doc **and** `fixb\ledger.txt` after **each** run).

---

## 1. Transcription of the existing evidence (NO compute)

### 1.1 The original mesh sweep at sink = 170 mm (NEW recipe), copied verbatim from `docs\S5_g7_verify_mesh.md`

Fixed across the sweep: indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`,
gap 0.5 mm, NEW recipe. Only the pad mesh changed. All runs `rc=0 NORMAL end_t=1`, no NAN.

| # | pad `(nx,ny,nz)` | pad nodes/tets | indentation (mm) | ξ | penetration (mm) | `F_n` (N) | σ_contact (MPa) | σ_law (MPa) | **σ/σ_law** | negJac |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | (3,3,2) coarse | 48 / 108 | 153.965 | 0.76983 | 15.535 | 34404.2 | 0.573404 | 0.291311 | **1.96836** | 11 |
| 2 | (6,6,4) medium | 245 / 864 | 148.736 | 0.74368 | 20.764 | 34601.7 | 0.576695 | 0.261581 | **2.20465** | 0 |
| 3 | (12,12,8) fine | 1521 / 6912 | 141.156 | 0.70578 | 28.344 | 34395.2 | 0.573254 | 0.222970 | **2.57099** | 51 |

**The drift observed there:** `σ/σ_law` moved **+12.00 %** (coarse→medium) then **+16.62 %**
(medium→fine) = net **+30.62 %**, *while the successive step grew*. **BUT** ξ also drifted
(0.76983 → 0.74368 → 0.70578) because the soft penalty (`penalty=0.1`) makes penetration
mesh-dependent (15.535 → 20.764 → 28.344 mm). Meanwhile `σ_contact = F_n/60000 mm²` was
**mesh-converged to <0.6 %** (0.5734 → 0.5767 → 0.5733). So **only the denominator
`σ_1D-law(ξ)` moved** — it is a strong function of ξ.

### 1.2 The hypothesis under test

> Because `σ_contact` is already convergence-grade (<0.6 %) and `σ_1D-law` depends **only on ξ**,
> the entire +30.6 % drift may be a **ξ artefact**, not a property of the contact solution.
> **Direct test:** *hold ξ fixed and let the sink vary per mesh.* Pick, per mesh, the sink that
> lands the **endpoint** at ξ = 0.70 ± 0.005, then compare `σ/σ_law` across meshes.

### 1.3 Cross-check evidence (from `docs\S5_g7_verify_element.md`, tet4 `(6,6,4)` only)

C's all-state read-back of the `(6,6,4)` tet4 sink-170 run gives the ξ→ratio curve:
`2.0573 @ ξ=0.7093`, `1.9190 @ ξ=0.6501`. So on the **medium** mesh, ratio rises steeply with ξ
in this band (≈+0.019 ratio per +0.01 ξ ≈ +1 %/0.01 ξ near ξ=0.70). If the coarse and fine meshes
are compared at *their own* ξ (0.770 and 0.706) the ratio difference is dominated by this
ξ-sensitivity — exactly H1. C did **not** put the meshes at a common ξ; this experiment does.

`σ_law(ξ)` local slope (from the same curve): `σ_law(0.7093)=0.2264`, `σ_law(0.6501)=0.1751`
→ `dσ_law/dξ ≈ 0.868 MPa/ξ`. Hence **matching ξ to ±0.005 leaves a residual σ_law (and ratio)
uncertainty of ≈ ±2.0 %** — this floor is stated up-front so the PASS tolerance is honest.

---

## 2. Protocol (pre-stated BEFORE any run)

### 2.1 Exact recipe (from `docs\S5_contact_udg.md` §7 — verbatim, not re-derived)

```xml
<contact name="indenter_pad" surface_pair="indenter_pad_pair" type="sliding-elastic">
  <laugon>PENALTY</laugon>      <penalty>0.1</penalty>      <auto_penalty>0</auto_penalty>
  <two_pass>0</two_pass>        <node_reloc>0</node_reloc>
  <symmetric_stiffness>0</symmetric_stiffness>
  <tolerance>0.005</tolerance>  <search_radius>20</search_radius>
  <fric_coeff>0.6</fric_coeff>
</contact>
```
Time stepper: `time_steps=2400`, `step_size = dtmax = 1/2400`, `cutback=0.125`,
`max_retries=20`; solver BFGS `max_ups=10`, `reform_each_time_step=1`, PARDISO, `max_refs=25`.
Load: `displacement`, `lateral_mm=(0,0)`; only `sinkage_mm` varies per run.

### 2.2 Held constant (identical to the sweep) — only the pad mesh + sink move

* pad sizes: `pad_block_tet(3,3,2,(300,300,200))`, `(6,6,4,(300,300,200))`, `(12,12,8,(300,300,200))`;
* indenter **fixed**: `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`, gap **0.5 mm**;
* pad material: fitted compressible Ogden (`k=1`, `pressure_model=2`) via `build_pad_contact_feb` default;
* symmetry BCs, `mu=0.6`, penalty 0.1, PARDISO.

### 2.3 Target and tolerance (pre-stated)

* **Target ξ = 0.70 ± 0.005** at the **ramp endpoint** (so `end_t=1` is a real converged state).
* **PASS criterion:** net `σ/σ_law` drift **coarse→fine ≤ 3.0 %** (and no growth in the successive
  step), i.e. well inside the ±2 % ξ-matching floor plus a little round-off, and far below the
  **+30.6 %** free-ξ drift. If the fixed-ξ drift is ≤3 %, the earlier +30 % is a ξ artefact.
* **FAIL criterion:** fixed-ξ net drift > 3 % (in particular if it stays ≈ +25–30 %) → the ratio
  is intrinsically mesh-dependent and the paper must quote absolute quantities only.

### 2.4 Mesh list + seeding guesses (pre-stated)

First-order guess `sink ≈ 170 × 0.70 / ξ_obs(at 170)`; one correction run from the run's own
ξ(travel) curve. **Budget ≤ 3 runs per mesh.**

| # | pad `(nx,ny,nz)` | ξ_obs @170 | first-order seed sink | note (from §1.3) |
|---|---|---|---|---|
| 1 | (3,3,2) | 0.76983 | **155** | coarse; expect ξ overshoot → lower sink |
| 2 | (6,6,4) | 0.74368 | **160** | C's curve cross-checks ≈159 |
| 3 | (12,12,8) | 0.70578 | **169** | fine; already near ξ=0.70 at 170 |

**Homing method (recorded per run):** each run writes **all** converged states to `.xplt`; the
harness reads `ξ(travel)` for every state. Because the ramp is linear and the response is
rate-independent, the **endpoint ξ at sink `S` equals ξ(travel=S) in any run**. So from the seed
run the harness interpolates `ξ(travel)=0.70` → predicted sink `T*`, which the correction run uses.
The monotonicity of ξ in travel is checked explicitly on the first iteration (pre-stated
assumption: ξ increases monotonically with sink at fixed mesh).

**Failure is data:** a mesh that cannot reach ξ=0.70 at a reasonable sink tops out; that is
recorded and **not** silently substituted with a different ξ.

**Per run recorded fields (all read from the real `.log` tail / last `.xplt` state):** `end_t`,
termination/NAN, `negJac`, nodes/tets, sink, indentation, ξ, penetration, `F_n`,
force balance (indenter-BC vs pad-top reaction), σ_contact, σ_1D-law, **σ/σ_law**.

---

## 3. Experiment log

_(appended one block per run, reading the real `.log`/`.xplt` — nothing held only in memory)_

<!-- harness appends run blocks below -->

### Run `fixb_coarse_s155` — pad `pad_block_tet(3,3,2,(300,300,200))`, sink = 155.0 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **21 s**.

| metric | endpoint value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2406 |
| first failed t | 0.00361111 |
| negJac | 58 |
| NAN | no |
| n_states | 2407 |
| indentation | 142.29 mm |
| ξ (actual) | **0.711458** |
| penetration | 12.208 mm |
| F_n | 24346.1 N |
| force balance (rel) | 4.02e-07 |
| σ_contact | 0.405769 MPa |
| σ_1D-law | 0.228433 MPa |
| **σ/σ_law** | **1.77632** |

**Homing read-back.** ξ is monotone non-decreasing in travel (0 decreasing steps of 2406).
Interpolated state at ξ=0.7: travel **152.15 mm** (indent 140, pen 11.654, F_n 22911.7, σ/σ_law 1.75546, fbal 9.162e-07).

> **Predicted sink for endpoint ξ=0.7: 152.15 mm** (endpoint ξ at sink S equals ξ(travel=S)).

### Run `fixb_coarse_s152p15` — pad `pad_block_tet(3,3,2,(300,300,200))`, sink = 152.15359 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **29 s**.

| metric | endpoint value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2406 |
| first failed t | 0.00387578 |
| negJac | 40 |
| NAN | no |
| n_states | 2407 |
| indentation | 140.01 mm |
| ξ (actual) | **0.700064** |
| penetration | 11.641 mm |
| F_n | 22877.3 N |
| force balance (rel) | 6.403e-08 |
| σ_contact | 0.381288 MPa |
| σ_1D-law | 0.217587 MPa |
| **σ/σ_law** | **1.75235** |

**Homing read-back.** ξ is monotone non-decreasing in travel (0 decreasing steps of 2406).
Interpolated state at ξ=0.7: travel **152.14 mm** (indent 140, pen 11.643, F_n 22883.4, σ/σ_law 1.7533, fbal 4.587e-07).

> **Predicted sink for endpoint ξ=0.7: 152.14 mm** (endpoint ξ at sink S equals ξ(travel=S)).

### Run `fixb_medium_s159` — pad `pad_block_tet(6,6,4,(300,300,200))`, sink = 159.0 mm

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **75 s**.

| metric | endpoint value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation | 140.34 mm |
| ξ (actual) | **0.70172** |
| penetration | 18.156 mm |
| F_n | 26646.4 N |
| force balance (rel) | 6.297e-07 |
| σ_contact | 0.444107 MPa |
| σ_1D-law | 0.219136 MPa |
| **σ/σ_law** | **2.02663** |

**Homing read-back.** ξ is monotone non-decreasing in travel (0 decreasing steps of 2400).
Interpolated state at ξ=0.7: travel **158.55 mm** (indent 140, pen 18.047, F_n 26381.6, σ/σ_law 2.02133, fbal 6.485e-07).

> **Predicted sink for endpoint ξ=0.7: 158.55 mm** (endpoint ξ at sink S equals ξ(travel=S)).

### Run `fixb_fine_s168` — pad `pad_block_tet(12,12,8,(300,300,200))`, sink = 168.0 mm

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **257 s**.

| metric | endpoint value |
|---|---|
| rc | 0 |
| termination | NORMAL |
| end_t | 1 |
| converged steps (nconv) | 2402 |
| first failed t | 0.983691 |
| negJac | 2 |
| NAN | no |
| n_states | 2403 |
| indentation | 139.75 mm |
| ξ (actual) | **0.698771** |
| penetration | 27.746 mm |
| F_n | 32563.7 N |
| force balance (rel) | 7.498e-09 |
| σ_contact | 0.542728 MPa |
| σ_1D-law | 0.216384 MPa |
| **σ/σ_law** | **2.50817** |

**Homing read-back.** ξ is monotone non-decreasing in travel (0 decreasing steps of 2402).
> max ξ at sink 168 is **0.698771** — i.e. 0.00123 **below** the target 0.70, which is **inside** the
> ±0.005 band, so this endpoint **is** the accepted ξ≈0.70 row (a correction run was unnecessary;
> a sink of ≈168.3 would hit 0.700 exactly). This is **not** a top-out — the distinction matters and
> is stated explicitly so the row is not mis-read as "fine could not reach ξ=0.70".

---

### 3.5 Homing steps actually used (these runs are data, not scaffolding)

| mesh | run | sink (mm) | endpoint ξ | role |
|---|---|---|---|---|
| (3,3,2) | `fixb_coarse_s155` | 155.000 | 0.71146 | seed — **overshoot** by +0.0115 → correction needed |
| (3,3,2) | `fixb_coarse_s152p15` | 152.154 | **0.70006** | **accepted** (Δξ = +0.00006) |
| (6,6,4) | `fixb_medium_s159` | 159.000 | **0.70172** | **accepted** (Δξ = +0.00172), no correction |
| (12,12,8) | `fixb_fine_s168` | 168.000 | **0.69877** | **accepted** (Δξ = −0.00123), no correction |

**Monotonicity assumption verified:** on every run ξ is monotone non-decreasing in travel (0 decreasing
steps), so the secant/interpolation prediction is well-behaved. **All four runs: `rc=0`, NORMAL,
`end_t=1`, no NAN.** Coarse seed 155 (ξ=0.71146, F_n 24346, ratio 1.7763, fbal 4.02e-07, negJac 58)
and its correction 152.154 (ξ=0.70006, F_n 22877, ratio 1.7523, fbal 6.40e-08, negJac 40) interpolate
to the same ξ=0.70 state (1.7555 vs 1.7533; F_n 22911.7 vs 22883.4) — the homing is self-consistent.

---

## 4. Consolidated result — `σ/σ_law` at FIXED ξ ≈ 0.70

### 4.1 The fixed-ξ grid

Indenter, gap (0.5 mm), material, BCs and the §7 recipe are **identical** to the sink-170 sweep; only
the pad mesh and the sink change, and the sink is chosen per mesh so the **endpoint** lands at
ξ = 0.70 ± 0.005.

| # | pad `(nx,ny,nz)` | pad nodes/tets | sink (mm) | `end_t` | negJac | NAN | indentation (mm) | ξ (actual) | penetration (mm) | `F_n` (N) | force-balance | σ_contact (MPa) | σ_law (MPa) | **σ/σ_law** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | (3,3,2) coarse | 48 / 108 | 152.154 | 1 | 40 | no | 140.013 | **0.70006** | 11.641 | 22 877.3 | 6.40e-08 | 0.38129 | 0.21759 | **1.75235** |
| 2 | (6,6,4) medium | 245 / 864 | 159.000 | 1 | 0 | no | 140.344 | **0.70172** | 18.156 | 26 646.4 | 6.30e-07 | 0.44411 | 0.21914 | **2.02663** |
| 3 | (12,12,8) fine | 1521 / 6912 | 168.000 | 1 | 2 | no | 139.754 | **0.69877** | 27.746 | 32 563.7 | 7.50e-09 | 0.54273 | 0.21638 | **2.50817** |

(σ_contact = `F_n` / 60 000 mm²; σ_law = `pad_stress_pa(ξ)/1e6`.)

### 4.2 Drift at fixed ξ vs the free-ξ sweep

| quantity | coarse→medium | medium→fine | **net coarse→fine** | free-ξ sweep (sink 170) net |
|---|---|---|---|---|
| **σ/σ_law** | **+15.65 %** | **+23.76 %** | **+43.13 %** | +30.62 % |
| `F_n` (= σ_contact) | +16.47 % | +22.21 % | **+42.34 %** | −0.03 % (✅ <0.6 %) |
| σ_1D-law(ξ) | +0.71 % | −1.26 % | −0.56 % | −23.46 % |
| ξ | +0.24 % | −0.42 % | −0.18 % | −8.32 % |
| penetration | +55.97 % | +52.82 % | +138.35 % | +82.46 % |

**Key contrast.** Holding ξ fixed **removed** the σ_law drift (from −23.5 % to −0.6 % ✅) but the
**numerator `F_n`/σ_contact — which was mesh-converged to <0.6 % at fixed sink — now drifts +42.3 %**.
The ratio therefore drifts **+43.1 %**, *more* than the +30.6 % free-ξ sweep, and the successive step
**grows** (+15.7 % → +23.8 %). The `σ/σ_law` sequence is not approaching a limit.

### 4.3 Why (mechanism — not over-claim)

* Under the soft penalty (`penalty=0.1`), **penetration is mesh-dependent** and does not shrink on
  refinement: at matched ξ it is 11.64 → 18.16 → 27.75 mm. Matching ξ pins the *average pad-top
  indentation* but **cannot** pin penetration, because the sink must rise (152.2 → 159.0 → 168.0 mm)
  exactly to absorb the extra penetration.
* The penalty contact force is carried by that penetration, so at fixed ξ `F_n` is **not** converged:
  22.88 → 26.65 → 32.56 kN. At fixed **sink** it *is* converged (<0.6 %), because global equilibrium
  fixes the resultant — but there the indent/pen split (hence ξ) drifts.
* So `F_n` and ξ are **never simultaneously converged**: fix the sink → `F_n` converges, ξ does not;
  fix ξ → ξ is pinned, `F_n` does not. Since `ratio = F_n / (60 000 · σ_law(ξ))`, neither control yields
  a converged ratio. Physically this is the flat-punch edge singularity: the resultant at a given
  *travel* is far-field-dominated (mesh-insensitive), but the force required for a given *indentation*
  is set by the near-edge traction, which no feasible mesh resolves.
* Consequence: the +30.6 % in `docs\S5_g7_verify_mesh.md` was **not** a pure ξ artefact. Making ξ equal
  across meshes does not collapse the ratio — the meshes carry different force at the same indentation.

---

## 5. Verdict — **FAIL** (the ratio is intrinsically mesh-dependent)

* **Pre-stated PASS criterion** (§2.3): net `σ/σ_law` drift coarse→fine **≤ 3.0 %**. **Observed: +43.13 %**
  — and the successive step *grows* (+15.65 % → +23.76 %). **FAIL by ~14× the tolerance.**
* **The decisive contrast:** free-ξ sweep drift = **+30.62 %**; fixed-ξ drift = **+43.13 %**. Holding
  ξ fixed did **not** remove the drift — it *increased* it, because the soft-penalty penetration
  transferred its mesh dependence into `F_n`. Therefore the earlier +30 % is **not** a ξ artefact, and
  the ratio is **not** quotable.
* **What the paper may quote at the 2 m landing — absolute quantities only:**
  * `F_n ≈ 34.4 kN` (mesh-converged to **< 0.6 %** across a 64× element-count change at the landing sink);
  * `σ_contact = F_n / 60 000 mm² ≈ 0.574 MPa` (same < 0.6 % convergence — 0.5734 / 0.5767 / 0.5733 MPa
    for (3,3,2)/(6,6,4)/(12,12,8)).
  * **Do not quote a stress-enhancement ratio `σ/σ_law` for this result.** Any single number (e.g. the
    single-`(6,6,4)` 2.205, or a fixed-ξ 2.03) is a mesh artefact that moves +15–24 % per refinement and
    has no converged limit in the feasible mesh range.
* **No single quotable `σ/σ_law`, with ξ and recipe, can be given** — that would require a PASS. This doc
  is therefore the record that the ratio is **retracted as a mesh-converged G7 deliverable**; only the
  absolute pair above survives.

**Recipe used (verbatim, `docs\S5_contact_udg.md` §7):** `laugon=PENALTY`, `penalty=0.1`,
`auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`,
`search_radius=20`, `fric_coeff=0.6`, `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`,
`max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`; `pad_block_tet(nx,ny,nz,(300,300,200))`;
`indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`, gap 0.5 mm; pad material = fitted
compressible Ogden; symmetry BCs; PARDISO.

### 5.1 Artifacts

* Harness `temp\opensim_fe\s5_t9\fixb\run.py`; ledger `fixb\ledger.txt`; per-run `*.summary.json`,
  `*.feb`, `*.log`. Large `*.xplt`/`*.hdf5` deleted after their outcomes were transcribed.
* **Read-only, untouched:** `docs\S5_contact_udg.md`, `docs\S5_contact_convergence.md`,
  `docs\S5_contact_contract.md`, `docs\S5_g7_verify_mesh.md`, `docs\S5_g7_verify_element.md` (and all
  sibling `S5_g7_verify_*`). No source defaults changed.

