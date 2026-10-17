# S5-contact hybrid/u-p formulation probe (2026-10-06)

> Live experiment log for the **mixed/hybrid u-p (nearly-incompressible) element
> formulation** as a way to let the soft-foam crash-pad FE survive large compression
> to nominal strain **xi = 0.85** (where fully-integrated tet4/hex8 both invert —
> see `docs/S5_contact_convergence.md` §7).
>
> Written incrementally as the experiments run, not at the end.  Previous agent on
> this problem timed out at 30 min with everything still in its head.  This file is
> the durable record.
>
> **Gate** (per S5 contract §G7 / task spec §2): the *uniaxial* test must
> converge at xi = 0.85 before any contact claim is made.  Contact is only
> attempted if the gate passes.

## Environment

* FEBio 4.13.0 (`D:\Program\FEBioStudio\bin\febio4.exe`)
* Pad material: 6-term compressible `Ogden` from `pad_foam.fit_pad_material()`
  (`<k>1</k>`, `<pressure_model>2</pressure_model>`).
* Uniaxial specimen: 200 mm cube, **free lateral** (x=0 / y=0 symmetry planes),
  bottom z fixed, top z prescribed `-xi * 200` mm.
* Test pad material written through `pad_contact._emit_pad_material` (the
  contact-deck writer) so the deck units are identical to the G7 contact deck
  — i.e. `c / 1e6 → MPa` and `k = 1` (MPa).  Same fitted coefficients either way.

## Experiment log

| # | Setting | rc | xi reached | negJac | FE/law | Notes |
|---|---------|----|-----------:|--------|--------|-------|
| 0a | hex8 cube, 200mm, nx=ny=nz=6, **standard elements**, xi=0.10 | 0 | 0.10 | 0 | (gate only) | baseline ✅ |
| 0b | hex8 cube, 200mm, nx=ny=nz=6, standard, xi=0.25 | 0 | 0.25 | 0 | (gate only) | baseline ✅ |
| 0c | hex8 cube, 200mm, nx=ny=nz=6, standard, xi=0.50 | FebioRunError | <0.50 | **24** | — | baseline FAIL — confirms S5_contact_convergence.md §7 |
| 0d | hex8 cube, 200mm, nx=ny=nz=6, standard, xi=0.85 | FebioRunError | <0.85 | **8**  | — | baseline FAIL |
| 0e | hex8 cube, 200mm, nx=ny=nz=3, standard, xi=0.85 | FebioRunError | <0.85 | **9**  | — | baseline FAIL (also coarse grid) |
| 0f | tet4 cube, 200mm, nx=ny=nz=4, standard, xi=0.10 | 0 | 0.10 | 0 | (gate only) | baseline ✅ |
| 0g | tet4 cube, 200mm, nx=ny=nz=4, standard, xi=0.25 | 0 | 0.25 | 0 | (gate only) | baseline ✅ |
| 0h | tet4 cube, 200mm, nx=ny=nz=4, standard, xi=0.50 | FebioRunError | <0.50 | **2** | — | baseline FAIL (worse than S5_contact_convergence.md §7 stated; the "ξ≈0.80" was for *hex→tet4 split* with `pad_foam` material — the contact-deck material via `_emit_pad_material` fails earlier) |
| 0i | tet4 cube, 200mm, nx=ny=nz=4, standard, xi=0.80 | FebioRunError | <0.80 | **32** | — | baseline FAIL |
| 0j | tet4 cube, 200mm, nx=ny=nz=4, standard, xi=0.85 | FebioRunError | <0.85 | **138** | — | baseline FAIL |
| 0k | tet4 cube, 200mm, nx=ny=nz=3, standard, xi=0.85 | FebioRunError | <0.85 | **40** | — | baseline FAIL |
| 0l | tet4 cube, 200mm, nx=ny=nz=2, standard, xi=0.85 | FebioRunError | <0.85 | **20** | — | baseline FAIL |

> Note on 0h vs S5_contact_convergence.md §7: that doc said "tet4 reaches ξ≈0.80
> (`pad_foam` already verified 2.54 %)".  That claim was for the
> `pad_foam.write_uniaxial_feb` writer (squat 20mm specimen, n_side=2, 6-term
> Ogden with `<k>k_mpa*1e6</k>` in Pa-units).  The contact-deck writer
> `_emit_pad_material` uses `<k>1</k>` in MPa — **the same numerical k** but in a
> different scale.  Result is the same model but **a much smaller relative
> bulk modulus**: at the fitted σ≈1 MPa the contact-deck k=1 MPa means
> k/σ≈1, whereas the unit-deck k=1e6 Pa = 1 MPa means k/σ≈1 — so the issue
> isn't the bulk.  The more likely reason tet4 inverts earlier here is the
> displacement target is **larger per element** (200mm specimen vs 20mm), which
> means more strain per element for the same xi.  Either way: baseline 0h/0i/0j
> confirms the contact path fails no later than xi=0.5.

## FEBio 4.13 mixed/hybrid u-p — what is actually available (verified 2026-10-06)

Sources for the spec assertions below:
- `D:\Program\FEBioStudio\sdk\include\FEBioMech\FE3FieldElasticSolidDomain.h` —
  the **three-field u-p** domain (Simo & Taylor 1991; trilinear u, piecewise-const
  p and Ĵ, augmented-Lagrangian pressure DOF; parameters `laugon`/`atol`/`minaug`/`maxaug`).
- `D:\Program\FEBioStudio\sdk\include\FEBioMech\FEUT4Domain.h` — **UT4** (Unified
  Tetrahedral) = nodally-integrated tet4 with isochoric stabilisation (Gee
  et al. 2009; parameters `alpha`/`iso_stab`).
- `D:\Program\FEBioStudio\sdk\include\FECore\fecore_enum.h` — element-type &
  integration-rule enum.
- `D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf` (also vendored as
  `temp\thums\febio_manual.txt`) — `<SolidDomain type="...">` enumeration.
- `temp\pyfebio_demo\_man4.txt` (FEBio User Manual LyX excerpt) — "By default,
  all of these materials make use of three-field elements (when available for a
  particular element type), as described by Simo and Taylor."

**Key facts** (verbatim):
- **There are NO element-type strings** like `uphex8` / `uptet4` / `up4` / `up8`
  in FEBio 4.13. The "mixed u-p" / "three-field" formulation is activated by
  `<SolidDomain type="three-field-solid">` (or it is default-on for any
  **uncoupled** material — `Ogden` qualifies — when the element shape
  supports it).
- **Three-field is supported on**: hex8, tet10, tet15, tet20, penta6, pyra5
  (NOT tet4 — see `FEBioMech\FE3FieldElasticSolidDomain.h` factory dispatch).
- **Three-field is DEFAULT-ON** for uncoupled materials when the element shape
  supports it. → Our current contact-deck hex8 path is *already* three-field;
  the failure to converge is the *remaining* Q-rule volumetric locking on
  squashed hex8 cells, not the lack of three-field.
- **tet4 has no three-field**: for tet4 the alternative anti-locking formulations
  are (a) **UT4** (`<Elements type="ut4">` + `type="ut4-solid"` SolidDomain),
  (b) upgrade to **tet10/tet15** (then three-field becomes active).
- **Q/U-rule alternative**: `elem_type="HEX8G1"` (UDG, U-rule 1-pt) is the
  "uniform-deformation-gradient" hex element — most robust for nearly-incompressible
  large-strain. Other choices: `HEX8G6` (= selective RI), `HEX8G8` (Q2 default).
- **`pressure_model` interaction**: independent. `pressure_model=2` is the
  Abaqus penalty `U(J) = k/2 (J−1)²`; three-field uses an augmented-Lagrangian
  pressure DOF *in addition to* the material's U(J) — no incompatibility, no
  change needed.

## Experiment log (hybrid)

| # | Setting | rc | xi reached | negJac | FE/law | Notes |
|---|---------|----|-----------:|--------|--------|-------|
| 0a..0l | (baselines) | … | … | … | … | (see §"Baselines") |
| **1a** | hex8 cube 200mm, 6×6×6, **`<SolidDomain type="three-field-solid">`**, xi=0.25 | 0 | 0.25 | 0 | — | ✅ |
| **1b** | hex8 cube 200mm, 6×6×6, **three-field**, xi=0.50 | FebioRunError | <0.50 | **24** | — | **FAIL** — three-field is default-on anyway; explicit type makes no difference |
| **1c** | hex8 cube 200mm, 6×6×6, **three-field**, xi=0.85 | FebioRunError | <0.85 | **8** | — | **FAIL** |
| **2a** | hex8 cube 200mm, 6×6×6, **`<SolidDomain type="udg-hex">`** (HEX8G1), xi=0.25 | 0 | 0.25 | 0 | — | ✅ |
| **2b** | hex8 cube 200mm, 6×6×6, **udg-hex**, xi=0.50 | 0 | 0.50 | 0 | — | **✅ — first formulation that survives past xi=0.5!** |
| **2c** | hex8 cube 200mm, 6×6×6, **udg-hex**, xi=0.85 | FebioRunError | <0.85 | (NAN) | — | **FAIL** — Newton line-search blows up (NAN), not element inversion. With more time-steps this should be recoverable. |
| **3a** | hex8 cube 200mm, 6×6×6, **`<SolidDomain type="sri-solid">`** (selective RI), xi=0.25 | 0 | 0.25 | 0 | — | ✅ |
| **3b** | hex8 cube 200mm, 6×6×6, **sri-solid**, xi=0.50 | FebioRunError | <0.50 | **4** | — | **FAIL** |
| **3c** | hex8 cube 200mm, 6×6×6, **sri-solid**, xi=0.85 | FebioRunError | <0.85 | **24** | — | **FAIL** |
| **4a** | tet4 cube 200mm, 4×4×4, **`<SolidDomain type="ut4-solid">`** (UT4 default), xi=0.25 | 0 | 0.25 | 0 | — | ✅ |
| **4b** | tet4 cube 200mm, 4×4×4, **ut4-solid**, xi=0.50 | FebioRunError | <0.50 | **8** | — | **FAIL** |
| **4c** | tet4 cube 200mm, 4×4×4, **ut4-solid**, xi=0.85 | FebioRunError | <0.85 | **80** | — | **FAIL** |
| **4d** | tet4 cube 200mm, 4×4×4, **ut4-solid `<iso_stab>1</iso_stab>`** (isochoric-only stab), xi=0.50 | FebioRunError | <0.50 | **14** | — | FAIL |
| **4e** | tet4 cube 200mm, 4×4×4, **ut4-solid `<iso_stab>1</iso_stab>`**, xi=0.85 | FebioRunError | <0.85 | **94** | — | FAIL |
| **4f** | tet4 cube 200mm, 4×4×4, **ut4-solid `<alpha>0.5</alpha>`** (more regular tet, less nodal), xi=0.50 | FebioRunError | <0.50 | **2** | — | FAIL |
| **4g** | tet4 cube 200mm, 4×4×4, **ut4-solid `<alpha>0.5</alpha>`**, xi=0.85 | FebioRunError | <0.85 | **108** | — | FAIL |

### Key takeaways so far

- **Three-field (the "mixed u-p") does not save us on hex8**: it's default-on
  for `Ogden` anyway, and even forcing `type="three-field-solid"` explicitly
  (1b) re-creates the same negJac=24 failure that `elastic-solid` exhibits
  at xi=0.50. This is consistent with the FEM literature: three-field solves
  the u-p constraint but the Q2 deviatoric part of HEX8G8 still locks when the
  cell aspect ratio gets extreme.
- **`HEX8G1` (UDG) is the only formulation that survives xi=0.50** (2b ✅);
  it uses an U-rule 1-pt volumetric + hourglass control instead of Q-rule +
  penalty.
- **`HEX8G1` fails at xi=0.85 with NAN, not negJac** (2c) — the Newton
  line-search diverges after the residual gets large. **This is a stepper /
  convergence problem, NOT an element-inversion problem** — meaning the
  physics can be made to converge with a better stepper configuration (more
  steps, tighter `dtmax`, Broyden QN).  This is **the first candidate that
  might pass the gate**.

### Next step: harden HEX8G1 at xi=0.85

Plans: (a) more time steps (400 → 2400), (b) shrink `dtmax` from `1/nsteps`
to `step_size` (cap at 1 step per solve), (c) switch QN to **Broyden** (often
more stable than BFGS at large strain), (d) add line-search,
(e) try fine grid (10×10×10).

---

## ✅ OUTCOME — the gate PASSES with UDG (HEX8G1)

`temp/opensim_fe/s5_t9/udg_harden.py` + its logs. **Runs converged to the end of
the ramp (`end_t = 1`) at xi = 0.85:**

| config | mesh | steps | QN | end_t |
|---|---|---|---|---|
| udg 6×6×6 ξ=0.5 (sanity) | 6³ | 2400 | BFGS | **1.0** ✅ |
| udg 6×6×6 ξ=0.85 | 6³ | 2400 | BFGS | **1.0** ✅ |
| udg 6×6×6 ξ=0.85 | 6³ | 2400 | BFGS | **1.0** ✅ |
| udg 10×10×10 ξ=0.85 | 10³ | 400 | BFGS | **1.0** ✅ |
| udg 10×10×10 ξ=0.85 | 10³ | 2400 | BFGS | **1.0** ✅ |
| udg 10×10×10 ξ=0.85 | 10³ | 2400 | Broyden | **1.0** ✅ |
| udg 12×12×12 ξ=0.85 | 12³ | 400 | Broyden | 0.667 (agent timed out) |

Working recipe (uniaxial): `<SolidDomain type="udg-hex">` + `hex8` elements +
**2400 steps**, `dtmax = 1/nsteps`, `cutback = 0.125`, `max_retries = 20`,
`<qn_method type="BFGS"><max_ups>10</max_ups></qn_method>`, **and a load-curve
ramp** on the prescribed displacement.

**Interpretation**: the large-strain limit was the **element formulation**, and
UDG's uniform-deformation-gradient / 1-pt volumetric treatment removes it. This
is the first thing in the whole S5 investigation that passes xi = 0.85.

## ⚠️ But UDG does NOT (yet) fix the CONTACT deck

Follow-up after the agent timed out (`temp/opensim_fe/s5_t9/udg_contact.py`;
`pad_domain_type="udg-hex"` is now an opt-in param of `build_pad_contact_feb` /
`run_pad_contact`):

| config | result |
|---|---|
| udg-hex, default `<Control>`, sink=170 | ❌ FAIL at t≈0.003 |
| udg-hex, hardened `<Control>`, sink=170 | ❌ **"NAN detected"**, residual 2.1e8, at t≈0.003 |
| udg-hex, hardened, sink=50 (tet4 handles this) | ❌ FAIL at t≈0.01 (≈0.5 mm indentation) |

Not element inversion — a **numerical blow-up at contact initiation**, i.e. a
UDG-element ↔ contact-formulation incompatibility (or a deck issue). A
`node_reloc=0` variant did not even parse. **Open thread.**

### Next steps for the contact
1. FEBio's constraints on **UDG + sliding contact** (does `HEX8G1` expose usable
   facets? interaction with `search_radius`/`node_reloc`?).
2. Indenter as UDG/rigid too, or a **tied** contact.
3. **tet10/tet15** — three-field becomes active for them; if they pass the
   uniaxial gate at ξ=0.85, use them for contact (tet contact facets are
   well-trodden).

**Bottom line**: the *material* problem is solved (UDG); the *contact* problem is
not. The trustworthy S5 numbers remain `xi <= 0.25`, FE/law ≈ 1.6; the G7 sink
is still unavailable **for the contact deck**.