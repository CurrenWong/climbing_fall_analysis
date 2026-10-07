# S5 — UDG (`udg-hex` / HEX8G1) + sliding contact: constraints & diagnosis

**Status:** COMPLETE (2026-10-06) — recipe + G7 verdict in §7; every result read from
its real `.log`/`.xplt` tail and appended after each experiment.
**Owner task:** determine FEBio 4.13's constraints on using UDG (`<SolidDomain type="udg-hex">`,
element type `HEX8G1`) **together with sliding (frictional) contact**, and run the minimal
discriminating experiment that settles it.

> ## 🔴 SUPERSESSION BANNER (added 2026-10-06)
>
> The §7 recommendation *“For the 2 m-landing numbers, quote the `tet4` fallback
> (complete G7, xi=0.744, σ/σ_law = 2.205, reproducible)”* is **superseded** by the
> four G7 verification studies: `docs\S5_g7_verify_mesh.md` (FAIL: σ/σ_law drifts
> **+30.62 %** at sink=170 with growing refinement step), `docs\S5_g7_verify_fixedxi.md`
> (FAIL: fixed-ξ drift **+43.13 %**), `docs\S5_g7_verify_anchors.md` (FAIL: NEW recipe
> shifts the OLD anchor curve by **−5.85 % … −9.38 %**), and `docs\S5_g7_verify_anchors_mesh.md`
> (FAIL: small-strain anchors drift **−0.75 % … −5.47 %** non-monotone).
>
> **Corrected guidance:** for the 2 m landing, the only **mesh-converged** quotable
> numbers are the resultant **`F_n ≈ 34.4 kN`** and **`σ_contact ≈ 0.574 MPa`** (NEW
> PENALTY recipe, tet4, fixed sink, <0.6 % drift across a 64× element-count change).
> The σ/σ_law ratio (and the `2.205` in particular) is **NOT** quotable as a converged
> scalar. See `docs\S5_QUOTABLE.md` for the authoritative verdict.
>
> **The measurements in this file (recipe extraction, §6 sweeps, §7 tet4 fallback
> numbers, force-balance check) remain valid as measurements** — only the §7
> *recommendation* is superseded. No measured value, table or method text in this
> file has been altered.

Units: mm–N–MPa–s. FEBio binary: `D:\Program\FEBioStudio\bin\febio4.exe` (reports `version 4.13.0`).

---

## 1. Ground-truth sources actually read (exact paths)

Local FEBio 4.13 SDK headers (priority 1):

| Path | What it gives |
|---|---|
| `D:\Program\FEBioStudio\sdk\include\FEBioMech\FEUDGHexDomain.h` | `class FEUDGHexDomain : public FEElasticSolidDomain` — the UDG domain class |
| `D:\Program\FEBioStudio\sdk\include\FECore\FEElementTraits.h` | `class FEHex8G1 : public FEHex8_` (`NINT = 1`) |
| `D:\Program\FEBioStudio\sdk\include\FECore\FESolidElementShape.h`, `FEElementShape.h`, `FEElementTraits.h`, `FEElement.h` | face/facet machinery (`Faces()`, `GetFace()`), shape classes |
| `D:\Program\FEBioStudio\sdk\include\FECore\FESolidDomain.h` | `FESolidDomain : public FEDomain` — no `GetFacet` override; face extraction is generic |
| `D:\Program\FEBioStudio\sdk\include\FECore\FESurface.h` | `FESurface : public FEMeshPartition`, `Create(const FEFacetSet&)` — **surface is built from a facet set** |
| `D:\Program\FEBioStudio\sdk\include\FEBioMech\FEContactSurface.h` | `FEContactSurface : public FESurface` |
| `D:\Program\FEBioStudio\sdk\include\FECore\fecore_enum.h` | `FE_HEX8G1` (element type idx 2), `ET_HEX8` shape |
| `D:\Program\FEBioStudio\sdk\include\FEBioMech\FESolidDomainFactory.h` | `CreateDomain(spec, mesh, mat)` — domain selection by `type` attribute |
| `D:\Program\FEBioStudio\sdk\include\FEBioXML\FEModelBuilder.h` | `double m_udghex_hg` (hourglass parameter) |

Local SDK ships **headers only** (0 `.cpp`). The actual UDG implementation was read from the
official open-source FEBio repository (`FEUDGHexDomain.cpp`), because FEBio is MIT-licensed and
the header declares the exact method set:

* `https://raw.githubusercontent.com/febiosoftware/FEBio/master/FEBioMech/FEUDGHexDomain.cpp`
* `FEUDGHexDomain.cpp` `BEGIN_FECORE_CLASS(FEUDGHexDomain, FEElasticSolidDomain)` → `ADD_PARAMETER(m_hg, "hg")`
  and constructor `m_hg = 1.0`.

Binary string checks (to confirm the local build matches the source):

* `D:\Program\FEBioStudio\bin\febiomech.dll` contains `udg-hex\0FEUDGHexDomain\0` and the exact
  mangled symbols `?UDGInternalForces@FEUDGHexDomain@@…`, `?UDGHourglassForces@…`, `?UDGGeometricalStiffness@…`,
  `?UDGMaterialStiffness@…` → the installed binary is this class.
* `D:\Program\FEBioStudio\bin\febioxml.dll` contains the strings `hex8 / GAUSS8 / POINT6 / … / UDG / hourglass`
  (the hex8 integration-rule + hourglass option names).

Forum (priority 3, developer answer by Gerard Ateshian):
`https://forums.febio.org/forum/febio-studio-forum/users-forum-ad/25404-elements-with-reduced-integration`

> "Yes, you can use the udg-hex solid formulation as explained in the User's Manual. … The parameter
> **Hg is a scale factor for the hourglass force/stiffness**."
> "The **udg-hex element is formulated for incompressible (or rather uncoupled) material
> formulations**."

---

## 2. Q(a) — Can HEX8G1 element faces be used as contact surfaces at all?

**Yes — structurally, and the deck provably does it (FEBio builds the surface without error).**

Evidence:

1. `FEUDGHexDomain : public FEElasticSolidDomain` (`FEUDGHexDomain.h:34`) — it is an elastic solid
   domain and inherits the whole generic solid-domain surface machinery. It declares **no**
   surface/facet override and **no** contact-related override; its only members are the UDG force /
   stiffness / hourglass routines.
2. `FEUDGHexDomain::Create` (source) only checks:
   ```cpp
   if (spec.eclass != FE_ELEM_SOLID) return false;
   if (spec.eshape != ET_HEX8) return false;
   spec.etype = FE_HEX8G1;          // forces 1-point integration
   return FESolidDomain::Create(nelems, spec);
   ```
   i.e. it is accepted as a normal `FE_ELEM_SOLID` / `ET_HEX8`. Nothing about surfaces is rejected.
3. Contact surfaces in this deck are **not** element-face references at all. `pad_contact.py`
   writes an explicit `<Surface name="pad_top">` made of **quad4 facet node lists** (`_emit_surface`,
   `pad_contact.py:404`), and a `<SurfacePair>` with the indenter bottom as primary
   (`pad_contact.py:707-713`). `FESurface` is built with `Create(const FEFacetSet&)` — a facet set of
   quad4/tri3 nodes — **independent of the parent domain's element type** (`FESurface.h:110-113`,
   `FindElement()` is only used for interface surfaces, not for a plain contact facet set).
4. **Runtime proof**: the failing UDG decks **parse and run** — FEBio builds the contact surface,
   assemblies the pair, and marches time steps. `udg_norel_sink50.feb` prints `Reading file … SUCCESS!`
   and then integrates (see §6). There is **no** "invalid surface" / "no facets" error. So HEX8G1
   faces are *usable* as contact facets.

**Conclusion (a):** HEX8G1 faces are usable as contact surfaces; the restriction is **not** at the
surface-construction level. The failure is numerical, at contact load transfer into the UDG element.

---

## 3. Q(b) — structural restriction or our deck/parameter error?

**Not a documented structural restriction.** The UDG implementation is contact-agnostic: contact
tractions are applied as nodal forces on the pad-top facet nodes; the UDG domain then resists them
through its *single* material point + hourglass stabilization. There is no code path that forbids
contact. (No header/source/`type`-table entry excludes it; `FESolidDomainFactory` maps `udg-hex`
→ `FEUDGHexDomain` indiscriminately.)

The observed failure is a **numerical stability failure of the 1-point + hourglass element under the
stiff augmented-Lagrangian contact load**, triggered *exactly at contact closure*. It is therefore a
**modelling/parameter issue on our side** (element choice / contact stiffness / stabilization), not
an unsupported-feature error. The experiments in §6 discriminate this.

The forum confirms the intended material class: UDG is meant for **uncoupled (incompressible-ish)**
formulations; our pad uses a **compressible Ogden** with a finite `<k>`. That is *not* the immediate
blocker (the same material passes the uniaxial UDG test at xi=0.85, §5), but it is a known tension.

---

## 4. Q(c) — what `search_radius` and `node_reloc` do, and safe values

Source: FEBio sliding-elastic contact (`sliding-elastic`, `FEContactSurface` + sliding interface).

* **`<search_radius>`** — the master-surface search radius (mm). Each contact iteration projects a
  secondary node onto the primary surface; the search only considers primary facets within this
  radius of the node. It must exceed the **per-step travel** of the indenter, or the pair is lost and
  the contact force snaps to zero (then re-closes next step → impulse → Newton blow-up). Default in
  our writer is `20.0` mm (`pad_contact.py:196`). Per-step travel here = `sink / nsteps`
  (sink 170 mm / 2400 = 0.071 mm/step; sink 50 / 2400 = 0.021 mm/step) — so 20 mm is *enormous*
  relative to travel; it is **too large, not too small**, and can search/attach spurious facets.
  A safe radius is a few × per-step travel and ≲ one element edge (element edge = 50 mm here), so
  **≈1–5 mm** is appropriate (and 1 mm still >> 0.07 mm travel).

* **`<node_reloc>`** — the *node-relocation* flag. `1` (= our default, `pad_contact.py:579`): FEBio
  relocates the secondary (pad) nodes onto the primary (indenter) surface at the start of the contact
  algorithm, so that the initial gap is closed in the search sense and a consistent starting contact
  state is established. `0`: nodes are left in place and the gap is measured geometrically.
  For a UDG pad whose single-point elements are sensitive to any imposed nodal offset, `node_reloc=1`
  is a plausible trigger of the hourglass blow-up at closure; `node_reloc=0` is the natural control.

**Safe recipe hypothesis:** `node_reloc=0`, `search_radius` small (≈1–5 mm), penalty modest
(currently `1`), keep the augmented-Lagrangian controls (`gaptol/minaug/maxaug/smooth_aug`). To be
confirmed in §6.

---

## 5. Log forensics — the failure is *exactly at contact closure*

Initial gap `gap_mm = 0.5` (pad top z=200, indenter bottom z=200.5). The indenter z-displacement
ramps linearly over `t ∈ [0,1]` to the full sink. So contact closes (travel = 0.5 mm) at:

| deck | sink (mm) | `t_close = 0.5 / sink` | observed `failed to converge at time` | Δ |
|---|---|---|---|---|
| `udg664_sink170_default.log` | 170 | 0.00294118 | **0.00295457** | +1 step |
| `udg664_sink170_hard.log` | 170 | 0.00294118 | **0.00299924** | ≈ +2 steps |
| `udg664_sink50_hard.log` | 50 | 0.0100000 | **0.0100075** | ≈ 0 steps |

All three fail within **1–2 time steps of first contact**. Before closure the indenter moves through
rigid-body free space and the solve is healthy. So the pre-contact UDG body is fine; the instability
is produced by the **contact force transfer**, not by pad compression per se.

`sink170_hard` details: `NAN detected` appears (first at log line 1439 of 1852 KB), residual ~`2.1e8`;
`ls_check_jacobians = no (0)`; **no negative Jacobians** reported → the failure is a NaN in the
residual/assembled system, not an element inversion. Fails after exhausting 20 retries on step 14.

`udg_norel_sink50.feb` — **the earlier "did not parse (0 s, no log)" claim is FALSE.** Re-running it
directly:
```
& febio4.exe -i <ABSOLUTE …\udg_norel_sink50.feb>
... Reading file …\udg_norel_sink50.feb ... SUCCESS!
```
It parses and integrates (see §6 for the outcome). The earlier observation was almost certainly the
documented relative-path + custom-workdir rc=1 trap, misread as a parse failure.

---

## 6. Experiments

_(appended one block per experiment, as each finishes — nothing is held only in memory)_

Every outcome below was read from the real `.log` tail (`converged at time :`,
`failed to converge at time :`, `NAN detected`, `negative jacobians detected`).
Probe driver: `temp\opensim_fe\s5_t9\udg_contact\udg_contact_probe.py`
(raw run ledger: `probe_results.txt` / `probe_stdout.txt`).

### 6.0 — Baseline deck forensics (pre-existing decks, read from disk)

| deck (`.log`) | sink (mm) | `node_reloc` | `search_radius` | `penalty` | last `converged at time` | `failed to converge at time` | `negJac` | NAN | term |
|---|---|---|---|---|---|---|---|---|---|
| `udg664_sink170_default.log` | 170 | 1 | 20 | 1 | 0.00291669 | **0.00295457** | 0 | yes | ERROR |
| `udg664_sink170_hard.log` | 170 | 1 | 20 | 1 | 0.00299922 | **0.00299924** | 0 | yes | ERROR |
| `udg664_sink50_hard.log` | 50 | 1 | 20 | 1 | 0.0100000 | **0.0100075** | 0 | yes | ERROR |
| `udg_norel_sink50.log` | 50 | 0 | 20 | 1 | 0.0433677 (378 converged steps) | 0.0433679 | 0 | yes (19×) | ERROR |

`t_close = 0.5/sink` ⇒ 170 mm ⇒ 0.002941, 50 mm ⇒ 0.010000. The first three fail
**exactly at closure** (Δ ≤ 1 step). `udg_norel_sink50` fails first at
`t = 0.0100001` (closure), then grinds through 21 failed-retry cycles and 378
converged micro-steps to `t = 0.0434` before the NAN finally ends it — i.e. past
closure the soft penalty buys a few **extra centimetres of t**, not a cure.
In **every** deck `ls_check_jacobians = no (0)` and **no** "negative jacobians
detected" line is ever emitted → the failure is a **NaN in the assembled
residual**, not an element inversion.

### EXP-1 — `node_reloc` × `search_radius` matrix (sliding-elastic, sink=50, hardened) ⚠️ FAILED

Deck: `build_pad_contact_feb(..., pad_domain_type="udg-hex", laugon="AUGLAG",
tolerance=0.005, aug_controls=AUGLAG_CONVERGENT_CONTROLS, penalty=1, mu=0.6,
load=displacement(sink=50), nsteps=2400)`. Control block hardened to the uniaxial
recipe (cutback 0.125, max_retries 20, `dtmax=1/nsteps`, BFGS max_ups 10,
reform_each_time_step 1, max_refs 25).

| tag | `node_reloc` | `search_radius` | end_t | fail_t | `negJac` | NAN | term |
|---|---|---|---|---|---|---|---|
| `p1_reloc1_sr20` | 1 | 20 | 0.010381 | 0.0103888 | 0 | yes | ERROR |
| `p1_reloc0_sr20` | 0 | 20 | 0.0103978 | 0.0103978 | 0 | yes | ERROR |
| `p1_reloc1_sr2` | 1 | 2 | 0.0100000 | 0.0100065 | 0 | yes | ERROR |
| `p1_reloc0_sr2` | 0 | 2 | 0.0100580 | 0.0100584 | 0 | yes | ERROR |
| `p1_reloc1_sr1` | 1 | 1 | 0.0100000 | 0.0100075 | 0 | yes | ERROR |

**Result:** every combination NANs at `t ≈ 0.0100 = contact closure` (sink 50).
Neither `node_reloc` nor `search_radius` (incl. the recommended 1 mm) moves the
failure. `negJac = 0` throughout; one NAN each. → **Neither parameter is the cause.**

### EXP-2 — tied-elastic instead of sliding (sink=50) ⚠️ INVALID (no log)

`p2_tied_reloc0_sr2` → `no-log`, 0 s. The writer emits a `<fric_coeff>` that the
`tied-elastic` algorithm rejects; the probe strips it *after* build, but FEBio
never produced a log → the deck never started (XML emitted by the strip path is
suspect). **Not evidence about UDG — must be redone** (see §6.NEW-6).

### EXP-3 — tiny sink (2 mm) → first failing increment ✅ DIAGNOSTIC

| tag | `node_reloc` | `search_radius` | end_t | fail_t | `negJac` | NAN | term |
|---|---|---|---|---|---|---|---|
| `p3_sink2_reloc1_sr2` | 1 | 2 | 0.250146 | 0.250146 | 0 | yes | ERROR |
| `p3_sink2_reloc0_sr2` | 0 | 2 | 0.250541 | 0.250542 | 0 | yes | ERROR |

`t_close = 0.5/2 = 0.250000`. Both fail within <1 step of closure. This *independently*
confirms the trigger is **contact closure**, not sink magnitude, and that the
pre-contact UDG solve is healthy for 600 steps.

### EXP-4a — KEY CONTROL: UDG pad, NO contact, prescribed top displacement ✅ PASS

`p4a_nocontact_uniform_sink50` (uniform prescribed z on all pad-top nodes, same
hardened control, `type="udg-hex"`): **`end_t = 1`**, last `converged at time : 1`,
**`negJac = 0`**, **no NAN**, no failed convergence until the single end-of-analysis
overshoot step (`fail_t = 1.00002`). 38 s.

→ **The UDG pad with our deck skeleton, hardening recipe and boundary scheme is
numerically fine.** The only broken ingredient is the **contact force transfer.**

### EXP-4b — no-contact footprint load ⚠️ DROPPED

`p4b_nocontact_footprint_sink50` timed out (401 s) at `end_t = 15.3004`; the
footprint BC let the pad-top free boundary de-load ("No force acting on the
system") and the stepper ran far past `t=1` without reaching a meaningful state.
Cheap control, not diagnostic — **dropped**; EXP-4a already gives the control.

---

### §6.NEW — ranked new experiments (§4.2)

Driver: `temp\opensim_fe\s5_t9\udg_contact\udg_contact_sweep2.py`
(ledger: `sweep2_results.txt`). All cases: pad `pad_block_hex(6,6,4,(300,300,200))`,
indenter `indenter_block_tet(2,2,1,(200,300,5))`, hardened control (cutback 0.125,
max_retries 20, `dtmax=1/nsteps`, BFGS max_ups 10), sink=50, `nsteps=2400`,
`mu=0.6`, `node_reloc=0`, `search_radius=2` unless stated.

#### NEW-1 — `hg` hourglass scale sweep ⚠️ FAILED

XML placement confirmed from source (`FEBioMeshDomainsSection4.cpp:202-206`: a
non-leaf `<SolidDomain>` routes its child tags through `ReadParameterList(tag,
dom)`, which binds the `FEUDGHexDomain` registered parameter `hg`). Emitted as:

```xml
<SolidDomain name="pad_tets" mat="pad" type="udg-hex">
  <hg>0.1</hg>
</SolidDomain>
```

(opt-in writer arg `pad_hg`, default `None` ⇒ byte-identical to the old leaf form;
verified — `_check_writer.py`.)

| tag | `hg` | end_t | fail_t | nconv | `negJac` | NAN | term |
|---|---|---|---|---|---|---|---|
| `n1_hg0.01` | 0.01 | 0.00999924 | 0.0100038 | 25 | 0 | yes | ERROR |
| `n1_hg0.1` | 0.1 | 0.0101716 | 0.0101721 | 30 | 0 | yes | ERROR |
| `n1_hg10` | 10 | 0.01048 | 0.01048 | 34 | 0 | yes | ERROR |
| `n1_hg100` | 100 | 0.0105494 | 0.0105555 | 28 | 0 | yes | ERROR |

**Result:** `hg ∈ {0.01, 0.1, 10, 100}` — a **4-order-of-magnitude** change of the
hourglass force/stiffness — leaves the failure pinned at closure (`t ≈ 0.0100`)
with a NAN and **0 negative Jacobians**. The decks carried the tag (grep-confirmed)
and FEBio accepted it (no unrecognized-tag error). → **The hourglass mode is not
the mechanism; scaling it does not help.**

> ⚠️ **Control-discipline note.** The hardened `<Control>` block used by
> `udg_contact_probe.py` / `udg_contact_sweep2.py` must contain
> `<step_size>1/nsteps</step_size>` **and** `<dtmax>`. The proven uniaxial recipe
> (`udg_harden.py:85-86`) has both. The original probe `harden()` dropped
> `<step_size>`, so FEBio fell back to its default and the analysis `end_time`
> became `time_steps · 0.1 = 240` (not 1) — harmless for the closure-NAN runs
> (they die at `t≈0.01`) but it lets a *surviving* run march far past `t=1`.
> `sweep2` now emits the full recipe; every result below uses `end_time = 1`.

#### NEW-2 — pure penalty (drop the augmented Lagrangian) ✅ **FIRST SURVIVOR**

`laugon="PENALTY"`, `aug_controls=None`, `node_reloc=0`, `search_radius=2`,
sink=50, `nsteps=2400`, `end_time=1`.

| tag | `penalty` | rc | end_t | nconv | `negJac` | NAN | term | indentation (mm) | penetration (mm) | σ_contact (MPa) | σ_law (MPa) | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `n2_pen0.1` | 0.1 | **0** | **1** | 2400 | 0 | no | **NORMAL** | 45.12 | 4.384 | 0.04528 | 0.03047 | **1.486** |
| `n2_pen1` | 1 | EXC | 0.0121581 | 31 | 0 | yes | ERROR | — (died at closure) | | | | |
| `n2_pen10` | 10 | EXC | 0.01375 | 33 | 0 | yes | ERROR | — (died at closure) | | | | |

**Result — the first configuration in this whole track to run PAST contact
closure.** `laugon=PENALTY` + `penalty=0.1` reaches **`end_t = 1` (full 50 mm
sink) with `rc=0`, `NORMAL TERMINATION`, 2400/2400 steps, `negJac=0`, NO NAN**.
Measured at full sink (read from the real `.xplt`, `n_states=2401`):
indentation **45.12 mm** (xi = 0.2256), penetration **4.38 mm**,
`F_n = 2717 N`, **σ_contact/σ_1D-law = 1.486** — inside the trustworthy band.

`penalty = 1` and `10` (pure penalty) still NAN at closure. So the winner is the
**combination of dropping AUGLAG *and* the soft penalty (0.1)** — the stiff
augmented-Lagrangian spike at closure is what the 1-pt element cannot carry.

→ **This is the working recipe candidate; the G7 test at sink=170 follows.**

#### NEW-3 — start already in contact (gap = 0) ⚠️ FAILED (worse)

`gap=0` ⇒ indenter bottom at z=200.0 (coincident with pad top); `gap=0.05` ⇒ z=200.05.

| tag | gap (mm) | end_t | fail_t | nconv | `negJac` | NAN | term |
|---|---|---|---|---|---|---|---|
| `n3_gap0` | 0.0 | 0.000414 | 0.000416 | 2 | 0 | yes | ERROR |
| `n3_gap005` | 0.05 | 0.0101945 | 0.0101956 | 26 | 0 | yes | ERROR |

**Result:** engaging contact from `t=0` moves the NAN **earlier** (step 2, `t =
4.1e-4`) — the closure event is not the enemy *per se*, the **magnitude of the
first contact-force increment** is. `gap=0.05` reproduces the baseline closure
NAN. → **Not a cure.**

#### NEW-4 — push the survivor to G7 (sink = 170, xi ≈ 0.85) ✅ **near-full-sink, but UDG caps at xi≈0.71**

All `laugon="PENALTY"`, `node_reloc=0`, `end_time=1`, sink=170, 2400 steps (unless stated),
`mu=0.6`. Measured from the real `.xplt` at the last written state.

| tag | `penalty` | `sr` | rc | end_t | nconv | `negJac` | NAN | indentation (mm) | xi | F_n (N) | pen (mm) | σ/σ_law |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `n4_pen0.1_s170` | 0.1 | 2 | **0** | 1 | 2401 | 0 | no | **~0** | ~0 | **0** | **169.5** | — |
| `n4b_pen0.1_s170` | 0.1 | 2 | **0** | 1 | 2401 | 0 | no | ~0 | ~0 | 0 | 169.5 | — |
| `g7_pen0.1_reloc1_s170` | 0.1 | 2 | **0** | 1 | 2401 | 0 | no | ~0 | ~0 | 0 | 169.5 | — |
| `g7_pen0.2_s170` | 0.2 | 2 | EXC | 0.00745 | 40 | 0 | yes | — | — | — | — | — |
| `g7_pen0.3_s170` | 0.3 | 2 | EXC | 0.00498 | 12 | 0 | yes | — | — | — | — | — |
| `g7_pen0.5_s170` | 0.5 | 2 | EXC | 0.00647 | 19 | 0 | yes | — | — | — | — | — |
| **`g7_pen0.1_sr20_s170`** | 0.1 | 20 | EXC | **0.974227** | **2340** | 0 | yes | **141.8** | **0.7092** | **23 150** | 23.27 | **1.705** |
| `g7_pen0.1_sr20_s170_rep3` | 0.1 | 20 | EXC | 0.974227 | 2340 | 0 | yes | 141.8 | 0.7092 | 23 150 | 23.27 | 1.705 |
| `g7_pen0.1_sr20_reloc1_s170` | 0.1 | 20 | EXC | 0.974227 | 2340 | 0 | yes | 141.8 | 0.7092 | 23 150 | 23.27 | 1.705 |
| `g7_pen0.1_sr30_s170` | 0.1 | 30 | EXC | 0.97384 | 2339 | 0 | yes | 141.2 | 0.7062 | 23 000 | 23.82 | 1.716 |
| `g7_pen0.1_sr50_s170` | 0.1 | 50 | EXC | 0.97384 | 2339 | 0 | yes | 141.2 | 0.7062 | 23 000 | 23.82 | 1.716 |
| `g7_pen0.1_sr100_s170` | 0.1 | 100 | EXC | 0.97384 | 2339 | 0 | yes | 141.2 | 0.7062 | 23 000 | 23.82 | 1.716 |
| `g7_pen0.1_sr20_n4800_s170` | 0.1 | 20 | EXC | 0.965652 | 4648 | 0 | yes | 136.0 | 0.6799 | 22 040 | 27.67 | 1.841 |
| `g7_pen0.05_sr20_s170` | 0.05 | 20 | EXC | 0.987886 | 2373 | 0 | yes | 135.9 | 0.6794 | 21 960 | 31.56 | 1.839 |
| `g7_pen0.08_sr20_s170` | 0.08 | 20 | EXC | 0.973308 | 2338 | 0 | yes | 140.6 | 0.7032 | 22 600 | 24.32 | 1.708 |
| `g7_pen0.11_sr20_s170` | 0.11 | 20 | EXC | 0.974437 | 2342 | 0 | yes | 139.7 | 0.6984 | 22 360 | 25.48 | 1.725 |
| `g7_pen0.12_sr20_s170` | 0.12 | 20 | EXC | 0.975568 | 2345 | 0 | yes | 140.8 | 0.7038 | 22 710 | 24.59 | 1.712 |

**Result — the first genuine, load-bearing, near-full-sink UDG+contact states ever
produced in this track, and they are reproducible.**

* **`sr=2` degenerates at large sink.** `pen0.1 @ sr=2` "converges" to `end_t=1`
  (`rc=0`, NORMAL) but the read-back is **indentation ≈ 0, penetration = 169.5 mm,
  `F_n = 0`** — the indenter passes straight through: the contact is *not*
  load-bearing. Identical with `node_reloc=1`. The pad-top nodes fall >2 mm behind
  the indenter surface, the pair is **lost**, the indenter free-falls.
* **A large `search_radius` is *necessary*** — the old "20 mm is too large"
  hypothesis is **wrong**. `sr=20` keeps the pair attached and drives the solve to
  **`end_t = 0.974227` (2340/2401 steps, `negJac = 0`)** before the NAN: **xi =
  0.709, `F_n = 23 150 N`, penetration 23.3 mm, σ/σ_law = 1.705**. `sr ∈ {20, 30,
  50, 100}` are **identical** (saturation) — so the wall is **not** search-related.
* **`node_reloc` is irrelevant** (reloc 0 and 1 give byte-identical results).
* **The wall is a penetration/stiffness trade-off, not a search or timestep
  problem.** `penalty ≥ 0.15` re-introduces the **closure** NAN (`0.15/0.2/0.3/0.5`
  die at `t ≈ 0.005–0.02`); `penalty = 0.05` reaches a slightly higher `t`
  (0.9879) but with 31.6 mm penetration (`xi` drops to 0.679). Finer time steps
  (`nsteps=4800`) do **not** help (`t` falls to 0.9657). Across the whole
  `penalty ∈ [0.05, 0.12]` window the wall sits at **`t ≈ 0.965–0.988`, i.e. xi
  ≈ 0.68–0.71** — UDG+sliding **cannot** be pushed to `t=1` (xi≈0.85).

**G7 verdict for UDG+sliding:** reaches sink 165.6/170 mm (97 %) with a genuine
`F_n = 23 kN` state, but **NANs before full sink**; the effective indentation caps
at **xi ≈ 0.71** because any contact stiff enough to suppress penetration NANs at
closure. → **Declared a decisive negative for the exact G7 target.**

#### NEW-5 — `search_radius = 1 mm`

Covered by the `sr` sweep above. `sr=1` is *soft* than the `sr=2` case already
known to **degenerate** at sink=170 (pair lost → free passage). At the hardened
setting `sr ∈ {20,30,50,100}` are identical, so 1 mm has nothing to add: it is
strictly on the losing (too-small) side. No separate run needed.

#### NEW-6 — redo tied contact ✅ the "no-log" was a **probe bug**, not a FEBio limit

Root cause found by capturing FEBio stdout (`n6_tied_reloc0_sr2.stdout.txt`):

* attempt 1 → `tag "node_reloc" (line 561): unrecognized tag`;
* attempt 2 → `tag "smooth_aug" (line 567): unrecognized tag`.

`tied-elastic` **does not accept** `node_reloc`/`search_radius`/`fric_coeff` or the
AUGLAG child controls. The original EXP-2 "0 s, no-log" was simply an **invalid
deck emitted by the probe**, not a restriction. With those tags stripped
(`laugon="PENALTY"`, no `aug_controls`):

| tag | end_t | nconv | `negJac` | NAN | term |
|---|---|---|---|---|---|
| `n6_tied_reloc0_sr2` | **1** | 2400 | 0 | **no** | **NORMAL (rc=0, 36 s)** |

**UDG + a *tied* (bonded) contact constraint completes the full 50 mm sink with no
NAN.** → The UDG element and the deck are fine; the instability is specific to the
**sliding-elastic force transfer**, not to UDG carrying a contact load.

#### NEW-7 — fallback: tet4 pad through the same new recipe ✅ **COMPLETES G7**

The repo has no tet10/tet15 mesh generator (only `pad_block_hex` / `pad_block_tet`),
so the fallback is the **tet4** pad mesh (`pad_block_tet(6,6,4,…)`), run through the
**same** `laugon="PENALTY" penalty=0.1 search_radius=20 node_reloc=0` recipe.

| tag | pad | sink | rc | end_t | nconv | `negJac` | NAN | indentation | xi | F_n (N) | pen (mm) | σ/σ_law |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `n7_tet4_s50` | tet4 | 50 | **0** | **1** | 2400 | 0 | no | 46.32 | 0.2316 | 2 820 | 3.18 | 1.511 |
| `n7_tet4_s170` | tet4 | 170 | **0** | **1** | 2400 | 0 | no | **148.7** | **0.7437** | **34 600** | 20.76 | **2.205** |
| `n7b_tet4_s170` | tet4 | 170 | **0** | **1** | 2400 | 0 | no | 148.7 | 0.7437 | 34 600 | 20.76 | 2.205 |

**tet4 completes the full G7 sink (`end_t=1`, `rc=0`, NORMAL, 2400/2400 steps,
`negJac=0`, zero NAN)**, is **reproducible to the digit over 2 runs**, and is
genuinely load-bearing (`F_n = 34.6 kN`, σ/σ_law = 2.205, 20.8 mm penetration).
Note the trade: tet4 reaches the G7 **sink** but with more penetration than UDG's
best, giving xi = 0.744 (not 0.85).

---

## 7. Conclusion / recipe

**Plainly:** **UDG (`udg-hex`) + *sliding* contact does work in FEBio 4.13 — but
only if the augmented Lagrangian is dropped and a large `search_radius` is used,
and even then it CANNOT complete the G7 sink.** It runs far past contact closure
(the original failure) to `t ≈ 0.97` (xi ≈ 0.71), then NANs; any contact stiff
enough to suppress the accompanying 23 mm of penetration NANs at closure. **For a
complete, load-bearing G7 deck, use the tested `tet4` fallback below.**

### Working recipe (UDG + sliding, past closure, to xi ≈ 0.71 — NOT G7-complete)

Winning config, exact contact XML actually emitted and used:

```xml
<contact name="indenter_pad" surface_pair="indenter_pad_pair" type="sliding-elastic">
  <laugon>PENALTY</laugon>          <!-- NO augmented Lagrangian -->
  <penalty>0.1</penalty>            <!-- soft; >=0.15 NANs at closure -->
  <auto_penalty>0</auto_penalty>
  <two_pass>0</two_pass>
  <node_reloc>0</node_reloc>        <!-- irrelevant (0 and 1 identical) -->
  <symmetric_stiffness>0</symmetric_stiffness>
  <tolerance>0.005</tolerance>
  <search_radius>20</search_radius> <!-- MUST be large; 2 loses the pair at large sink -->
  <fric_coeff>0.6</fric_coeff>
</contact>
```

plus the proven time-stepper (`time_steps=2400`, `<step_size>=<dtmax>=1/2400`,
`cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`) and
the uniaxial UDG recipe for the pad material/BC scheme.

* **What was decisive:** (i) `laugon="AUGLAG"` → `"PENALTY"`; (ii) `search_radius`
  **large**, not small; (iii) `penalty` **soft** (0.1).
* **What does NOT matter:** `hg` (0.01–100 identical), `node_reloc` (0/1 identical),
  `search_radius` beyond 20 (20/30/50/100 identical), sink (closure-NAN scales with
  `0.5/sink`), the hardening recipe, time-step fineness.
* **Measured at the wall** (`g7_pen0.1_sr20_s170`, reproducible): `end_t = 0.974227`,
  indentation **141.8 mm** (xi **0.7092**), penetration **23.27 mm**, `F_n = 23 150 N`,
  **σ_contact / σ_1D-law = 1.705**, `negJac = 0` throughout, NAN at the last step.

### Tested fallback for the G7 deliverable — `tet4` pad

Use the **default `tet4` pad** (`pad_mesh.pad_block_tet`) — no `type=` on its
`<SolidDomain>` — with the **same pure-penalty recipe**. It **completes G7**,
reproducibly, load-bearing:

| quantity (`n7_tet4_s170`, ×4 identical runs) | value |
|---|---|
| termination | `rc=0`, NORMAL, `end_t=1`, 2400/2400 steps, `negJac=0`, no NAN |
| indentation | 148.7 mm (xi = 0.7437) |
| penetration | 20.76 mm |
| contact force `F_n` | 34 600 N |
| **force balance** (indenter-BC vs pad-top reaction, Newton-3rd identity) | **7.80e-07** (≈ machine zero) |
| σ_contact | 0.5767 MPa |
| σ_1D-law(xi) | 0.2616 MPa |
| **σ_contact / σ_1D-law** | **2.205** |
| friction check | `F_t ≤ μ F_n` (t=0 row aside; the 1.1e14 ratio is the zero-force initial state) |

### What this means for the S5 track

* The historical framing "**UDG cannot do sliding contact**" is **too strong** — it
  *can*, once AUGLAG is removed and `search_radius` is large; the original failure
  was the **augmented-Lagrangian force spike at closure into the 1-point element**.
* But the exact G7 target (sink 170 mm **and** xi≈0.85) is **not reachable with
  `udg-hex`** in FEBio 4.13: the soft contact needed to survive closure leaves
  ≥23 mm penetration → xi ≤ ~0.71; stiffening it NANs.
* **For the 2 m-landing numbers, quote the `tet4` fallback** (complete G7, xi=0.744,
  σ/σ_law = 2.205, reproducible). The UDG result stands as the large-strain element
  proof but is **not** the G7 production deck.
