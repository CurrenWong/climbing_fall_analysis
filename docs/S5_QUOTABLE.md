# S5 — What may be quoted (authoritative verdict)

**Status:** COMPLETE (2026-10-06) — the single entry point for any reader who needs a number out of S5-contact. **Start here.** This doc is the authoritative reconciliation of the **six** G7 verification studies against the older "1.56/1.59/1.62/1.64" and "2.205" recommendations that still appear in the downstream docs. **§8 (2026-10-07) adds the wave5 bone-integrated contact results — a DIFFERENT quantity set (real THUMS mesh + a modeling hold spring); never mix a §8 number with a §0–§7 number.**

**The contradiction it resolves.** Four independent verification studies (`S5_g7_verify_anchors.md`, `S5_g7_verify_mesh.md`, `S5_g7_verify_fixedxi.md`, `S5_g7_verify_anchors_mesh.md`) have now returned FAIL for every σ÷σ_law ratio. The fifth (`S5_g7_verify_element.md`) returns PASS only on its narrow matched-ξ question; its §5 recommendation to "quote 2.205" is itself superseded here. The sixth (`S5_g7_verify_metric.md`) tested whether an alternative stress metric would converge and returned **FAIL** — closing the last open thread (§5). Yet two docs still recommend quoting 2.205, and several downstream docs still cite 1.56 / 1.59 / 1.62 / 1.64 / 2.205 / 2.70 / 15.14 with no recipe or mesh qualifier. **Anyone reading the repo today can walk out with a wrong number.** This doc fixes that.

---

## 0. Ten-second scan — the verdict table

| Quantity | Value | QUOTABLE? | Evidence |
|---|---|:---:|---|
| `F_n` at 2 m landing (sink 170, NEW PENALTY recipe, mesh-converged at fixed sink) | **≈ 34.4 kN** (34 404 / 34 602 / 34 395 N across (3,3,2)/(6,6,4)/(12,12,8) — drift **−0.026 %**) | **✅ QUOTABLE** | `S5_g7_verify_mesh.md` §3; reproducible ×4 |
| `σ_contact = F_n / 60 000 mm²` (same basis) | **≈ 0.574 MPa** (0.5734 / 0.5767 / 0.5733 MPa — drift **−0.026 %**) | **✅ QUOTABLE** | `S5_g7_verify_mesh.md` §3; same provenance |
| `σ_1D-law(ξ)` | a `pad_stress_pa(ξ)` material function (Pa → /1e6 = MPa) | **✅ QUOTABLE** | mesh-independent by construction; `pad_foam.pad_stress_pa` |
| Element agreement `udg-hex` vs `tet4` at matched ξ ≤ 0.65 | **≤ ±2.5 %** in σ/σ_law | **✅ QUOTABLE** (method validation only) | `S5_g7_verify_element.md` §3.3 / §5 |
| OLD-recipe anchor reproducibility | **≤ 0.1 %** vs the historical 1.56/1.59/1.62/1.64 (AUGLAG pen=1, n=1200, default `<Control>`) | **✅ QUOTABLE** (credibility of the harness under the OLD recipe; label the recipe) | `S5_g7_verify_anchors.md` §4.3 |
| NEW-recipe anchors on coarse (3,3,2) tet4, sink 20/30/40/50 mm | 1.4670 / 1.4767 / 1.4818 / 1.4863 | ⚠️ **only as a coarse-mesh reproduction** (cite the mesh and the recipe together) | `S5_g7_verify_anchors.md` §3 |
| NEW-recipe small-strain anchors as a mesh-converged σ/σ_law curve | any single number, or the bare coarse curve | **❌ NOT-QUOTABLE** | `S5_g7_verify_anchors_mesh.md` §4 (FAIL, drift −0.75 %…−5.47 %, refinement step growing/non-monotone) |
| Any σ/σ_law at sink=170 (G7), single mesh | `2.205` (tet4 (6,6,4)) | **❌ NOT-QUOTABLE as "converged"** | `S5_g7_verify_mesh.md` §3 (FAIL, drift +30.62 % with growing step) |
| Any σ/σ_law at sink=170 (G7), fixed ξ ≈ 0.70 across meshes | 1.752 / 2.027 / 2.508 across coarse/medium/fine | **❌ NOT-QUOTABLE** | `S5_g7_verify_fixedxi.md` §4–§5 (FAIL, drift +43.13 %, worse than free-ξ) |
| σ/σ_law at matched ξ = 0.7092 across elements | tet4 2.0573 vs udg-hex 1.7047 (pre-NAN) | **❌ NOT-QUOTABLE as a G7 number** | `S5_g7_verify_element.md` §3.4 / §5 (only honest at ξ ≤ 0.65, then ≤±2.5 %) |
| Cross-recipe comparison (OLD AUGLAG 1.56–1.64 vs NEW PENALTY 1.47–1.49 vs G7 2.205) | the comparison itself | **❌ NOT-QUOTABLE** (recipes are different curves: NEW − OLD = **−5.85 % … −9.38 %** monotone, ~60–100× the ≤0.1 % noise floor; OLD-on-NEW-stepper dies at closure) | `S5_g7_verify_anchors.md` §6 |
| `15.14 cm` / `2.70×` (large-strain "2.70× = flat-punch physics" claim) | — | **❌ NOT-QUOTABLE — retracted** | `S5_contact_convergence.md` §3.2 / `S5_contact_contract.md` §G7; outside the convergence envelope |
| `F_n` as mesh-converged **at small strain (sink ≤ 50 mm)** | — | **❌ NOT-QUOTABLE as converged** (drift −8.4 … −12.5 % across the 3 meshes) | `S5_g7_verify_anchors_mesh.md` §3 |
| Element formulation driving the 1.705 vs 2.205 gap | — | **❌ NOT-QUOTABLE as "element effect"** (≈⅓ ξ-mismatch + ≈⅔ UDG terminal degradation, not formulation) | `S5_g7_verify_element.md` §4 / §5 |

---

## 1. Evidence — the five studies, file by file

The five verification docs are the evidence base; this section points to the exact section/line, not retells them.

### 1.1 `docs\S5_g7_verify_anchors.md` — verdict **FAIL** (§6)

Cross-recipe shift at the anchor mesh `(3,3,2)` tet4. The OLD-recipe (`AUGLAG pen=1`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20`, `n=1200`) anchors reproduce the historical 1.56 / 1.59 / 1.62 / 1.64 to **≤ 0.1 %** (§4.3). The NEW recipe (`PENALTY 0.1`, `tolerance=0.005`, `node_reloc=0`, `search_radius=20`, `n=2400` + hardened `<Control>`) on the same mesh gives **1.4670 / 1.4767 / 1.4818 / 1.4863** (§3) — a **−5.85 % … −9.38 % monotone cross-recipe shift** (§4.1). `F_n` shifts **−6.72 % … −10.45 %** (§4.2). Penetration sign also flips (OLD: sub-mm micro-gap; NEW: +0.08…+0.62 mm). The OLD contact recipe cannot run on the NEW stepper (dies at closure, `negJac=55`, §5) — the recipes are not interchangeable.

**Consequence.** The two curves are not the same; quoting 2.205 alongside the 1.56–1.64 anchors is a **cross-recipe comparison** and is not valid as-is.

### 1.2 `docs\S5_g7_verify_mesh.md` — verdict **FAIL** (§3)

Free-ξ sweep at sink = 170 mm (NEW recipe) across pad meshes `(3,3,2)` / `(6,6,4)` / `(12,12,8)`. `F_n` and `σ_contact` are mesh-converged at **< 0.6 %** across a 64× element-count change (§3 table). The σ/σ_law ratio drifts **+12.00 % → +16.62 %**, net **+30.62 %** (1.968 → 2.205 → 2.571) and **the step grows with refinement** — it is not approaching a limit. Indentation drifts **−8.32 %**; penetration **+82.46 %**; σ_1D-law(ξ) **−23.46 %**. `negJac` go `11 → 0 → 51` — refinement does not stabilise the element response.

**Consequence.** The single-mesh G7 number `2.205` **must not** be quoted as a mesh-converged result. Only the resultant `F_n ≈ 34.4 kN` is mesh-converged at fixed sink.

### 1.3 `docs\S5_g7_verify_fixedxi.md` — verdict **FAIL** (§5)

The decisive experiment that separates "ξ artefact" (H1) from "intrinsic mesh dependence" (H2). Holding ξ ≈ 0.70 fixed per mesh (sink chosen per mesh: 152.15 / 159.0 / 168.0 mm), the ratio drifts **+15.65 % → +23.76 %**, net **+43.13 %** (§4.2) — *more* than the +30.62 % free-ξ sweep, with the step *growing*. The σ_1D-law drift collapses from −23.5 % → −0.6 % ✅ (so H1 was partly right for the denominator), but `F_n` now drifts **+42.3 %** (so the soft-penalty penetration moved the mesh dependence from the denominator into the numerator). **Conclusion:** `F_n` and ξ are **never simultaneously converged**; the ratio is intrinsically mesh-dependent.

**Consequence.** No single quotable σ/σ_law with ξ and recipe can be given. Quote `F_n` (fixed sink) or σ_contact — not the ratio.

### 1.4 `docs\S5_g7_verify_anchors_mesh.md` — verdict **FAIL** (§4)

The NEW-recipe small-strain anchors (1.4670 / 1.4767 / 1.4818 / 1.4863 from `S5_g7_verify_anchors.md` §3) re-run on `(3,3,2)` / `(6,6,4)` / `(12,12,8)` at sink 20/30/40/50 mm — 12 runs, all `rc=0 NORMAL end_t=1 no NAN`. Per-sink σ/σ_law drift coarse→fine is **−5.47 % / −4.64 % / −3.03 % / −0.75 %** with a **non-monotone, growing refinement step** (`|c→m| = 0.2–1.7 %`, `|m→f| = 2.4–5.3 %`); F_n also drifts **−8.4 … −12.5 %** at small strain (contrast G7, where F_n was converged to <0.6 %). Mechanism is identical to G7 (penetration grows ×9–×24 with refinement, ξ drops, σ_law(ξ) drops), just weaker.

**Consequence.** The NEW-recipe anchors are mesh-sensitive by construction; no single number may be quoted as converged. The maximum defensible statement is the **mesh-bounded band per anchor** (§5, table).

### 1.5 `docs\S5_g7_verify_element.md` — verdict **PASS on its own question**, but its §5 recommendation is **superseded** here

At matched ξ ≤ 0.65 the `udg-hex` (HEX8G1) and `tet4` elements agree to **≤ ±2.5 %** in σ/σ_law on the same mesh geometry and recipe (§3.3 / §5). That is a valid method validation — the element formulation itself is *not* the driver of the gap in the load-bearing range. The headline 1.705-vs-2.205 gap is **≈⅓ ξ-mismatch (29.3 % → 20.7 %) + ≈⅔ UDG terminal degradation** (UDG's 1.7047 is its pre-NAN last state, dipped off its own smooth curve from 1.874 at ξ=0.65) — **not** an element discrepancy (§4).

**However.** Its §5 recommendation to "quote the `tet4` completion: σ_contact/σ_1D-law = 2.205 @ ξ = 0.744" is **superseded** by `S5_g7_verify_mesh.md` §3 (2.205 is not mesh-converged at the G7 sink) and `S5_g7_verify_fixedxi.md` §5 (no single quotable σ/σ_law with ξ and recipe can be given). The §5 lines **260–263** carry the supersession banner in-place; the measurements in that file (the matched-ξ table, the element-formulation A/B, the per-ξ run blocks) remain valid as measurements and may be quoted for the matched-ξ question itself.

---

## 2. Root cause (one paragraph)

The σ/σ_law ratio is **mesh-broken, not solver-broken**. The cause is the recipe: `penalty=0.1` (soft) is forced by the contact-closure constraint — it is the only `penalty` value that survives the FOAM-on-FOAM closure at sink 170 mm without NANing (`S5_contact_udg.md` §7; `S5_g7_verify_anchors.md` §5 shows the OLD recipe dies on the NEW stepper). Soft penalty ⇒ mesh-dependent penetration (verified ×24 at anchor sink 20 mm, ×9 at sink 50 mm, ×1.8 at sink 170 mm across the (3,3,2)→(12,12,8) sweep). With mesh-dependent penetration, "fixed travel" and "fixed indentation" states cannot be aligned: fixing the sink pins `F_n` (the global equilibrium resultant) but lets ξ drift; fixing ξ pins the pad-top indentation but lets `F_n` drift (§1.3 §4.3). On top of that, the flat-punch **edge is a stress singularity**, so any footprint-averaged σ_contact is mesh-dependent (a mesh artefact of the singular traction profile at the punch rim, not a resolved field). **Therefore no σ/σ_law ratio is quotable as a converged scalar at any sink in this mesh family**; the metric is broken, not the solver, the material, or the harness.

---

## 3. What the paper may now quote (precise wording)

The numbers below are the only ones currently defensible for citation, with the explicit recipe/mesh/ξ stated.

1. **`F_n ≈ 34.4 kN`** at the 2 m landing (NEW PENALTY recipe, tet4 pad, sink 170 mm, indenter `(2,2,1,(200,300,5))` offset `(50,0,200.5)`, gap 0.5 mm, μ = 0.6). Mesh-converged to **< 0.6 %** across `(3,3,2)`→`(6,6,4)`→`(12,12,8)` (a 64× element-count change); reproducible ×4; force-balance ≈ 1e-7 (Newton-3rd identity check holds). Source: `S5_g7_verify_mesh.md` §2.1 / §3.

2. **`σ_contact = F_n / 60 000 mm² ≈ 0.574 MPa`** at the 2 m landing — same basis, same convergence, same reproducibility as `F_n`. Source: `S5_g7_verify_mesh.md` §3.

3. **`σ_1D-law(ξ)`** is `pad_stress_pa(ξ)/1e6` (MPa). Mesh-independent by construction; it is a fitted material function, not an FE result. Source: `src/climbing/coupling/pad_foam.py::pad_stress_pa`.

4. **Element agreement `≤ ±2.5 %`** between `udg-hex` (HEX8G1) and `tet4` **at matched ξ ≤ 0.65**, same recipe, same mesh geometry. This is a **method validation** (it shows the formulation choice does not drive the gap in the load-bearing range), not a G7 result. Source: `S5_g7_verify_element.md` §3.3 / §5.

5. **OLD-recipe anchor reproducibility `≤ 0.1 %`** (AUGLAG pen=1, n=1200, default `<Control>`, tet4 (3,3,2)). This establishes the **numerical noise floor of the harness** under the OLD recipe. Always cite the recipe together with the number. Source: `S5_g7_verify_anchors.md` §4.3.

6. **NEW-recipe coarse-mesh anchors 1.4670 / 1.4767 / 1.4818 / 1.4863 at sink 20/30/40/50 mm** — only as a **coarse-mesh reproduction of the diagnostic**; the paper may cite these together with the mesh-bounded band from `S5_g7_verify_anchors_mesh.md` §5 (1.39–1.47 / 1.41–1.48 / 1.44–1.49 / 1.47–1.51 across the 3-mesh sweep) and the explicit "not mesh-converged, refinement step non-decreasing" qualifier.

7. **Flat-punch enhancement `≈ 1.6`** in the trusted small-strain range `ξ ≤ 0.25` (sink ≤ 50 mm, OLD AUGLAG recipe, penetration ≈ 0). This is the `1.56 → 1.64` smooth rise. Source: `S5_contact_convergence.md` §3.1 (verbatim transcription in `S5_g7_verify_anchors.md` §1.1).

---

## 4. What the paper must NOT quote (with reasons)

1. **Any σ/σ_law ratio as a "converged" G7 number.** Four independent FAILs: free-ξ drift +30.62 % with growing step (`S5_g7_verify_mesh.md` §3); fixed-ξ drift +43.13 % with growing step (`S5_g7_verify_fixedxi.md` §4); small-strain drift −0.75 %…−5.47 % with non-monotone step (`S5_g7_verify_anchors_mesh.md` §3); cross-recipe shift −5.85 %…−9.38 % (`S5_g7_verify_anchors.md` §4). Including the `2.205` itself, including any single fixed-ξ pick (1.752 / 2.027 / 2.508), including any matched-ξ pick at ξ > 0.65.

2. **`2.205` as "the converged G7 σ/σ_law."** It is a single-mesh measurement; the sweep moves +30.6 % away from it on refinement, with the step growing. It is **the most misleading possible single number** to publish.

3. **Cross-recipe comparison of the historical `1.56/1.59/1.62/1.64` anchors with any G7 number.** The OLD AUGLAG pen=1 curve and the NEW PENALTY 0.1 curve are different curves, separated by **−5.85 % … −9.38 %** monotone (60–100× the ≤0.1 % noise floor). They are not the same delivery and cannot be plotted on the same y-axis without the recipe label.

4. **`15.14 cm` and `2.70×`** as a large-strain "flat-punch physics" number. Retracted (`S5_contact_convergence.md` §3.2; `S5_contact_contract.md` §G7). The recipe state that produced it (`PENALTY 0.1` at sink 170 mm) sits outside the convergence envelope — three identical runs gave three different negJac counts and all failed. The 2.70× figure came from a single lucky run, and the framing "2.70× = large-strain flat-punch physics" is therefore wrong on two counts: (a) the value is not reproducible, and (b) even if it were, the value is mesh-dependent in the measured range.

5. **`F_n` as "mesh-converged" at small strain** (sink ≤ 50 mm). It is *not* at small strain — it drifts −8.4 % … −12.5 % across the 3-mesh sweep (`S5_g7_verify_anchors_mesh.md` §3). It *is* converged at G7 (sink 170, <0.6 %, `S5_g7_verify_mesh.md` §3). The two statements are not contradictory — small-strain response is dominated by local flat-punch contact compliance (mesh-sensitive under soft penalty); G7 is dominated by the bulk material response (mesh-insensitive) — but the paper must distinguish them, not let a "mesh-converged F_n" claim leak between regimes.

7. **The 1.705 vs 2.205 pair as an "element-formulation difference."** It is **≈⅓ ξ-mismatch (29.3 % → 20.7 %) + ≈⅔ UDG terminal degradation** (`S5_g7_verify_element.md` §4). Cite at matched ξ ≤ 0.65 if a matched-ξ number is wanted; never cite 1.705 next to 2.205 as if they were the same point.

---

## 5. Alternative metrics — `docs\S5_g7_verify_metric.md` — **COMPLETE, verdict FAIL**

The sibling study asked whether a *different* stress metric would be mesh-independent: average the contact traction over the footprint **interior only**, excluding a margin `m` from the punch boundary (footprint `[50,250]×[0,300]` mm, 60 000 mm²). **It does not converge either.**

* Pre-stated grid: **5 margins (m = 0 / 5 / 10 / 20 / 30 mm) × 3 pad meshes** at sink 170 under the same recipe, ≤3 % PASS tolerance. All three runs reproduce `S5_g7_verify_mesh.md` to 3–4 sig figs — the harness self-check passed before any conclusion was drawn.
* **σ_interior drifts +25.9 % at m = 0 → +32.3 % at m = 20**, and the successive refinement step **grows** (+9.5 % → +14.4 %) — FAIL by ~8–10× the tolerance. The interior **ratio** drifts +64 % → +73 %, i.e. **worse** than the full-footprint average (+30.62 %) and worse than fixed-ξ (+43.13 %).
* **Mechanism — this corrects the flat-punch-edge-singularity hypothesis** proposed in §2: the dominant contamination is *not* the punch rim but the **soft-penalty pressure spill**. Under `penalty=0.1` the contact pressure spreads beyond the rigid footprint by roughly one cell width; coarse meshes spread far, fine meshes barely — so *removing* the rim makes the interior average grow with refinement instead of settling.
* **Fallback also tested:** the mid-depth full-cross-section traction `σ_mid = F_n / 90 000 mm²` is mesh-converged at **−0.13 %** — but it is `F_n` rescaled, so it carries no enhancement information and cannot serve as a σ÷σ_law substitute.

**Consequence.** The §0 verdict table stands unchanged. **No footprint-based stress metric converges at this penalty stiffness**; the paper must quote absolute quantities only. The one open thread this doc flagged is now closed, and nothing moves from NOT-QUOTABLE to QUOTABLE.

---

## 6. Pointer index — which doc proves which claim

| Claim in this verdict | Doc + section |
|---|---|
| F_n ≈ 34.4 kN, σ_contact ≈ 0.574 MPa, mesh-converged <0.6 % at G7 | `S5_g7_verify_mesh.md` §2.1 / §3 |
| σ/σ_law at G7 drifts +30.62 % with growing step | `S5_g7_verify_mesh.md` §3 |
| σ/σ_law at fixed ξ ≈ 0.70 drifts +43.13 % (worse than free-ξ) | `S5_g7_verify_fixedxi.md` §4.2 / §5 |
| NEW-recipe small-strain anchors drift −0.75 %…−5.47 % with non-monotone step | `S5_g7_verify_anchors_mesh.md` §3 / §4 |
| Mesh-bounded band per anchor (the only honest single anchor) | `S5_g7_verify_anchors_mesh.md` §5 |
| Cross-recipe shift −5.85 %…−9.38 %, OLD on NEW-stepper dies | `S5_g7_verify_anchors.md` §4.1 / §5 / §6 |
| OLD-recipe reproducibility ≤0.1 % (noise floor) | `S5_g7_verify_anchors.md` §4.3 |
| Element agreement ≤ ±2.5 % at ξ ≤ 0.65 (matched-ξ) | `S5_g7_verify_element.md` §3.3 / §5 |
| 1.705-vs-2.205 ≈⅓ ξ + ≈⅔ UDG terminal degradation (not formulation) | `S5_g7_verify_element.md` §4 |
| 15.14 cm / 2.70× retracted | `S5_contact_convergence.md` §3.2 / `S5_contact_contract.md` §G7 |
| Soft-penalty mechanism (penalty=0.1 forced by closure; penetration mesh-dependent) | `S5_contact_udg.md` §7 + `S5_g7_verify_anchors.md` §5 + `S5_g7_verify_fixedxi.md` §4.3 |
| σ/σ_law at sink=170 mesh-converged "as converged": **NO** (this is the headline of this doc) | this file §0–§4 |
| Interior-margin stress metric FAIL (+25.9 % → +32.3 %, step growing; cause is penalty spill, not the punch rim) | `S5_g7_verify_metric.md` §5 / §6 |
| `σ_mid = F_n/90 000 mm²` mesh-converged at −0.13 %, but it is `F_n` rescaled (no enhancement content) | `S5_g7_verify_metric.md` §6 |
| No footprint-based stress metric converges at `penalty=0.1` | `S5_g7_verify_metric.md` §6 (final FAIL) |

---

## 7. Reproduction pointers

All five studies have on-disk ledgers and per-run `.feb` / `.log` / `.xplt`. No new experiment was performed to author this verdict — it is a reconciliation of the five prior studies against the existing downstream text.

```powershell
$env:PYTHONPATH="src"
# the five verification studies
ls docs\S5_g7_verify_anchors.md       # §1 OLD/§3 NEW at anchor mesh, §4 A/B, §5 contact-vs-stepper isolation, §6 FAIL
ls docs\S5_g7_verify_mesh.md          # §2 sweep at sink 170, §3 FAIL
ls docs\S5_g7_verify_fixedxi.md       # §3 fixed-ξ grid, §4.2 drift table, §5 FAIL
ls docs\S5_g7_verify_anchors_mesh.md  # §2.1 3×4 grid, §4 verdict, §5 mesh-bounded band
ls docs\S5_g7_verify_element.md       # §3 matched-ξ, §4 decomposition, §5 PASS-on-question / §5-recommendation-superseded
ls docs\S5_g7_verify_metric.md        # §4 margin×mesh matrix, §5 drift vs margin, §6 FAIL (penalty spill)

# existing downstream text this doc reconciles
ls docs\S5_contact_udg.md             # §7 tet4 fallback (supersession banner added)
ls docs\S5_contact_convergence.md     # §3.1 trusted anchors, §3.2 G7 failure
ls docs\S5_contact_contract.md        # §G7 G7' correction

# the cross-recipe isolation in temp\ (the recipe/stepper cross at sink 50)
ls temp\opensim_fe\s5_t9\g7verify_a\  # anc_old_s50, anc_old_s20/30/40, anc_new_s20/30/40/50, anc_old_s50_newstep, anc_new_s50_oldstep
ls temp\opensim_fe\s5_t9\g7verify_b\  # sink=170 sweep (3,3,2)/(6,6,4)/(12,12,8)
ls temp\opensim_fe\s5_t9\fixb\        # fixed-ξ grid
ls temp\opensim_fe\s5_t9\fixa\        # NEW-anchor mesh x-grid
ls temp\opensim_fe\s5_t9\g7verify_c\  # matched-ξ element runs
ls temp\opensim_fe\s5_t9\fixmetric\   # interior-margin metric (3 runs @ sink 170, 5 margins each)
```

---

## 8. Wave5 — bone-integrated contact on the real THUMS calcaneus (2026-10-07)

**§0–§7 are the standalone-pad G7 verdict** (tet4 pad, prescribed sink, NEW PENALTY). **§8 is a different quantity set**: the **bone-integrated** model (real THUMS `calcaneus_r_anatframe.npz`, 137 plantar quad4) driven by the subtalar-joint load, and it carries a **modeling hold spring**. **Never place a §8 number next to a §0–§7 number.** Full report: `results/opensim_fe/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md`; plan/log: `docs/S5_wave5_matrix.md`.

### 8.1 Quotable ✅

| Quantity | Value | Condition | Evidence |
|---|---|---|---|
| **Bone-integrated contact is usable (method)** | real mesh 10/10 contact rows ran; **8 CONVERGED+CARRYING** (`rc=0`, `end_t=1.0`); wave4 was **0/7** | recipe = plantar **quad4** + G7 `<Control>` (2400 steps, `dtmax`=1÷2400, `cutback`=0.125, `max_retries`=20, BFGS) + opt-in hold spring + **`use_rigid=True`** | `..._WAVE5_REPORT.md` §0/§2 |
| **Load-path finding (method)** | `PressureLoad` (`use_rigid=False`) diverges at contact closure **t≈0.419**; rigid path converges | same mesh/recipe | §1/§2 |
| **Hold spring required (method)** | remove hold ⇒ **DIVERGED** (2–3 s, pre-closure rigid modes) | — | §3 |
| **Bone-integrated contact force `F_n` (pen100)** | **≈ 20.4 kN** @ h = 2.0 m (applied 21.086 kN; hold spring contributes only **3.4 %**; penetration **2.19 mm**) | `penalty=100`, `use_rigid=True`, hold 1 N/mm, μ=0.6, anatframe 137q | §2 (`h2.0_contact_rigid_mu0.6_pad_pen100`) |
| **Penetration magnitude (cite penalty with it)** | pen0.1 → **40–48 mm** (h=2.0–4.5); pen100 → **≈ 2.2 mm** | penalty label mandatory | §2/§4 |

### 8.2 NOT quotable ❌

| Quantity | Reason |
|---|---|
| **pen0.1 contact `F_n`** (e.g. 13.8 kN @ h2.0) as a physical force | it is `applied − hold-spring force`, and the hold spring (1 N/mm/node × 175 = **175 N/mm**) is an **un-calibrated numerical device**; at pen0.1 it carries **7310/21086 = 35 %** ⇒ `F_n` depends on the arbitrary spring stiffness. Near-robust only at pen100 (3.4 %). |
| **bone-integrated `F_n` next to the standalone-pad `F_n ≈ 34.4 kN`** | different mesh, different load path (subtalar load vs prescribed sink), and the integrated model includes a modeling spring |
| **any wave5 σ/σ_law or footprint-averaged stress** | same soft-penalty spill as §2; wave5 did no mesh-convergence study and resolved no contact-pressure field |
| **bone-integrated `contact_area`** (e.g. 2953.8 mm² at pen0.1) as a contact patch | it is a geometric count ("triangles with all 3 nodes below the plate plane"), not the contact-pressure patch |