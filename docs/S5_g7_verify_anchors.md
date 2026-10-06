# S5 — G7 anchor-continuity check: does the NEW contact recipe sit on the SAME curve as the trusted anchors?

**Status:** COMPLETE (2026-10-06) — **verdict: FAIL** (§6). Written incrementally — one
block per run, appended to disk the moment each run finished (no result ever held only in
memory). Every value read from the real `.log` tail / `.xplt` read-back.

**Question.** The trusted historical small-strain anchors
(FE σ / 1D-law σ ≈ **1.56 / 1.59 / 1.62 / 1.64** at sink **20/30/40/50 mm**) were produced
with the **OLD** contact recipe (`laugon="AUGLAG"`, `penalty=1.0`, `tolerance=0.005`,
`aug_controls=AUGLAG_CONVERGENT_CONTROLS`). The **NEW** G7 result uses a **different**
recipe (`laugon="PENALTY"`, `penalty=0.1`, `node_reloc=0`, `search_radius=20`).
**A number is only quotable alongside the anchors if it lies on the SAME curve.** This doc
decides whether the recipe change moves the curve.

Units: mm–N–MPa–s. FEBio binary: `D:\Program\FEBioStudio\bin\febio4.exe` (v4.13.0).

Reproducibility dir (NEW, collides with nobody): `temp\opensim_fe\s5_t9\g7verify_a\`.

---

## 1. Historical anchors — EXACT transcription

Transcribed verbatim from `docs\S5_contact_convergence.md`.

### 1.1 The historical anchor table (§3.1 of `S5_contact_convergence.md`)

> ### 3.1 可靠区间扫描（AUGLAG pen=1 + 配方，n=1200）
>
> | sink (mm) | ξ | F_n (N) | **FE/law** | 穿透 (mm) | 结果 |
> |---|---|---|---|---|---|
> | 20 | 0.098 | 1589 | 1.56 | −0.164 | ✅ |
> | 30 | 0.149 | 2163 | 1.59 | −0.228 | ✅ |
> | 40 | 0.199 | 2700 | 1.62 | −0.235 | ✅ |
> | 50 | 0.248 | 3232 | 1.64 | −0.086 | ✅ |
> | 60 | — | — | — | — | ❌ negJac=53 |
> | 70 | — | — | — | — | ❌ negJac=20 |
>
> 穿透为负 = 接触面有**亚毫米微间隙**（健康）。FE/law 从 1.56 平滑升到 1.64 ⇒ 这是**真实的平冲头
> 几何效应**，但**量级只有 ~1.6，不是 2.7**。

(Section heading quoted verbatim: `### 3.1 可靠区间扫描（AUGLAG pen=1 + 配方，n=1200）`.)

### 1.2 The OLD recipe (§2 of `S5_contact_convergence.md`)

> ```python
> build_pad_contact_feb(..., laugon="AUGLAG", penalty=1.0,
>                       tolerance=0.005,
>                       aug_controls=AUGLAG_CONVERGENT_CONTROLS)  # gaptol=0.001, minaug=3, maxaug=200, smooth_aug=1
> ```

`AUGLAG_CONVERGENT_CONTROLS` (from `src/climbing/coupling/pad_contact.py:208`):
`gaptol=0.001, minaug=3, maxaug=200, smooth_aug=1`.

### 1.3 Which element / mesh / steps produced the anchors

The §3.1 table is the output of `temp\opensim_fe\s5_t9\envelope.py`
(evidence: §3.1 lists `sink=40` and the `sink=60/70 → negJac=53/20` failures; those
exact negJac values appear only in `envelope\env_sink60.log` / `env_sink70.log`, and
`envelope.py` is the only 20–70 sweep). `envelope.py` fixes:

| parameter | value (from `envelope.py`) |
|---|---|
| **pad element** | **tet4** (`pad_mesh.pad_block_tet`) — no `type=` on the pad `<SolidDomain>` |
| **pad mesh** | **`pad_block_tet(3, 3, 2, (300, 300, 200))`** → 100 mm nominal elements |
| indenter | `indenter_block_tet(2, 2, 1, (200, 300, 5), offset_mm=(50, 0, 200.5))` |
| gap | 0.5 mm (indenter bottom z = 200.5, pad top z = 200.0) |
| μ | 0.6 |
| **contact** | `laugon="AUGLAG"`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS` |
| `node_reloc` | `1` (library default) |
| `search_radius` | `20.0` (library default) |
| **time steps** | **`n_solver_steps=1200`** |
| load | `{"mode": "displacement", "sinkage_mm": sink, "lateral_mm": (0,0)}` |

`sink=20/30/50` were also produced by `verify_lib.py` with the identical mesh/recipe
(`verify_lib.py:28-36,49-55`), so the whole 20–50 band has a single, consistent provenance.
The `finish\A_sink30_n1200.feb` deck is the same writer too.

### 1.4 The NEW recipe being tested (from `docs\S5_contact_udg.md` §7)

> ```xml
> <contact name="indenter_pad" surface_pair="indenter_pad_pair" type="sliding-elastic">
>   <laugon>PENALTY</laugon>
>   <penalty>0.1</penalty>
>   <auto_penalty>0</auto_penalty>
>   <two_pass>0</two_pass>
>   <node_reloc>0</node_reloc>
>   <symmetric_stiffness>0</symmetric_stiffness>
>   <tolerance>0.005</tolerance>
>   <search_radius>20</search_radius>
>   <fric_coeff>0.6</fric_coeff>
> </contact>
> ```
> plus the proven time-stepper (`time_steps=2400`, `<step_size>=<dtmax>=1/2400`,
> `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`).

The G7 point to be placed on the curve (`S5_contact_udg.md` §7, tet4 fallback,
`pad_block_tet(6,6,4)`, sink 170): `end_t=1`, `F_n=34 600 N`, indentation 148.7 mm
(ξ=0.7437), penetration 20.76 mm, force balance 7.80e-07, **σ/σ_law = 2.205**.

> **Like-for-like rule.** The anchors are tet4 at **3×3×2**; the G7 point is tet4 at
> **6×6×4**. The recipe effect must therefore be measured at the **anchor mesh (3×3×2)**
> for the 4 anchor sinks — *not* by reading the G7 6×6×4 run.

---

## 2. Experiment log

_(appended one block per run, from the real `.log` tail / `.xplt` read-back.)_

Driver: `temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py`
(ledger: `temp\opensim_fe\s5_t9\g7verify_a\results.txt`).

### anc_new_s20 — recipe **NEW**, sink = 20 mm

recipe: `laugon=PENALTY`, `penalty=0.1`, `tolerance=0.005`, `aug_controls=None`, `node_reloc=0`, `search_radius=20`, `n_solver_steps=2400` + hardened `<Control>`
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.025377 |
| nconv | 2406 |
| negJac | 28 |
| NAN | no |
| indentation (mm) | 19.42 |
| xi | 0.0971 |
| penetration (mm) | 0.08226 |
| F_n (N) | 1482.26 |
| sigma_contact (MPa) | 0.0247043 |
| sigma_1D-law (MPa) | 0.0168396 |
| **sigma/sigma_law** | **1.4670** |
| force_balance_rel | 3.243e-07 |
| n_states | 2407 |

_(24s)_

---

### anc_new_s30 — recipe **NEW**, sink = 30 mm

recipe: `laugon=PENALTY`, `penalty=0.1`, `tolerance=0.005`, `aug_controls=None`, `node_reloc=0`, `search_radius=20`, `n_solver_steps=2400` + hardened `<Control>`
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2406 |
| negJac | 0 |
| NAN | no |
| indentation (mm) | 29.35 |
| xi | 0.1467 |
| penetration (mm) | 0.1512 |
| F_n (N) | 1997.49 |
| sigma_contact (MPa) | 0.0332915 |
| sigma_1D-law (MPa) | 0.0225439 |
| **sigma/sigma_law** | **1.4767** |
| force_balance_rel | 1.222e-07 |
| n_states | 2407 |

_(27s)_

---

### anc_new_s40 — recipe **NEW**, sink = 40 mm

recipe: `laugon=PENALTY`, `penalty=0.1`, `tolerance=0.005`, `aug_controls=None`, `node_reloc=0`, `search_radius=20`, `n_solver_steps=2400` + hardened `<Control>`
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0127778 |
| nconv | 2406 |
| negJac | 7 |
| NAN | no |
| indentation (mm) | 39.19 |
| xi | 0.1960 |
| penetration (mm) | 0.3069 |
| F_n (N) | 2444.06 |
| sigma_contact (MPa) | 0.0407343 |
| sigma_1D-law (MPa) | 0.0274892 |
| **sigma/sigma_law** | **1.4818** |
| force_balance_rel | 4.839e-08 |
| n_states | 2407 |

_(27s)_

---

### anc_new_s50 — recipe **NEW**, sink = 50 mm

recipe: `laugon=PENALTY`, `penalty=0.1`, `tolerance=0.005`, `aug_controls=None`, `node_reloc=0`, `search_radius=20`, `n_solver_steps=2400` + hardened `<Control>`
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0107533 |
| nconv | 2424 |
| negJac | 76 |
| NAN | no |
| indentation (mm) | 48.88 |
| xi | 0.2444 |
| penetration (mm) | 0.6241 |
| F_n (N) | 2894.07 |
| sigma_contact (MPa) | 0.0482346 |
| sigma_1D-law (MPa) | 0.0324529 |
| **sigma/sigma_law** | **1.4863** |
| force_balance_rel | 3.869e-07 |
| n_states | 2425 |

_(49s)_

---

### anc_old_s50 — recipe **OLD**, sink = 50 mm

recipe: `laugon=AUGLAG`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20`, `n_solver_steps=1200` (no hardening — verbatim historical)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 0.999996 |
| fail_t | 0.0108333 |
| nconv | 1478 |
| negJac | 32 |
| NAN | no |
| indentation (mm) | 49.59 |
| xi | 0.2479 |
| penetration (mm) | -0.0863 |
| F_n (N) | 3231.76 |
| sigma_contact (MPa) | 0.0538627 |
| sigma_1D-law (MPa) | 0.0328409 |
| **sigma/sigma_law** | **1.6401** |
| force_balance_rel | 7.673e-09 |
| n_states | 1479 |

_(69s)_

---

### anc_old_s20 — recipe **OLD**, sink = 20 mm

recipe: `laugon=AUGLAG`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20`, `n_solver_steps=1200` (no hardening — verbatim historical)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 0.999996 |
| fail_t | 0.0257575 |
| nconv | 1299 |
| negJac | 5 |
| NAN | no |
| indentation (mm) | 19.66 |
| xi | 0.0983 |
| penetration (mm) | -0.1637 |
| F_n (N) | 1588.97 |
| sigma_contact (MPa) | 0.0264829 |
| sigma_1D-law (MPa) | 0.0169966 |
| **sigma/sigma_law** | **1.5581** |
| force_balance_rel | 1.921e-08 |
| n_states | 1300 |

_(25s)_

---

### anc_old_s30 — recipe **OLD**, sink = 30 mm

recipe: `laugon=AUGLAG`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20`, `n_solver_steps=1200` (no hardening — verbatim historical)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 0.999996 |
| fail_t | 0.0171969 |
| nconv | 1321 |
| negJac | 3 |
| NAN | no |
| indentation (mm) | 29.73 |
| xi | 0.1486 |
| penetration (mm) | -0.2276 |
| F_n (N) | 2163.46 |
| sigma_contact (MPa) | 0.0360577 |
| sigma_1D-law (MPa) | 0.0227418 |
| **sigma/sigma_law** | **1.5855** |
| force_balance_rel | 6.7e-08 |
| n_states | 1322 |

_(20s)_

---

### anc_old_s40 — recipe **OLD**, sink = 40 mm

recipe: `laugon=AUGLAG`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20`, `n_solver_steps=1200` (no hardening — verbatim historical)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 0.999996 |
| fail_t | - |
| nconv | 1303 |
| negJac | 0 |
| NAN | no |
| indentation (mm) | 39.73 |
| xi | 0.1987 |
| penetration (mm) | -0.2349 |
| F_n (N) | 2698.5 |
| sigma_contact (MPa) | 0.0449751 |
| sigma_1D-law (MPa) | 0.0277573 |
| **sigma/sigma_law** | **1.6203** |
| force_balance_rel | 1.555e-08 |
| n_states | 1304 |

_(16s)_

---

### anc_new_s50_oldstep — contact **NEW** / stepper **OLD**, sink = 50 mm

recipe: contact=`laugon=PENALTY`, `penalty=0.1`, `tolerance=0.005`, `aug_controls=None`, `node_reloc=0`, `search_radius=20` · stepper=`n_solver_steps=1200`, default `<Control>` (verbatim historical)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 0.999996 |
| fail_t | - |
| nconv | 1200 |
| negJac | 0 |
| NAN | no |
| indentation (mm) | 48.88 |
| xi | 0.2444 |
| penetration (mm) | 0.6181 |
| F_n (N) | 2928.55 |
| sigma_contact (MPa) | 0.0488092 |
| sigma_1D-law (MPa) | 0.032456 |
| **sigma/sigma_law** | **1.5039** |
| force_balance_rel | 3.891e-07 |
| n_states | 1201 |

_(10s)_

---

### anc_old_s50_newstep — contact **OLD** / stepper **NEW**, sink = 50 mm

recipe: contact=`laugon=AUGLAG`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`, `node_reloc=1`, `search_radius=20` · stepper=`n_solver_steps=2400` + hardened `<Control>` (step_size=dtmax=1/2400, cutback 0.125, retries 20, BFGS max_ups 10)
mesh: `pad_block_tet(3,3,2,(300,300,200))` tet4 · indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · mu 0.6

| field | value |
|---|---|
| termination | ERROR_TERM |
| end_t | 0.010039 |
| fail_t | 0.010039 |
| nconv | 31 |
| negJac | 55 |
| NAN | no |
| indentation (mm) | 2.198 |
| xi | 0.0110 |
| penetration (mm) | -2.196 |
| F_n (N) | 1561.77 |
| sigma_contact (MPa) | 0.0260294 |
| sigma_1D-law (MPa) | 0.00247243 |
| **sigma/sigma_law** | **10.5279** |
| force_balance_rel | 0.8557 |
| n_states | 32 |

_(3s)_

---

## 3. The SAME 4 anchor sinks re-run under the NEW recipe

Anchor mesh (`pad_block_tet(3,3,2)` tet4) held fixed. NEW contact = `laugon=PENALTY`
`penalty=0.1` `tolerance=0.005` `node_reloc=0` `search_radius=20`; stepper = 2400 steps
+ hardened `<Control>`. σ_contact = F_n / 60 000 mm² (indenter top area); σ_1D-law =
`pad_stress_pa(xi)/1e6`.

| sink (mm) | ξ | indentation (mm) | penetration (mm) | F_n (N) | σ_contact (MPa) | σ_1D-law (MPa) | **σ/σ_law** | force-balance | end_t | negJac | NAN |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.0971 | 19.42 | +0.0823 | 1482.26 | 0.024704 | 0.016840 | **1.4670** | 3.24e-07 | 1 | 28 | no |
| 30 | 0.1467 | 29.35 | +0.1512 | 1997.49 | 0.033292 | 0.022545 | **1.4767** | 1.22e-07 | 1 | 0 | no |
| 40 | 0.1960 | 39.19 | +0.3069 | 2444.06 | 0.040734 | 0.027489 | **1.4818** | 4.84e-08 | 1 | 7 | no |
| 50 | 0.2444 | 48.88 | +0.6241 | 2894.07 | 0.048235 | 0.032452 | **1.4863** | 3.87e-07 | 1 | 76 | no |

All four NEW runs: `rc=0`, **NORMAL TERMINATION, `end_t=1`, 2400/2400 steps, no NAN**;
transient `negJac` events occurred but were recovered by cutback (final states healthy).
`force-balance` is ≈1e-7 (machine-zero) on every run. Penetration is **positive** under
the NEW soft penalty (pad surface lags the indenter by 0.08–0.62 mm), whereas the OLD
AUGLAG recipe gave a **sub-millimetre micro-gap** (negative values, §1.1).

## 4. Direct OLD-vs-NEW A/B at the anchor mesh

The OLD-recipe re-runs reproduce the historical anchors **to the digit** — proving the
harness is faithful and fixing the numerical noise floor at **≤0.1 %**.

### 4.1 σ/σ_law

| sink (mm) | historical quoted | OLD re-run | NEW | Δ (NEW − OLD) | **% change** |
|---|---|---|---|---|---|
| 20 | 1.56 | 1.5581 | 1.4670 | −0.0911 | **−5.85 %** |
| 30 | 1.59 | 1.5855 | 1.4767 | −0.1088 | **−6.86 %** |
| 40 | 1.62 | 1.6203 | 1.4818 | −0.1385 | **−8.55 %** |
| 50 | 1.64 | 1.6401 | 1.4863 | −0.1538 | **−9.38 %** |

### 4.2 F_n (N) and penetration (mm) — same A/B

| sink | historical F_n | OLD F_n | NEW F_n | ΔF_n | % change | OLD pen | NEW pen |
|---|---|---|---|---|---|---|---|
| 20 | 1589 | 1588.97 | 1482.26 | −106.71 | −6.72 % | −0.164 | +0.0823 |
| 30 | 2163 | 2163.46 | 1997.49 | −165.97 | −7.67 % | −0.228 | +0.1512 |
| 40 | 2700 | 2698.50 | 2444.06 | −254.44 | −9.43 % | −0.235 | +0.3069 |
| 50 | 3232 | 3231.76 | 2894.07 | −337.69 | −10.45 % | −0.086 | +0.6241 |

The NEW curve is **systematically below** the OLD curve, and the gap **grows monotonically
with sink** (5.9 % → 9.4 % on σ/σ_law) — so it is a genuine curve shift, not a constant offset.

### 4.3 OLD-recipe fidelity check (why the noise floor is ≤0.1 %)

| sink | historical σ/σ_law | OLD re-run | historical F_n | OLD F_n | historical pen | OLD pen |
|---|---|---|---|---|---|---|
| 20 | 1.56 | 1.5581 | 1589 | 1588.97 | −0.164 | −0.1637 |
| 30 | 1.59 | 1.5855 | 2163 | 2163.46 | −0.228 | −0.2276 |
| 40 | 1.62 | 1.6203 | 2700 | 2698.50 | −0.235 | −0.2349 |
| 50 | 1.64 | 1.6401 | 3232 | 3231.76 | −0.086 | −0.0863 |

`sink=50` OLD also re-emits the historical log signature exactly: `nconv=1478`,
`negJac=32`, `end_t=0.999996` (identical to `envelope\env_sink50.log`).

## 5. Isolating the CONTACT recipe from the time-stepper

The historical runs used `n=1200` default control; the NEW G7 deck uses `n=2400`
hardened. Two orthogonal controls at sink=50 separate the two effects:

| case | contact | stepper | σ/σ_law | F_n (N) | pen (mm) | term |
|---|---|---|---|---|---|---|
| `anc_old_s50` | AUGLAG pen=1 | n=1200 default | 1.6401 | 3231.76 | −0.0863 | NORMAL |
| `anc_new_s50_oldstep` | PENALTY 0.1 | n=1200 default | **1.5039** | 2928.55 | +0.6181 | NORMAL |
| `anc_new_s50` | PENALTY 0.1 | n=2400 hardened | **1.4863** | 2894.07 | +0.6241 | NORMAL |
| `anc_old_s50_newstep` | AUGLAG pen=1 | n=2400 hardened | — | — | — | **FAIL** (`negJac=55`, closure `t=0.0100`) |

* **Contact-only swap** (hold stepper = n=1200): `1.6401 → 1.5039` = **−0.1362 (−8.30 %)**.
* **Stepper-only swap** (hold contact = PENALTY 0.1): `1.5039 → 1.4863` = **−0.0176 (−1.17 %)**.
* So the **contact recipe is ~7× the stepper effect**, and both push the curve **down**.
* The OLD contact recipe **cannot even run** on the NEW stepper (dies at closure) — the two
  recipes are not interchangeable; the NEW stepper only works *because* the soft penalty
  it is paired with.

## 6. Verdict — **FAIL**

**The recipe change moves the curve. The historical anchors and the NEW G7 point are NOT on
the same curve.**

| decision input | value |
|---|---|
| σ/σ_law deviation across the 4 anchors (NEW − OLD) | **−5.85 % … −9.38 %** (monotone ↑ with sink) |
| deviation at the A/B anchor (sink 50) | **−0.1538 ratio, i.e. −9.38 %** |
| contact-recipe-only effect at fixed stepper (sink 50) | **−8.30 %** |
| stepper-only effect at fixed contact (sink 50) | −1.17 % |
| numerical noise floor (OLD re-run vs historical) | **≤0.1 %** |
| significance | recipe effect is **~60–100× the noise floor** |
| stated drift tolerance the doc already uses (`FORCE_BALANCE_TOL = 0.02`, 2 %) | exceeded by **3–5×** |

**Numerical margin:** at the recommended A/B anchor the recipe shift is
Δ(σ/σ_law) = **−0.154** (−9.4 %); the smallest shift anywhere in the band is
**−0.091** (−5.9 %). Both are far above the ≤0.1 % reproduction noise floor.

**Consequence for the paper.** The G7 number (σ/σ_law = **2.205**, `n7_tet4_s170`,
NEW recipe) sits on the **NEW curve**, which is ~6–9 % below the validated-anchor
**OLD curve** at ξ ≤ 0.25 (and the soft-penalty penetration mechanism shifts it further
at high strain). Quoting 2.205 alongside the ~1.6 anchors is a **cross-recipe comparison
and is not valid as-is**.

**To make a quotable number:** regenerate the small-strain anchors under the **NEW**
recipe (expected ≈ 1.47 / 1.48 / 1.48 / 1.49 at sink 20/30/40/50 — i.e. the curve drops
by ~6–9 %), then quote G7 against *that* curve. Alternatively, state the G7 value only as a
NEW-recipe result and explicitly retract cross-comparison with the OLD anchors.

> Caveat / scope: this test isolates the **recipe** at the fixed anchor mesh (3×3×2). The
> G7 point additionally uses a **finer mesh** (6×6×4). Mesh refinement is a separate
> variable not tested here; the recipe effect alone is already large enough to disqualify
> the cross-quote.

## 7. Reproduce

```powershell
$env:PYTHONPATH="src"
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py anc_new_s20
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py anc_new_s30
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py anc_new_s40
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py anc_new_s50
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\g7verify_a\run_recipe_ab.py anc_old_s50
# OLD at the other sinks + the stepper/contact isolation controls:
.venv\Scripts\python.exe -c "import sys; sys.path.insert(0,'temp/opensim_fe/s5_t9/g7verify_a'); import run_recipe_ab as r; [r.run_case(t, recipe='old', sink=s, nsteps=1200, timeout=1200) for t,s in [('anc_old_s20',20.0),('anc_old_s30',30.0),('anc_old_s40',40.0)]]"
.venv\Scripts\python.exe -c "import sys; sys.path.insert(0,'temp/opensim_fe/s5_t9/g7verify_a'); import run_recipe_ab as r; r.run_case('anc_new_s50_oldstep', recipe='new', stepper='old', sink=50.0, nsteps=1200, timeout=1200); r.run_case('anc_old_s50_newstep', recipe='old', stepper='new', sink=50.0, nsteps=2400, timeout=1200)"
```

Ledger: `temp\opensim_fe\s5_t9\g7verify_a\results.txt`. Decks `*.feb` + logs kept.

