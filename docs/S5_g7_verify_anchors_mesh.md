# S5 — G7 anchor set: mesh convergence of the NEW-recipe small-strain anchors

**Status:** COMPLETE (2026-10-06) — **verdict: FAIL** (§4). All 12 runs (3 pad meshes × 4 sinks)
finished `rc=0 / NORMAL / end_t=1 / no NAN`. Written **incrementally** — §1 (transcriptions) and
§2 (the planned grid) were written *before any compute*; the per-run blocks in §2 were appended
to disk the moment each single FE run finished. No result was ever held only in memory.
**Answer:** the NEW-recipe small-strain anchors are **NOT mesh-converged** — the ratio is
mesh-sensitive; no clean σ/σ_law curve may be quoted (only mesh-bounded bands, §5).

**Two questions this doc answers.**

1. **(i) Canonical anchor set.** The NEW-recipe anchors
   (sink 20/30/40/50 → σ/σ_law = 1.4670 / 1.4767 / 1.4818 / 1.4863) were run on **one mesh
   only** (`pad_block_tet(3,3,2)`) as a diagnostic. They can only become the *canonical*
   anchor set if they are **mesh-converged in the trusted small-strain range**.
2. **(ii) Is the small-strain ratio mesh-sensitive at all?** It **is** mesh-sensitive at G7
   (sink 170: σ/σ_law drifts **+12.0 % → +16.6 %** coarse→medium→fine, per
   `docs\S5_g7_verify_mesh.md`). The open question: does that instability also bite at
   **ξ ≤ 0.25** (sink ≤ 50)?

**Units:** mm–N–MPa–s. FEBio binary `D:\Program\FEBioStudio\bin\febio4.exe` (4.13.0).
**Run dir (new; collides with nobody):** `temp\opensim_fe\s5_t9\fixa\`.
**Harness:** `temp\opensim_fe\s5_t9\fixa\run_anchor_mesh.py` (extends the proven
`g7verify_a\run_recipe_ab.py`; one case per invocation; appends to this doc + `results.txt`
after **each** run).

---

## 1. Transcriptions (no compute)

### 1.1 OLD-recipe anchors (historical, trusted)

Transcribed from `docs\S5_g7_verify_anchors.md` §1.1 / §4.3 (which transcribed
`docs\S5_contact_convergence.md` §3.1 verbatim).

**Mesh:** tet4 `pad_block_tet(3,3,2,(300,300,200))` (48 nodes / 108 tets).
**Recipe:** `laugon="AUGLAG"`, `penalty=1.0`, `tolerance=0.005`,
`aug_controls=AUGLAG_CONVERGENT_CONTROLS` (gaptol=0.001, minaug=3, maxaug=200, smooth_aug=1),
`node_reloc=1`, `search_radius=20`, `n_solver_steps=1200`, default `<Control>`.
**Indenter:** `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`; gap 0.5 mm; μ 0.6.

| sink (mm) | ξ | F_n (N) | penetration (mm) | **σ/σ_law** | hist. quote | OLD re-run | provenance |
|---|---|---|---|---|---|---|---|
| 20 | 0.0983 | 1588.97 | −0.1637 | 1.5581 | 1.56 | 1.5581 | envelope.py / verify_lib.py |
| 30 | 0.1486 | 2163.46 | −0.2276 | 1.5855 | 1.59 | 1.5855 | envelope.py / verify_lib.py |
| 40 | 0.1987 | 2698.50 | −0.2349 | 1.6203 | 1.62 | 1.6203 | envelope.py |
| 50 | 0.2479 | 3231.76 | −0.0863 | 1.6401 | 1.64 | 1.6401 | envelope.py / verify_lib.py |

OLD-recipe re-runs reproduce the historical quotes to **≤ 0.1 %** — that is the numerical
noise floor of this harness. All OLD runs: `rc=0`, NORMAL termination, `end_t≈0.999996`,
no NAN. Penetration is **negative** (sub-mm micro-gap).

### 1.2 NEW-recipe anchors (the set under test for mesh convergence)

Transcribed from `docs\S5_g7_verify_anchors.md` §3. **Same mesh** as §1.1 (tet4
`(3,3,2)`), same indenter/gap/μ, **different contact recipe** (§1.3).

| sink (mm) | ξ | indentation (mm) | penetration (mm) | F_n (N) | σ_contact (MPa) | σ_1D-law (MPa) | **σ/σ_law** | force-balance | end_t | negJac | NAN |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 0.0971 | 19.42 | +0.0823 | 1482.26 | 0.024704 | 0.016840 | **1.4670** | 3.24e-07 | 1 | 28 | no |
| 30 | 0.1467 | 29.35 | +0.1512 | 1997.49 | 0.033292 | 0.022545 | **1.4767** | 1.22e-07 | 1 | 0 | no |
| 40 | 0.1960 | 39.19 | +0.3069 | 2444.06 | 0.040734 | 0.027489 | **1.4818** | 4.84e-08 | 1 | 7 | no |
| 50 | 0.2444 | 48.88 | +0.6241 | 2894.07 | 0.048235 | 0.032452 | **1.4863** | 3.87e-07 | 1 | 76 | no |

All four: `rc=0`, NORMAL termination, `end_t=1`, 2400/2400 steps, no NAN; force-balance
≈1e-7 (machine zero).

**Recipe shift (NEW − OLD), established, not re-proved here:** −5.85 % (s20) → −9.38 % (s50),
monotone; ~60–100× the ≤0.1 % noise floor. So the two recipes are **not** on the same curve
(`docs\S5_g7_verify_anchors.md` §6 verdict FAIL).

### 1.3 The NEW recipe under test (verbatim from `docs\S5_contact_udg.md` §7)

```xml
<contact name="indenter_pad" surface_pair="indenter_pad_pair" type="sliding-elastic">
  <laugon>PENALTY</laugon>          <!-- NO augmented Lagrangian -->
  <penalty>0.1</penalty>            <!-- soft; >=0.15 NANs at closure -->
  <auto_penalty>0</auto_penalty>
  <two_pass>0</two_pass>
  <node_reloc>0</node_reloc>
  <symmetric_stiffness>0</symmetric_stiffness>
  <tolerance>0.005</tolerance>
  <search_radius>20</search_radius>
  <fric_coeff>0.6</fric_coeff>
</contact>
```
plus `time_steps=2400`, `<step_size>=<dtmax>=1/2400`, `cutback=0.125`, `max_retries=20`,
`opt_iter=15`, BFGS `max_ups=10`, `reform_each_time_step=1`, PARDISO, `max_refs=25`.

### 1.4 G7 large-strain reference (mesh-divergence mechanism — the thing we test whether it bites)

From `docs\S5_g7_verify_mesh.md` §2, sink = 170 mm, identical NEW recipe:

| pad mesh | nodes/tets | indentation (mm) | ξ | penetration (mm) | F_n (N) | σ_contact | σ_law | **σ/σ_law** |
|---|---|---|---|---|---|---|---|---|
| (3,3,2) | 48/108 | 153.965 | 0.76983 | 15.535 | 34404.2 | 0.573404 | 0.291311 | 1.96836 |
| (6,6,4) | 245/864 | 148.736 | 0.74368 | 20.764 | 34601.7 | 0.576695 | 0.261581 | 2.20465 |
| (12,12,8) | 1521/6912 | 141.156 | 0.70578 | 28.344 | 34395.2 | 0.573254 | 0.222970 | 2.57099 |

Mechanism: soft `penalty=0.1` makes **penetration mesh-dependent** (+33.7 % → +36.5 %), which
drags **ξ down** and therefore drags the denominator `σ_1D-law(ξ)` down — while `F_n` and
σ_contact stay mesh-converged (< 0.6 %). The ratio therefore diverges. **This doc tests whether
that chain also runs at ξ ≤ 0.25.**

---

## 2. Mesh × sink grid (planned BEFORE any compute)

**Only the pad mesh varies.** Indenter, gap, recipe, load, BCs are held fixed at §1.2/§1.3.

| pad `(nx,ny,nz)` | expected pad nodes / tets | element size | label |
|---|---|---|---|
| `pad_block_tet(3,3,2,(300,300,200))` | 48 / 108 | 100 mm | coarse |
| `pad_block_tet(6,6,4,(300,300,200))` | 245 / 864 | 50 mm | medium |
| `pad_block_tet(12,12,8,(300,300,200))` | 1521 / 6912 | 25 mm | fine |

× sink ∈ {20, 30, 40, 50} mm → **12 cases**. Priority order if the full grid is too slow:
**sink=50 across all 3 meshes first** (that is the mesh check), then complete the grid on the
mesh declared canonical.

Per-run fields (all from the real `.log` tail + `.xplt`/`.hdf5` read-back): `end_t`, `negJac`,
NAN, indentation, ξ, penetration, `F_n`, force-balance, σ_contact, σ_1D-law, **σ/σ_law**.
σ_contact = F_n / (indenter top area). σ_1D-law = `pad_stress_pa(ξ)/1e6`.

<!-- runner appends one block per finished run below -->

### Run `fxa_3x3x2_s50` — pad `(3,3,2)` (coarse), sink = 50 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **23 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0103851 |
| nconv | 2407 |
| negJac | 3 |
| NAN | no |
| n_states | 2408 |
| indentation (mm) | 48.876 |
| xi | 0.2444 |
| penetration (mm) | 0.62432 |
| F_n (N) | 2892.59 |
| sigma_contact (MPa) | 0.0482099 |
| sigma_1D-law (MPa) | 0.0324527 |
| **sigma/sigma_law** | **1.48554** |
| force_balance_rel | 6.001e-08 |

---

### Run `fxa_6x6x4_s50` — pad `(6,6,4)` (medium), sink = 50 mm

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **55 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 46.322 |
| xi | 0.2316 |
| penetration (mm) | 3.1779 |
| F_n (N) | 2819.52 |
| sigma_contact (MPa) | 0.0469919 |
| sigma_1D-law (MPa) | 0.0310925 |
| **sigma/sigma_law** | **1.51136** |
| force_balance_rel | 2.353e-07 |

---

### Run `fxa_12x12x8_s50` — pad `(12,12,8)` (fine), sink = 50 mm

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **263 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 44.121 |
| xi | 0.2206 |
| penetration (mm) | 5.3785 |
| F_n (N) | 2650.14 |
| sigma_contact (MPa) | 0.044169 |
| sigma_1D-law (MPa) | 0.0299565 |
| **sigma/sigma_law** | **1.47444** |
| force_balance_rel | 1.833e-07 |

---

### Run `fxa_3x3x2_s20` — pad `(3,3,2)` (coarse), sink = 20 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **20 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0254167 |
| nconv | 2406 |
| negJac | 41 |
| NAN | no |
| n_states | 2407 |
| indentation (mm) | 19.418 |
| xi | 0.0971 |
| penetration (mm) | 0.082249 |
| F_n (N) | 1482.05 |
| sigma_contact (MPa) | 0.0247008 |
| sigma_1D-law (MPa) | 0.0168396 |
| **sigma/sigma_law** | **1.46683** |
| force_balance_rel | 8.365e-09 |

---

### Run `fxa_6x6x4_s20` — pad `(6,6,4)` (medium), sink = 20 mm

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **45 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 18.874 |
| xi | 0.0944 |
| penetration (mm) | 0.62577 |
| F_n (N) | 1447.96 |
| sigma_contact (MPa) | 0.0241327 |
| sigma_1D-law (MPa) | 0.0164891 |
| **sigma/sigma_law** | **1.46355** |
| force_balance_rel | 5.467e-08 |

---

### Run `fxa_12x12x8_s20` — pad `(12,12,8)` (fine), sink = 20 mm

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **262 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 17.524 |
| xi | 0.0876 |
| penetration (mm) | 1.9765 |
| F_n (N) | 1297.54 |
| sigma_contact (MPa) | 0.0216257 |
| sigma_1D-law (MPa) | 0.0155961 |
| **sigma/sigma_law** | **1.38662** |
| force_balance_rel | 3.335e-08 |

---

### Run `fxa_3x3x2_s30` — pad `(3,3,2)` (coarse), sink = 30 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **14 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0175662 |
| nconv | 2409 |
| negJac | 13 |
| NAN | no |
| n_states | 2410 |
| indentation (mm) | 29.349 |
| xi | 0.1467 |
| penetration (mm) | 0.15123 |
| F_n (N) | 1997.47 |
| sigma_contact (MPa) | 0.0332911 |
| sigma_1D-law (MPa) | 0.0225439 |
| **sigma/sigma_law** | **1.47672** |
| force_balance_rel | 2.063e-07 |

---

### Run `fxa_6x6x4_s30` — pad `(6,6,4)` (medium), sink = 30 mm

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **41 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 28.228 |
| xi | 0.1411 |
| penetration (mm) | 1.2724 |
| F_n (N) | 1948.6 |
| sigma_contact (MPa) | 0.0324766 |
| sigma_1D-law (MPa) | 0.0219519 |
| **sigma/sigma_law** | **1.47944** |
| force_balance_rel | 4.356e-08 |

---

### Run `fxa_12x12x8_s30` — pad `(12,12,8)` (fine), sink = 30 mm

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **195 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 26.502 |
| xi | 0.1325 |
| penetration (mm) | 2.9975 |
| F_n (N) | 1776.07 |
| sigma_contact (MPa) | 0.0296012 |
| sigma_1D-law (MPa) | 0.0210205 |
| **sigma/sigma_law** | **1.40821** |
| force_balance_rel | 1.997e-08 |

---

### Run `fxa_3x3x2_s40` — pad `(3,3,2)` (coarse), sink = 40 mm

pad: **48 nodes / 108 tet4**; combined deck: 66 nodes / 132 elems.  wall-clock **13 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | 0.0127778 |
| nconv | 2406 |
| negJac | 31 |
| NAN | no |
| n_states | 2407 |
| indentation (mm) | 39.193 |
| xi | 0.1960 |
| penetration (mm) | 0.30694 |
| F_n (N) | 2444.07 |
| sigma_contact (MPa) | 0.0407345 |
| sigma_1D-law (MPa) | 0.0274892 |
| **sigma/sigma_law** | **1.48184** |
| force_balance_rel | 6.321e-08 |

---

### Run `fxa_6x6x4_s40` — pad `(6,6,4)` (medium), sink = 40 mm

pad: **245 nodes / 864 tet4**; combined deck: 263 nodes / 888 elems.  wall-clock **42 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 37.237 |
| xi | 0.1862 |
| penetration (mm) | 2.2626 |
| F_n (N) | 2377.51 |
| sigma_contact (MPa) | 0.0396252 |
| sigma_1D-law (MPa) | 0.0265232 |
| **sigma/sigma_law** | **1.49398** |
| force_balance_rel | 1.16e-07 |

---

### Run `fxa_12x12x8_s40` — pad `(12,12,8)` (fine), sink = 40 mm

pad: **1521 nodes / 6912 tet4**; combined deck: 1539 nodes / 6936 elems.  wall-clock **367 s**; rc=0.

recipe: NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`, `two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`, `fric_coeff=0.6`; stepper `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`.

indenter `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))` · gap 0.5 mm · **only the pad mesh varies**.

| field | value |
|---|---|
| termination | NORMAL |
| end_t | 1 |
| fail_t | - |
| nconv | 2400 |
| negJac | 0 |
| NAN | no |
| n_states | 2401 |
| indentation (mm) | 35.307 |
| xi | 0.1765 |
| penetration (mm) | 4.1934 |
| F_n (N) | 2204.39 |
| sigma_contact (MPa) | 0.0367399 |
| sigma_1D-law (MPa) | 0.0255677 |
| **sigma/sigma_law** | **1.43697** |
| force_balance_rel | 1.474e-08 |

---

## 2.1 Consolidated anchor table — NEW recipe, 3 pad meshes × 4 sinks

All 12 runs: `rc=0`, **NORMAL termination, `end_t=1`, no NAN**; transient `negJac` events only on
the coarse mesh were recovered by cutback. Indenter / gap (0.5 mm) / recipe / load held fixed;
**only the pad mesh varies**. σ_contact = `F_n / 60 000 mm²` (indenter top 200×300). σ_1D-law =
`pad_stress_pa(ξ)/1e6`. Coarse column re-emits the diagnostic anchors of §1.2 to ≤0.05 %.

| sink | pad mesh | nodes/tets | end_t | negJac | NAN | indentation (mm) | ξ | penetration (mm) | F_n (N) | force-balance | σ_contact (MPa) | σ_1D-law (MPa) | **σ/σ_law** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | coarse (3,3,2) | 48/108 | 1 | 41 | no | 19.418 | 0.0971 | +0.0822 | 1482.05 | 8.37e-09 | 0.024701 | 0.016840 | **1.46683** |
| 20 | medium (6,6,4) | 245/864 | 1 | 0 | no | 18.874 | 0.0944 | +0.6258 | 1447.96 | 5.47e-08 | 0.024133 | 0.016489 | **1.46355** |
| 20 | fine (12,12,8) | 1521/6912 | 1 | 0 | no | 17.524 | 0.0876 | +1.9765 | 1297.54 | 3.34e-08 | 0.021626 | 0.015596 | **1.38662** |
| 30 | coarse (3,3,2) | 48/108 | 1 | 13 | no | 29.349 | 0.1467 | +0.1512 | 1997.47 | 2.06e-07 | 0.033291 | 0.022544 | **1.47672** |
| 30 | medium (6,6,4) | 245/864 | 1 | 0 | no | 28.228 | 0.1411 | +1.2724 | 1948.60 | 4.36e-08 | 0.032477 | 0.021952 | **1.47944** |
| 30 | fine (12,12,8) | 1521/6912 | 1 | 0 | no | 26.502 | 0.1325 | +2.9975 | 1776.07 | 2.00e-08 | 0.029601 | 0.021021 | **1.40821** |
| 40 | coarse (3,3,2) | 48/108 | 1 | 31 | no | 39.193 | 0.1960 | +0.3069 | 2444.07 | 6.32e-08 | 0.040734 | 0.027489 | **1.48184** |
| 40 | medium (6,6,4) | 245/864 | 1 | 0 | no | 37.237 | 0.1862 | +2.2626 | 2377.51 | 1.16e-07 | 0.039625 | 0.026523 | **1.49398** |
| 40 | fine (12,12,8) | 1521/6912 | 1 | 0 | no | 35.307 | 0.1765 | +4.1934 | 2204.39 | 1.47e-08 | 0.036740 | 0.025568 | **1.43697** |
| 50 | coarse (3,3,2) | 48/108 | 1 | 3 | no | 48.876 | 0.2444 | +0.6243 | 2892.59 | 6.00e-08 | 0.048210 | 0.032453 | **1.48554** |
| 50 | medium (6,6,4) | 245/864 | 1 | 0 | no | 46.322 | 0.2316 | +3.1779 | 2819.52 | 2.35e-07 | 0.046992 | 0.031093 | **1.51136** |
| 50 | fine (12,12,8) | 1521/6912 | 1 | 0 | no | 44.121 | 0.2206 | +5.3785 | 2650.14 | 1.83e-07 | 0.044169 | 0.029957 | **1.47444** |

Wall-clock: coarse 13–23 s, medium 41–55 s, fine 195–367 s.

---

## 3. Mesh drift — per sink

Drift relative to the **coarser** member of each pair (`c→m`, `m→f`) and net coarse→fine
(`c→f`); spread = (max−min)/coarse over the 3 meshes. All in %.

| sink | σ/σ_law coarse | medium | fine | **c→m** | **m→f** | **net c→f** | **spread** | F_n net | ξ net | indentation net | penetration c→f |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 1.46683 | 1.46355 | 1.38662 | −0.22 | **−5.26** | **−5.47** | **5.47** | −12.45 | −9.78 | −9.75 | 0.082 → 1.98 mm (×24) |
| 30 | 1.47672 | 1.47944 | 1.40821 | +0.18 | **−4.81** | **−4.64** | **4.82** | −11.08 | −9.68 | −9.70 | 0.151 → 3.00 mm (×20) |
| 40 | 1.48184 | 1.49398 | 1.43697 | +0.82 | **−3.82** | **−3.03** | **3.85** | −9.81 | −9.95 | −9.92 | 0.307 → 4.19 mm (×14) |
| 50 | 1.48554 | 1.51136 | 1.47444 | +1.74 | **−2.44** | **−0.75** | **2.49** | −8.38 | −9.74 | −9.73 | 0.624 → 5.38 mm (×9) |

**Reading of the table.**

* **The ratio is mesh-sensitive, and NON-MONOTONE.** coarse↔medium agree to 0.2–1.7 %, but
  fine sits **2.4–5.3 % below** medium. The sequence **rises then falls** (s30/40/50) or falls
  with an accelerating step (s20) — it is not approaching a limit.
* **The refinement step does not shrink — it grows.** `|c→m|` is 0.2–1.7 %; `|m→f|` is
  2.4–5.3 %. A quantity whose successive refinement step *increases* is not in its asymptotic
  range, so no mesh in this family can be declared "the converged one."
* **The ratio's ingredients drift ~10 % net each — the ratio only survives by cancellation.**
  `F_n` (and therefore σ_contact) drifts **−8.4 … −12.5 %**; `ξ`/indentation drift
  **−9.7 … −10.0 %**; σ_1D-law(ξ) drifts −6.8 … −7.7 %. The ratio is the quotient of two
  quantities moving the *same* direction, so its ~9–12 % ingredient drift partly cancels to
  0.75–5.5 %.
* **F_n is NOT mesh-converged at small strain** (−8.4 … −12.5 % net) — the opposite of G7,
  where F_n was converged to <0.6 %. At ξ ≤ 0.25 the response is dominated by the local
  flat-punch contact compliance (mesh-sensitive under soft penalty); at G7 it is dominated by
  the bulk material response (mesh-insensitive).
* **Penetration is the driving artifact.** Under `penalty=0.1` the surface penetration grows
  monotonically with refinement — **×24 (s20), ×20 (s30), ×14 (s40), ×9 (s50)** — exactly the
  G7 mechanism, but starting from a much smaller base.

---

## 4. Verdict — **FAIL** (the small-strain σ/σ_law is NOT mesh-converged)

**Tolerance judged.** Mesh convergence is declared PASS only if **both** hold:
(a) the **net coarse→fine drift ≤ 2 %** (the repo's own `FORCE_BALANCE_TOL = 0.02` convention,
used as the G7 mesh bar), **and** (b) the **successive refinement step is non-increasing**
(standard asymptotic-range requirement). Neither the anchors as a set, nor any individual sink,
satisfies (b); s20/s30/s40 also fail (a).

| sink | net c→f | |m→f| vs |c→m| | passes (a) ≤2 %? | passes (b)? | verdict |
|---|---|---|---|---|---|
| 20 | −5.47 % | 5.26 > 0.22 (grows) | ❌ | ❌ | **FAIL** |
| 30 | −4.64 % | 4.81 > 0.18 (grows) | ❌ | ❌ | **FAIL** |
| 40 | −3.03 % | 3.82 > 0.82 (grows) | ❌ | ❌ | **FAIL** |
| 50 | −0.75 % | 2.44 > 1.74 (grows, non-monotone) | ✅ | ❌ | **FAIL** |

* **Answers question (i):** the NEW-recipe diagnostic anchors (1.4670 / 1.4767 / 1.4818 /
  1.4863) are the **coarse `(3,3,2)`** values and reproduce this sweep to ≤0.05 %. They **cannot**
  be promoted to a canonical *mesh-converged* anchor set: refining the pad moves the ratio by
  **−5.5 % / −4.6 % / −3.0 % / −0.75 %** (coarse→fine), with a **growing, non-monotone**
  refinement step (fine is 2.4–5.3 % below medium). The set is mesh-sensitive.
* **Answers question (ii):** **yes**, the small-strain ratio **is** mesh-sensitive — a
  **2.5–5.5 % spread** across the 3 meshes (down from **+30.6 %** at G7, but the same
  *direction of cause*: soft-penalty penetration grows with refinement → ξ falls → σ_law(ξ)
  falls). The chain that made G7 FAIL **is active at ξ ≤ 0.25**, just weaker.
* **Mechanism (not over-claim).** The drift is a **soft-penalty contact artifact**, not a
  material/convergence property: `penetration` grows monotonically with refinement
  (×9–×24), which steals travel from the pad, lowers `ξ` and `F_n`, and perturbs
  σ/σ_law. The exact recipe is fixed by `docs\S5_contact_udg.md` §7 and may not be changed, so
  within the canonical recipe the ratios are mesh-sensitive by construction.

### 4.1 Contrast with the G7 mesh study (`docs\S5_g7_verify_mesh.md`)

| | G7 (sink 170) | anchors (sink 20–50) |
|---|---|---|
| σ/σ_law net c→f | **+30.6 %** (rising) | **−0.75 … −5.47 %** (falling/non-monotone) |
| refinement step | +12.0 % → +16.6 % (grows) | 0.2–1.7 % → 2.4–5.3 % (grows) |
| F_n net c→f | −0.03 % ✅ converged | **−8.4 … −12.5 %** ❌ not converged |
| σ_contact net | −0.03 % ✅ | −8.4 … −12.5 % ❌ |
| penetration c→f | 15.5 → 28.3 mm (×1.8) | 0.08 → 5.4 mm (×9–×24) |
| verdict | FAIL | **FAIL** |

So the mesh instability is **weaker in magnitude and opposite in sign** at small strain, but it
is **not absent** and **not convergent** on this mesh family.

---

## 5. Which numbers the paper may cite

**Recipe for every citable number:** NEW — `laugon=PENALTY`, `penalty=0.1`, `auto_penalty=0`,
`two_pass=0`, `node_reloc=0`, `symmetric_stiffness=0`, `tolerance=0.005`, `search_radius=20`,
`fric_coeff=0.6`, `time_steps=2400`, `step_size=dtmax=1/2400`, `cutback=0.125`,
`max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`; tet4 pad; indenter
`(2,2,1,(200,300,5))` at offset `(50,0,200.5)`; gap 0.5 mm; μ 0.6.

**May cite — with the stated qualification:**

1. **σ/σ_law small-strain curve — only as a bounded band, never as a single converged value.**
   State the mesh explicitly. The only defensible anchor numbers are the **coarse `(3,3,2)`**
   values (they reproduce the existing diagnostic to ≤0.05 %), quoted with the coarse→fine
   mesh band as the uncertainty:

   | sink (mm) | ξ | σ/σ_law `(3,3,2)` `(6,6,4)` `(12,12,8)` (coarse/medium/fine) | mesh band (c→f) | defensible band |
   |---|---|---|---|---|
   | 20 | 0.097 | 1.467 / 1.464 / 1.387 | −5.47 % | **1.39 – 1.47** |
   | 30 | 0.147 | 1.477 / 1.479 / 1.408 | −4.64 % | **1.41 – 1.48** |
   | 40 | 0.196 | 1.482 / 1.494 / 1.437 | −3.03 % | **1.44 – 1.49** |
   | 50 | 0.244 | 1.486 / 1.511 / 1.474 | −0.75 % | **1.47 – 1.51** |

   Any quote must say: *"NEW recipe, tet4; value varies by X % over a 64× element-count change
   with a non-decreasing refinement step, so it is mesh-sensitive, not mesh-converged."*

2. **σ_contact = F_n/60 000 mm²** — cite with the **same −8 … −12.5 % mesh band**.
3. **F_n** — cite with the **−8 … −12.5 % mesh band** (contrast G7, where F_n was converged).
4. **ξ / indentation** — cite with the **−9.7 … −10.0 % mesh band**.
5. **σ_1D-law(ξ)** — this is a **material function** (`pad_stress_pa`), mesh-independent; the
   observed −6.8 … −7.7 % drift is inherited solely from the measured ξ, not a material effect.

**May NOT cite (this doc's FAIL result):**

* Any single σ/σ_law anchor (including 1.4670 / 1.4767 / 1.4818 / 1.4863) as a
  **mesh-converged** value, or the coarse-only curve without the mesh band above.
* The **retracted** "15.14 cm / 2.70×".
* Any **cross-recipe** comparison (OLD AUGLAG 1.56–1.64 vs NEW PENALTY 1.47–1.49) — the two
  recipes are different curves (`docs\S5_g7_verify_anchors.md` §6), and the mesh band is
  **comparable in size** to that recipe shift (5.5 % vs 5.9–9.4 %), so the two effects are
  entangled.
* The G7 2.205 as lying on the same converged curve as these anchors (G7 is a separate,
  itself-non-converged large-strain point).
* `F_n` as a mesh-converged resultant **at these small-strain anchor levels** (it is converged
  only at G7).

**Bottom line.** Under the canonical NEW recipe, **no σ/σ_law value is mesh-converged in the
ξ ≤ 0.25 range, so no clean ratio curve may be quoted.** The maximum defensible statement is a
**mesh-bounded band per anchor** (table above) with the mesh and the non-decreasing refinement
step both disclosed. The 1D-law σ is the only mesh-independent quantity here.

---

## 6. Reproduce

```powershell
$env:PYTHONPATH="src"
# whole grid (crash-resumable; skips tags already in results.txt, appends doc after each run)
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\fixa\run_anchor_mesh.py
# single case:
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\fixa\run_anchor_mesh.py fxa_12x12x8_s40
```

Ledger: `temp\opensim_fe\s5_t9\fixa\results.txt`. Decks `*.feb` + logs + state files kept in
`temp\opensim_fe\s5_t9\fixa\`. Coarse column reproduces `docs\S5_g7_verify_anchors.md` §3 to
≤0.05 % (the ≤0.1 % noise floor established there).

