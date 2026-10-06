# S5 G7 — element-formulation attribution: `udg-hex` vs `tet4` at MATCHED ξ

**Status:** COMPLETE (2026-10-06) — written incrementally; every number below was
appended only after its run finished and its real `.log`/`.xplt` tail was read.
**Verdict in §5: the 1.705 vs 2.205 pair is apples-to-oranges; at matched ξ the two
elements agree to ≤ ±2.5 % for ξ ≤ 0.65.**

> ## 🔴 SUPERSESSION BANNER — §5 recommendation only (added 2026-10-06)
>
> The §5 recommendation *“Quote the `tet4` completion: σ_contact/σ_1D-law = 2.205 @ ξ = 0.744”*
> is **superseded** by `docs\S5_g7_verify_mesh.md` (the `2.205` is not mesh-converged at the
> G7 sink — σ/σ_law drifts **+30.62 %** across the 3-mesh sweep with the refinement step
> growing) and `docs\S5_g7_verify_fixedxi.md` (no single quotable σ/σ_law with ξ and recipe
> can be given — at fixed ξ ≈ 0.70 the drift is **+43.13 %**).
>
> **Corrected guidance:** the only **mesh-converged** numbers in this file that may be
> quoted for the 2 m landing are the **`F_n ≈ 34.4 kN`** and **`σ_contact ≈ 0.574 MPa`**
> at fixed sink (NEW PENALTY recipe, tet4, <0.6 % drift). The σ/σ_law ratio is **NOT**
> quotable as a converged scalar. The matched-ξ element agreement (≤ ±2.5 % at ξ ≤ 0.65)
> remains valid as **method validation** — see `docs\S5_QUOTABLE.md`.
>
> **The measurements in this file (the matched-ξ table in §3.1 / §3.3, the
> endpoint-sink table in §3.2, the decomposition in §4) remain valid as measurements.**
> Only the §5 *recommendation* is superseded. No measured value, table or method text
> in this file has been altered.

**Owner task:** the write-up currently uses **two inconsistent σ_contact/σ_1D-law
numbers** for the same G7 contact result:

| element (pad) | end state | ξ | σ_contact (MPa) | σ_1D-law (MPa) | **σ/σ_law** |
|---|---|---:|---:|---:|---:|
| `udg-hex` (HEX8G1) | dies at ξ=0.709 (end_t=0.9742, NAN) | **0.7092** | 0.3858 | 0.2263 | **1.705** |
| `tet4` | completes G7 (end_t=1, NORMAL) | **0.7437** | 0.5767 | 0.2616 | **2.205** |

They differ by **+29.3 %** in the ratio, **but they are quoted at different ξ**
(0.7092 vs 0.7437). The question this doc settles:

> **How much of the 1.705 → 2.205 gap is genuine element formulation, and how much
> is just the different ξ at which each was measured?**

Method: compare the two elements **at the same ξ, on the same mesh sizes, with the
same recipe** (the NEW §7 recipe of `docs/S5_contact_udg.md`, used verbatim — not
re-derived). Only the element type differs.

---

## 1. Ground truth copied from the owned docs (NOT re-derived)

Source: `docs/S5_contact_udg.md` §6.NEW-4 / §6.NEW-7 / §7
(this doc does NOT modify that file).

* `udg-hex` G7 run `g7_pen0.1_sr20_s170` (reproducible ×3): `end_t=0.974227`,
  **indentation 141.8 mm**, **ξ=0.7092**, penetration 23.27 mm, **F_n=23 150 N**,
  σ_contact **0.3858 MPa**, σ_1D-law **0.2263 MPa**, **σ/σ_law=1.705**, `negJac=0`,
  NAN at the last step.
* `tet4` G7 run `n7_tet4_s170` (reproducible ×2, force-balance 7.80e-07):
  `end_t=1` NORMAL, **indentation 148.7 mm**, **ξ=0.7437**, penetration 20.76 mm,
  **F_n=34 600 N**, σ_contact **0.5767 MPa**, σ_1D-law **0.2616 MPa**,
  **σ/σ_law=2.205**.
* Both use pad `6×6×4` over `(300,300,200) mm` (UDG = `pad_block_hex`,
  tet4 = `pad_block_tet`); indenter `indenter_block_tet(2,2,1,(200,300,5),
  offset_mm=(50,0,200.5))`; gap 0.5 mm; same material; same recipe.

**Naive ratio 2.205 / 1.705 − 1 = +29.3 %** — this is the number the paper must
NOT present as an element effect until ξ is matched.

---

## 2. Protocol (pre-stated BEFORE any run)

**Recipe (exact, from `docs/S5_contact_udg.md` §7 — verbatim):**

```xml
<laugon>PENALTY</laugon>
<penalty>0.1</penalty>
<auto_penalty>0</auto_penalty>
<two_pass>0</two_pass>
<node_reloc>0</node_reloc>
<symmetric_stiffness>0</symmetric_stiffness>
<tolerance>0.005</tolerance>
<search_radius>20</search_radius>
<fric_coeff>0.6</fric_coeff>
```
Time stepper: `time_steps=2400`, `<step_size>=<dtmax>=1/2400`, `cutback=0.125`,
`max_retries=20`; solver BFGS `max_ups=10`, `reform_each_time_step=1`.

Held constant across both elements: same indenter, same gap (0.5 mm), same material,
same recipe, same mesh **sizes** `6×6×4`. Mesh geometry differs only in
hex-vs-6-tet cell split (the tet counterpart of the same structured grid).

**Matched-ξ extraction.** The indenter ramp is linear; every run writes all
converged states to the `.xplt`. So for each run I read **every** state
(`fe_post.read_displacement(..., last=False)`), compute ξ(state) = indentation/200,
and report the state **closest to each target ξ**. This gives an exact matched-ξ
comparison even when the ramp endpoints differ. Independently, I also run
per-element **endpoint** sinks chosen so the *final* state lands near the target ξ,
so the table contains genuine converged endpoints too.

**Planned sink values (pre-stated; targets ξ ≈ 0.35, 0.50, 0.65):**

| target ξ | indent (mm) | `udg-hex` sink (mm) | `tet4` sink (mm) |
|---|---:|---:|---:|
| 0.35 | 70 | 85 | 80 |
| 0.50 | 100 | 120 | 113 |
| 0.65 | 130 | 158 | 148 |

Sinks differ per element because the two elements penetrate differently at the same
travel (UDG penetrates more at large sink → needs more travel for the same
indentation). Plus the **G7 anchors** at sink=170 for both, and the tet4 state read
at the **common ξ = 0.7092** (the largest ξ both reach) — this is the cleanest
matched-ξ pair.

**Reported per state:** `end_t`, `negJac`, NAN, indentation, penetration, `F_n`,
force-balance (indenter-BC vs pad-top reaction), σ_contact, σ_1D-law, σ/σ_law, and
the **Δ between elements (%)**.

**Failure is data:** if `udg-hex` NANs before a target ξ, I record where it dies and
read the highest ξ it did reach.

**Run dir:** `temp\opensim_fe\s5_t9\g7verify_c\` (new; collides with nobody).
Harness reuses `udg_contact_sweep2.py` machinery (the previous agent's harness that
already knows `pad_domain_type="udg-hex"`).

---

## 3. Experiment log

_(appended one block per run, reading the real `.log`/`.xplt` tail — nothing held only in memory)_

### 3.1 — G7 anchors (sink 170) reproduce the two quoted numbers, and expose the matched-ξ states

Harness: `temp\opensim_fe\s5_t9\g7verify_c\run.py` (ledger `g7verify_results.txt`,
per-run `<tag>.summary.json`). Both runs use `search_radius=20`, `penalty=0.1`,
`laugon=PENALTY`, `node_reloc=0`, 2400 steps — the exact §7 recipe.

| run | pad | sink | rc | term | end_t | nconv | negJac | NAN | n_states |
|---|---|---:|---|---:|---:|---:|---:|---|---:|
| `g7_udg_s170` | udg-hex | 170 | EXC | ERROR | 0.974227 | 2340 | 0 | yes | 2341 |
| `g7_tet_s170` | tet4 | 170 | 0 | **NORMAL** | 1 | 2400 | 0 | no | 2401 |

**Reproduction is exact to the established numbers.** UDG endpoint: ξ=0.70924,
indent 141.85 mm, pen 23.27 mm, F_n 23 146 N, σ_c 0.38576, σ_law 0.22629,
**ratio 1.7047**, force-balance **3.08e-07**. tet4 endpoint: ξ=0.74368,
indent 148.74 mm, pen 20.76 mm, F_n 34 602 N, σ_c 0.57669, σ_law 0.26158,
**ratio 2.2047**, force-balance **7.80e-07**.

**Matched-ξ extraction** (read EVERY converged state from each `.xplt`; report the
state closest to each target ξ). This is the core result:

| ξ (matched) | pad | state | σ_contact (MPa) | σ_1D-law (MPa) | **σ/σ_law** | F_n (N) | pen (mm) | fbal |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.3499 | udg-hex | 1122 | 0.07746 | 0.04710 | **1.6446** | 4 647 | 8.99 | 1.17e-07 |
| 0.3501 | tet4 | 1074 | 0.07611 | 0.04713 | **1.6149** | 4 567 | 5.55 | 5.71e-05 |
| 0.5000 | udg-hex | 1634 | 0.15804 | 0.08893 | **1.7770** | 9 482 | 15.25 | 2.26e-07 |
| 0.4999 | tet4 | 1570 | 0.15537 | 0.08891 | **1.7475** | 9 322 | 10.73 | 1.61e-07 |
| 0.6500 | udg-hex | 2163 | 0.32789 | 0.17498 | **1.8739** | 19 674 | 22.72 | 1.08e-08 |
| 0.6501 | tet4 | 2066 | 0.33599 | 0.17508 | **1.9190** | 20 159 | 15.82 | 6.18e-07 |
| 0.7092 | udg-hex | 2340 | 0.38576 | 0.22629 | **1.7047** | 23 146 | 23.27 | 3.08e-07 |
| 0.7093 | tet4 | 2274 | 0.46570 | 0.22636 | **2.0573** | 27 942 | 18.71 | 5.95e-07 |

**Δ(tet4 − udg-hex)/udg-hex in σ/σ_law:**

| ξ | Δ σ/σ_law |
|---:|---:|
| 0.35 | **−1.8 %** |
| 0.50 | **−1.7 %** |
| 0.65 | **+2.4 %** |
| 0.709 | **+20.7 %** |

**Read:** the σ/σ_law ratio is a strong function of ξ for **both** elements (≈1.6 at
ξ=0.35 → ≈1.9 at ξ=0.65). At matched ξ ≤ 0.65 the two elements agree to **within
±2.5 %**. The divergence only opens up at the very last UDG state (ξ=0.709), where
UDG's σ_contact is ~20 % *below* tet4 at the same ξ — and that state is the
NAN-precursor endpoint. Note the UDG ratio even **dips** from 1.874 (ξ=0.65) to
1.705 (ξ=0.709), i.e. the UDG endpoint is off its own smooth curve.

### 3.2 — Endpoint runs at the target ξ (both elements terminate NORMAL)

Same recipe, same mesh **sizes**, sinks chosen per element so the *final* state lands
near the target ξ. **All six runs: `rc=0`, `NORMAL`, `end_t=1`, 2400/2400 steps,
`negJac=0`, no NAN.** (UDG reaches ξ≈0.71 at sink 170 and NANs *there*; at every
sink ≤158 it completes cleanly — the wall is a high-ξ phenomenon, not a sink one.)

| pad | sink (mm) | endpoint ξ | indent | pen (mm) | F_n (N) | σ_c (MPa) | σ_law | **σ/σ_law** | fbal | term |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| udg-hex | 85 | 0.3756 | 75.12 | 9.38 | 5 108 | 0.08513 | 0.05207 | 1.6349 | 5.1e-08 | NORMAL |
| tet4 | 80 | 0.3674 | 73.47 | 6.03 | 4 821 | 0.08035 | 0.05040 | 1.5942 | 7.4e-08 | NORMAL |
| udg-hex | 120 | 0.5153 | 103.06 | 16.44 | 10 273 | 0.17121 | 0.09533 | 1.7961 | 2.5e-08 | NORMAL |
| tet4 | 113 | 0.5075 | 101.50 | 11.00 | 9 599 | 0.15999 | 0.09202 | 1.7386 | 1.1e-07 | NORMAL |
| udg-hex | 158 | 0.6651 | 133.01 | 24.49 | 20 679 | 0.34466 | 0.18697 | 1.8434 | 7.2e-07 | NORMAL |
| tet4 | 148 | 0.6571 | 131.42 | 16.08 | 20 838 | 0.34730 | 0.18055 | 1.9236 | 1.6e-06 | NORMAL |

### 3.3 — Combined matched-ξ table (two independent measurements per ξ)

Each ξ is measured **twice**: once by state-extraction inside the sink-170 anchor
run (A), once by state-extraction inside the per-element endpoint run (E). The two
agree to ≲2.6 % (small load-path/increment dependence), so the result is robust.

| ξ (udg) | ξ (tet4) | source | udg σ/σ_law | tet4 σ/σ_law | **Δ = (tet−udg)/udg** |
|---:|---:|:--:|---:|---:|---:|
| 0.3499 | 0.3501 | A @170 | 1.6446 | 1.6149 | **−1.8 %** |
| 0.3500 | 0.3501 | E @85/80 | 1.6197 | 1.5833 | **−2.2 %** |
| 0.4999 | 0.4999 | A @170 | 1.7770 | 1.7475 | **−1.7 %** |
| 0.4999 | 0.5001 | E @120/113 | 1.7746 | 1.7314 | **−2.4 %** |
| 0.6500 | 0.6501 | A @170 | 1.8739 | 1.9190 | **+2.4 %** |
| 0.6499 | 0.6500 | E @158/148 | 1.8742 | 1.9109 | **+2.0 %** |
| **0.7092** | **0.7093** | A @170 | **1.7047** | **2.0573** | **+20.7 %** |

Full per-state detail for the matched picks (indentation, penetration, F_n, σ_c,
σ_law, force-balance) is in `g7verify_results.txt` and `*.summary.json`.

> **Δ in σ_contact is the same as Δ in the ratio** at matched ξ (the σ_1D-law term
> is a function of ξ only, so it cancels): at ξ=0.35 Δσ_c = −1.7 %, ξ=0.50 −1.7 %,
> ξ=0.65 +2.5 %, ξ=0.709 **+20.7 %**. At ξ=0.709, udg F_n=23 146 N vs tet4
> F_n=27 942 N — UDG carries ~17 % less contact force at the same indentation.

### 3.4 — The ratio is a strong function of ξ for BOTH elements

σ/σ_law rises monotonically with ξ from ≈1.6 (ξ=0.35) to ≈1.9 (ξ=0.65) for **either**
element. It is *not* a material constant; it is a measure of how the FE contact
force departs from the 1-D pad law as the pad densifies. Any two numbers quoted at
different ξ are not comparable — which is exactly the defect in the current pair.

Note also the **UDG endpoint is off its own curve**: UDG's ratio is 1.874 at
ξ=0.65 but only **1.705 at ξ=0.709** — it *dips* over the last 0.06 of ξ. That dip
is the signature of the contact degradation immediately preceding the NAN (the
final stored state is state 2340, the last converged step before the blow-up).

---

## 4. Decomposition of the 1.705 vs 2.205 gap

Original pair: `udg = 1.7047 @ ξ=0.7092` vs `tet4 = 2.2047 @ ξ=0.7437`
→ ratio ratio `2.2047 / 1.7047 = 1.293`, i.e. **+29.3 %**.

Decompose by holding ξ:

| step | value | vs udg endpoint |
|---|---:|---:|
| udg, as quoted | 1.7047 @ 0.7092 | — |
| tet4, as quoted | 2.2047 @ 0.7437 | **+29.3 %** |
| tet4 brought to ξ=0.7092 | 2.0573 | **+20.7 %** |
| tet4 brought to ξ=0.6500 | 1.9190 | +12.6 % |
| tet4 brought to ξ=0.5000 | 1.7475 | +2.5 % |

* **The different ξ accounts for ~8.6 pp of the 29.3 pp gap** (29.3 % → 20.7 %).
  That is a **~30 % share**, not the whole thing.
* The residual **+20.7 % at matched ξ=0.709** is *not* a stable element-material
  difference: it exists only because the UDG value is taken at its terminal
  (pre-NAN) state. Compare at ξ=0.65, where both elements are healthy and still
  on-curve: **Δ = +2.0 – +2.4 %**. Compare at ξ=0.50: **Δ = −1.7 – −2.4 %**. At
  ξ=0.35: **Δ = −1.8 – −2.2 %**.

**Therefore the 1.705 vs 2.205 pair is an apples-to-oranges comparison**: it
conflates (i) two different ξ **and** (ii) a healthy tet4 endpoint with a degraded
UDG pre-failure endpoint. The element formulation itself is *not* the driver of the
gap in the load-bearing range.

---

## 5. Verdict

**Not a genuine element discrepancy at matched ξ.**

* In the stable range **ξ ≤ 0.65**, `udg-hex` and `tet4` agree to **≤ ±2.5 %** in
  σ_contact/σ_1D-law on the same mesh geometry and recipe. The σ/σ_law ratio for
  *both* elements rises with ξ (≈1.6 → ≈1.9), so the ratio is a function of
  indentation, essentially independent of element type.
* The single point where the elements diverge (**+20.7 % at ξ=0.709**) is UDG's
  **last converged state before it NANs** — its σ_contact is ~17 % below tet4 and
  its own ratio has already dipped below its ξ=0.65 value. It is a failure
  signature, not a formulation property.
* **Matched-ξ gap, stated numerically:** at ξ=0.7092 the gap is **+20.7 %**; at
  ξ=0.65 it is **+2.0–2.4 %**; at ξ≤0.50 it is **≈ −2 %**. The headline 29.3 % is
  ~⅓ ξ-mismatch and ~⅔ UDG's terminal degradation.

**Which number is defensible for the paper at G7?**

* **Quote the `tet4` completion: σ_contact/σ_1D-law = 2.205 @ ξ = 0.744**
  (completed run, `rc=0` NORMAL, force-balance 7.8e-07, reproducible) — **provided
  the paper states ξ = 0.744 explicitly** and does not present it as if it were at
  UDG's ξ.
* **The UDG number `1.705` is NOT defensible as a G7 figure.** It is quoted at a
  *lower* ξ (0.709) **and** at a pre-failure state that is off UDG's own trend; it
  under-reports and coincides with the NAN.
* If a **matched-ξ** number is wanted anywhere, use **2.06 at ξ=0.709** (tet4), or
  **1.91 at ξ=0.65** (tet4), or quote the pair at ξ=0.65 where both are stable
  (**1.87 vs 1.91, Δ=+2.0 %**). Never present 1.705 and 2.205 side by side without
  their ξ.

---

## 6. Method / geometry / artifacts

* **Same geometry, only element type differs.** `pad_block_hex(6,6,4,(300,300,200))`
  and `pad_block_tet(6,6,4,(300,300,200))` share **byte-identical node coordinates
  and ordering** (verified: 245 nodes, `np.allclose` = True; the tet mesh is the
  exact 6-tet subdivision of the hex cells). The pad-top node set is identical.
  So the A/B isolates element formulation alone.
* **Held constant:** pad size/shape, indenter
  `indenter_block_tet(2,2,1,(200,300,5), offset_mm=(50,0,200.5))`, gap 0.5 mm,
  fitted Ogden material, symmetry BCs, and the exact §7 recipe
  (`laugon=PENALTY`, `penalty=0.1`, `search_radius=20`, `node_reloc=0`,
  `tolerance=0.005`, `mu=0.6`, 2400 steps, `dtmax=step_size=1/2400`,
  `cutback=0.125`, `max_retries=20`, BFGS `max_ups=10`, `reform_each_time_step=1`).
* **Force balance** was computed for every run (indenter-BC z-reaction vs pad-top
  z-reaction): all **≤ 1.6e-06** relative — consistent with the prior agent's
  7.80e-07 on the tet4 G7 run.
* **Artifacts:** harness `temp/opensim_fe/s5_t9/g7verify_c/run.py`; ledger
  `g7verify_results.txt`; per-run `*.summary.json`, `*.feb`, `*.log`, `*.xplt`.
* **Not used:** `docs/S5_contact_udg.md` was read only (never modified);
  `docs/S5_g7_verify_anchors.md` / `S5_g7_verify_mesh.md` untouched.

