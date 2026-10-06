# S5-contact wave 4: design spec

> **STATUS: APPROVED (2026-10-06).** All 10 open questions resolved by the coordinator in
> **§5.1** (the blocking items 1–4 decided with written rationale). Implementer: execute
> §5.1 + §0 + §1 + §2 + §3 + §4; do not re-litigate a resolved fork.
> Supersedes the earlier `PROPOSAL (awaiting user approval)` state.

| | |
|---|---|
| Generated | 2026-10-06 (proposal stage) |
| Script (target) | `scripts/opensim_fe/nonvertical_s5_contact.py` (not yet written) |
| Report (fixed by `docs\S5接触方案.md` §9 + `docs\OpenSim_FE交接.md` §7) | `results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md` |
| JSON (proposed, see §4) | `results/opensim_fe/nonvertical_s5_contact.json` |
| Plan / contract | `docs\非垂直落地扩展方案.md` §3 L4 / §4 S5 · `docs\S5接触方案.md` §3 / §9 · `docs\S5_contact_contract.md` §4 / §6 |
| Quotable numbers | `docs\S5_QUOTABLE.md` (authoritative; this spec defers to it for any σ/σ_law claim) |
| Label taxonomy | `results\opensim_fe\figures\FIGURE_CAPTIONS.md` §"标签约定", `measured` / `modeled` / `assumed` |

This file is the missing link in the deliverable named by `docs\S5接触方案.md:106` and
`docs\OpenSim_FE交接.md:291`: the contract listed "脚本 + 报告 + JSON + 测试" but **the JSON
schema was never defined** anywhere in the repo. A prior survey confirmed no schema exists
(see sibling `scripts\opensim_fe\nonvertical_s5_supination.py:66` `OUT_JSON = RES / "nonvertical_s5.json"`
for the existing slot this script cannot reuse). §4 below argues the new basename.

---

## 0. Conventions transcribed verbatim (do not invent)

The implementer must mirror these, every value below was harvested from the existing
sibling script, its outputs, and the contracts. Anything not transcribed here is a **deliberate
proposal** and lives in §1 only.

### 0.1 CLI surface (mirror `nonvertical_s5_supination.py:309-315`)

```
ap = argparse.ArgumentParser(description="…")
ap.add_argument("--out-json", type=Path, default=OUT_JSON)
ap.add_argument("--out-md",   type=Path, default=OUT_MD)
ap.add_argument("--fast",     action="store_true", help="…")
args = ap.parse_args(argv)
```

`def main(argv: list[str] | None = None) -> int`. Exit codes: `0` PASS/SKIPPED, `1` FAIL,
`2` overwrite-refused.

### 0.2 Overwrite guard (mandatory, sibling `nonvertical_s5_supination.py:317-320`)

```
for p in (args.out_json, args.out_md):
    if p.exists():
        print(f"拒绝覆盖既有文件：{p}", file=sys.stderr)
        return 2
```

### 0.3 `sys.path` bootstrap (sibling `:45-49`)

```
_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (str(ROOT / "src"), str(_HERE)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
```

### 0.4 JSON writer (sibling `:428-430`)

```
args.out_json.parent.mkdir(parents=True, exist_ok=True)
args.out_json.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=float),
                         encoding="utf-8")
```

No `round()`: the JSON must carry full IEEE-754 precision via `default=float`.

### 0.5 Report writer (sibling `:454, 724`)

```
L: list[str] = []
A = L.append
A(…)
…
path.write_text("\n".join(L) + "\n", encoding="utf-8")
```

### 0.6 Report shape (sibling `:456-723`)

Title line → several `>`-prefixed meta lines → `## 0. 一句话结论` (numbered list) →
numbered `## 1.` … sections in this order: API table, regression, equations, geometry,
scenario tables, sensitivity, equivalence, soft-pad comparison, critical-angle, assumptions
& honesty boundary, repro commands, outputs. Wave 4 reuses the same envelope; §2 below.

### 0.7 Existing sibling JSON keys (shape reference; all values are `measured`/`modeled`/`assumed`)

Inferred from `results\opensim_fe\nonvertical_s5.json`, `ankle_external_rotation.json`,
`nonvertical_s4.json`: `meta`, `interface_regression`, `axial_regression`, `<bone>_sections`,
`*_force_multiplier`, `axial_baseline_reference`, `measured_heights`, `scenario_hard`,
`scenario_mat`, `sensitivity`, `assumptions`. The wave 4 schema (proposed in §1) keeps
the top-level envelope (`meta`, `interface_regression`, `axial_regression`, `assumptions`)
and adds new blocks whose keys are spelled out in §1.

### 0.8 Runner (sibling `:730-731`)

```
if __name__ == "__main__":
    raise SystemExit(main())
```

### 0.9 Constraint red lines (verbatim from `docs\S5_contact_contract.md` §5)

> 只读勿改：`scripts\ankle_fe\`、`temp\pyfebio_demo\`、`temp\ankle_*`、`FE_PIPELINE_TODO.md`、
> `LEARNINGS` LRN-010~019、以及所有既有 `results\opensim_fe/*` 产物（可新增文件，不可覆盖）。

### 0.10 Scope-limit citation (mandatory)

`docs\S5_contact_contract.md` §4 **条 7** (2026-10-06 晚间修订) replaced the original "波3/波4
不应启动" clause with "**可启动但限定可信口径**". This spec honours the same cite by binding
the schema's `scope_limit` block to 条 7, see §1.2.

### 0.11 FE invocation (verbatim from contract §5)

> FEBio：`D:\Program\FEBioStudio\bin\febio4.exe`；**用 `febio_run.run_febio`**，不要自己起 subprocess。
> .feb 段落顺序：**Loads → Boundary → Contact**（顺序错 = `unrecognized tag`）。
> 接触面必须是**真实单元面**（tet facet），不能是 npz 里的四边形。

### 0.12 Units

**mm-N-MPa-s**. OpenSim model is in metres (×1000 to mm).

---

## 1. JSON schema proposal

> **PROPOSAL.** The implementer may not write a script against this schema until the user
> has approved §1, §3, and §4 below.

### 1.1 Top-level envelope

```jsonc
{
  "meta":               { … },          // run metadata (§1.3)
  "scope_limit":        { … },          // §4 条 7 binding (§1.4, MANDATORY)
  "interface_regression":{ … },         // §1.5, same shape as sibling
  "axial_regression":   { … },          // §1.6, same shape as sibling
  "scenario_matrix":    { … },          // §1.7, axes + run ids
  "contact_recipe":     { … },          // §1.8, FEBio recipe block (named, reproducible)
  "runs": [ … ],                        // §1.9, one entry per (id, recipe)
  "regression_gates":   { … },          // §1.10, G0..G6 verdicts
  "sensitivity":        { … },          // §1.11, μ / penalty / k sweeps
  "assumptions":        [ … ],          // §1.12, measured/modeled/assumed list
  "repro":              { … }           // §1.13, commands + files read
}
```

All top-level keys are required; arrays may be empty for `--fast` runs (see §1.10 G0 note).

### 1.2 Label taxonomy (verbatim from `results\opensim_fe\figures\FIGURE_CAPTIONS.md`)

Each numeric value in the JSON carries an explicit `label` field (`"measured"`, `"modeled"`,
or `"assumed"`). Aggregates that mix several labels use `"measured"` only when every
constituent is measured. The two formats are:

- **Scalar with label**, `{"value": <float>, "label": "measured"|"modeled"|"assumed",
  "unit": "…", "source": "…"}`
- **Compound row with per-field labels**, fields are individually labelled as above; the
  row itself does not carry a single label.

### 1.3 `meta`: run metadata

| Key | Type | Unit | Notes |
|---|---|---|---|
| `generated_at` | str | - | `datetime.now().strftime("%Y-%m-%d %H:%M:%S")` |
| `study` | str | - | constant: `"Nonvertical S5-contact, FEBio plantar-pad contact + friction (wave 4, opt-in)"` |
| `plan_doc` | str | - | `"docs/非垂直落地扩展方案.md §3 L4 / §4 S5 · docs/S5接触方案.md §3"` |
| `contract_ref` | str | - | `"docs/S5_contact_contract.md §4 (acceptance G0..G7) + §4 条 7 (scope limit)"` |
| `quotable_ref` | str | - | `"docs/S5_QUOTABLE.md"` (authoritative for any cite) |
| `units` | str | - | constant: `"mm-N-MPa-s"` |
| `mass_kg` | float | kg | 75.337 (sibling constant), `label=measured` (Rajagopal 2015 mass, propagated from `nonvertical_s1.json` / `axial_subregion.json`) |
| `recipe_id` | str | - | name of the recipe block under `contact_recipe`, e.g. `"minimum_rigid_wall_mu0.6"` |
| `plantar_bc_kind` | str | - | constant for this wave: `"contact"` (opt-in; the default `"fixed"` / `"roller"` / `"spring"` are unchanged) |
| `febio4_exe` | str | - | `"D:\\Program\\FEBioStudio\\bin\\febio4.exe"` (contract §5) |
| `commands` | list[str] | - | one-line reproduction command(s) |
| `measured_vs_modeled` | object | - | two lists: `"measured"` and `"modeled_assumed"`, mirroring sibling `:386-399` |

### 1.4 `scope_limit`: MANDATORY binding to contract §4 条 7

```jsonc
{
  "contract_clause":  "docs/S5_contact_contract.md §4 条 7 (2026-10-06 晚间修订)",
  "scope_summary":    "Wave 4 may start, but every claim must stay inside the trusted
                       footprint of `docs/S5_QUOTABLE.md`.",
  "no_quotable_sigma_ratio": {
    "value":           true,
    "label":           "measured",
    "rule":            "No σ/σ_law ratio may be reported as a converged scalar in this script,
                       including any single-mesh pick and any fixed-ξ pick. Wave 4 outputs
                       σ_contact (absolute MPa, at fixed sink) only, never σ/σ_law.",
    "evidence":        ["docs/S5_QUOTABLE.md §4 (NOT-QUOTABLE #1, #2)",
                        "docs/S5_g7_verify_mesh.md §3 (FAIL, free-ξ)",
                        "docs/S5_g7_verify_fixedxi.md §4 (FAIL, fixed-ξ)",
                        "docs/S5_g7_verify_anchors_mesh.md §4 (FAIL, small-strain)"]
  },
  "citation_policy": {
    "value":           "QUOTABLE-only",
    "label":           "assumed",
    "rule":            "Any number cited in the report body or in the `meta.measured_vs_modeled`
                       list must be a row of `docs/S5_QUOTABLE.md` §0 / §3, with the row's
                       recipe + mesh + ξ quoted together. Cross-recipe comparisons and the
                       retracted large-strain figures are forbidden (QUOTABLE §4 #3, #4).",
    "evidence":        ["docs/S5_QUOTABLE.md §3 (allow-list)", "§4 (deny-list)"]
  },
  "recipe_constraint": {
    "value":           "minimum_version_only",
    "label":           "assumed",
    "rule":            "Wave 4 implements the minimum version (S5接触方案 §3.1, rigid wall +
                       Coulomb μ, no foam solid). The FOAM version (§3.2) is out of scope
                       and is blocked by G7 / G7' (S5_contact_contract.md §4 G7' + 条 6).
                       The minimum version is the only one inside the wave-4 scope limit.",
    "evidence":        ["docs/S5_contact_contract.md §3.1 / §4 条 8",
                        "docs/S5接触方案.md §3.1"]
  },
  "fe_output_constraint": {
    "value":           "absolute_quantities_only",
    "label":           "assumed",
    "rule":            "FE outputs are F_n (kN, contact resultant), σ_contact_max / p95 / p99
                       / mean / median (MPa, footprint-averaged) and penetration (mm).
                       No footprint-based σ/σ_law substitute is computed (S5_QUOTABLE.md §5,
                       S5_g7_verify_metric.md §6 FAIL). The mid-depth σ_mid = F_n / 90 000 mm²
                       may be reported as a fourth column with explicit F_n-rescale note.",
    "evidence":        ["docs/S5_QUOTABLE.md §5", "docs/S5_g7_verify_metric.md §6"]
  }
}
```

The four nested `value` / `rule` / `evidence` triples make the clause machine-checkable: a
downstream linter may reject the run if any of the four `value`s flips to `false`.

### 1.5 `interface_regression`: same shape as sibling `:112-131`

| Key | Type | Label | Unit | Source |
|---|---|---|---|---|
| `default_normal` | list[float, float, float] | measured | - | OpenSim `ground_reaction().ground_normal` |
| `default_is_vertical` | bool | measured | - | sibling `:125` |
| `roll_deg_0_bitwise_equals_default` | bool | measured | - | sibling `:126` |
| `ground_normal_explicit_bitwise_equals_default` | bool | measured | - | sibling `:127` |
| `supination_beta0_is_zero` | bool | measured | - | sibling `:128` |
| `supination_beta0_f_lat_n` | float | measured | N | sibling `:129` |

For the **minimum-version contact** script there is no `supination_load` call (β is out of
scope; see §3). The relevant interface checks are therefore reduced to:

- `plantar_bc_default_unchanged`, bool: `plantar_bc("fixed")` returns the same BC node-set
  as the pre-wave-4 sibling run (sibling `nonvertical_s4.json::cases[axial_baseline_fixed]`,
  baseline `max=213.50 MPa`, `gauge_max=186.76 MPa`).
- `contact_optin_loads_when_default`, bool: `plantar_bc("contact")` produces a `<Contact`
  section in the .feb but `plantar_bc("fixed")` does not (gate G0 isolation).

The reduced fields are listed as `None` if the script never instantiates the contact branch
(see §3 scope: wave 4 exercises the contact branch only when `--plantar-bc=contact` is
explicit; the default path is unchanged).

### 1.6 `axial_regression`: same shape as sibling `:134-163`

```jsonc
{
  "cached_file":        "results\\opensim_fe\\s1_h5_opensim.json",
  "measured_default_run": { /* five scalars, sibling :144-151 */ },
  "cached":             { /* mirror of the above read from disk */ },
  "rel_err":            { /* five relative errors, sibling :152-155 */ },
  "max_rel_err":        float,    // max of rel_err dict, label=measured
  "verdict":            "PASS"|"FAIL"|"SKIPPED"   // threshold 1e-6 (sibling :162)
}
```

This is the **hard gate** called out in `docs\S5_contact_contract.md` §4 G0. The script
re-runs the default axial link (`ground_reaction` → `run_dead_drop` → `subtalar_reaction`)
and compares against the pre-wave-4 cache. Fail ⇒ exit code 1, no JSON / report written
(beyond the overwrite guard, see §0.2).

### 1.7 `scenario_matrix`: axes + run ids (grounded in `docs\非垂直落地扩展方案.md` §2)

```jsonc
{
  "axes": {
    "height_m":     [2.0, 3.0, 3.5, 4.5],                // label=measured; sibling constant HEIGHTS :52
    "plantar_bc":    ["fixed", "contact"],                // fixed = hard gate; contact = the wave
    "contact_kind": ["rigid_wall"],                      // §3.1 only in this wave; FOAM out of scope (§1.4)
    "mu":           [0.0, 0.3, 0.6],                     // Coulomb friction coefficient
    "pad_tilt_deg": [0.0],                               // ② held at axial baseline (§3)
    "foot_mode":    ["two-foot"],                        // ④ held at axial baseline (§3)
    "posture":      ["stiff"],                           // ③ held at axial baseline (§3)
    "cop_offset_mm":[0.0],                               // ⑤ held at axial baseline (§3)
    "angular_momentum": ["zero"],                        // ⑥ held at axial baseline (§3)
    "velocity_dir": ["vertical"],                        // ① held at axial baseline (§3)
    "pad_scenario": ["on_pad"],                          // on_pad=True; rigid baseline = G0 hard gate
    "spring_k_N_per_mm": [null, 100.0]                   // contact vs spring-k reference
  },
  "scope_notes":  { /* per-axis: "in scope for wave 4" / "held at axial baseline (rationale)" */ },
  "runs": [ /* list of run ids; one entry per (height × plantar_bc × contact_kind × mu × pad_tilt × foot_mode × posture × cop_offset × velocity_dir × angular_momentum × pad_scenario × spring_k) */ ]
}
```

§3 below states the per-axis scope explicitly. The `scope_notes` object is required
(proposal) so a downstream linter can audit "wave 4 scope" without re-reading the spec.

### 1.8 `contact_recipe`: the FEBio recipe block (named, reproducible)

Each run under `runs` references one named recipe here. The minimum version has **one**
recipe:

```jsonc
{
  "minimum_rigid_wall_mu0.6": {
    "feb_section_order":  ["Loads", "Boundary", "Contact"],   // contract §5 (Loads → Boundary → Contact)
    "contact_block": {
      "febio_type":       "sliding-frictional",               // FEBio contact formulation
      "master_surface":   "<pad top facet set>",              // real tet facet, not npz quad (contract §5)
      "slave_surface":    "<calcaneus plantar facet set>",
      "mu":               0.6,                                // Coulomb coefficient (QUOTABLE §3 #1)
      "penalty":          1.0,                                // OLD AUGLAG pen=1 (QUOTABLE §3 #5)
      "laugon":           "AUGLAG",
      "tolerance":        0.005,
      "aug_controls":     "AUGLAG_CONVERGENT_CONTROLS",       // contract §4 G7
      "search_radius":    20,
      "n":                1200
    },
    "rigid_wall_block": {
      "febio_type":       "rigid wall",                       // §3.1 (no foam solid)
      "normal_unit":      [0, -1, 0],                         // OpenSim ground frame
      "offset_mm":        0.5,                                // gap, sibling :76
      "blocked_dof":      ["tx", "ty", "tz"]                  // 3-DOF pinned
    }
  }
}
```

The exact recipe is **placeholder** in this design (the implementer fills it from the FE
recipe extracted in `docs\S5_contact_notes.md` + the QUOTABLE §3 #1 row). What the schema
**does** pin down is that the recipe lives in its own named block, has a fixed shape, and
is referenced by run id, so a downstream audit can compare recipes without parsing every
run row.

### 1.9 `runs`: one entry per (id, recipe)

| Field | Type | Unit | Label | Source / formula |
|---|---|---|---|---|
| `id` | str | - | - | synthetic; e.g. `"h2.0_contact_rigid_mu0.6_pad"` |
| `recipe_id` | str | - | - | references `contact_recipe.<name>` |
| `height_m` | float | m | measured | sibling `HEIGHTS` :52 |
| `plantar_bc` | str | - | - | `"fixed"` (hard gate) or `"contact"` (wave) |
| `pad_scenario` | str | - | - | `"on_pad"` / `"hard_surface"` |
| `spring_k_N_per_mm` | float \| null | N/mm | assumed | `null` for `plantar_bc="contact"` and `plantar_bc="fixed"`; `100.0` for `plantar_bc="spring"` reference |
| `febio_rc` | int | - | measured | FEBio 4.13 `rc` (`0` = NORMAL) |
| `febio_end_t` | float | - | measured | last time-step reached |
| `wall_s` | float | s | measured | wall time of the FEBio call |
| `feb_path` | str | - | - | absolute path to `.feb` |
| `xplt_path` | str | - | - | absolute path to `.xplt` |
| `F_n_n` | float \| null | N | measured \| null | contact resultant normal force at last step (= QUOTABLE §3 #1 / #2 source) |
| `F_t_n` | float \| null | N | measured \| null | tangential contact force at last step |
| `friction_cone_ratio_max` | float \| null | - | measured | max over all time-steps of `|F_t| / (μ · F_n)`; gate G4 (must ≤ 1.0) |
| `force_balance_relative_err` | float \| null | - | measured | max over time of `|Σ contact force − applied force| / |applied force|`; gate G3 (must ≤ 0.02) |
| `sigma_contact_max_mpa` | float | MPa | measured | `max(σ_contact)` over all slave-surface elements |
| `sigma_contact_p99_mpa` | float | MPa | measured | 99th-percentile |
| `sigma_contact_p95_mpa` | float | MPa | measured | 95th-percentile |
| `sigma_contact_mean_mpa` | float | MPa | measured | arithmetic mean |
| `sigma_contact_median_mpa` | float | MPa | measured | 50th-percentile |
| `sigma_mid_mpa` | float | MPa | measured | `F_n_n / 90 000 mm²` (QUOTABLE §5 / `S5_g7_verify_metric.md` §6; explicitly noted as `F_n` rescale) |
| `penetration_mm` | float | mm | measured | max penetration into the rigid wall over all time-steps. The trusted small-strain anchor range has penetration ≈ 0 (QUOTABLE §2, ≤ 0.25 mm); do **not** claim a bound at the G7 sink |
| `contact_area_mm2` | float | mm² | measured | footprint of slave elements in contact (last step) |
| `footprint_mm2` | float | mm² | measured | full plantar surface area (for context) |
| `gauge_max_mpa` | float | MPa | measured | gauge-regularised max stress on the calcaneus (sibling `nonvertical_s4.json` shape) |
| `gauge_p95_mpa` | float | MPa | measured | same, p95 |
| `peak_elem` | int | - | measured | element index hosting `gauge_max` |
| `peak_centroid_mm` | list[float, float, float] | mm | measured | centroid of `peak_elem` |
| `at_plantar_rim` | bool | - | measured | sibling `nonvertical_s4.json::cases[*].at_plantar_rim` (A7 indicator) |
| `run_verdict` | `PASS` \| `FAIL` \| `SKIPPED` | (n/a) | (n/a) | computed per-run from G4 (`friction_cone_ratio_max ≤ 1.0`) and G3 (`force_balance_relative_err ≤ 0.02`) |

The fields `sigma_contact_*` and `F_n_n` are explicit **absolute quantities** (QUOTABLE §3
#1 / #2 / §5); **no** `sigma_contact_over_sigma_law` field exists in the schema. This is the
mechanical enforcement of `scope_limit.no_quotable_sigma_ratio` (§1.4). The linter can
grep the JSON for the forbidden key string.

### 1.10 `regression_gates`: G0..G6 verdicts (G7 dropped, out of scope per §1.4)

| Gate | Verdict key | Rule (from `docs\S5_contact_contract.md` §4) | Source in JSON |
|---|---|---|---|
| **G0** | `g0_default_unchanged` | `pytest tests\test_opensim_fe.py tests\test_febio_run.py -q → 19 passed / 1 skipped`; default axial regression rel_err ≤ 1e-6 vs `s1_h5_opensim.json` cache (sibling `:162`) | `axial_regression.verdict` + an external `pytest` exit code (reported separately) |
| **G1** | `g1_pad_stress_reproduced` | `ξ ∈ [0, 0.85]` σ-curve reproduced by FE **to max_rel_err ≤ 5 %** (FOAM only) | **SKIPPED** for the minimum version (no foam solid; cite scope_limit.recipe_constraint). Reported as `"N/A_minimum_version"` with reason |
| **G2** | `g2_limit_self_check` | `μ=0 + pad extreme-stiff → spring k→∞ / rigid-support result` | `runs` filtered by `μ=0` and `pad_scenario=hard_surface` vs `plantar_bc=fixed` baseline; ratio within tolerance |
| **G3** | `g3_contact_conservation` | contact force sum = applied force, ≤ 2 %; penetration bounded | `runs[*].force_balance_relative_err` (script computes from reaction forces vs applied) |
| **G4** | `g4_friction_cone` | any time-step `|F_t| ≤ μ · F_n` | `runs[*].friction_cone_ratio_max` ≤ 1.0 |
| **G5** | `g5_a7_relief` | contact gauge_max ≪ fixed gauge_max | ratio of `runs` filtered by `plantar_bc=contact` vs `plantar_bc=fixed` at the same height |
| **G6** | `g6_sensitivity_monotone_bounded` | μ / penalty / k sweep output monotone/bounded | `sensitivity.*` verdicts (see §1.11) |

Each verdict is `{"verdict": "PASS"|"FAIL"|"SKIPPED"|"N/A_minimum_version",
"max_rel_err": float | null, "evidence": "<which run ids>"}`.

G1 is **`N/A_minimum_version`** because the FOAM version (§3.2 of `S5接触方案.md`) is
out of scope per `scope_limit.recipe_constraint`. G7 is **omitted** because wave 4 only
implements the minimum version; G7' is omitted for the same reason (contract §4 G7' is a
FOAM-version correction). The linter must reject a JSON containing a `g7` key.

### 1.11 `sensitivity`: μ / penalty / k sweeps

```jsonc
{
  "mu_sweep": [ /* one row per filter id, each:
                 {"mu": 0.0..0.9 step,
                  "monotone_in_F_t": bool,         // G6
                  "monotone_in_penetration": bool, // G6
                  "friction_cone_holds": bool,     // G4 across the sweep
                  "evidence_run_ids": [...]} */ ],
  "penalty_sweep": [ /* one row per penalty value, each:
                     {"penalty": 1..50,            // log-spaced
                      "febio_rc_ok": bool,
                      "force_balance_relative_err": float,
                      "evidence_run_ids": [...]} */ ],
  "k_sweep":       [ /* one row per k, only for plantar_bc="spring" reference:
                     {"spring_k_N_per_mm": 10..1e5,
                      "monotone_in_gauge_max": bool,  // comparison with rigid limit
                      "evidence_run_ids": [...]} */ ],
  "gate_g6_verdict": "PASS"|"FAIL"|"SKIPPED",
  "gate_g6_max_rel_err": float | null
}
```

These rows mirror the sibling `nonvertical_s5_supination.py::sensitivity` block.

### 1.12 `assumptions`: measured/modeled/assumed list

A `list[str]` of plain-text assumptions (sibling `nonvertical_s5_supination.py:412-425`),
one per line. Each entry begins with the label tag: `[measured]`, `[modeled]`, or `[assumed]`.
The list must include (at minimum) the following (proposal; the implementer may add
project-specific ones but cannot remove these):

1. `[assumed]` **Minimum version (rigid wall + μ, no foam)** is the only wave-4 configuration
   per `docs\S5接触方案.md` §3.1 + `docs\S5_contact_contract.md` §4 条 8. G1 / G7 / G7'
   are out of scope.
2. `[assumed]` Recipe parameters (`μ`, `penalty`, `laugon`, `aug_controls`, `search_radius`,
   `n`) come from `docs\S5_QUOTABLE.md` §3 #5 (OLD AUGLAG pen=1 reproducibly reproduces the
   historical supination-family small-strain anchors on tet4 (3,3,2) to ≤ 0.1 %; the recipe
   label is carried with any cited number).
3. `[measured]` μ = 0.6 is the QUOTABLE §3 #1 row's setup; carrying it does **not** mean
   it has been independently measured for foot-pad friction. Treat as cited, not re-derived.
4. `[measured]` The two anchor quantities (F_n ≈ 34.4 kN at fixed sink 170 mm,
   σ_contact ≈ 0.574 MPa at the same sink, both mesh-converged < 0.6 %) are reproducible
   under this recipe (`S5_g7_verify_mesh.md` §3, `S5_g7_verify_anchors.md` §4.3).
5. `[assumed]` **No σ/σ_law ratio is quoted.** The script computes σ_contact (absolute
   MPa) only; σ/σ_law would have to come from `pad_stress_pa(ξ)` and the result would
   inherit the mesh-convergence FAILs of §1.4 (cite rows).
6. `[assumed]` Citations in the report body come from `docs\S5_QUOTABLE.md` QUOTABLE rows
   only; cross-recipe comparisons are forbidden (QUOTABLE §4 #3).
7. `[assumed]` Pad geometry: rigid wall only; no foam solid; `pad_stress_pa(ξ)` is a
   material function (QUOTABLE §3 #3), not a FE result.
8. `[measured]` Bone + cartesian + plantar facet come from the THUMS calcaneus mesh used
   by sibling `nonvertical_s4.json` (path: `temp\opensim_fe\thums_calcaneus\calcaneus_r_anatframe.npz`).
9. `[modeled]` The default-regression baseline (`plantar_bc="fixed"`) reproduces
   `nonvertical_s4.json::cases[axial_baseline_fixed]` (`max=213.50 MPa`,
   `gauge_max=186.76 MPa`) bit-for-bit; this is G0.
10. `[measured]` All FEBio calls go through `febio_run.run_febio` (contract §5; never
    `subprocess`). `.feb` section order is `Loads → Boundary → Contact` (contract §5).
    The contact surface is a real tet facet, not an npz quadrilateral (contract §5).
11. `[assumed]` Hard gate G0 is the only verification that the default path is
    byte-for-byte preserved. Anything that touches a non-`"fixed"` branch is opt-in and
    never runs by default.

The list above is the proposal; the implementer may add entries, and may split compound
assumptions into multiple, but cannot drop entries 1, 5, 6, or 11 (those are bound to
contract clauses and to §1.4 of this spec).

### 1.13 `repro`: commands + files read

```jsonc
{
  "commands":     [ /* single-line shell strings */ ],
  "files_read":   [ /* absolute paths, sibling pattern */ ],
  "files_written":[ /* absolute paths of JSON + MD; both must match args.out_* */ ],
  "venv":         ".venv\\Scripts\\python.exe",
  "env":          {"PYTHONPATH": "src"}
}
```

`files_read` includes the cached axial (`s1_h5_opensim.json`), the THUMS calcaneus mesh,
and any other inputs (e.g. `nonvertical_s4.json` for the G0 baseline), but **not** the
JSON this script writes, which goes to `files_written`.

### 1.14 Complete example (one run only, for shape)

```jsonc
{
  "meta": {
    "generated_at":   "2026-10-06 12:00:00",
    "study":          "Nonvertical S5-contact, FEBio plantar-pad contact + friction (wave 4, opt-in)",
    "plan_doc":       "docs/非垂直落地扩展方案.md §3 L4 / §4 S5 · docs/S5接触方案.md §3",
    "contract_ref":   "docs/S5_contact_contract.md §4 (acceptance G0..G7) + §4 条 7 (scope limit)",
    "quotable_ref":   "docs/S5_QUOTABLE.md",
    "units":          "mm-N-MPa-s",
    "mass_kg":        {"value": 75.337, "label": "measured", "unit": "kg",
                       "source": "Rajagopal 2015 (propagated from nonvertical_s1.json)"},
    "recipe_id":      "minimum_rigid_wall_mu0.6",
    "plantar_bc_kind": "contact",
    "febio4_exe":     "D:\\Program\\FEBioStudio\\bin\\febio4.exe",
    "commands":       [".venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact.py"],
    "measured_vs_modeled": {
      "measured": [
        "OpenSim GRF + subtalar reaction at 5 m (default link, hard gate)",
        "FEBio 4.13 STATIC solution under script-instantiated run_febio (recipe above)",
        "contact resultant F_n, F_t, friction_cone_ratio_max (read from .xplt post-processing)",
        "calcaneus mesh + plantar facet set from THUMS AM50 (sibling path)",
        "OpenSim ground frame normal unit"
      ],
      "modeled_assumed": [
        "Coulomb friction cone |F_t| ≤ μ·F_n (FEBio default sliding-frictional)",
        "Rigid wall = analytic flat surface, no foam solid (S5接触方案 §3.1)",
        "μ = 0.6 carried from QUOTABLE §3 #1 (treated as cited setup, not re-derived)",
        "G1/G7/G7' out of scope (scope_limit.recipe_constraint)"
      ]
    }
  },
  "scope_limit": {
    "contract_clause":        "docs/S5_contact_contract.md §4 条 7 (2026-10-06 晚间修订)",
    "scope_summary":          "Wave 4 may start, but every claim must stay inside the trusted footprint of docs/S5_QUOTABLE.md.",
    "no_quotable_sigma_ratio":{"value": true, "label": "measured",
                               "rule": "No σ/σ_law ratio is reported; σ_contact (absolute MPa) only.",
                               "evidence": ["docs/S5_QUOTABLE.md §4 #1"]},
    "citation_policy":        {"value": "QUOTABLE-only", "label": "assumed",
                               "rule": "Only QUOTABLE rows may be cited; cross-recipe forbidden.",
                               "evidence": ["docs/S5_QUOTABLE.md §3 / §4"]},
    "recipe_constraint":      {"value": "minimum_version_only", "label": "assumed",
                               "rule": "Rigid wall + μ only; no foam solid.",
                               "evidence": ["docs/S5_contact_contract.md §3.1 / §4 条 8"]},
    "fe_output_constraint":   {"value": "absolute_quantities_only", "label": "assumed",
                               "rule": "σ_contact / F_n / penetration only; no σ/σ_law.",
                               "evidence": ["docs/S5_QUOTABLE.md §5"]}
  },
  "interface_regression": {
    "default_normal":                          [0.0, 1.0, 0.0],
    "default_is_vertical":                     true,
    "roll_deg_0_bitwise_equals_default":       true,
    "ground_normal_explicit_bitwise_equals_default": true,
    "supination_beta0_is_zero":                 null,
    "supination_beta0_f_lat_n":                null,
    "plantar_bc_default_unchanged":            true,
    "contact_optin_loads_when_default":        false
  },
  "axial_regression": {
    "cached_file":           "results\\opensim_fe\\s1_h5_opensim.json",
    "measured_default_run":  {"impulse_ns": 789.328, "impulse_rel_err": 0.058,
                              "grf_peak_n": 53065.0, "subtalar_peak_vertical_n": 26383.1,
                              "subtalar_peak_force_n": 26383.1, "subtalar_peak_moment_nm": 20.03},
    "cached":                {"impulse_ns": 789.328, "impulse_rel_err": 0.058,
                              "grf_peak_n": 53065.0, "subtalar_peak_vertical_n": 26383.1,
                              "subtalar_peak_force_n": 26383.1, "subtalar_peak_moment_nm": 20.03},
    "rel_err":               {"impulse_ns": 0.0, "impulse_rel_err": 0.0, "grf_peak_n": 0.0,
                              "subtalar_peak_vertical_n": 0.0, "subtalar_peak_force_n": 0.0,
                              "subtalar_peak_moment_nm": 0.0},
    "max_rel_err":           0.0,
    "verdict":               "PASS"
  },
  "scenario_matrix": {
    "axes": {
      "height_m":          [2.0, 3.0, 3.5, 4.5],
      "plantar_bc":         ["fixed", "contact"],
      "contact_kind":      ["rigid_wall"],
      "mu":                [0.0, 0.3, 0.6],
      "pad_tilt_deg":      [0.0],
      "foot_mode":         ["two-foot"],
      "posture":           ["stiff"],
      "cop_offset_mm":     [0.0],
      "angular_momentum":  ["zero"],
      "velocity_dir":      ["vertical"],
      "pad_scenario":      ["on_pad"],
      "spring_k_N_per_mm": [null, 100.0]
    },
    "scope_notes": {
      "plantar_bc.contact": "in scope for wave 4 (the new opt-in capability)",
      "plantar_bc.fixed":   "in scope (G0 hard gate + G5 A7 reference)",
      "plantar_bc.spring":  "in scope as reference (G2 limit self-check: μ=0 + rigid wall → spring k→∞)",
      "mu":                 "in scope (G4 + G6 μ_sweep)",
      "pad_tilt_deg":       "held at 0.0 (② ground tilt is upstream of wave 4; see §3)",
      "foot_mode":          "held at two-foot (④ single-foot is upstream; see §3)",
      "posture":            "held at stiff (③ posture is upstream; see §3)",
      "cop_offset_mm":      "held at 0.0 (⑤ CoP / moment is upstream; see §3)",
      "angular_momentum":   "held at zero (⑥ ω is upstream; see §3)",
      "velocity_dir":       "held at vertical (① velocity direction is upstream; see §3)",
      "contact_kind.foam":  "out of scope (FOAM version blocked by G7 / G7')"
    },
    "runs": [
      "h2.0_fixed_pad_spring100",                  // G0 + G5 reference
      "h3.0_fixed_pad_spring100",
      "h3.5_fixed_pad_spring100",
      "h4.5_fixed_pad_spring100",
      "h2.0_contact_rigid_mu0.0_pad",              // G2 (μ=0) + G6
      "h2.0_contact_rigid_mu0.3_pad",              // G4 + G6
      "h2.0_contact_rigid_mu0.6_pad",              // G4 + G5 + G3 + QUOTABLE anchor
      "h3.0_contact_rigid_mu0.6_pad",
      "h3.5_contact_rigid_mu0.6_pad",
      "h4.5_contact_rigid_mu0.6_pad"
    ]
  },
  "contact_recipe": {
    "minimum_rigid_wall_mu0.6": {
      "feb_section_order": ["Loads", "Boundary", "Contact"],
      "contact_block": {
        "febio_type":      "sliding-frictional",
        "master_surface":  "<pad top facet set, real tet facet>",
        "slave_surface":   "<calcaneus plantar facet set>",
        "mu":              {"value": 0.6, "label": "measured", "unit": "dimensionless",
                            "source": "docs/S5_QUOTABLE.md §3 #1 (carried setup, not re-derived)"},
        "penalty":         {"value": 1.0, "label": "assumed", "unit": "dimensionless",
                            "source": "QUOTABLE §3 #5 (OLD AUGLAG pen=1 recipe)"},
        "laugon":          "AUGLAG",
        "tolerance":       0.005,
        "aug_controls":    "AUGLAG_CONVERGENT_CONTROLS",
        "search_radius":   20,
        "n":               1200
      },
      "rigid_wall_block": {
        "febio_type":     "rigid wall",
        "normal_unit":    [0, -1, 0],
        "offset_mm":      {"value": 0.5, "label": "assumed", "unit": "mm",
                           "source": "sibling QUOTABLE §3 #1 setup"},
        "blocked_dof":    ["tx", "ty", "tz"]
      }
    }
  },
  "runs": [
    {
      "id": "h2.0_contact_rigid_mu0.6_pad",
      "recipe_id":          "minimum_rigid_wall_mu0.6",
      "height_m":           {"value": 2.0, "label": "measured", "unit": "m"},
      "plantar_bc":         "contact",
      "pad_scenario":       "on_pad",
      "spring_k_N_per_mm":  null,
      "febio_rc":           0,
      "febio_end_t":       1.0,
      "wall_s":            {"value": 0.32, "label": "measured", "unit": "s"},
      "feb_path":          "D:\\Project\\climbing_fall_analysis\\temp\\opensim_fe\\s5_contact\\h2.0_contact_mu0.6.feb",
      "xplt_path":         "D:\\Project\\climbing_fall_analysis\\temp\\opensim_fe\\s5_contact\\h2.0_contact_mu0.6.xplt",
      "F_n_n":             {"value": 34395.0, "label": "measured", "unit": "N",
                            "source": "FEBio .xplt contact pressure integral (QUOTABLE §3 #1 row)"},
      "F_t_n":             {"value": 4120.0, "label": "measured", "unit": "N"},
      "friction_cone_ratio_max":{"value": 0.20, "label": "measured", "unit": "dimensionless"},
      "force_balance_relative_err": {"value": 0.004, "label": "measured", "unit": "dimensionless"},
      "sigma_contact_max_mpa": {"value": 0.59, "label": "measured", "unit": "MPa"},
      "sigma_contact_p99_mpa":  {"value": 0.58, "label": "measured", "unit": "MPa"},
      "sigma_contact_p95_mpa":  {"value": 0.58, "label": "measured", "unit": "MPa"},
      "sigma_contact_mean_mpa": {"value": 0.57, "label": "measured", "unit": "MPa"},
      "sigma_contact_median_mpa":{"value": 0.57,"label": "measured", "unit": "MPa"},
      "sigma_mid_mpa":          {"value": 0.573, "label": "measured", "unit": "MPa",
                                 "source": "F_n / 90 000 mm² (QUOTABLE §5 / S5_g7_verify_metric.md §6, F_n rescale)"},
      "penetration_mm":         {"value": 0.10, "label": "measured", "unit": "mm"},
      "contact_area_mm2":       {"value": 60000.0, "label": "measured", "unit": "mm²"},
      "footprint_mm2":          {"value": 60000.0, "label": "measured", "unit": "mm²"},
      "gauge_max_mpa":          {"value": 175.0, "label": "measured", "unit": "MPa"},
      "gauge_p95_mpa":          {"value": 122.0, "label": "measured", "unit": "MPa"},
      "peak_elem":              198,
      "peak_centroid_mm":       [12.5, 7.5, 16.8],
      "at_plantar_rim":         false,
      "run_verdict":            "PASS"
    }
  ],
  "regression_gates": {
    "g0_default_unchanged":  {"verdict": "PASS", "max_rel_err": 0.0,
                              "evidence": "axial_regression.verdict == PASS"},
    "g1_pad_stress_reproduced":{"verdict": "N/A_minimum_version", "max_rel_err": null,
                                "evidence": "scope_limit.recipe_constraint"},
    "g2_limit_self_check":   {"verdict": "PASS", "max_rel_err": 0.012,
                              "evidence": ["runs[μ=0]"]},
    "g3_contact_conservation":{"verdict": "PASS", "max_rel_err": 0.005,
                               "evidence": ["runs[*]"]},
    "g4_friction_cone":      {"verdict": "PASS", "max_rel_err": null,
                              "evidence": ["runs[*].friction_cone_ratio_max ≤ 1.0"]},
    "g5_a7_relief":          {"verdict": "PASS", "max_rel_err": null,
                              "evidence": ["runs[plantar_bc=contact] vs [plantar_bc=fixed]"]},
    "g6_sensitivity_monotone_bounded": {"verdict": "PASS", "max_rel_err": null,
                                        "evidence": ["sensitivity.*"]}
  },
  "sensitivity": {
    "mu_sweep": [
      {"mu": 0.0, "monotone_in_F_t": true,  "monotone_in_penetration": true,  "friction_cone_holds": true,  "evidence_run_ids": ["h2.0_contact_rigid_mu0.0_pad"]},
      {"mu": 0.3, "monotone_in_F_t": true,  "monotone_in_penetration": true,  "friction_cone_holds": true,  "evidence_run_ids": ["h2.0_contact_rigid_mu0.3_pad"]},
      {"mu": 0.6, "monotone_in_F_t": true,  "monotone_in_penetration": true,  "friction_cone_holds": true,  "evidence_run_ids": ["h2.0_contact_rigid_mu0.6_pad"]}
    ],
    "penalty_sweep": [],
    "k_sweep":       [],
    "gate_g6_verdict": "PASS",
    "gate_g6_max_rel_err": null
  },
  "assumptions": [
    "[assumed] Minimum version (rigid wall + μ, no foam) per docs/S5接触方案.md §3.1 + docs/S5_contact_contract.md §4 条 8. G1 / G7 / G7' out of scope.",
    "[assumed] Recipe parameters from docs/S5_QUOTABLE.md §3 #5 (OLD AUGLAG pen=1).",
    "[measured] μ = 0.6 carried from QUOTABLE §3 #1 (cited setup, not re-derived).",
    "[measured] F_n ≈ 34.4 kN at fixed sink 170 mm under this recipe (S5_g7_verify_mesh.md §3).",
    "[assumed] No σ/σ_law ratio is quoted; σ_contact (absolute MPa) only.",
    "[assumed] Citations from QUOTABLE only; cross-recipe forbidden.",
    "[assumed] Pad geometry = rigid wall only; pad_stress_pa(ξ) is a material function.",
    "[measured] Bone + plantar facet from THUMS calcaneus mesh (sibling path).",
    "[modeled] Default (plantar_bc=fixed) reproduces nonvertical_s4.json::cases[axial_baseline_fixed] bit-for-bit (G0).",
    "[measured] All FEBio calls via febio_run.run_febio (contract §5).",
    "[assumed] Hard gate G0 is the only default-path verification."
  ],
  "repro": {
    "commands":      [".venv\\Scripts\\python.exe scripts\\opensim_fe\\nonvertical_s5_contact.py"],
    "files_read":    ["results\\opensim_fe\\s1_h5_opensim.json",
                      "results\\opensim_fe\\nonvertical_s4.json",
                      "temp\\opensim_fe\\thums_calcaneus\\calcaneus_r_anatframe.npz"],
    "files_written": ["results\\opensim_fe\\nonvertical_s5_contact.json",
                      "results\\opensim_fe\\NONVERTICAL_S5_CONTACT_REPORT.md"],
    "venv":          ".venv\\Scripts\\python.exe",
    "env":           {"PYTHONPATH": "src"}
  }
}
```

(The example carries one `runs[*]` row for compactness. Production runs contain the full
10-row matrix under §1.7, see the example ids in `scenario_matrix.runs`.)

---

## 2. Report structure: `NONVERTICAL_S5_CONTACT_REPORT.md`

The shape mirrors `NONVERTICAL_S5_SUPINATION_REPORT.md` (sibling `:456-723`) but
reorganises the section order around FE contact rather than 1D analytic. The report has
thirteen numbered sections (`## 0.` through `## 12.`):

- **`## 0. 一句话结论`**: three to five numbered bullets, in the sibling style
  (`:468-486`). Bullet 1 = G0 hard gate verdict + max_rel_err. Bullet 2 = scope summary
  (the §1.4 `scope_limit` block in one sentence). Bullet 3 = G5 A7 reduction (the headline
  numerical claim). Bullet 4 = G4 friction cone compliance across the sweep. Bullet 5 =
  "see §10 for honesty boundary".

- **`## 1. API（opt-in contact kind）`**: table mirroring sibling `:499-507`. Adds one
  row: `febio_plantar_bc_contact`, new `plantar_bc="contact"` opt-in that triggers a
  `<Contact>` block in the `.feb` and the rigid-wall minimum-version contact recipe. The
  default `plantar_bc="fixed"` is unchanged (G0 hard gate).

- **`## 2. 默认回归（硬门槛）`**: sibling `:515-538`. Two sub-tables:
  `interface_regression` table (sibling format) + `axial_regression` table with the
  `s1_h5_opensim.json` cache comparison (5 rows, max_rel_err, verdict).

- **`## 3. 接触配方（FEBio recipe block）`**: names the recipe block from §1.8, prints
  the `feb_section_order`, the `contact_block` fields with units + labels, and the
  `rigid_wall_block`. Cites `docs\S5_QUOTABLE.md` §3 #1 / #5 row provenance. Explicit
  reminder that the FOAM recipe (G1 / G7) is **not** in scope.

- **`## 4. FE 接触输出（per-run 摘要）`**: one row per run id from
  `scenario_matrix.runs`. Columns: id, `F_n` (kN), `F_t` (N), `friction_cone_ratio_max`,
  `σ_contact_max/p95/median` (MPa), `σ_mid` (MPa, with `F_n` rescale note), penetration
  (mm), `contact_area_mm²`, `gauge_max/p95` (MPa), `febio_rc`, `wall_s`. Every cell
  carries its label tag inline (e.g. `[measured]`) so a reader can audit.

- **`## 5. 场景表（hard gate + sweep）`**: one table per `pad_scenario` (`on_pad`,
  `hard_surface`). Heights × `plantar_bc` × `μ`. Each row links to the corresponding
  `runs[*].id`.

- **`## 6. 验收 G0..G6 逐条`**: verbatim from `docs\S5_contact_contract.md` §4 table,
  one row per gate. Each row reads from `regression_gates` (§1.10), includes the verdict
  (PASS / FAIL / SKIPPED / N/A_minimum_version), the max_rel_err, and the evidence run
  ids. Footnotes: G1 is `N/A_minimum_version` per §1.4; G7 is omitted.

- **`## 7. 摩擦锥 + 接触守恒`**: G3 + G4 detail. G3: print the per-run `force_balance_relative_err`
  against applied load. G4: print the max `friction_cone_ratio_max` per run (must ≤ 1.0).

- **`## 8. 敏感性（μ / penalty / k）`**: G6 detail. Three sub-tables (one per sweep)
  from `sensitivity` (§1.11). Each cell shows the per-filter verdict, with the evidence
  run ids. Explicit "monotone/bounded" definitions.

- **`## 9. A7 消失验证（fixed vs contact gauge 对照）`**: G5 detail. Two columns:
  `plantar_bc="fixed"` baseline vs `plantar_bc="contact"` minimum-version. Compute the
  gauge_max drop %; cite `nonvertical_s4.json::cases[axial_baseline_fixed]`
  (`max=213.50 MPa`, `gauge_max=186.76 MPa`) as the baseline; the contact row is a
  measured value from `runs`. **Do not** mix in any FOAM-version or σ/σ_law ratio.

- **`## 10. 假设（ASSUMPTIONS）与诚实边界`**: numbered list from `assumptions`
  (§1.12), plus a short **"未解决的风险与下一步"** subsection. Add a "structural boundary"
  paragraph (sibling `:692-699` style) noting that wave 4 is the **minimum version**, that
  G1 / G7 / G7' are out of scope, and that this script is **only** the FE-contact
  capability: the contact-pressure distribution is FE-computed, but the **absolute value**
  of σ_contact inherits the QUOTABLE §3 #1 / #2 mesh-converged bounds (no further
  refinement).

- **`## 11. 复现命令`**: sibling `:705-717`. Two-line invocation: one for default (G0
  only), one for the full wave (contact branch).

- **`## 12. 产物`**: sibling `:719-725`. Two paths: the JSON and the MD, both pointing
  to the §4 basename slots.

Mapping to `docs\S5_contact_contract.md` §6 (the mandatory per-subagent report format).
The report must cover all five items of §6:

| contract §6 item | report section |
|---|---|
| 1. files added / changed (absolute paths) | `## 12. 产物` |
| 2. commands run, exit codes, key output | `## 11. 复现命令` plus `runs[*].febio_rc` / `febio_end_t` / `wall_s` in `## 4` |
| 3. G0..G6 per-gate verdict + numbers | `## 6. 验收 G0..G6 逐条` |
| 4. unresolved risks and next steps | `## 10.` §"未解决的风险与下一步" |
| 5. one reproducible shortest command | `## 11. 复现命令` (default one-liner) |

The `>`-prefixed meta lines at the top of the report follow the sibling style (`:458-463`):
generated_at, repo path, contract refs, and one explicit line: "本脚本只报绝对量，**任何
σ/σ_law 比值一律不得作为'收敛'单一值引用**（`docs/S5_QUOTABLE.md` §4 #1、#2）。

---

## 3. Scenario matrix: grounded in the 6-DoF decomposition

The 6-DoF decomposition (`docs\非垂直落地扩展方案.md` §2) gives the axes; wave 4 covers
a small, opt-in subset and holds the rest at the axial baseline. The **scope_notes**
block in `scenario_matrix` (§1.7) is the machine-checkable record of this table.

| # | DoF | In scope for wave 4? | Why |
|---|---|---|---|
| ① | **撞击速度方向** | ❌ held at axial baseline (`velocity_dir=["vertical"]`) | S1 (非垂直落地扩展方案 §4) covers ①. Re-running them in wave 4 would mean re-running the OpenSim chain, which is **out of scope** for an FE-contact script. `axial_regression` (G0) is the only place the velocity-direction link is touched, and only at the default |
| ② | **地面 / 垫面倾角** | ❌ held at axial baseline (`pad_tilt_deg=[0.0]`) | ② is L1 (pad contact geometry) and is upstream of the FE contact. Wave 4 puts the pad tilt at 0 and verifies the contact recipe at the axial baseline. The sensitivity sweep on `pad_tilt_deg` is **out of scope**; if needed later, a follow-up wave is preferred to keep this one small |
| ③ | **姿势** | ❌ held at axial baseline (`posture=["stiff"]`) | S3 (方案 §4) covers ③. Stiff legs (锁死腿) are the G0 baseline and the only posture the FE contact runs against in this wave |
| ④ | **左右不对称（单脚）** | ❌ held at axial baseline (`foot_mode=["two-foot"]`) | S2 (方案 §4) covers ④. The FE contact in wave 4 keeps the split=(0.5, 0.5) default. Single-foot is a follow-up if the multi-foot contact recipe proves non-trivial |
| ⑤ | **着力点 / CoP 偏移** | ❌ held at axial baseline (`cop_offset_mm=[0.0]`) | S4 (方案 §4) covers ⑤ via the 3D wrench (already in scope for FE via `use_rigid=True` + `RigidForceLoad` + `RigidMomentLoad`, sibling `nonvertical_s4.json::cases[rigid_tilt_*]`). Wave 4 does **not** exercise CoP offset; the FE contact recipe applies to the full-axial case. Adding CoP would re-test the force-balance and force-torque recipe already validated in s4 |
| ⑥ | **撞击角动量** | ❌ held at axial baseline (`angular_momentum=["zero"]`) | ⑥ is upstream (initial conditions in `run_dead_drop`). Same reason as ①: re-running OpenSim with non-zero ω is out of scope |

What **is** in scope:

- **`plantar_bc="contact"`** (the new opt-in capability wave 3 implements). Wave 4 runs
  the FE contact through scenario variations on top of this capability. Wave 4
  treats `plantar_bc="contact"` as opt-in (default `plantar_bc="fixed"` unchanged, G0 hard
  gate). The contact branch is a new `.feb` `<Contact>` section (sibling
  `nonvertical_s4.json::cases[*].plantar_bc` already enumerates `"fixed"` and `"spring"`;
  `"contact"` is the new value).
- **Coulomb μ sweep**: μ ∈ {0.0, 0.3, 0.6}. This is the G4 (friction cone) and G6
  (sensitivity) coverage. μ = 0.0 is the G2 limit self-check (no friction ⇒ the contact
  reduces to a rigid wall + zero shear).
- **Height sweep**: h ∈ {2.0, 3.0, 3.5, 4.5} m (sibling `HEIGHTS` constant, :52). The
  height sweep is the same as the sibling for cross-report comparability.
- **`plantar_bc="spring"` reference** at k = 100 N/mm. Used as the G2 limit self-check
  baseline (sibling `nonvertical_s4.json::cases[axial_baseline_spring]`,
  `max=151.90 MPa`, `gauge_max=151.90 MPa`). The script compares `μ=0 + rigid wall` vs
  `plantar_bc=fixed` to verify the limit.

**Axial regression hard gate (required row, contract §4 G0)**: `plantar_bc="fixed"` at
h ∈ {2.0, 3.0, 3.5, 4.5} m, `pad_scenario="on_pad"`, `spring_k_N_per_mm=null`. This row
**must** reproduce `nonvertical_s4.json::cases[axial_baseline_fixed]` (`max=213.50 MPa`,
`gauge_max=186.76 MPa`) to the sibling's 1e-6 threshold (sibling `:162`). Any other row
in the matrix is irrelevant if this row fails.

The total run count from the axis products in §1.7 is `4 × (1 fixed + 1 spring + 3 contact
μs) = 20 runs` for the minimum matrix. (The 1 foam `contact_kind` branch is **out of
scope** per `scope_limit.recipe_constraint` and not enumerated.)

---

## 4. Filename / collision decision

### 4.1 The collision

| Slot | Owner | Status |
|---|---|---|
| `results\opensim_fe\nonvertical_s5.json` | `scripts\opensim_fe\nonvertical_s5_supination.py` (sibling, 731 lines) | **occupied 2026-10-04; read-only in derivative only** |
| `results\opensim_fe\nonvertical_s5_*.json` | (vacant) | the wave-4 slot |
| `results\opensim_fe\nonvertical_s6_*.json` | (vacant) | future |
| `results\opensim_fe\ankle_external_rotation.json` | `scripts\opensim_fe\nonvertical_s6_external_rotation.py` (s6) | occupied 2026-10-05 |

The s5-supination script wrote to `nonvertical_s5.json` (line `:66`) and to
`NONVERTICAL_S5_SUPINATION_REPORT.md` (line `:67`). The s6 script (external rotation)
**split off** to `ankle_external_rotation.json` + `ANKLE_EXTERNAL_ROTATION_REPORT.md`
rather than reusing `nonvertical_s6.json`, on exactly this collision principle: the
s-family scripts share the `nonvertical_s*.json` slot but the moment a variant introduces
a different study name, it gets its own descriptive name. (See
`results\opensim_fe\ankle_external_rotation.json:4` `"study": "Phase-S6 踝外旋 / SER …"`,
which is a different study from `results\opensim_fe\nonvertical_s5.json:4` `"study":
"Phase-S5 踝旋后 / 内翻（冠状面）"`.)

### 4.2 Proposed basename

**JSON**: `results\opensim_fe\nonvertical_s5_contact.json`

**Report**: `results\opensim_fe\NONVERTICAL_S5_CONTACT_REPORT.md` (fixed, no choice).

### 4.3 Justification

1. **Family slot**: the proposed basename sits inside the `nonvertical_s5_*` family. The
   supination script currently owns the bare `nonvertical_s5.json`; if it were renamed
   for consistency it would be `nonvertical_s5_supination.json`, which is exactly the
   sibling of the proposed `nonvertical_s5_contact.json`. The `s5_contact` suffix keeps
   the "this is also an s5 variant" semantic.
2. **No collision**: `nonvertical_s5_contact.json` is not present in
   `results\opensim_fe\` (a prior ls confirmed only `nonvertical_s5.json` exists in the
   `nonvertical_s5_*` family). The overwrite guard (§0.2) would reject any collision
   anyway, but the basename avoids it.
3. **Symmetric with the report name**: the report name is fixed by the contract /
   spec as `NONVERTICAL_S5_CONTACT_REPORT.md`. The JSON basename mirrors that name;
   `nonvertical_s5_contact.json` is the lower-case / dot-ext form of the same token.
4. **Consistent with s6 precedent**: s6 split off to `ankle_external_rotation.json` (a
   descriptive basename) rather than `nonvertical_s6.json`. Wave 4 follows the same
   precedent by using `nonvertical_s5_contact.json` (a descriptive basename within the
   s5 family) rather than the more generic `nonvertical_s5b.json` (which would lose the
   study-identification).
5. **Linter-friendly**: the basename encodes the variant (`contact`), so a downstream
   linter can match `nonvertical_s5_contact.json` against the `NONVERTICAL_S5_CONTACT_REPORT.md`
   report name and against the `nonvertical_s5_*.json` family glob.

### 4.4 Alternative basenames considered and rejected

| Alt basename | Why rejected |
|---|---|
| `nonvertical_s5b.json` | loses study identification; not consistent with s6 precedent |
| `contact_fe.json` | breaks the `nonvertical_s*.json` family glob that other tools may rely on |
| `febio_contact.json` | even more generic; not in any family |
| `s5_contact_wave4.json` | over-encodes the wave; the waveless baseline lives in the QUOTABLE §0 / §3 scope table |
| `nonvertical_s5.json` (overwrite) | **forbidden**, occupied by supination; the overwrite guard (§0.2) returns exit code 2 |

---

## 5. Open questions for the user

The implementer cannot decide these alone. Each is a fork the user must resolve before
`scripts\opensim_fe\nonvertical_s5_contact.py` can be written.

1. **Recipe choice**, `docs\S5_QUOTABLE.md` §0 #1 cites the NEW PENALTY 0.1 recipe for
   `F_n ≈ 34.4 kN` (mesh-converged < 0.6 % at fixed sink 170 mm), but §3 #5 cites the
   OLD AUGLAG pen=1 recipe for the historical supination-family anchor reproduction (≤ 0.1 %
   on tet4 (3,3,2)). The two recipes are **not interchangeable** (QUOTABLE §1.1 / §4 #3). Which
   recipe does wave 4 use? (The §1.3 example above picks OLD AUGLAG pen=1 for G6 μ
   sweep comparability with `nonvertical_s5_supination.py`'s audit history. The user
   may prefer NEW PENALTY 0.1 to align with QUOTABLE §0 #1.)

2. **Call FEBio or wave 3 walk-through?**, the §1.9 `runs[*]` schema assumes the wave-4
   script **calls** FEBio via `febio_run.run_febio` (contract §5). Alternatively, wave 4
   could be a **non-FEBio** "FE-contact-prepared" script that the wave-3-capable
   `plantar_bc="contact"` opt-in invokes (i.e. wave 4 is a scenario-shape audit, not a
   runner). The latter avoids re-running FEBio in this wave and keeps the run time small
   (which matters given the 30-minute ceiling). Which is the intent?

3. **μ sweep granularity**, §1.7 proposes μ ∈ {0.0, 0.3, 0.6}. Should it be a single μ
   value (one row) or a finer grid (e.g. μ ∈ {0.0, 0.2, 0.4, 0.6, 0.8})? The QUOTABLE
   §3 #1 row uses μ = 0.6 only; the G4 + G6 sweep intent is broader. Three coarser points
   or five finer points?

4. **G1 / G7 placement in the report**, §2 §6 reports G1 as `N/A_minimum_version` per
   §1.4. Should G7 be reported as "out of scope (FOAM version only)" or omitted entirely?
   The contract §4 lists G7 in the acceptance table; omitting it loses a paper-reference
   trail. Reporting it as `N/A_minimum_version` keeps the trail. Which?

5. **9-bone coverage marker**, the FE contact is on the calcaneus only. How do the
   other 8 bones get marked in the report? `nonvertical_s5_supination.py` reports on the
   fibula only (sibling `NONVERTICAL_S5_SUPINATION_REPORT.md` is a fibula-only
   single-bone study). Should wave 4 follow the fibula-only pattern (calcaneus-only) or
   attempt a 9-bone sweep? The latter would significantly expand the run matrix.

6. **Where does the .feb / .xplt land?**, §1.9 `feb_path` / `xplt_path` point to
   `temp\opensim_fe\s5_contact\…`. `nonvertical_s4.py` puts its data under
   `temp\opensim_fe\s4_nonvertical\…`. The s5 sibling (supination) writes nothing to disk
   because it's a non-FEBio script. Wave 4 is the first s5-family script that touches
   FEBio, so the directory naming convention has no prior in the s5 family. The §1.9
   proposal is `temp\opensim_fe\s5_contact\…`. Acceptable, or should it be `s5_wave4`,
   `s5c`, or another tag?

7. **Hard-surface label**, `pad_scenario="hard_surface"` is in the §1.7 axis but not
   enumerated in the example `runs[*]` (which only lists `pad_scenario="on_pad"`). The
   hard-surface runs are the G2 limit self-check + the G5 A7 reference; they are
   in-scope. Confirm?

8. **σ_contact_max vs σ_contact_p99 / p95 / median**: §1.9 lists five σ_contact
   percentiles (max / p99 / p95 / mean / median). The QUOTABLE §0 / §3 cite only the
   full-footprint average (`σ_contact = F_n / 60 000 mm² ≈ 0.574 MPa`); the finer
   percentiles are diagnostic-only. Should the JSON carry all five, or only `max / p95 /
   F_n`-derived average (closer to QUOTABLE)? More columns = more diagnostic; fewer =
   closer to QUOTABLE citation discipline.

9. **`friction_cone_ratio_max` over time or peak?**, §1.9 lists `friction_cone_ratio_max` as a
   scalar (the max over time). The G4 gate is "any time-step |F_t| ≤ μ·F_n", so the
   correct diagnostic is the **time-maximum**. The JSON carries the max. Acceptable, or
   should it carry a time series (much larger)?

10. **Hard gate on `temperature`?**, the QUOTABLE §3 #1 row assumes standard conditions;
    no thermal coupling. Wave 4 should follow the same. Confirm? (Default = no thermal.)

The implementer must not start work until items 1-4 are resolved; items 5-10 are
scope-shaping and the user may delegate them.

### 5.1 Resolutions — coordinator decisions, 2026-10-06 (ALL 10 RESOLVED)

> 这些是**技术分叉**（配方 / 是否跑 FEBio / 扫描密度 / 门禁记法），由协调层拍板并记录理由；
> 非产品级决策。实现者**按下表执行，不得自行改判**。

| # | 决定 | 理由 |
|---|---|---|
| **1** | **用 NEW `PENALTY 0.1`**（`laugon="PENALTY"`, `penalty=0.1`, `search_radius=20`, `node_reloc=0`, `tolerance=0.005`, `fric_coeff=0.6`） | ① OLD AUGLAG 在**接触闭合处必死**，波4 要跑真实接触就过不去；② QUOTABLE §0 #1 的 `F_n ≈ 34.4 kN` **就是这个配方**；③ 跨配方对比是禁令 —— 全脚本**只允许一个配方**。原文举例选 OLD 是为了跟旋后线的审计史对比，但旋后线的数字**不在 QUOTABLE 集内**，不构成理由。 |
| **2** | **真跑 FEBio**（`climbing.coupling.febio_run.run_febio`，契约 §5） | 契约 §6 要求"命令 + 退出码 + 关键输出"与"G0…G6 逐条带数字判定"。**非 FEBio 的场景审计给不出 G3（力守恒）/ G4（摩擦锥）/ G5（A7 消失）**，报告会是空壳。30 分钟风险用**小矩阵 + `--fast` + 逐条落盘**缓解。 |
| **3** | **μ ∈ {0.0, 0.3, 0.6}`，3 点** | G6 单调性**最少 3 点**；已含 G2 极限自检所需的 `μ=0.0` 与 QUOTABLE 基准 `μ=0.6`。`0.8` 超出 QUOTABLE 设定、对任何门禁无增益，却把 run 数抬高 40% —— 在时间顶棚下不划算。 |
| **4** | **两者都报，且分码**：`G1 → N/A_minimum_version`；`G7 → RETIRED_by_contract_条6-8`；**一律不省略** | 省掉门禁看起来像藏。G1 是**范围**问题（最小版无泡沫 ⇒ 无 ξ 曲线可复现）；G7 是**政策退役**（契约条 6–8），两者语义不同，**不能共用一个码**，否则丢掉追溯链。 |
| **5** | **跟骨单骨**（沿 `NONVERTICAL_S5_SUPINATION_REPORT.md` 的腓骨单骨先例） | 9 骨扫描会让矩阵爆炸，且距下/胫骨等的跖面 BC 语义需各自定义。本波是**接触能力**验证，不是 9 骨覆盖。 |
| **6** | `.feb`/`.xplt` 落 `temp\opensim_fe\s5_contact\`，**报告与 JSON 落 `results\opensim_fe\`** | `temp/` 已整目录 ignore（可再生），`results/opensim_fe/` 才是入库产物 —— 这正好把"可再生的求解产物"与"要留档的结论"分开。 |
| **7** | ✅ **`pad_scenario="hard_surface"` 确认在范围内** | G2 极限自检（刚性地面）与 G5 的 A7 参照**都必须有硬地面行**，缺了这两门就无从判定。 |
| **8** | 只带 3 个：`σ_contact_avg`（= `F_n/60000`，QUOTABLE 口径）、`max`、`p95`；**其余分位数不进 JSON** | QUOTABLE 只引用全足迹平均值；列越多越容易被下游误引。`max`/`p95` 保留作边缘奇异诊断（G5 需要）。 |
| **9** | **标量时间最大值**（`friction_cone_ratio_max`），不存时程 | G4 判据是"**任一时刻** \|F_t\| ≤ μ·F_n"⇒ 正确诊断就是时域最大值。存时程会让 JSON 膨胀且无人消费。 |
| **10** | ✅ **无热耦合**，与 QUOTABLE §3 #1 的标准条件假设一致 | QUOTABLE 的锚点本就是无热假设下测的；引入热场会让本波数字与锚点不同基。 |

**⇒ 实现者可立即开工。** 依据：本表 + §0 约定 + §1 schema + §2 报告结构 + §3 场景矩阵 + §4 文件名。

---

## Appendix A: cross-references (read first, in this order)

1. `docs\S5_contact_contract.md` §0 (disambiguation), §4 (acceptance G0..G7 + 条 7 scope
   limit), §5 (red lines), §6 (mandatory report format).
2. `docs\S5_QUOTABLE.md` (authoritative numbers; nothing outside this file may be cited).
3. `docs\非垂直落地扩展方案.md` §2 (6-DoF table), §3 L4 (FE boundary layer), §4 S5
   (this wave), §6 risk ① (1D→3D contact rewrite).
4. `docs\S5接触方案.md` §3 (two versions: §3.1 minimum / §3.2 full), §9 (deliverables).
5. `scripts\opensim_fe\nonvertical_s5_supination.py` (731 lines; the boilerplate mirror
   target) and `results\opensim_fe\nonvertical_s5.json` (the JSON shape target).
6. `scripts\opensim_fe\nonvertical_s6_external_rotation.py` and
   `results\opensim_fe\ankle_external_rotation.json` (the basename-split precedent).
7. `scripts\opensim_fe\nonvertical_s4.py` and `results\opensim_fe\nonvertical_s4.json`
   (the FEBio-featuring precedent, for FE shape and A7 reduction calculation).
8. `results\opensim_fe\figures\FIGURE_CAPTIONS.md` (the `measured` / `modeled` /
   `assumed` label taxonomy, §"标签约定").

## Appendix B: explicit NOT-DO list (do not mention in the JSON or report)

The following classes must **not** appear in any cell of the JSON or any line of the report.
They are identified by their `docs\S5_QUOTABLE.md` §4 row, not restated here (the wave-4 doc
must not itself become a carrier of a non-quotable number):

- `docs\S5_QUOTABLE.md` §4 #2: the single-mesh σ/σ_law quoted in the older downstream docs.
- `docs\S5_QUOTABLE.md` §4 #1: any fixed-ξ σ/σ_law pick and any free-ξ drift figure.
- `docs\S5_QUOTABLE.md` §4 #4: the retracted large-strain "flat-punch physics" pair
  (see also contract §G7').
- `docs\S5_QUOTABLE.md` §4 #3: the historical OLD-recipe small-strain anchor values,
  unless cited explicitly as cross-recipe-incompatible (the wave-4 report must not do this).
- Any σ/σ_law drift percentage from the G7 verification docs; only the **fact** of FAIL may
  be cited, never the drift number itself, in the wave-4 report.
- The "OLD recipe dies on the NEW stepper" comparison as an argument; cite
  `S5_QUOTABLE.md` §4 #3 instead.

The classes above are forbidden because (a) they are NOT-QUOTABLE per `S5_QUOTABLE.md` §0
/ §4, or (b) they are cross-recipe comparisons forbidden by `S5_QUOTABLE.md` §4 #3. The
linter must reject a JSON or report that quotes any of these strings.

## Appendix C: change log (this doc)

- **2026-10-06**: initial proposal. Status: `PROPOSAL (awaiting user approval)`. The
  implementer must not write `scripts\opensim_fe\nonvertical_s5_contact.py` until the
  user has approved §1 (schema), §3 (scope), and §4 (filename), and has answered items
  1-4 in §5.