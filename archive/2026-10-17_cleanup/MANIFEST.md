# 归档清单 2026-10-17_cleanup

> 仅 **move**（可恢复）。事由：清理 2026-10-16~17 返工事故后的过期产物与草稿；
> 约定见 `CANONICAL.md` §3。

## 已移动

- `results/b2_report.txt` → `archive/2026-10-17_cleanup/results_invalid/b2_report.txt`
- `results/smoke.txt` → `archive/2026-10-17_cleanup/results_invalid/smoke.txt`
- `results/calibration_report.txt` → `archive/2026-10-17_cleanup/results_invalid/calibration_report.txt`
- `results/rope_report.txt` → `archive/2026-10-17_cleanup/results_invalid/rope_report.txt`
- `scripts/verify_flex_observable.py` → `archive/2026-10-17_cleanup/scripts_invalid/verify_flex_observable.py`
- `SECOND_CONTACT_TODO.md` → `archive/2026-10-17_cleanup/docs_stale/SECOND_CONTACT_TODO.md`

## 拒移（仍被引用，需人工处置）

- `results/UNLOAD_PROCESS_REPORT.md` ← docs/methods/P3_卸力过程最小物理化_预注册.md
- `results/second_contact_probe.json` ← docs/methods/SECOND_CONTACT_预注册.md, results/SECOND_CONTACT_REPORT.md, scripts/second_contact_probe.py
- `results/calibration_ff_scan.csv` ← scripts/calibrate_rope.py
- `results/calibration_rope.csv` ← scripts/calibrate_rope.py
- `docs/notes/ANKLE_600N_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/reports/踝子线_scope决策_2026-10-07.md
- `docs/notes/ANKLE_FE_TODO.md` ← docs/methods/OpenSim_FE复现方案.md, docs/methods/踝关节FE建模交接.md, scripts/ankle_fe/geom.py
- `docs/notes/ANKLE_LIG_TODO.md` ← docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/CONV_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/PYFEBIO_TODO.md` ← docs/methods/OpenSim_FE复现方案.md
- `temp/_inventory.py` ← results/opensim_fe/P2_BC_SCOPE_INVENTORY.md
- `temp/check_nplantar.py` ← results/opensim_fe/V1_CORRECTION_PER_NODE_K.md
- `temp/g5_visco.py` ← docs/methods/G5_QLV粘弹性_预注册.md, results/G5_TISSUE_REPORT.md
- `temp/loadchain_impact.py` ← results/opensim_fe/LOADCHAIN_RESIDUAL_CLOSURE.md
- `temp/meshseq2.py` ← results/opensim_fe/S10_MESH_CONVERGENCE.md
- `temp/parts_sum.py` ← results/opensim_fe/FOOT_MODELS_INVENTORY.md
- `temp/s10_graded.py` ← results/opensim_fe/V1_ABSOLUTE_VALUE_UNDERDETERMINED.md
- `temp/s10_graded_fold.py` ← results/opensim_fe/V1_ABSOLUTE_VALUE_UNDERDETERMINED.md
- `temp/s10_graded_low.py` ← results/opensim_fe/V1_ABSOLUTE_VALUE_UNDERDETERMINED.md
- `temp/s10_seq2.py` ← results/opensim_fe/S10_MESH_CONVERGENCE.md
- `temp/s10_seq_run.py` ← results/opensim_fe/S10_MESH_CONVERGENCE.md
- `temp/s10_tet2domain.py` ← results/opensim_fe/S10_MESH_CONVERGENCE.md
- `temp/s10_two_knob.py` ← results/opensim_fe/V1_ABSOLUTE_VALUE_UNDERDETERMINED.md
- `temp/v1_1d_force2.py` ← results/opensim_fe/V1_1D_APPROXIMATION.md, results/opensim_fe/V1_CALCANEUS_THRESHOLD_SOURCING.md
- `temp/v1_1d_low.py` ← results/opensim_fe/V1_1D_APPROXIMATION.md, results/opensim_fe/V1_CALCANEUS_THRESHOLD_SOURCING.md
- `temp/v1_arrdump.py` ← scripts/pad_energy_check.py
- `temp/v1_b2a_smoke.py` ← results/LEG2SEG_REPORT.md
- `temp/v1_b2a_timing.py` ← results/LEG2SEG_REPORT.md
- `temp/v1_c1_wrapup.py` ← results/ECC_MUSCLE_REPORT.md
- `temp/v1_cap_compliance.py` ← results/opensim_fe/V1_LOADSHARE_VERIFICATION.md
- `temp/v1_cmp.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_IS_INERT.md
- `temp/v1_contactA_probe.py` ← docs/methods/P3_足部接触柔度_预注册.md
- `temp/v1_contactA_scan.py` ← docs/methods/P3_足部接触柔度_预注册.md
- `temp/v1_ctrl.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_IS_INERT.md
- `temp/v1_dyn2.py` ← results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md
- `temp/v1_dyn_scan.py` ← results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md
- `temp/v1_family_full.py` ← results/opensim_fe/V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md
- `temp/v1_fixed_family.py` ← V1_OPT_A_TODO.md, results/opensim_fe/V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md
- `temp/v1_flip2.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_WINDING_FIX.md
- `temp/v1_footpart.py` ← results/opensim_fe/FOOT_MODELS_INVENTORY.md
- `temp/v1_fpanchor_predry.py` ← LEG_K_ANCHOR_TODO.md, docs/methods/P3_腿刚度力台锚定_预注册.md
- `temp/v1_freeflight.py` ← results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md
- `temp/v1_gap_scan.py` ← results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md
- `temp/v1_holdk_big.py` ← results/opensim_fe/V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md
- `temp/v1_hotspot_gauge.py` ← results/opensim_fe/V1_FAILURE_LOCUS.md
- `temp/v1_hotspot_move.py` ← results/opensim_fe/V1_FAILURE_LOCUS.md
- `temp/v1_inspect.py` ← results/opensim_fe/V1_LOADSHARE_VERIFICATION.md
- `temp/v1_iso.py` ← results/opensim_fe/V1_OPTIMIZATION_OPTIONS.md, results/opensim_fe/V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md
- `temp/v1_jam.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_IS_INERT.md
- `temp/v1_kflex.py` ← results/opensim_fe/V1_WHY_20KN_AND_GAP_LOCUS.md
- `temp/v1_legtravel.py` ← results/opensim_fe/V1_WHY_20KN_AND_GAP_LOCUS.md
- `temp/v1_loadshare_p1.py` ← results/opensim_fe/V1_LOADSHARE_VERIFICATION.md
- `temp/v1_matched.py` ← results/opensim_fe/V1_CORRECTION_PER_NODE_K.md
- `temp/v1_mre.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_WINDING_FIX.md
- `temp/v1_mre2.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_WINDING_FIX.md
- `temp/v1_mre_read.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_WINDING_FIX.md
- `temp/v1_newcond.py` ← results/opensim_fe/V1_OPTIMIZATION_OPTIONS.md
- `temp/v1_normal.py` ← results/opensim_fe/V1_CONTACT_UNCONVERGED_ADDENDUM.md, results/opensim_fe/V1_PLANTAR_CONTACT_IS_INERT.md
- `temp/v1_pad_k.py` ← results/opensim_fe/V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md
- `temp/v1_padjac.py` ← results/opensim_fe/V1_CONTACT_UNCONVERGED_ADDENDUM.md, results/opensim_fe/V1_PLANTAR_CONTACT_IS_INERT.md
- `temp/v1_padk_run.py` ← results/opensim_fe/V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md
- `temp/v1_path_geom.py` ← results/opensim_fe/V1_LOAD_PATH_FACTS.md
- `temp/v1_payoff.py` ← results/opensim_fe/V1_PLANTAR_CONTACT_WINDING_FIX.md
- `temp/v1_phys_corner.py` ← results/opensim_fe/V1_CORRECTION_PER_NODE_K.md, results/opensim_fe/V1_PHYSICAL_ANCHOR_AND_FINAL_VERDICT.md
- `temp/v1_phys_spring.py` ← results/opensim_fe/V1_PHYSICAL_ANCHOR_AND_FINAL_VERDICT.md
- `temp/v1_plantar_split.py` ← results/opensim_fe/V1_LOAD_PATH_FACTS.md
- `temp/v1_sym_scan.py` ← results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md
- `temp/v1_ucurve.py` ← results/opensim_fe/V1_PLANTAR_BC_FAMILY_ASYMPTOTE.md
- `temp/v1_where.py` ← results/opensim_fe/V1_FAILURE_LOCUS.md
- `temp/v1_wholefoot.py` ← results/opensim_fe/FOOT_MODELS_INVENTORY.md

## 未列入本批（仍为现行，勿动）

- `results/fp_anchor_probe.json`（7 个脚本的输入锚）
- `results/*_REPORT.md` 中被 `pad.py` / `tests` 引用的 3 个
- `docs/reports/最终结论_2026-10-07.md`（项目规则：旧版结论件保留不改）
