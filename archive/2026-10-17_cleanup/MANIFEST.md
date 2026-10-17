# 归档清单 2026-10-17_cleanup

> 仅 **move**（可恢复）。事由：清理 2026-10-16~17 返工事故后的过期产物与草稿；
> 约定见 `CANONICAL.md` §3。

## 已移动

- `docs/s5-contact/S5_contact_hybrid.md` → `archive/2026-10-17_cleanup/docs_old/S5_contact_hybrid.md`
- `temp/_fix_ignore.py` → `archive/2026-10-17_cleanup/temp_drafts/_fix_ignore.py`
- `temp/_push_learnings.py` → `archive/2026-10-17_cleanup/temp_drafts/_push_learnings.py`

## 拒移（仍被引用，需人工处置）

- `results/UNLOAD_PROCESS_REPORT.md` ← docs/methods/P3_卸力过程最小物理化_预注册.md
- `results/second_contact_probe.json` ← docs/methods/SECOND_CONTACT_预注册.md, results/SECOND_CONTACT_REPORT.md, scripts/second_contact_probe.py
- `results/calibration_ff_scan.csv` ← scripts/calibrate_rope.py
- `results/calibration_rope.csv` ← scripts/calibrate_rope.py
- `docs/notes/ANKLE_600N_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/reports/踝子线_scope决策_2026-10-07.md
- `docs/notes/ANKLE_FE_TODO.md` ← docs/methods/OpenSim_FE复现方案.md, scripts/ankle_fe/geom.py
- `docs/notes/ANKLE_LIG_TODO.md` ← docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/CONV_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/PYFEBIO_TODO.md` ← docs/methods/OpenSim_FE复现方案.md
- `docs/methods/论文对齐基线与场景探索计划.md` ← docs/paper/汇总报告_骨与韧带损伤建模.md, docs/reports/OpenSim_FE交接.md, docs/reports/对齐_Y25_假设与结论.md, scripts/opensim_fe/pad_stiffness.py
- `docs/methods/踝FE可用模型使用说明.md` ← docs/README.md, docs/notes/ANKLE_FIX_TODO.md, docs/reports/踝子线_scope决策_2026-10-07.md
- `docs/methods/踝关节FE建模交接.md` ← docs/notes/ANKLE_ROADMAP_TODO.md, docs/notes/FE_PIPELINE_TODO.md
- `docs/methods/踝韧带判据方案.md` ← docs/paper/汇总报告_骨与韧带损伤建模.md, docs/reports/OpenSim_FE交接.md, results/opensim_fe/ankle_ligament_criterion.json, results/opensim_fe/ANKLE_LIGAMENT_CRITERION_REPORT.md, scripts/opensim_fe/ankle_ligament_criterion.py
- `docs/methods/非垂直落地扩展方案.md` ← docs/paper/汇总报告_骨与韧带损伤建模.md, docs/reports/OpenSim_FE交接.md, docs/reports/场景探索结论.md, docs/reports/对齐_Y25_假设与结论.md, results/opensim_fe/ankle_external_rotation.json, results/opensim_fe/ANKLE_EXTERNAL_ROTATION_REPORT.md, results/opensim_fe/nonvertical_s1.json, results/opensim_fe/NONVERTICAL_S1_REPORT.md, results/opensim_fe/nonvertical_s2.json, results/opensim_fe/NONVERTICAL_S2_REPORT.md, results/opensim_fe/nonvertical_s3.json, results/opensim_fe/NONVERTICAL_S3_POSTURE_REPORT.md, results/opensim_fe/nonvertical_s4.json, results/opensim_fe/NONVERTICAL_S4_REPORT.md, results/opensim_fe/nonvertical_s5.json, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/NONVERTICAL_S5_SUPINATION_REPORT.md, results/opensim_fe/figures/FIGURE_CAPTIONS.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, scripts/opensim_fe/nonvertical_s1.py, scripts/opensim_fe/nonvertical_s2.py, scripts/opensim_fe/nonvertical_s3_posture.py, scripts/opensim_fe/nonvertical_s4.py, scripts/opensim_fe/nonvertical_s5_contact.py, scripts/opensim_fe/nonvertical_s5_supination.py, scripts/opensim_fe/nonvertical_s6_external_rotation.py
- `docs/notes/ANKLE_600N_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/reports/踝子线_scope决策_2026-10-07.md
- `docs/notes/ANKLE_FE_TODO.md` ← docs/methods/OpenSim_FE复现方案.md, scripts/ankle_fe/geom.py
- `docs/notes/ANKLE_LIG_TODO.md` ← docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/CONV_TODO.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/notes/PYFEBIO_TODO.md` ← docs/methods/OpenSim_FE复现方案.md
- `docs/reports/A1_距下关节_交接报告.md` ← docs/notes/ANKLE_ROADMAP_TODO.md, docs/reports/踝子线_scope决策_2026-10-07.md, scripts/ankle_fe/ankle_lig_feb.py
- `docs/reports/Anderson双论文_关键信息.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/reports/两篇临床论文启示与建议.md` ← docs/paper/汇总报告_骨与韧带损伤建模.md, docs/reports/OpenSim_FE交接.md, results/opensim_fe/ankle_external_rotation.json, results/opensim_fe/ANKLE_EXTERNAL_ROTATION_REPORT.md, results/opensim_fe/ankle_ligament_criterion.json, results/opensim_fe/ANKLE_LIGAMENT_CRITERION_REPORT.md, results/opensim_fe/nonvertical_s5.json, results/opensim_fe/NONVERTICAL_S5_SUPINATION_REPORT.md, results/opensim_fe/pad_supination_tradeoff.json, scripts/opensim_fe/ankle_ligament_criterion.py, scripts/opensim_fe/nonvertical_s5_supination.py, scripts/opensim_fe/nonvertical_s6_external_rotation.py, scripts/opensim_fe/pad_supination_tradeoff.py
- `docs/reports/吴恺论文_重解析关键信息.md` ← docs/notes/ANKLE_FIX_TODO.md, docs/notes/FE_PIPELINE_TODO.md, docs/reports/总体进度.md
- `docs/s5-contact/S5_contact_contract.md` ← docs/README.md, docs/reports/OpenSim_FE交接.md, docs/reports/S5_QUOTABLE.md, docs/reports/总体进度.md, docs/reports/阶段二推断_输出清单.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md, results/opensim_fe/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md, results/opensim_fe/sink130_140_150_closure.json, results/opensim_fe/SINK_130_140_150_CLOSURE_REPORT.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_REPORT.md, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md, scripts/opensim_fe/nonvertical_s5_contact.py, scripts/opensim_fe/nonvertical_s5_contact_wave5.py, scripts/opensim_fe/run_sink_closure_20261010.py, src/climbing/coupling/pad_contact.py, src/climbing/coupling/pad_foam.py, tests/test_pad_foam.py, tests/test_pad_mesh.py
- `docs/s5-contact/S5_contact_convergence.md` ← docs/reports/OpenSim_FE交接.md, docs/reports/S5_QUOTABLE.md, results/opensim_fe/figures/FIGURE_CAPTIONS.md, src/climbing/coupling/pad_contact.py, src/climbing/coupling/pad_mesh.py, tests/test_pad_contact.py, tests/test_pad_mesh.py
- `docs/s5-contact/S5_contact_facet_fix.md` ← docs/notes/FE_BC_GOTCHAS.md, scripts/opensim_fe/proof_s5_facet_fix.py, scripts/opensim_fe/thums_feb.py, src/climbing/coupling/joint_contact.py, src/climbing/coupling/plantar_bc.py, tests/test_s5_contact_facet_fix.py
- `docs/s5-contact/S5_contact_notes.md` ← docs/reports/OpenSim_FE交接.md, src/climbing/coupling/pad_contact.py
- `docs/s5-contact/S5_contact_udg.md` ← docs/notes/PENALTY_FIX_RESEARCH.md, docs/reports/OpenSim_FE交接.md, docs/reports/S5_QUOTABLE.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_REPORT.md, scripts/opensim_fe/g4_fe_plate.py, scripts/opensim_fe/nonvertical_s5_contact.py, scripts/opensim_fe/thums_feb.py, src/climbing/coupling/plantar_bc.py, tests/test_s5_contact_facet_fix.py
- `docs/s5-contact/S5_g7_verify_anchors.md` ← docs/reports/OpenSim_FE交接.md, docs/reports/S5_QUOTABLE.md, results/opensim_fe/figures/FIGURE_CAPTIONS.md, src/climbing/coupling/pad_contact.py
- `docs/s5-contact/S5_g7_verify_anchors_mesh.md` ← docs/reports/S5_QUOTABLE.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, scripts/opensim_fe/nonvertical_s5_contact.py, src/climbing/coupling/pad_contact.py
- `docs/s5-contact/S5_g7_verify_element.md` ← docs/reports/S5_QUOTABLE.md
- `docs/s5-contact/S5_g7_verify_fixedxi.md` ← docs/reports/S5_QUOTABLE.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, scripts/opensim_fe/nonvertical_s5_contact.py, src/climbing/coupling/pad_contact.py
- `docs/s5-contact/S5_g7_verify_mesh.md` ← docs/reports/OpenSim_FE交接.md, docs/reports/S5_QUOTABLE.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/figures/FIGURE_CAPTIONS.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, scripts/opensim_fe/nonvertical_s5_contact.py, src/climbing/coupling/pad_contact.py
- `docs/s5-contact/S5_g7_verify_metric.md` ← docs/reports/S5_QUOTABLE.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, scripts/opensim_fe/nonvertical_s5_contact.py
- `docs/s5-contact/S5_wave3_contact.md` ← docs/reports/OpenSim_FE交接.md, docs/reports/总体进度.md, scripts/opensim_fe/nonvertical_s5_contact.py
- `docs/s5-contact/S5_wave4_design.md` ← docs/reports/总体进度.md, results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_REPORT.md, scripts/opensim_fe/nonvertical_s5_contact.py
- `docs/s5-contact/S5_wave4_impl.md` ← docs/reports/OpenSim_FE交接.md, docs/reports/总体进度.md
- `docs/s5-contact/S5_wave5_matrix.md` ← docs/reports/S5_QUOTABLE.md, docs/reports/总体进度.md, results/opensim_fe/nonvertical_s5_contact_wave5.json, results/opensim_fe/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact_wave5.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact_wave5.partial6.json, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_WAVE5_REPORT.md, scripts/opensim_fe/nonvertical_s5_contact.py, scripts/opensim_fe/nonvertical_s5_contact_wave5.py
- `docs/s5-contact/S5接触方案.md` ← docs/paper/汇总报告_骨与韧带损伤建模.md, docs/reports/OpenSim_FE交接.md, docs/reports/总体进度.md, results/opensim_fe/nonvertical_s5_contact.json, results/opensim_fe/NONVERTICAL_S5_CONTACT_REPORT.md, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.json, results/opensim_fe/rerun_20261008/nonvertical_s5_contact.v1.json, results/opensim_fe/rerun_20261008/NONVERTICAL_S5_CONTACT_REPORT.md, scripts/opensim_fe/nonvertical_s5_contact.py, src/climbing/coupling/plantar_bc.py
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
