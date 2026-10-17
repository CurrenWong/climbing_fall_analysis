==============================================================================================
① 目录结构（文件数 / 功能归属）
==============================================================================================

目录                                    文件数   功能（按内容判读）
src/climbing                           40   **核心库**（1D 两大分支 + 耦合层）
src/climbing/coupling                  22   **FE 耦合层**（OpenSim→FEBio 链、接触、材料）
scripts                               196   ⚠️ **根目录是"平的"**：1D 落地线探针 + 杂项工具混放
scripts/opensim_fe                    133   **FEBio/OpenSim 探针与驱动**（62）
scripts/ankle_fe                       17   **踝子线专用**（建模/网格/验收）
scripts/p5_cond                        11   **论文数字校验**（5 个，成套）
tests                                  34   **回归测试**（按模块一一对应，命名规约清晰）
docs/methods                           32   预注册/设计（按 G*/P* 编号，清晰）
docs/notes                             17   踩坑手册 + 旧 TODO（手册 3 本，活跃）
docs/reports                           20   结论/交接/进度（含枢纽 `总体进度.md`）
docs/paper                              3   论文正文（v2）
docs/s5-contact                        16   ⚠️ **S5 阶段遗留**（16 个，10-07）
docs/archive                            6   已归档文档
results                               443   产物（报告 + json）
results/opensim_fe                    391   ⚠️ **最大杂堆**（FE 产物 300+）
temp                                14492   ⚠️ **草稿区**，但 `temp/opensim_fe` 有 **12,593** 个 FE 产物
temp/opensim_fe                     12583   
model                                 100   FE 模型输入
paper                                  29   论文源（LaTeX/图表）
video                                  33   短视频线（独立）
legacy                                 11   历史遗留
archive                              2065   本次整理归档

==============================================================================================
② 跨文件重复函数名（同一天里「复制公式」的温床）
==============================================================================================
  main                             150 文件  例: scripts/ankle_fe/force_root.py, scripts/ankle_fe/lig_strain.py, scripts/ankle_fe/penalty_auto0_probe.py
  build                             12 文件  例: scripts/ankle_fe/check_c2.py, scripts/opensim_fe/axial_subregion.py, scripts/opensim_fe/fibula_identifiability_analysis.py
  rel                                8 文件  例: scripts/archive_cleanup.py, scripts/archive_doc_clusters.py, scripts/audit_2026_10_17.py
  run                                8 文件  例: scripts/chain_soft_probe.py, scripts/ecc_muscle_probe.py, scripts/opensim_fe/g9_fibula_fe_support.py
  write_report                       8 文件  例: scripts/opensim_fe/fibula_fe_audit.py, scripts/opensim_fe/fracture_matrix.py, scripts/opensim_fe/g3_subject_variability.py
  force_n                            5 文件  例: src/climbing/bodies.py, src/climbing/chain2d.py, src/climbing/chain2d_crush.py
  settled                            5 文件  例: src/climbing/bodies.py, src/climbing/chain2d.py, src/climbing/chain2d_crush.py
  statistic                          5 文件  例: scripts/opensim_fe/fe_bc_linearity.py, scripts/opensim_fe/fe_bc_perheight.py, scripts/opensim_fe/fe_bc_probe.py
  peak_force_kn                      4 文件  例: src/climbing/bodies.py, src/climbing/bone.py, src/climbing/chain2d.py
  rhs                                4 文件  例: src/climbing/bodies.py, src/climbing/head_branch.py, src/climbing/pad.py
  positive_control                   4 文件  例: scripts/fp_anchor_probe.py, scripts/fp_contact_probe.py, scripts/fp_stroke_probe.py
  scaled_pad                         4 文件  例: scripts/mc_stat_robustness.py, scripts/p5_cond/pad_stiffness_h20.py, scripts/p5_cond/pad_stiffness_sweep.py
  parse_feb                          4 文件  例: scripts/ankle_fe/acceptance.py, scripts/ankle_fe/check_c2.py, scripts/ankle_fe/lig_strain.py
  load                               4 文件  例: scripts/opensim_fe/fe_bc_agg_sens.py, scripts/opensim_fe/make_model_figures.py, scripts/opensim_fe/v1_max_caliber_check.py

==============================================================================================
③ ★ 全局可变「门」（gate dict）—— 隐式状态，脚本间会互相影响
==============================================================================================
  RIGID_STANDING_CHAIN       定义于 src/climbing/pad.py
  RIGID_CONTACT              定义于 src/climbing/pad.py
  STROKE_ANCHOR              定义于 src/climbing/pad.py
  TWO_SEG_LEG                定义于 src/climbing/pad.py

  谁在改它们（脚本侧 write 次数）：
    scripts/ecc_muscle_probe.py                          RIGID_CONTACT, TWO_SEG_LEG
    scripts/fp_contact_probe.py                          RIGID_CONTACT
    scripts/fp_stroke_probe.py                           RIGID_CONTACT, STROKE_ANCHOR
    scripts/leg2seg_probe.py                             RIGID_CONTACT, TWO_SEG_LEG
    scripts/passive_path_probe.py                        RIGID_CONTACT, TWO_SEG_LEG
    scripts/recompute_hstar.py                           RIGID_CONTACT, TWO_SEG_LEG
    scripts/second_contact_probe.py                      RIGID_CONTACT, TWO_SEG_LEG
    scripts/unload_process_probe.py                      RIGID_CONTACT, TWO_SEG_LEG

==============================================================================================
④ 已知坑类的复发检测（模式扫描）
==============================================================================================

  ▸ 对采样轨迹做数值微分（采样依赖）：**25 处**
      scripts/ecc_muscle_probe.py:66
      scripts/pad_energy_check.py:102
      scripts/unload_process_probe.py:77
      src/climbing/bodies.py:337
      src/climbing/bodies.py:341
      src/climbing/bodies.py:534
      …另 19 处

  ▸ 三层目录脚本（与 tests 的 parents[1] 不同，易混）：**74 处**
      scripts/ankle_fe/geom.py:42
      scripts/opensim_fe/a1_anchor_diag.py:16
      scripts/opensim_fe/a1_boundary_convention.py:13
      scripts/opensim_fe/a1_boundary_window.py:45
      scripts/opensim_fe/a1_calcaneus_stat.py:14
      scripts/opensim_fe/a1_gauge_check.py:15
      …另 68 处

  ▸ 脚本自注入 sys.path（须与本文件层级匹配）：**197 处**
      scripts/ankle_fe/geom.py:40
      scripts/calibrate_rope.py:28
      scripts/chain_crush_probe.py:35
      scripts/chain_soft_probe.py:14
      scripts/chain_tissue_probe.py:16
      scripts/chain_vtarget_probe.py:31
      …另 191 处

  ▸ M(θ) 表接口（契约见 CANONICAL.md §0）：**20 处**
      scripts/recompute_hstar.py:3
      scripts/recompute_hstar.py:4
      scripts/recompute_hstar.py:31
      scripts/recompute_hstar.py:49
      scripts/second_contact_probe.py:30
      scripts/second_contact_probe.py:46
      …另 14 处

  ▸ 腿力观测量（须与 rhs 同源）：**32 处**
      scripts/ecc_muscle_probe.py:70
      scripts/fp_anchor_probe.py:68
      scripts/fp_anchor_probe.py:70
      scripts/fp_anchor_probe.py:82
      scripts/fp_anchor_probe.py:90
      scripts/fp_contact_probe.py:66
      …另 26 处

  ▸ 按「键」判定门的姿势（传对象会静默绕门）：**52 处**
      scripts/fp_anchor_probe.py:65
      scripts/fp_anchor_probe.py:88
      scripts/fp_anchor_probe.py:197
      scripts/phase_p1_montecarlo.py:779
      scripts/phase_p1_montecarlo.py:874
      scripts/second_contact_probe.py:119
      …另 46 处

  ▸ pyfebio 接触 .xplt 无名域缺陷：**41 处**
      scripts/ankle_fe/build_feb.py:705
      scripts/opensim_fe/bone_feb.py:13
      scripts/opensim_fe/bone_feb.py:14
      scripts/opensim_fe/bone_feb.py:1207
      scripts/opensim_fe/bone_feb.py:1208
      scripts/opensim_fe/bone_feb.py:1224
      …另 35 处
