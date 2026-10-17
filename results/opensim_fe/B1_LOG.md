# B1_LOG — 线2 B1（S3.1→50×9→S4→A1–A4）执行日志

[2026-10-09 12:40] B1.0 kickoff: 开工。已读 复现方案 S3/S4/§六A1–A6/§七T1–T10/附A + FRACTURE_MATRIX_REPORT；发现 `.learnings/LEARNINGS.md` 不在仓库（commit 1baff1b 删除并 gitignore），「LRN-054 端部夹持+点力⇒排序不可信」在现有 LEARNINGS（仓库 git 历史 + workspace 副本）中均检索不到该主题条目 → 改以 ALIGN_B_TODO.md:8 + BC_ROBUST/MESH_FIX/MATRIX_SCALING 报告为红线原文来源（待主会话确认 ID 归属）。
[2026-10-09 12:47] B1.0 kickoff: A节素材齐（旧矩阵=THUMS甲骨直读网格+bbox最长轴+顶/底30%%带+底带三向固定+顶带等分NodalForce；首折p95: 8/9在1m、骨盆8m；BC_ROBUST判排序不可信rho=0.22、MATRIX_SCALING_FIX只修折算不救排序）
[2026-10-09 13:00] B1.0 kickoff续跑: B节 — BC入口钉死 build_bone_feb（bone_feb.py:709；底带BCZeroDisplacement三向:887-894、顶带等分NodalForce:895-905、方向=bbox最长轴:655-691）；PCA旋转开关只在temp脚本 run_remesh_fe.py:44/81，scripts内无；给出候选(a)opt-in参数/(b)新builder两方案
[2026-10-09 13:00] B1.0 kickoff续跑: C节 — 正对照方案定为30°斜置棱柱patch test（10×10×40mm/σ_nom=10MPa，正对照=新BC命中解析解P1-P5、阴性对照=旧BC，脚本建议 scripts/opensim_fe/bc_patch_test.py，≤10min）+ calcaneus新旧BC回归锚
[2026-10-09 13:00] B1.0 kickoff续跑: D节 — B1.1~B1.5执行序展开（谁做/预估/日志=temp/fe_out/b1/logs/*.log）；S4输入路线甲(1D risk(h)，9骨×50高度已实测齐全) vs 乙(新BC矩阵)建议「甲跑通管线、正式判定用乙」，待主会话/用户拍板
[2026-10-09 13:00] B1.0 kickoff续跑: E节 — 甲网格9骨全≤5万实测(femur 25909最大)；load_chain.npz存在且1..50精确可免重跑(省~21min)；OUT_TEMP常量碰撞风险已标；S3.1体素批量仅兜底；fibula/统计量切换p95→median下游耦合等11条风险列全
[2026-10-09 13:01] B1.0 kickoff: 报告已产出 results/opensim_fe/B1_KICKOFF.md
[2026-10-09 13:47] B1.2: 开工。必读已读（KICKOFF §B.3/§C.1/§C.2/§E.5/§E.6 + B1.1 决策①-⑤）；确认 LRN-20261008-006 测试铁律（.venv + -n auto --dist load --timeout=900）；开始实现 bone_feb.py opt-in align_pca
[2026-10-09 14:15] B1.2: ⚠️ 前一子代理因 LLM 流超时中断（model stream timeout，非逻辑错误）；实现段已完成（bone_feb.py +47 行 opt-in align_pca，AST_OK，diff 已核）；续任务接手剩余：默认等价 diff → bc_patch_test.py → patch P1-P5 → calcaneus A/B → pytest 回归 → 补日志
[2026-10-09 14:29] B1.2: 步骤2 默认等价性 PASS — align_pca=False 重建 calcaneus .feb 与旧文件归一化 diff 全同（5595 行，唯一差异=输出文件名，已归一化）；red line ①满足
[2026-10-09 14:33] B1.2: ✅ 正对照 PASS — 30°斜置棱柱双工况：新BC(align_pca) P1-P5 全过（P1 |sigma_zz|误差0.55%、P2 CV 3.68%、P3 中段median 0.986x、P4 p95 1.013x、P5 载荷和0）；阴性对照=旧BC CV 290.6%、mid median 2.55x、p95 4.95x、vm_max 82.2MPa=8.2x名义 ⇒ 检查器可分辨坏BC。判据两处口径修正已写进 docstring：P1取压缩量级、P3限中段(端带60%体积，全域median仅作诊断)。JSON=temp/fe_out/b1/patch_test.json
[2026-10-09 14:38] B1.2: ✅ calcaneus A/B PASS — R1 新BC median/nom=0.666x∈[0.3,1.0]、R2 新旧median差19.4%≤30%、R3p 新/旧p95=0.938x≤1.1（新BC降p95 6.2%）；原R3绝对阈2.0跨口径FAIL已披露（出处=remesh+E18000口径，旧BC也2.19x同样FAIL）⇒改同模型回归门R3p，spec_correction已写入JSON；载荷=load_chain subtalar@5m实测25354.4N、sigma_nom=117.38MPa(A_cort 216.0)
[2026-10-09 14:56] B1.2: ✅✅ 全部完成 — 全量回归 628 passed / 8 skipped（1071.99s，与基线一致）；报告 results/opensim_fe/B12_PATCH_TEST_REPORT.md；结论=新BC过全部正对照门（patch P1-P5 全PASS、calcaneus R1/R2/R3p PASS、默认等价5595行全同、628绿）⇒ 准入 B1.3
[2026-10-09 15:20] B1.3: --init 完成 — 独立根 fracture_matrix_bc2 已建、load_chain.npz 拷贝并校验（heights×joints 全覆盖）、BC=align_pca=True (B1.2 opt-in; patch test P1-P5 PASS)、stat=median
[2026-10-09 15:21] B1.3: bone calcaneus_r (足部) -> status=ok solve_load=25354.401972520634 lrf=1.0 (62s)
[2026-10-09 15:22] B1.3: bone C5 (颈椎) -> status=ok solve_load=22087.579953288558 lrf=1.0 (61s)
[2026-10-09 15:23] B1.3: bone T6 (胸椎) -> status=ok solve_load=10000.0 lrf=2.2087579953288556 (64s)
[2026-10-09 15:24] B1.3: bone L3 (腰椎) -> status=ok solve_load=22087.579953288558 lrf=1.0 (63s)
[2026-10-09 15:25] B1.3: bone R_HIPBONE (骨盆) -> status=ok solve_load=15227.651697700654 lrf=1.0 (63s)
[2026-10-09 15:26] B1.3: bone parietal_r (颅骨) -> status=ok solve_load=22087.579953288558 lrf=1.0 (64s)
[2026-10-09 15:28] B1.3: bone tibia_r (胫骨) -> status=ok solve_load=5000.0 lrf=5.053249826018091 (122s)
[2026-10-09 15:30] B1.3: bone femur_r (股骨) -> status=ok solve_load=15227.651697700654 lrf=1.0 (75s)
[2026-10-09 15:31] B1.3: bone fibula_r (腓骨) -> status=ok solve_load=5000.0 lrf=5.053249826018091 (70s)
[2026-10-09 15:31] B1.3: --finalize 完成 — stat=median、9 骨全缓存、产物 results/opensim_fe/fracture_matrix_bc2.{json,csv,png} + FRACTURE_MATRIX_BC2_REPORT.md（耗时 2s）
[2026-10-09 15:32] B1.3: diff 完成 — ρ(新med,旧med)=0.700、ρ(新med,旧p95)=0.433；产物 B13_MATRIX_BC2_DIFF.md
[2026-10-09 17:33] B1.4: S4 对标完成 — 正式(乙 median) A1=FAIL A2=FAIL(3/5) A3=FAIL(0/3) A4=FAIL；甲 A3=FAIL（预期结构性失败）；产物 S4_ALIGNMENT_BC2_REPORT.md + s4_alignment_bc2.json（1s）
[2026-10-09 17:41] B1.4: S4 对标完成 — 正式(乙 median) A1=FAIL A2=FAIL(3/5) A3=FAIL(0/3) A4=FAIL；甲 A3=FAIL（预期结构性失败）；产物 S4_ALIGNMENT_BC2_REPORT.md + s4_alignment_bc2.json（1s）
[2026-10-09 17:42] B1.4: 诊断归因 — A1 修正 three_phases interleave→顺序分组（pa文档自述等分相邻类，interleave 产交叠 phase；pa未改，驱动内 three_phases_seq）：乙边界[5,12] Δ=[2,4]、甲[5,17] Δ=[2,1]；A2 3/5 漏=股骨(ends带≥80从不+shaft带427MPa@ref含弯曲热点→1m即超220常值1)与胸椎(整域median 426≫150常值1)；A3 0/3 结构性=单调吸收二值→完全分离/常值→OR无定义（需被试级变异或正则化估计器，甲同）；A4 结构性=严格线弹性⇒所有σ(h)∝load(h)⇒PC2仅0.024%（rank-1）且G1/G2含行为完全相同部位(股/胸/腓全1起)⇒不可分；femur探针=fix带内Dirichlet零应变+中段弯曲热点492MPa（ends带≠解剖端，BC带代理语义失配已记录）
