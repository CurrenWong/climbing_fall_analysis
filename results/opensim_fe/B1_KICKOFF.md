# B1_KICKOFF — B1 线2 kickoff（续跑）：现状判定 · BC 代码入口 · 正对照 · 执行序

> 产出：B1.0 子代理续跑（A 节素材采集已完成于 `results/opensim_fe/B1_LOG.md`，本文直接引用其结论不重查）。
> 范围：只读分析 + 本报告（未改任何 src/scripts/测试/已有文档）。
> **红线出处**（LRN-054 经 git 历史 + workspace 双查证伪、已作废）：
> `results/opensim_fe/BC_ROBUST_METRIC_REPORT.md`（§5/§6/§7）+ `MESH_FIX_REPORT.md`（§5/§6/§8）+
> `MATRIX_SCALING_FIX_REPORT.md` + `ALIGN_B_TODO.md`（BC 约定节）。
> 事实类论断带 file:line/节号；推断标「未验证」。

---

## 0. 结论先行（TL;DR）

1. **旧矩阵为什么不可信（已钉死，见 A 节 + 红线报告）**：BC = 斜置骨沿 bbox 最长轴加载 + 顶/底 30% 坐标带 + 底带三向固定 + 顶带等分 NodalForce + 主统计量 p95 ⇒ 排序不可信（ρ=0.22，`BC_ROBUST_METRIC_REPORT.md` §6.1）、二值矩阵 8/9@1m；端带剔除 p95 反升 ⇒ **中段弯曲主导、端夹持假设证伪**（§6.2）。折算 bug 已修（现代码 `fracture_matrix.py:575-576` 已乘 `load_reduction_factor`，修复说明见 `regenerate_fracture_matrix.py:1-15`），**方向与统计量未修**。
2. **BC 代码入口唯一**：`scripts/opensim_fe/bone_feb.py::build_bone_feb`（:709）——固定 = 底带节点集 `bottom_fix` 上 `BCZeroDisplacement(x,y,z dof=1)` 三向（:887-894）；载荷 = 顶带节点集 `top_load` 上**等分** `NodalForce`（:879、:895-905）；方向 = bbox 最长轴或显式 `axis∈{0,1,2}`（`_select_top_bottom` :655-691）。**无面压力、无反力输出**（plotfile 仅 displacement+stress，:912-918）。
3. **「旋转网格 → 骨轴」已有可运行实现，但只在 temp 脚本**：`temp/opensim_fe/fracture_matrix_remesh/run_remesh_fe.py:44`（`align_pca`）+ `:81`（`--align-pca`）；`scripts/` 内 grep 无任何 axis/rotate/PCA 开关，`fracture_matrix.py:360-361` 仍是 `axis=None` 裸调。MESH_FIX §6 实测该轴向对齐把胫骨 p95 94.6×→26.2×、股骨 59.8×→16.9×、跟骨 1.5×（**但端夹持奇异性仍在**）。
4. **最小侵入两候选（详见 §B.3）**：(a) `build_bone_feb` 加 opt-in 参数（`align_pca=False` + `load_mode="nodal_equal"` 默认不变）——改动仅 `bone_feb.py`，推荐；(b) 同目录新 builder 函数（复用 helpers）——零触碰现函数，但有 writer 代码分叉风险。**决策留 B1.1（待主会话/用户拍板）**。统计量用现成 `--stat median`（`fracture_matrix.py:1021`），不改 `PRIMARY_STAT="p95"` 默认（:93）。
5. **正对照方案已设计（§C）**：30° 斜置方形棱柱（10×10×40 mm，E=1000、ν=0.3、F=1000 N → σ_nom=F/A=10.0 MPa）双工况——新 BC 必须命中解析解（中段 |σ̄_zz−σ_nom|/σ_nom ≤5%），旧 BC 作为阴性对照记录；脚本建议 `scripts/opensim_fe/bc_patch_test.py`，全程 ≤10 min。**PASS 前不准上矩阵**（ALIGN_B 红线）。
6. **S4 输入路线（待拍板，详见 §D.2）**：甲（1D risk(h)）现成——`risk_1d_loadshare.py` 与 `risk_1d_bending.py` 均已覆盖 9 部位 × 50 高度（`risk_by_height` 实测 50 元素），BC_ROBUST §6.3 明确推荐；但 `PAPER_ALIGNMENT_REPORT.md` §4.1/4.2 实测 1D 低高度饱和 ⇒ 完全分离无 OR（A3 结构性不可算）、聚类仅 3 类。**建议：甲先跑通管线、正式 A1–A4 用乙（B1.3 新 BC 矩阵）** —— 待拍板。
7. **E 节关键事实（已验证）**：甲原生网格 9 骨全部 ≤5 万（最大 femur 25,909 solid 单元，见 §E.1 表）；`temp/opensim_fe/fracture_matrix/load_chain.npz` 存在且 heights=1..50 精确 + 4 关节齐全 ⇒ **可免重跑**（省实测 ~21 min：fd 617 s + 关节反力 637 s）；但 `fracture_matrix.OUT_TEMP` 固定（:83），新 BC 重跑**必须换独立 temp 根**，否则 `--force-fe` 会覆盖旧 BC 的 `<b>_fe.json`（:390-394）。

---

## B. BC 代码入口与最小侵入修改方案

### B.1 `build_bone_feb` 语义（逐参数，file:line 全部实读）

调用链：`fracture_matrix._build_run`（fracture_matrix.py:354-362）→ `bone_feb.build_bone_feb`（bone_feb.py:709）→（可选）`run_bone_feb`（:963，FEBio 执行，返回 `{rc, elapsed_s,...}`）→ `post_process`（:1097，xplt→hdf5→逐域 σ_vm 统计，`_stats_from_voigt` :1085-1093 给 max/p99/p95/mean/median）。

| 项 | 语义 | 出处 |
|---|---|---|
| `top_frac` / `bot_frac`（默认 0.30/0.30） | 沿选定轴的 **bbox 坐标带**：顶带 = `coord ≥ a_max − top_frac·span`，底带 = `coord ≤ a_min + bot_frac·span`。作用对象是**全部节点**（含壳节点、内部节点、侧面节点），**不是端面/面片**——顶带里每个节点（不论在不在横截面上）都拿到等分力 | bone_feb.py:682-687；节点集写出 :840-846 |
| `axis`（默认 `None`） | `None` = bbox 三向跨度最大者（`argmax(hi−lo)`）；显式时仅接受全局 0/1/2。**没有任意方向（骨轴）支持** | bone_feb.py:671-676 |
| 固定方式 | 节点集 `bottom_fix` 上 `BCZeroDisplacement(x_dof=1, y_dof=1, z_dof=1)` —— **三向全固定**（非 roller、非仅法向） | bone_feb.py:887-894（BCZeroDisplacement 在 :890），节点集 :843-845 |
| 载荷方式 | 节点集 `top_load` 上**单个** `NodalForce` 卡，向量 = `-load_n · ê_axis / N_top`（等分）；**无面压力/Traction、无按分摊面积加权**。方向恒为全局 −axis（压缩） | bone_feb.py:879（注释 "equal-share top NodalForce"）、:895-905（load_vec 在 :897） |
| 载荷施加节奏 | `time_steps` 个等步长 STATIC；⚠️ `fracture_matrix` 的 help 明写 **`>1` 会让这套 NodalForce 配置误 ramp 并发散** | bone_feb.py:879-885；fracture_matrix.py:1023-1025 |
| 输出 | plotfile 仅 `displacement` + `stress`——**反力未导出**（平衡判据需另行加 plot 变量） | bone_feb.py:912-918 |
| 返回 meta | `axis_idx/bbox/n_top_nodes/n_bot_nodes/load_per_node_n/load_desc...`（供日志与溯源） | bone_feb.py:930-959 |

辅助函数：`_select_top_bottom`（:655-691，带选择逻辑本体）、`_fix_febio413`（:694-707，FEBio 4.13 solver 块兼容）、`extract_mesh`（:420，THUMS 直读，依赖 `kmesh_io`/`shell_io`，:75-76）、`read_materials`（:210）。

**结构性问题清单（对红线的代码级落点）**：
- ① 方向：`axis` 只会选全局轴 ⇒ 斜置骨（胫/腓 vs bbox 轴 42.4°/42.8°，`MESH_FIX_REPORT.md` §5）必然被打成弯曲；
- ② 等分 NodalForce 打在**坐标带全部节点**上（含侧面/内部），带边缘载荷突变 ⇒ 应力奇异源之一；
- ③ 底带三向固定 ⇒ 端部过约束（`PLANTAR_BC_REPORT.md` 已在单骨跟骨上量化过同类效应：改弹性地基降 gauge 16–24%，见 `scripts/opensim_fe/plantar_bc_probe.py:520-523` 的机制结论：**法向刚度主导、切向不主导**——可作为底带改 roller/spring 设计的参考，但它走的是另一条单骨管线的 patch 机制 `plantar_bc_probe.py:104-115`，与通用 builder 无关）；
- ④ 主统计量 p95 是 BC 奇异敏感指标（`MESH_FIX_REPORT.md` §8.2：网格越干净 p95 可能越高）。

### B.2 `fracture_matrix.py` 调用形态 & 旋转开关现状

- 调用（已核实与任务描述一致）：`fracture_matrix.py:360-361`
  `build_bone_feb(part.spec, feb, load_n=load, top_frac=0.30, bot_frac=0.30, axis=None, with_shell=True, time_steps=time_steps, mesh=mesh)` —— **全程 `axis=None`，无任何对齐参数**。
- 周边事实：`OUT_TEMP = temp/opensim_fe/fracture_matrix`（:83，**常量**）；`HEIGHTS=1..50`（:89）；`STATS=("max","p95","mean","median")`、`PRIMARY_STAT="p95"`（:92-93）；载荷链缓存复用判断（:280-284）；逐骨 FE 缓存 `<b>_fe.json`（:390-394，`reuse=not --force-fe`）；CLI（:1017-1025，含 `--stat`、`--force-fe`、`--time-steps` 及其发散警告）。
- **旋转/PCA 开关是否已进代码：没有。** `scripts/opensim_fe/` 下 grep `align-pca|align_pca|run_remesh_fe` 零命中；`bone_feb.py` 内仅 `_select_top_bottom` 的 `axis` 全局轴选择；PCA 相关代码全在**配准/1D 分析**（`registration*.py`、`risk_1d_bending.py`）与**图件旋转**（`make_fracture_figures.py:139`），与 BC 无关。
  - 唯一可运行的「PCA 长轴→+Z 再 build」实现：`temp/opensim_fe/fracture_matrix_remesh/run_remesh_fe.py:44-64`（`align_pca`，绕质心旋转）+ `:81`（`--align-pca` 开关）+ `:92-93`（先转 mesh 再走原 builder）+ `:130-133`（build 仍 `axis=None`——转完后 bbox 最长轴=骨轴）。**它在 temp/（易失、非 scripts/），MESH_FIX 报告 §9 自述为"新驱动/诊断脚本，不改既有模块"。**
  - 结果基线（MESH_FIX §6 表）：轴向后 p95/名义 = 跟骨 1.5×、胫骨 26.2×、股骨 16.9× —— 方向修正**有效但不够**（端夹持奇异性残留），median 始终物理（≈0.3–0.7×）。
- 折算状态：BC_ROBUST §5 指出的「少乘 lrf」**已在现代码修复**（fracture_matrix.py:575-576 与 :581-583 均乘 `lrf`；`regenerate_fracture_matrix.py` 即为该修复后的缓存重放驱动）。**引用红线时注意：§5 描述的是修复前状态。**

### B.3 两个最小侵入候选

共同前提（红线）：不改 `fracture_matrix.py:360-361` 现调用的语义、不改任何默认值、不动现有测试；新 BC 的所有效果必须 opt-in。

#### 候选 (a)：`build_bone_feb` 加 opt-in 参数【推荐，暂定】

```python
def build_bone_feb(spec, out_feb, *, load_n, top_frac=0.30, bot_frac=0.30,
                   axis=None, with_shell=True, time_steps=1, mesh=None,
                   align_pca: bool = False,          # 新增，默认 False=现状
                   load_mode: str = "nodal_equal"):  # 新增，默认=现状
```

- `align_pca=True`：build 前把 `mesh["nodes"]` 副本绕质心旋转使 PCA 长轴→+Z（**移植** `run_remesh_fe.py:44-64` 既有已验证实现），随后 `axis=None` 即取骨轴——方向问题一步解决，带选择/固定/载荷逻辑零改动。⚠️ 必须 **copy 后旋转**，`fracture_matrix` 只提取一次 mesh（fracture_matrix.py:408）并在降载重试循环里反复传同一 dict（:428-433）。
- `load_mode`：
  - `"nodal_equal"`（现状默认）；
  - `"nodal_tributary"`：等分 → 按节点分摊面积/关联体积加权（消除侧面节点吃等分力的非物理分布；分摊面积算法需 patch test 校准，「未验证」）；
  - （可选第三档）`"prescribed"`：顶带 `BCPrescribedDisplacement` 位移加载 + 反力校核——需同时给 plotfile 加 reaction 变量（bone_feb.py:916 改 `all_vars`，仅在 opt-in 分支生效）。
- **改动文件清单**：仅 `scripts/opensim_fe/bone_feb.py`（新 kwarg + 一个私有 `_align_pca_nodes()` helper + 分支）。
- **风险**：① 旋转改变节点坐标 ⇒ 下游基于坐标/ bbox 的约定受影响——`risk_1d_bending.py:94`（`THETA_LIMB_DEG` 取自 MESH_FIX 旧网格）与 `paper_alignment.py:33`（载荷轴=bbox 最长轴）在**新矩阵重跑后**需要重新核对（旋转后 θ→0 是构造使然，不是发现）；② `align_pca` 若做成 builder 内部行为，mesh 元数据（`meta.bbox_*`）变成旋转后坐标系，报告里必须写明（返回 dict 已有 `load_desc` 可带标注）。
- **默认行为**：`align_pca=False, load_mode="nodal_equal"` ⇒ 与现状逐字节等价；现有调用点（fracture_matrix.py:360-361、bone_feb_validate.py、temp 驱动）不改一行也不受影响。现有测试零命中 `build_bone_feb`（grep tests/ 验证），无测试破坏面。

#### 候选 (b)：新 builder 函数【备选】

- 新文件 `scripts/opensim_fe/bone_feb_bc2.py`，函数 `build_bone_feb_bc2(...)`：整段 writer 自包含（或 `from bone_feb import _select_top_bottom, _fix_febio413, run_bone_feb, post_process` 复用私有 helper——同目录脚本可直接 import，`run_remesh_fe.py` 已示范该模式）。
- **改动文件清单**：仅新增文件（+B1.3 的新驱动指向它）；`bone_feb.py` **零触碰**。
- **风险**：writer 主体（模型组装/材料/域/4.13 修复/全精度坐标，bone_feb.py:777-920 约 140 行）一旦复制即**双源分叉**——后续 FEBio 坑（4.13 solver 块 :694-707、nameless shell `:988-1041`、`%e` 坐标精度 :781-791）只修一处会静默漂移；缓解 = 只复制、不再改 `bone_feb.py`，且 patch test 对两条路径都跑。
- **默认行为**：天然不破坏（新函数无人调用即无影响）。

**选择建议（留 B1.1 拍板）**：若新 BC = 「对齐 + 载荷加权 + 统计量」三件套（都不改 writer 骨架）→ 选 (a)；若底带要换 roller/prescribed、顶带要上 PressureLoad（需面片 facet 集，pyfebio 有 `Surface`/`PressureLoad`/`TractionLoad`，见 `.venv/Lib/site-packages/pyfebio/mesh.py:221`、`loads.py:60-75`，**但 builder 现在没建 facet 集，属新代码路径**）→ 分支会把 `build_bone_feb` 撑成多模式大函数，倾向 (b)。**（待主会话/用户拍板）**

---

## C. 正对照方案（新 BC 的已知解验证）

> 约束：新 BC 必须带**已知正对照**才准上矩阵（ALIGN_B 红线）；检查器要有已知正对照（SOUL 红线）。方案按「有单一解释的解析解」设计，不依赖反力正负号判方向。

### C.1 主方案：斜置棱柱单轴压缩 patch test（解析对照，含阴性对照）

- **脚本**：新建 `scripts/opensim_fe/bc_patch_test.py`（新文件；骨架复制 `run_remesh_fe.py:96-133` 的「合成 mesh dict → 直接挂 `bf.Material` → `build_bone_feb` → `run_bone_feb` → `post_process`」路径，**不需要** `read_materials`/THUMS deck）。
- **几何/材料/载荷（已知解）**：方形棱柱 10×10×40 mm（A=100 mm²），结构化 hex8 网格 10×10×40=4000 单元（节点 11×11×41=4961），E=1000 MPa、ν=0.3，**整体绕 Y 轴旋转 30°**（复现斜置骨病理），F=1000 N → 解析 σ_nom = F/A = **10.00 MPa**（σ_vm 解析 = 10.00，中段单轴状态）。选方棱柱而非圆柱：端面分摊面积无离散化误差，解析值精确。
- **双工况（正/负对照配对）**：
  1. **阴性对照 = 旧 BC**（`align_pca=False`，`axis=None`）：斜置 + bbox 轴加载 ⇒ 预期中段弯曲（σ̄ 偏离 + 端部热点），**记录不设门槛**——它验证检查器能分辨「坏 BC」；
  2. **正对照 = 新 BC**（`align_pca=True` + 选定 `load_mode`）：对齐后端带=横截面 ⇒ 必须命中解析解。
- **判据（新 BC 全部 PASS 才算正对照成立）**：

| # | 判据 | 阈值 | 说明 |
|---|---|---|---|
| P1 | 中段（两端各 25% 长度以外）体积平均 \|σ̄_zz − σ_nom\|/σ_nom | **≤5%**（目标 ≤2%） | 解析真值对照，主判据；σ_zz 从 hdf5 voigt 直读（`post_process` 同源读法，bone_feb.py:1126-1147；σ_vm≠σ_zz 的判据不能用 σ_vm 对 F/A） |
| P2 | 中段 σ_zz 变异系数 CV | ≤10% | 均匀性，抓弯曲残留 |
| P3 | 全域 median(σ_vm)/σ_nom | ∈[0.8, 1.2] | 单轴状态 σ_vm=σ_nom；median 对端部奇异稳健（MESH_FIX §6「median 始终物理」） |
| P4 | 全域 p95(σ_vm)/σ_nom | ≤1.5 | 端部奇异护栏（非主判据；棱柱规则网格预期远低于此） |
| P5 | 构建期输入自检：Σ NodalForce 向量模 = F | =F±1e-6 | 载荷没写错（builder meta `load_per_node_n × n_top_nodes`） |
| P6（可选） | 反力平衡 \|ΣR − F\|/F | ≤1% | 仅当新 BC 同时 opt-in reaction 输出时启用（bone_feb.py:916 加 plot var）；否则跳过并注明 |

- **跑多久**：4000 hex 静力解 FEBio 秒级–1 分钟级（参照 BONE_FEB_BUILDER_REPORT 同量级网格、`run_bone_feb` 默认 timeout 900 s）；两工况 + 网格生成 + 出 JSON **合计 ≤10 min**（估算，未验证）。
- **落盘**：原始输出 tee → `temp/opensim_fe/b1/logs/patch_test_<variant>.log`；结果 JSON → `temp/opensim_fe/b1/patch_test.json`；PASS/FAIL 汇总行 → `results/opensim_fe/B1_LOG.md`。脚本 FAIL 时非零退出（将来可选挂 pytest，**新增**测试文件不违反「不改现有测试」红线——是否加留 B1.1）。
- **通过后**才能进入 B1.3 矩阵重跑。

### C.2 辅方案（回归锚，非已知解）：calcaneus 新旧 BC 对比

- **脚本**：并入 `bc_patch_test.py` 的第二子命令，或新建 `scripts/opensim_fe/bc_calcaneus_ab.py`。
- **对象**：跟骨 = 旧矩阵里唯一 p95 物理者（0.37×，`MESH_FIX_REPORT.md` §5）；甲原生网格现成（`temp/opensim_fe/fracture_matrix/calcaneus_r/`，1544 节点/3957 solid 单元，混拓扑，历史 2 attempts 收敛于 25354 N）。
- **做法**：同一 `extract_mesh` 结果跑两遍——旧 BC（现状参数）与新 BC（`align_pca`+`load_mode`），载荷取 `load_chain.npz` 中 subtalar_r @5 m（25354 N，历史 ref），**写入新 temp 根不覆盖旧缓存**。
- **判据（回归性质，无解析真值，标注为合理性锚）**：新 BC median/σ_nom ∈ [0.3, 1.0]（MESH_FIX §6 median≈0.3–0.7× 基线）；新旧 median 相对差 ≤30%；新 BC p95/σ_nom ≤2×。任一越界 ⇒ 新 BC 在紧凑骨上引入了新偏差，回 B1.1 复盘。
- **跑多久**：单骨双解，估 ≤15 min（历史 attempts 少；未验证）。

---

## D. 执行序建议（B1.1→B1.5 可执行展开）

> 日志约定（沿 ALIGN_B）：进度行追加 `results/opensim_fe/B1_LOG.md`；原始输出 tee 到 `temp/opensim_fe/b1/logs/*.log`（目录本次已建）。

### D.1 步骤表

| 步骤 | 可执行内容 | 谁做 | 预估 | 日志 |
|---|---|---|---|---|
| **B1.1 BC 设计评审（拍板会）** | 读本报告 §B/§C，四选一拍板：① 候选 (a) vs (b)；② `load_mode` 取值（tributary/prescribed/维持 equal）；③ 底带 fixed 是否换 roller/spring（参考 PLANTAR_BC 结论：法向刚度主导）；④ 统计量 = `--stat median`（既有 flag，不改默认）+ 新 BC 重跑的 temp/结果命名（建议 `temp/opensim_fe/fracture_matrix_bc2/` + `results/opensim_fe/fracture_matrix_bc2.*`） | **主会话 + 用户** | 0.5–1 h | 决策行 → `B1_LOG.md`；决策同时补进 `ALIGN_B_TODO.md` 决策记录（主会话执行） |
| **B1.2 BC 实现 + 正对照** | 按拍板实现 opt-in（或新 builder）→ 写 `bc_patch_test.py` → 跑 C.1 双工况判据 P1–P5（+P6 可选）→ 跑 C.2 calcaneus A/B → 全量回归 `tests/`（`.venv\Scripts\python.exe -n auto --dist load`，MEMORY §🧪 红律；现有测试必须全绿） | **子代理**（主会话复核 PASS 行） | 2–4 h | `temp/opensim_fe/b1/logs/bc_impl_*.log`、`patch_test_*.log`、`calcaneus_ab_*.log`；PASS/FAIL → `B1_LOG.md` |
| **B1.3 9×50 矩阵重跑（新 BC）** | 新驱动（建议 `scripts/opensim_fe/rerun_matrix_bc2.py`，模式参考 `run_remesh_fe.py`/`regenerate_fracture_matrix.py`，**不改** `fracture_matrix.py`）：设独立 `OUT_TEMP` 副本根 → **拷贝 `load_chain.npz` 进新根**（12.5 KB，heights/joints 校验通过即命中 :280-284 复用，免 ~21 min 载荷链）→ 9 骨 FE 全新跑（新根无旧缓存）→ `--stat median` 出 `fracture_matrix_bc2.{json,csv,png}` + 报告 → 与旧矩阵 diff（排序 Spearman、首折 p95/median 对照、逐骨 σ 比值表） | **子代理**（长任务，可后台） | 载荷链 0 min（缓存）；FE 估 0.5–2 h（**未验证**——旧 `_fe.json` 不记时长，fibula 历史 7 attempts 是尾部风险）；diff+报告 +0.5 h | `temp/opensim_fe/b1/logs/matrix_bc2_*.log`（驱动 stdout）；进度/完成行 → `B1_LOG.md` |
| **B1.4 S4 对标** | 输入按 §D.2 拍板结果 → **复用** `paper_alignment.py:161-229`（`logistic_fit` 含完全分离检测 + Jaccard 平均联结聚类）与 :499-507（A3 OR 排序判定）→ **新增 A4 PCA 双路径**（grep 证实 `paper_alignment.py` 无 PCA/协方差代码，S4.3 需补写，估 +2–4 h）→ 按 A1(±2 m)/A2(≥4/5)/A3(±30%)/A4(可分) 出判定表（判据原文：`docs/methods/OpenSim_FE复现方案.md` §六 表，A5/A6 参考位） | **子代理**写码 + **主会话**核判定 | 0.5–1 天（含 PCA 新增与试跑） | `temp/opensim_fe/b1/logs/s4_*.log`；判定行 → `B1_LOG.md` |
| **B1.5 收尾** | mesh 收敛声明（A2 移交项，若适用）+ 论文/交接文档/出处索引同步（`OpenSim_FE交接.md`、MEMORY 索引）+ 经验落盘（`.learnings/LEARNINGS.md`：BC 设计结论、patch test 判据） | **主会话** | 1–2 h | 完成行 → `B1_LOG.md` |

### D.2 S4 输入路线建议（甲 vs 乙）——**待主会话/用户拍板**

**甲 = 1D 名义 risk(h)（BC-free）现状**（两脚本均为只读追加型、不跑 FEBio、不改默认）：

| 脚本 | 输入 | 输出 | 9 部位？ | risk(h) |
|---|---|---|---|---|
| `risk_1d_loadshare.py`（Route-1 + 载荷分担：踝胫/腓并联 k=E·A/L 分流、轴向骨 mass-above 模型） | `fracture_matrix.json` + `bc_robust_metric_route1.json`（docstring :37-41；`FORCE_SOURCE` 表 :79-93） | `results/opensim_fe/risk_1d_loadshare.json` + `RISK_1D_LOADSHARE_REPORT.md` | ✅ 9/9（FORCE_SOURCE 全覆盖） | ✅ `rows_*[*].risk_by_height` 50 元素（实测） |
| `risk_1d_bending.py`（梁弯曲项 σ=F·cosθ/A + M·c/I，θ 取 MESH_FIX §5 实测/FE 网格 PCA） | 上者 + `fracture_matrix.json` + route1 + remesh 网格 + `<b>.hdf5`（docstring :46-56） | `results/opensim_fe/risk_1d_bending.json` + `RISK_1D_BENDING_REPORT.md`（2026-10-04 已生成） | ✅ 9/9（`theta_deg` 9 键实测） | ✅ `rows_primary_beta4[*].risk_by_height` 50 元素（实测）；β∈{2,4,8} 敏感性齐备 |

- **支持甲的证据**：`BC_ROBUST_METRIC_REPORT.md` §6.3 明确「S4 输入建议改用 Route 1 的 1D 名义 risk(h)（BC-free、9 骨区分度好）；当前矩阵 p95/二值结果不应直接进 S4」；产物现成，今天就能喂。
- **甲的硬伤（比 BC_ROBUST 更晚拿到的证据）**：`PAPER_ALIGNMENT_REPORT.md` §4.1/§4.2/判据 b——1D 矩阵**低高度即饱和**，可骨折部位全部**完全分离** ⇒ 腓骨/骨盆/颅骨**取不到 OR（A3 结构性失败）**；Jaccard 聚类只出 3 类（论文 5 类）。甲路线按现状**过不了 A3、A2 命中也存疑**（`PAPER_ALIGNMENT_REPORT.md:213`：真实不对齐，非计算失败）。
- **乙 = FE 矩阵重跑**：BC_ROBUST §6.4 要求「方向沿骨轴 + 修正折算（已做）+ BC-健壮统计量（避免 p95）」⇒ 即 B1.3；且 `paper_alignment.py` 的 Route-FE 分支**写死 p95**（:336-362、:425），换 median 需同步加 stat 开关（B1.4 内做，子代理顺手）。
- **建议（待拍板）**：**甲立即跑通 S4 管线**（验证 Logistic/聚类/判据代码路径，产物已有）；**正式 A1–A4 判定以乙（B1.3 新 BC 矩阵）为输入**，甲作为 BC-free 交叉参照并列报告。理由：A3 在甲上结构性失败、A1/A2 在饱和矩阵上不可判；而乙是三份红线报告共同指向的修复路径。折中变体：A1/A2 甲乙双跑对照、A3/A4 只认乙。
  **——以上为建议，最终路线「甲 / 乙 / 甲+乙」待主会话/用户拍板（B1.4 前必须定）。**

---

## E. 风险 / 不确定点

1. **甲原生网格 ≤5 万约束：已满足（实测）**。逐骨规模（读自 `temp/opensim_fe/fracture_matrix/<b>/<b>_fe.json` 的 mesh 块，9 骨 status 全 ok）：

   | 骨 | 节点 | solid 单元 | shell 单元 | 骨 | 节点 | solid | shell |
   |---|---:|---:|---:|---|---:|---:|---:|
   | femur_r | 8,536 | 25,909 | 0 | T6 | 933 | 3,186 | 712 |
   | tibia_r | 7,245 | 23,550 | 0 | parietal_r | 2,952 | 1,430 | 2,860 |
   | R_HIPBONE | 3,175 | 11,767 | 2,232 | C5 | 550 | 1,660 | 466 |
   | L3 | 2,203 | 8,792 | 1,340 | calcaneus_r | 1,544 | 3,957 | 0 |
   | fibula_r | 1,837 | 3,820 | 0 | | | | |

   上限来源：`extract_part.py:59-61`（`HDF5_TET_LIMIT = 50_000`，xplt→hdf5 attribute 上限）+ `S31_SCOUT_REPORT.md:146-147`（LRN-027 坑 4）。最大 25,909 ≪ 50,000，且 9 个 `.hdf5` 均已存在（写成功即证通过）——**满足**。
2. **S3.1 批量网格化（scout 体素路线）在新方案里：原则上不需要**。新方案仍是「甲优先」（ALIGN_B 约定），9 骨原生直读网格现成；S3.1 的批量前提 gmsh 尺寸控制失效**仍未解决**（`S31_SCOUT_REPORT.md:129-131`【阻塞级】：`Mesh.MeshSizeMax` 不生效，任何真骨 >5 万 tet）。乙路线降级为**按需单骨兜底**，且 MESH_FIX 的干净 remesh 网格可直接复用（`temp/opensim_fe/fracture_matrix_remesh/meshes/`）。触发兜底的条件：某骨甲网格在新 BC 下不收敛且旧 remesh 也不行（现实里腓骨是唯一高危者，见 6）。
3. **`extract_part.py` / `shell_io.py` 的角色**：`extract_part.py` = 乙路线 scout 工具（任意 PID → STL → tet，三条网格路径、体素 Delaunay 兜底 :180/:227-232，内建 5 万上限）——新方案中**仅兜底时用**；`shell_io.py` = THUMS `*ELEMENT_SHELL` 流式读取器（中轴骨/颅骨皮质壳，`shell_io.py:1-35`），是甲路线 `bone_feb.extract_mesh` 的直接依赖（`bone_feb.py:33,75-76`）——B1.3 重跑会继续经由它，**与 BC 改动正交，不需要动**（9 骨壳域读取历史成功，见上表 shell 列）。
4. **50 高度载荷链缓存：可复用（已实证）**。`temp/opensim_fe/fracture_matrix/load_chain.npz` 存在（2026-10-04，12,502 B），实测 `heights` = 1..50 精确、4 关节 `__vert/__abs/__t_ms` 全齐 ⇒ 命中 `fracture_matrix.py:280-284` 复用分支，免重跑省实测 617 s（FD）+ 637 s（关节反力）≈ **21 min**。⚠️ 换新 temp 根时须**手工拷贝**该文件进新根，否则新根 miss 会重跑整条载荷链（`--force-loads` 语义同）。
5. **FE 缓存碰撞（新 BC vs 旧产物）**：`OUT_TEMP` 是常量（`fracture_matrix.py:83`），FE 缓存路径 `<b>/<b>_fe.json`（:390-394）。若图省事在原根上 `--force-fe`，会**覆盖旧 BC 的 9 份 FE 记录**（MESH_FIX §8.6「未覆盖既有产物」的自我约束即破）。⇒ B1.3 必须走独立根/新驱动（D.1 已列入）；`--force-fe` 只在新根内使用。
6. **腓骨收敛/屈曲**：旧 BC 下 7 条 attempt 记录，主口径降载到 **1000 N** 才收敛（`fe_solve_load_n=1000`、降载因子 25.27；500 N 为线性自检工况，自检比值 1.71–1.89 ≠ 期望 2.0，实测非线性）；MESH_FIX §7 轴向对齐后**全候选载荷仍发散**（细杆屈曲 + 数值脆弱，同一 1 kN 两次结果不一致）。新 BC 若仍发散，回退策略 = `solve_bone_fe` 既有降载 + 线性折回（`fracture_matrix.py:377` 起的候选循环，折回现已乘对）或转乙 remesh；**腓骨在 B1.3 的结果状态要单独如实报告**（未验证能否过）。
7. **统计量切换的下游耦合**：主统计量 p95 是 `PRIMARY_STAT` 默认（:93），改默认被红线禁止 ⇒ 用 `--stat median`（:1021-1022，opt-in 合法）；但 `paper_alignment.py` Route-FE **写死 p95**（:336-362、:425、:453），`PAPER_ALIGNMENT` 判据表也基于 p95 ⇒ B1.4 需同步给它加 stat 参数（子代理小改，新参数不破坏默认）。
8. **`time_steps>1` 与新载荷方式的相性**：`fracture_matrix.py:1023-1025` help 明写 >1 会让现 NodalForce 配置误 ramp 发散（bone_feb docstring :737-744 又建议软骨用增量——两处并存的历史矛盾）。若新 `load_mode` 要走增量加载，必须先小样本验证 ramp 行为（未验证）。
9. **旋转坐标系对下游约定的影响面（未验证）**：`align_pca` 后节点坐标系变 ⇒ `risk_1d_bending` 的 θ 表（THETA_LIMB_DEG 取自 MESH_FIX §5 旧坐标系）、`paper_alignment` 的「bbox 最长轴=载荷轴」约定、以及任何读 `<b>.hdf5` 坐标的脚本（`make_fracture_figures` 等）在**新产物目录下**需重核。缓解：新根产物与旧产物目录隔离，旧链路完全不动。
10. **C.2 无解析真值**：calcaneus 新旧对比是回归锚不是已知解；其「0.37× 物理」本身是 MESH_FIX §5 的口径结论（A_cort 名义），跨口径引用时注意（MESH_FIX §8.3：A_fe 与 A_cort 两种名义并存）。
11. **未验证清单（诚实边界）**：新 BC 求解时长；`nodal_tributary` 分摊面积算法精度；B1.3 全程耗时（旧记录无 FE 时长字段）；A2 mesh 收敛声明对本线是否适用（取决于 B1.4 输入路线拍板）；A4 PCA 新增代码的数值口径。
