# OpenSim + FEBio 复现 [Y25] —— 进度交接（供上下文压缩后继续）

> 主方案：[`docs/OpenSim_FE复现方案.md`](OpenSim_FE复现方案.md)（本文件是它的**执行状态快照**）。
> 阶段性总结：[`docs/OpenSim_FE阶段性总结.md`](OpenSim_FE阶段性总结.md)（一页看"到哪了/什么可信/还没解决"）；**总报告：[`docs/汇总报告_骨与韧带损伤建模.md`](汇总报告_骨与韧带损伤建模.md)**。
> 更新于 2026-10-05。**进度：论文对齐（核心 ✓）+ 临床锚定 + 场景判据 R1–R6 + 非垂直 S1–S4 全部完成；S5（FE 脚-垫接触）计划。新增：论文 v2 稿 + 全套配图（图 1–9）+ 坠落仿真动画。** 细条目见 §7，**配图/动画见 §4b**。

---

## 1. 目标与架构（不变）

用 **OpenSim（多体动力学）+ Python/FEBio（有限元）** 复现 [Y25]《成人直立位双足坠地跨
高度骨骼损伤机制的有限元分析》的**方法链**，再把"刚性地面 1–50 m"移植到"室内抱石软垫
≤4.5 m"。

```
GRF(pad.py) → OpenSim 正动力学(关节反力) → load_transfer(关节 wrench→FE BC)
            → FEBio 单骨 von Mises → 50×9 骨折矩阵 → Logistic/聚类对标 [Y25]
```

**根本前提**：OpenSim 是多刚体，**算不出 von Mises**；必须显式拆成"多体出动载荷 + 单骨 FE 出应力"。

---

## 2. 环境（本机已核实）

| 项 | 值 |
|---|---|
| 解释器 | `D:\Project\climbing_fall_analysis\.venv\Scripts\python.exe`（Py 3.12） |
| 运行约定 | 仓库根、`$env:PYTHONPATH="src"`（脚本自带 sys.path 注入） |
| OpenSim | 4.6（pip wheel） |
| pyfebio | 0.3.0（**只能写不能读**，单向 Python→.feb→febio4） |
| FEBio | `D:\Program\FEBioStudio\bin\febio4.exe`（4.13.0，不在 PATH） |
| gmsh SDK | `D:\Program\gmsh-sdk\gmsh-4.15.2-Windows64-sdk`（PATH+sys.path+`gmsh.initialize(LIB)`） |
| 其他 | meshio 5.3.5 / pyvista 0.49.0 / trimesh 5.1.1 / scipy / numpy / matplotlib |

---

## 3. 已完成（S0 – S3.C）

- **S0 环境 ✅** — `febio_run.py` 的 `find_febio()/probe_febio()/run_febio()`；`tests/test_febio_run.py` **8/8**。
- **S1 多体半段（S1.1–S1.3）✅** — GRF→FD→距下关节反力。@5 m：T1 **+5.8 %**、T2 **0.00 N**、
  距下关节 **26.4 kN @ 20.1 ms**；`tests/test_opensim_fe.py` **11 passed / 1 skipped**。
- **S1 全链路（S1.4–S1.5）✅** — 解剖级跟骨 STL + 面压 `p=F/A` → FEBio 正常终止。
- **S2 全高度扫描 + 过程区正则化 ✅** — `s2_height_sweep.py` 扫 1…50 m：
  T9 单调 **18/18 ✅**、T1 **+5.6~+9.3 % ✅**。定义 **ℓ=4 mm**（RNG/RVE 依据：松质骨
  RVE≈5 mm，皮质骨过程区 1–5 mm）作操作判据；`fe_post.von_mises_gauge` 给体积加权下界。
  图 `results/opensim_fe/s2_sigma_h.png`。
- **S3.0 前处理（decimate + gmsh）✅** — `meshing.stl_to_tet_gmsh` 现为**主路径**：
  抽面 `reduction=0.80`（6718→1342 面，固定几何）→ `classifySurfaces(0.001 rad)` →
  `generate(2)` 重网格化 → `generate(3)`。质心 J<0.3 占比 **17.2 %→5–8 %**、体积误差 **0.00 %**。
  **gauge(ℓ) 在 cl=2.5–4 mm 收敛：48.7 / 48.7 / 53.4 MPa ✅**（正则化目标达成）。
- **S3.2 载荷传递 + S3.2b 跟腱力 ✅（结论：跟腱非主因）** ——
  - `load_transfer.py`（S3.2）：关节力 + 跟腱力 → `FeLoadSpec`；`F_ach = a·F_max(triceps)`，
    **F_max = 10885 N**（`probe_achilles.py` 实测 soleus+gasmed+gaslat）。
  - `meshing`：新增 `achilles` 面（后上结节，x<−25 & y>0；实测质心 (−39.8, 3.9, 1.6) mm）。
  - `febio_model`：新增 `achilles_n`/`achilles_dir` → `TractionLoad`，方向 ≈(0.15, 0.99, 0)。
  - **实测**（`s3_achilles_sweep.py` @5/10/20 m，a=0 vs 1）：**gauge_max 逐位相同**
    （26.31/42.58/68.39 MPa），仅体平均 +6~7%。诊断：后上区 σ_vm 随 a 升
    （a=0/1/5 → 0.2/20/101 MPa），但**峰值恒在距下关节棱边**。
  - **⇒ 跟腱在生理范围（a≤1, ≤10.9 kN）不改变跟骨峰值应力 → 不是 7–9 m 差距的原因。**
    真正的嫌疑转为**距下关节均匀压力+跖面全固定造成的边缘奇异**（§7 新增 S3.2c）。
  - ⚠️ 网格坑：gmsh 平滑路径单骨 **>5 万 tet**（cl=4→59k, cl=5→58k，`MeshSizeMax` 不被
    尊重）→ 触发 HDF5 attribute 上限；S3.2b 扫描用 **`vtp_to_tet` 体素路径**（6.2k tet）。
  - ⚠️ **几何来源**：`calcaneus_r.stl` 属于 **IITD 踝足模型**（`model/Model/`，
    `IITD_Ankle_foot_model_NLSRclass.OSIM`），**不是** Rajagopal calcn_r 帧；OpenSim
    的跟腱附着点 (−110,46,79) mm 落在 STL 之外。故跟腱**位置/方向按 IITD STL 解剖确定**，
    不映射 Rajagopal 路径。
- **S3.C（方向 C：借 THUMS AM50 跟骨几何+材料）✅ 完成，未复现绝对高度** ——
  - **抽取**：`kmesh_io.py` + `extract_thums_calcaneus.py` → `temp/opensim_fe/thums_calcaneus/`。
    THUMS 甲板 ASCII；右跟骨 = `R_CALCANEUS_SPON`（pid 81001200，3323 塌缩六面体记法
    四面体，`*MAT_ELASTIC` RO=1e-9、E=73.4 MPa、ν=0.45）+ `R_CALCANEUS_CORT`
    （pid 81001300，634 hex，`*MAT_PIECEWISE_LINEAR_PLASTICITY` E=15000 MPa、ν=0.3、
    σy=224 MPa）；体积 **96.2 cm³**。甲板内**无**针对跟骨的 `*MAT_ADD_EROSION`
    （仅有的侵蚀卡挂在 `*MAT_SIMPLIFIED_RUBBER` 韧带且全为 0）；骨折 override 只改股/胫/腓。
  - **配准（两版）**：(a) `registration.py` ICP 到 `r_foot.vtp` → **发现 `r_foot.vtp`
    就是 `calcn_r` 挂的整只脚网格（193.6 mm），模型无独立跟骨几何**；ICP RMS 5.65 mm，
    原点落骨后缘。(b) `registration_anatomical.py` 解剖放置：距下关节面质心→原点、
    跖面→−Y、前突→+X；R 正交 det=1；载荷块 **222 mm²(14 面) → 2495 mm²(119 面)**。
    calcn local→ground 旋转 = **单位阵**（X前/Y上/Z内外、原点=距下中心）。
  - **FE**：`thums_feb.py` 双材料域 cortical（hex8 E15000/ν0.3）+ trabecular
    （tet4 E73.4/ν0.45）；距下 `PressureLoad p=F/A`、跟腱 `TractionLoad`、跖面三向固定；
    FEBio 4.13 **rc=0**。判据 = 皮质 σ_vm ≥ 150 MPa。
  - **扫描**：`s2_thums_sweep.py --mesh anat` → `results/opensim_fe/s2_thums_height_sweep_anat.{json,csv,png}`、
    `S2_THUMS_REPORT_ANAT.md`；**38/38 rc=0**。**a=0 首次骨折 ≈2.0-3 m**
    （gauge 1m=114.9 / 2m=149.9 / 3m=156.6 MPa）；**a=1 gauge ≡ 裸 max**（跟腱边缘奇异）→"1 m"非物理。
  - **结论**：**未复现 [Y25] 足部骨折簇 7-9 m**（偏高 ~3-4×）；残余差距在 **BC 层面**
    （跖面全固定奇异 / 150 MPa 判据 / 跟腱 TractionLoad 边缘 / 多体 F_subt 量级），
    超出 C 的"不修动力学/BC"边界。
  - ⚠️ **caveat**：旧 S1-S3 `"calcaneus"` 域实为 `r_foot.vtp` **整足网格（121.7 cm³）**，
    新模型才是真实 AM50 跟骨（CORT 18279 + SPON 77794 mm³）；**数值不可直接比较**。

### 模型成分 / 假设清单（方向 C 跟骨模型 · 必读）

> **一句话**：方向 C 模型 = **AM50 跟骨几何 + 材料** × **Rajagopal2015 载荷** × **我们自定的 BC/判据**。
> 它**不是完整 AM50 仿真**，也**不是 [Y25] 复现** —— 只能看趋势/方向，**绝对高度不可移植**。

| 成分 | 值 / 来源 | 是否 AM50 |
|---|---|---|
| 跟骨**几何** | THUMS 甲板 PID 81001200(SPON) / 81001300(CORT) | ✅ AM50 |
| 材料**常数** | CORT E=15000 MPa/ν=0.3；SPON E=73.4 MPa/ν=0.45 | ✅ AM50（仅常数） |
| 材料**本构** | **线弹性**（只用 E/ν）；THUMS 的 `*MAT_PIECEWISE_LINEAR_PLASTICITY` σy=224 **未用** | ❌ 简化 |
| 松质单元 | 塌缩六面体 → **tet4** | ❌ 近似 |
| **动力学/载荷** `F_sub(h)` | **Rajagopal2015**（通用模型，75.34 kg）正动力学 | ❌ 非 AM50 |
| 跟腱力 | `a·F_MAX_TRICEPS_N`（=a·10885 N，Rajagopal 三头肌 F_max 上界） | ❌ 非 AM50 |
| **边界条件** | 跖面全固定 + 距下均匀压力 `p=F/A` + 跟腱 `TractionLoad`（见下） | ❌ **自定义** |
| **骨折判据** | 皮质 σ_vm ≥ **150 MPa**（[Y25] 表值），**后处理**、无单元删除 | ❌ 自定义 |
| 分析类型 | FEBio **STATIC** 单步；**无接触、无动力学** | ❌ 简化 |

**BC 用的是谁的**：**完全自建**（`thums_feb.py` / `febio_model.py`），与 THUMS 甲板、LS-DYNA 无关，
也不是照抄 [Y25]。THUMS 甲板里的 `*BOUNDARY_SPC` / `*CONTACT` / 初速度**一点没用**（只借几何+材料，
且**不跑 LS-DYNA**）。三条 BC 具体为：

- 跖面：`BCZeroDisplacement(node_set="plantar", x_dof=y_dof=z_dof=1)` —— 三向全固定（刚性地面理想化）；
- 距下关节：`PressureLoad p=F/A`（本次走 `use_rigid=False` 压力回退路径；备用 `use_rigid=True` 的幽灵刚体
  `RigidForceLoad`）；方向 `(0,-1,0)` 来自 `load_transfer.SUBTALAR_DIR`；
- 跟腱：`TractionLoad`，方向 `(0.15,0.985,0)`，大小 `a·10885 N`；
- OpenSim **只提供载荷大小** `F_sub(h)`，**不提供** FE 边界。

⇒ 因此方向 C 的 `gauge_max(h)` 只做**筛选/趋势**；`~2–3 m` 与 [Y25] `7–9 m` 的差是 **BC 层面**的
（跖面全固定奇异 / 均匀压力 / 跟腱边缘 / 判据），**不是**几何或材料的差。

### 场景 / 姿态假设（整条管线 · 必读）

> 与上面的"模型成分"无关：这是**载荷 / 姿态**层面的假设，**贯穿 S1 – S3.C 的全部结果**。
> 引用任何高度 / 应力数字前，先确认场景是否匹配。

| 假设 | 现状 | 代价 |
|---|---|---|
| **落地方向** | **纯竖直**（初速度只设 `pelvis_ty`） | 斜向 / 带水平速度的落地无法表达 |
| **落地姿势** | **直立、双足、腿绷直**（`LEG_JOINTS` 全锁死） | 单脚 / 屈膝 / 躯干倾斜无 |
| **GRF** | **1D 标量竖直力**（`pad.py` 双质点 → `ground_reaction`） | 无方向、无 CoP、无摩擦 |
| **左右分配** | **对称**（`split=(0.5,0.5)`） | 单脚先落无 |
| **关节 wrench** | 3D 力 + 力矩**被 scalarize** 成标量 × 固定方向 | 切向力、力矩丢失 |
| **FE 跖面** | **三向全固定**（`BCZeroDisplacement(plantar,…)`） | 固定端奇异（= `~2–3 m` 差距的 BC 根因） |

⚠️ **现有全部结果都只在"竖直 + 对称双足 + 绷直腿"下成立**；真实抱石坠落（单脚先落、斜落、旋转）
**表达不出来**。扩展到非垂直的路线图见 **`docs/非垂直落地扩展方案.md`**（**S1–S4 ✅ 已实现并回归 / S5 计划**）。
（注：中段 `joint_loads.subtalar_reaction` **已是 3D 力+力矩**，卡住的是两端的 1D 接口。）

### S3.D · BC 诊断与次级因素排除（✅ 完成 2026-10-04）

- **跖面 BC 实验**（`results/opensim_fe/PLANTAR_BC_REPORT.md`）：跖面 fixed 186.76 →
  **三向 Winkler 弹簧 k≈1e3 142.7 MPa（降 24%）**、峰值**离开跖面棱边**；首折外推
  **2.11 → 5.43 m**（**向 [Y25] 7–9 m 靠近但未到达**）。**机制 = 法向（竖直）支承刚度**，
  不是"三向/切向"（真放开切向反而 +14%）。默认 BC 仍 `fixed`。
- **③ 距下表示 + STATIC→DYNAMIC**（`JOINT_LOAD_REPORT.md`）：均匀/梯度/节点载荷 gauge **±10% 内**
  且**方向敏感**；DYNAMIC 瞬态仅 **+1.1%** ⇒ **STATIC 充分**（T8 通过）。峰值恒在跖面棱边。
- **④a 多体 77 kg 校准**（`MASS_CALIB_REPORT.md`）：F_subt(5 m) **+1.3%** → 首折 2.02→1.99 m（**可忽略**）。
- **④b 跖腱膜**（`FASCIA_REPORT.md`）：固定跖面下效应**精确为 0**（被反力吸收）；弹性支承下 +0.28%@5 kN（**小**）。
- **⇒ 结论**：残余差距主因**唯一指向跖面法向刚度固定端奇异**；即便修到最好仍 **5.43 m < 7–9 m**
  → 属**范式差**（D5）。四层对齐与清单见 `docs/对齐_Y25_假设与结论.md`。

### 关键数字（勿再引用错）

| 量 | 值 |
|---|---|
| 跟骨几何体积 | **57.7 cm³**（`model/Model/Geometry/calcaneus_r.stl`，非 r_foot.vtp 的 121.7 cm³"整足"） |
| @5 m 距下关节纵向力 | **26.4 kN**（随 h 增长**快于 √h**：50 m→121.4 kN，α≈0.66） |
| 力控线弹性线性度 | σ ∝ F，实测比 **0.99–1.00** ✅ |
| 首次骨折（**gauge ℓ=4 mm**，平滑网格，**实测**） | **≈26.5 m** |
| 首次骨折（旧体素网格裸峰值） | 7.0 m —— **伪影，勿用** |
| [Y25] 预期 | 7–9 m（**绝对高度不可移植**） |

---

## 4. 关键文件

| 文件 | 作用 |
|---|---|
| `src/climbing/coupling/meshing.py` | **主路径 `stl_to_tet_gmsh`**（decimate+gmsh）；兜底 `vtp_to_tet`（voxel repair） |
| `src/climbing/coupling/febio_model.py` | `build_calcaneus_feb(mesh, path, E_mpa, nu, load_n, use_rigid, time_steps)` |
| `src/climbing/coupling/fe_post.py` | `von_mises_stats`（max/p95/p99）、`von_mises_gauge`（正则化）、`peak_von_mises` |
| `src/climbing/coupling/load_transfer.py` | **S3.2** 关节力 + 跟腱力 → `FeLoadSpec`（`transfer()` / `.build()`）；`F_MAX_TRICEPS_N=10885` |
| `scripts/opensim_fe/s3_achilles_sweep.py` | **S3.2b** 高度 × 激活 a 扫描（gauge > 150 MPa 判首次骨折） |
| `src/climbing/coupling/febio_run.py` | `run_febio`（非 0 退出码抛错，绝不静默） |
| `src/climbing/coupling/opensim_grf.py` / `opensim_fall.py` / `joint_loads.py` | S1 多体半段 |
| `scripts/opensim_fe/s1_vertical_slice.py` | 单高度全链路 |
| `scripts/opensim_fe/s2_height_sweep.py` | 高度扫描（`REG_LEN_MM=4.0`、`--gauge-mm`） |
| `scripts/opensim_fe/s2_plot.py` | σ_vm(h) 图 |
| `scripts/opensim_fe/meshing_driver.py` | VTP/STL→tet→.feb→FEBio→σ_vm 单命令 |
| `scripts/opensim_fe/kmesh_io.py` | THUMS `.k` 甲板 ASCII 解析（节点/单元/材料） |
| `scripts/opensim_fe/extract_thums_calcaneus.py` | 抽右跟骨两 part（SPON/CORT）→ `temp/opensim_fe/thums_calcaneus/` |
| `scripts/opensim_fe/registration.py` | ICP 配准到 OpenSim `r_foot.vtp`（发现：该 VTP 是整足网格） |
| `scripts/opensim_fe/register_calcaneus.py` | 跟骨配准入口 |
| `scripts/opensim_fe/registration_anatomical.py` | 解剖配准：距下关节面质心→原点、跖面→−Y、前突→+X |
| `scripts/opensim_fe/thums_feb.py` | 双材料域 `.feb` 生成（cortical hex8 + trabecular tet4）；可选参数 `plantar_bc`/`joint_load`/`analysis`/`fascia_n`（**默认=现状**） |
| `scripts/opensim_fe/thums_feb_run.py` | 方向 C 单高度 FE 运行入口 |
| `scripts/opensim_fe/s2_thums_sweep.py` | 方向 C 全高度扫描（`--mesh anat`） |
| `src/climbing/coupling/plantar_bc.py` | 跖面 BC：`fixed` / `roller` / `spring`(Winkler)；`apply_plantar_bc()` |
| `scripts/opensim_fe/plantar_bc_probe.py` | 跖面 BC 敏感性探针（→ `PLANTAR_BC_REPORT.md`） |
| `scripts/opensim_fe/joint_load_probe.py` | 距下载荷表示 + STATIC/DYNAMIC 探针（→ `JOINT_LOAD_REPORT.md`） |
| `scripts/opensim_fe/mass_calib_probe.py` | 多体 77 kg 校准探针（→ `MASS_CALIB_REPORT.md`） |
| `scripts/opensim_fe/fascia_probe.py` | 跖腱膜载荷通路探针（→ `FASCIA_REPORT.md`） |
| `scripts/opensim_fe/scout_parts.py` | 甲板 `*PART` 全量索引（PID+名称）→ `parts_index.json` |
| `scripts/opensim_fe/extract_part.py` | 任意 PID → STL → tet（体素路径、坐标无关、`--max-tets` 强制 ≤5 万） |
| `scripts/opensim_fe/make_paper_figures.py` / `make_model_figures.py` | 图 1–7 数据图（管线 / 对齐 / 单脚 / 旋后 / 垫权衡 / 韧带 / 外旋） |
| `scripts/opensim_fe/make_fracture_figures.py` | **图 8/9**：真实骨几何 + 受伤体位渲染（图8 内翻28°+外旋20°、图9 背屈24°）+ 骨折线切片 + 3D 机制箭头 + 空白感知中文标注；输出 4 方位（`--` 无参默认二者都出） |
| `scripts/opensim_fe/make_fall_animation.py` | **坠落仿真动画**：`run_dead_drop` 正动力学 states → 逐刚体 `getTransformInGround` 摆姿 → 3D+GRF 曲线合成帧 → ffmpeg MP4/GIF；`--passive` 出 S3b 解锁对比版 |
| `scripts/opensim_fe/make_pad_paradox_figure.py` | **图 10「垫子悖论」**：`pad_supination_tradeoff.json` → `total_risk(k)` 曲线 + risk=1 阈值 + 刚性极限 |
| `scripts/opensim_fe/make_scenario_animations.py` | **危险姿势演示动画**：`--scenario ser/talus`（三踝 SER / 距骨背屈），垫 vs 刚性 GRF 对比 + 踝部特写内嵌 + 实时 `M_inv`/risk 读数 |

探测/复现脚本在 `temp/opensim_fe/probe_*.py`（`gauge_family` / `linearity` / `decimate` / `gmsh_remesh`）。

---

## 4b. 论文配图与坠落动画（2026-10-05 新增）

> 论文稿：`docs/论文初稿_抱石落地骨与韧带损伤建模.md`（**v2 问题驱动**：《软垫上为何仍会骨折？》；引言→方法→结果→讨论→结论；图 1–9 内嵌；参考文献 **[1]–[36]**）。
> 所有图注：`results/opensim_fe/figures/FIGURE_CAPTIONS.md`（含每图的**数据源 / 标签 / 诚实边界**）。
> 复现脚本：`scripts/opensim_fe/{make_paper_figures,make_model_figures,make_fracture_figures,make_fall_animation}.py`。

| 图 | 文件 | 内容 | 数据源 |
|---|---|---|---|
| 1–7 | `figures/fig1…fig7*.png` | 管线 / 对齐 / 单脚 / 旋后(R1) / 垫权衡 / 韧带 / 外旋(S6) | 各自 JSON（`axial_subregion`…`ankle_external_rotation`） |
| **8** | `figures/fig8_trimalleolar_ser.png` | **三踝骨折机制**：真实骨几何 + **SER 受伤体位（内翻28° + 外旋20°）**，**J 后方观** | 无（机制示意）+ 外踝臂定量见 `ankle_external_rotation.json` |
| **9** | `figures/fig9_talus_mechanism.png` | **距骨骨折机制**：真实骨几何 + **背屈24°** 体位，**G 内侧观** | 无（机制示意） |
| 动画 | `figures/fall_simulation.mp4` / `.gif` | **全身坠落仿真过程**（自由落体段 + 冲击段慢放 + GRF(t) 曲线面板） | `run_dead_drop` 正动力学 states |
| **10** | `figures/fig_pad_paradox.png` | **「垫子悖论」**：旋后着地时 `total_risk(k)` 随垫子刚度单调递减（越软越危险，无内部最优） | `pad_supination_tradeoff.json`（R4） |
| **动画 A** | `figures/fall_scenario_ser.mp4` / `.gif` | **三踝 / SER 姿势演示**：单脚 + 内翻30° + 外旋25°，垫 vs 刚性 + 踝部特写 | `ankle_supination` 口径（见 FIGURE_CAPTIONS 附 2） |
| **动画 B** | `figures/fall_scenario_talus.mp4` / `.gif` | **距骨 / 背屈楔入演示**：单脚 + 踝背屈25°（真实锁定姿势） | 同上（机制演示，无量化判据） |

> ⚠️ **载荷口径（务必注明）**：两个场景动画默认是**单脚** `split=(1,0)` + 77 kg（比 R1/S5 报告严格 ≈2×：518 vs 257 N·m）。要**精确复现 S5 报告数字**用 `--two-foot`（`split=(0.5,0.5)` + 75.337 kg → `fall_scenario_ser_twofoot.mp4`，M_inv ≈264 N·m）。
> ⚠️ **动画 ≠ 计算的部分**（画面已标注）：① 自由落体段是纯运动学插值，非 FD；② 垫子压缩是**渲染约定**（FD 无接触约束，16 cm 是自由行程）；③ SER 的「内翻/外旋」**不是模型自由度**，画面足部为中立位；④ 右侧刚性虚线仅作**载荷**对比，姿态仍是垫工况解。
> 完整审计表见 `results/opensim_fe/figures/FIGURE_CAPTIONS.md` §附 2。

**图 8/9 的做法（可复用）**：THUMS AM50 真实骨网格（`temp/opensim_fe/parts*/…/*_surface.stl`；**胫骨用 gmsh 流形面 + 多遍 Taubin**，其余不平滑以免足趾起刺）→ 用 `pyvista` 渲染；**骨折线 = 对骨网格做平面切片**（随视角正确遮挡）、向外偏 1.3–3.4 mm 避免共面打架；**机制箭头 = 真 3D 弧线**（Rodrigues 旋转，方向按物理判据）；**标注 = 空白感知自适应放置**（逐像素采样背景、避让已放标签、引线短）。四方位（A 外侧 / D 前方 / G 内侧 / J 后方）任选。

**受伤体位的临床依据（关键，勿再写反）**：
- 三踝/后踝骨折最常见机制 = Lauge-Hansen **旋后-外旋（SER）**；**旋后（supination）= 内翻 = 脚心向内**，**外旋**是其**致伤力矩** → **脚位是内翻，不是外翻**（论文 §3.4/§3.6 已改正）。
- 队列：后踝骨折中**三踝-SER 占 71.0%**（最常见，vs 三踝-旋前外旋 9.1%）— **参考文献 [36]** Li Y 等, *J Orthop Surg Res* 2023;18:507。

**坠落动画的做法**：`ground_reaction(height_m, mass_kg, on_pad)`（`pad.py` 1D 双质点）→ `run_dead_drop(grf, scale_mass_kg=77)` → `FallResult.states`（`TimeSeriesTable`，**列取值用 `getDependentColumn(lb).to_numpy()`**，`getMatrix()` 无 `.get()`）→ 逐帧 `setStateVariableValue` + `realizePosition` + `body.getTransformInGround(state)` 变换 `.vtp` 网格 → `pyvista` 渲染 + matplotlib 合成 GRF 曲线 → ffmpeg。
- 本次：2.0 m / 77 kg / 软垫 → v₀=6.26 m/s、**峰值 34.8 kN、34 ms**（对照刚性地面 42.4 kN / 20 ms）。
- ⚠️ **自由落体段是纯运动学插值（v=g·t），只有冲击段是 FD 解算**；这是 S1「竖直 + 对称双足 + 绷直腿」工况，**不声称是真实抱石坠落运动学**。

---

## 4c. S5（FE 脚-垫接触）状态（2026-10-06）

> 规格：`docs/S5接触方案.md`；**接口契约 + 验收门**：`docs/S5_contact_contract.md`（含「G7 修正说明」）；
> 接触配方（从唯一跑通的 `scripts/ankle_fe/build_feb.py` 抄出，62 处 `file:line`）：`docs/S5_contact_notes.md`。
> **注意三个 "S5" 别串**（契约 §0）：本文是 **S5-contact**（FE 接触），不是 S5-swap（换软垫）、不是 S5-supination（旋后 R1）。

**新增代码**（全部 opt-in，默认路径逐位不变）：
| 文件 | 作用 |
|---|---|
| `src/climbing/coupling/pad_foam.py` | 垫材料：**6 项 Ogden + `pressure_model=2` + `k=1 MPa`**，拟合 `pad.py` 的 σ(ξ)（自由侧面**单轴 FE 复现到 2.54%**） |
| `src/climbing/coupling/pad_mesh.py` | 垫块结构化网格（hex→6 tet），具名面 top/bottom/sides、零负 Jacobian |
| `src/climbing/coupling/pad_contact.py` | 接触 deck：两固体域 + 真实 tri3 facet + `SurfacePair` + sliding-elastic + μ；载荷模式 `displacement`/`pressure`/`force`/`traction` |
| `src/climbing/coupling/fe_post.py` | **追加** `build_hdf5` / `read_displacement` / `read_reaction_forces` / `contact_force_sum`（FEBio 4.13 无名 shell domain 修补） |
| `tests/test_{pad_foam,pad_mesh,pad_contact,fe_post_io}.py` | S5 全套 **86 passed**；G0 硬门槛 **19 passed / 1 skipped** |

**验收实况**（`python -m climbing.coupling.pad_contact`）：
- **G1 ✅** FE 单轴 vs pad 律 **max_rel = 2.54%**，单调。
- **G2 ✅** 刚度扫描 E=5/20/100/1000 MPa → 压入 0.49/0.45/0.30/0.065 mm，单调收敛向刚性极限。
- **G3 ⚠️** 量级已修正（2.4 GN → **2.0 kN**，单位问题已修）；**但位移驱动下本质是牛顿第三定律恒等式**，代码已如实自我标注。
- **G4 ✅** 摩擦锥 `max|F_t|/(μF_n) = 0.44 ≤ 1`。
- **G7 ❌（诚实 FAIL，原因已彻底查清）**：2 m 坠落要求垫压掉 ~85% 厚度（ξ≈0.85）—— 这**超出 `hex→tet4` 泡沫的稳健域**。载荷驱动（580 kPa）**第一步即反演**；位移驱动到 sink ≥ 60 mm 也反演。

**🔴 接触稳健化结论（2026-10-06，决定性）— 详见 `docs/S5_contact_convergence.md`**

- **可靠区间**：**ξ ≤ 0.25**（压头下沉 ≤ 50 mm）→ `rc=0`、**穿透 ≈ 0**（<0.25 mm）、**FE σ / 1D 律 σ = 1.56 → 1.64** 平滑单调。这是**可信的平冲头增强 ≈ 1.6**。
- **G7 工况不可解且不可复现**【**仅对旧 AUGLAG 配方成立 —— 后被新配方部分推翻，见本条下方修订**】：同一配置连跑 3 次 → **negJac 27 / 3 / 15，全部失败**（此前一次 rc=0 是运气）；细化 (6,6,4)/(8,8,6)/(10,10,8) 与提高 penalty 同样反演。原因是**极端压入 + tet4 泡沫畸变**，不是"某个参数没调对"。
  > 🔄 **同日晚间修订**：换成 `laugon="PENALTY"` + `penalty=0.1`（**保留 `search_radius=20`，大而非小**）后 G7 **可复现 ×4**、`rc=0 NORMAL / end_t=1 / negJac=0`、力平衡 ≈1e-7、`F_n ≈ 34.4 kN` 跨 3 网格漂移 **<0.6 %** ⇒ **「不可复现」已解除**；但 **σ/σ_law 仍不网格收敛**（4 项 FAIL）。现行口径见 `docs/S5_QUOTABLE.md` 与 `S5_contact_contract.md` §4 **条 6–8**。
- ⚠️ **因此 15.14 cm / 2.70× 一律不得引用**。此前「2.70× = 大应变平冲头物理」的结论**作废**（它取自收敛包线**之外**的解）；T8 的**材料验证**（同一 Ogden 的均质自由侧面 FE ≡ 1D 律 ≤0.5%）**仍成立**，作废的只是其对高 sink 接触比值的推断。
- **收敛配方已固化进库**：`laugon="AUGLAG"` + `penalty=1` + `tolerance=0.005` + `aug_controls=AUGLAG_CONVERGENT_CONTROLS`（FEBio 4.13 的 AUGLAG **默认子参数会发散**；不传 ⇒ 默认行为**逐字不变**）。
- **探针**：`temp/opensim_fe/s5_t9/{verify_lib,envelope,repro,finish,mesh_refine}.py`。

> **⚠️ 引用限定（2026-10-06，针对本节提及的 1.56/1.59/1.62/1.64 与 34.6 kN/2.205 对）**
>
> **(i) 配方**：1.56–1.64 是 **OLD AUGLAG pen=1** 配方在 tet4 (3,3,2) 单网格上重现历史值
> （≤0.1% 重现噪声）；2.205 / 34.6 kN 来自 **NEW PENALTY 0.1** 配方在 sink=170 mm tet4
> (6,6,4) 上（见 `S5_contact_udg.md` §7 tet4 fallback 表）。
>
> **(ii) 网格未收敛**：1.56–1.64 未做 mesh-convergence 扫描；NEW PENALTY 配方在 3 网格
> (3,3,2)→(6,6,4)→(12,12,8) 上 σ/σ_law 漂移 **+30.62 %**（`S5_g7_verify_mesh.md` §3），
> σ/σ_law **不是 mesh-converged**；固定 sink 时 `F_n ≈ 34.4 kN` 与 `σ_contact ≈ 0.574 MPa`
> 在 3 网格上漂移 **<0.6 %** ⇒ **这两个是当前唯一在 2 m 落地可引用且网格收敛的量**。
>
> **(iii) 不可跨配方对比**：NEW vs OLD σ/σ_law 在 4 个 anchor sinks 上单调偏移
> **−5.85 %…−9.38 %**（`S5_g7_verify_anchors.md` §4.1），约为 ≤0.1% 重现噪声的 60–100× ⇒
> 两条曲线不可叠加。OLD 配方装上 NEW stepper 在 closure `negJac=55` 立即失败 ⇒
> 两配方**不可互换**。
>
> **(iv) 权威裁决**：`docs\S5_QUOTABLE.md` —— `F_n ≈ 34.4 kN`、`σ_contact ≈ 0.574 MPa`、
> `σ_1D-law(ξ)`（材料函数）、element agreement ≤±2.5% @ ξ ≤ 0.65 是当前可引用集；任何 σ/σ_law
> 比值（含 2.205）一律不得作为「收敛」单一值引用；15.14 cm / 2.70× 已撤销。
> 本节保留为**事件时间线记录**，实际引用数字一律以 `S5_QUOTABLE.md` 为准。

**结论（2026-10-06 晚间修订；取代原「波3/波4 在那之前不启动」）**：接触在**中等压入区间稳健可信**（增强 **~1.6**，穿透≈0）。G7 工况的「**不可复现**」已**解除** —— `laugon="PENALTY"` + `penalty=0.1` ⇒ **复现 ×4**、`rc=0 NORMAL`、力平衡 ≈1e-7，`F_n ≈ 34.4 kN` / `σ_contact ≈ 0.574 MPa` 跨 3 网格漂移 **<0.6 %**；但 **G7 的 σ/σ_law 仍不网格收敛**（4 项独立 FAIL）。
⇒ **波3/波4 可以启动，但只准在可信口径内做**：引用一律以 `docs\S5_QUOTABLE.md` 的 **QUOTABLE** 行为准（绝对量 + ξ ≤ 0.25 包线），**任何 σ/σ_law 比值不得作为收敛值**。完整修订见 `S5_contact_contract.md` §4 结论条 6–8。

**未做**：波3（`plantar_bc="contact"` 接进 `febio_model.build_calcaneus_feb` / `thums_feb.py`）、波4（`scripts/opensim_fe/nonvertical_s5_contact.py` + `NONVERTICAL_S5_CONTACT_REPORT.md` + JSON）。

---

## 5. 不变量状态（§七）

| # | 不变量 | 状态 |
|---|---|---|
| T1 | 冲量守恒 `∫GRFdt≈m√(2gh)` <10% | ✅ +5.6~+9.3 % |
| T2 | 重力基线 | ✅ @5 m 残差 0.00 N |
| T5 | 网格质量 J>0.3 | ⚠️ 改善后仍 5–8 % 单元 J<0.3（voxel 17.2 %） |
| T7 | 体积守恒 | ✅ 0.00 %（vs 抽面后参考） |
| T9 | σ_vm(h) 单调 | ✅ 18/18 段非降 |
| T6/T8/T3/T4/T10 | 收敛/静动/平衡/分配/单位 | 部分未系统化（S3 补） |

---

## 6. 关键经验（`.learnings/LEARNINGS.md`，本项目条目）

- **LRN-020** 取消后台任务前先取结果——它可能已跑完。
- **LRN-021** OpenSim 算不出 von Mises——多体-FE 耦合的架构前提。
- **LRN-022** OpenSim/FEBio 工具链实操速查（本机已核实）。
- **LRN-023** 骨 FE "峰值应力"常是**边界奇异**、随网格发散——判据要用分布/正则化。
- **LRN-024** 选几何先验"**标本身份**"（体积/解剖量能识破张冠李戴）。
- **LRN-025** gmsh `classifySurfaces` 失败是**角度敏感**（0.001 rad 可解）。
- **LRN-026** 强度判据必须**正则化**：裸峰值发散、体平均过保守、用过程区尺度 ℓ。
- **LRN-027**（= `LRN-20261003-027`）用**粗网格"对上"文献数字，可能只是正则化副作用**——先做收敛验证再对标。
- **LRN-20261004-027** "**整只脚**"陷阱：`calcn_r` 挂的 `r_foot.vtp` 是**整足**、不是跟骨（用前先量体积辨身份）。
- **LRN-20261004-028** 配准**锚错靶**（整足/错位）→ 载荷面畸变、直接污染判据；配准必须锚**解剖标志**（距下关节面质心）。
- **LRN-20261004-029** `gauge ≡ 裸 max` = **加载边界角点奇异**信号（未正则化），非骨内过程区峰值——报数前先识别并标注。

> ⚠️ 编号为 `LRN-YYYYMMDD-NNN`：`LRN-027`（无日期）= `20261003` 家族；方向 C 三条属 `20261004` 家族，**勿混**。

---

## 7. 下一步

> 🎯 **当前下一步 = 清单⑤：S3.1 九骨批量网格 + S4 Logistic/聚类对标 [Y25]**（`docs/对齐_Y25_假设与结论.md` §4 ⑤）。
> 这是"对结论"的**唯一路径**；绝对高度不可移植（D5）。子步骤：
> ① 定 **9 部位**（足/胫/腓/股/骨盆/腰/胸/颈/颅）几何来源（优先 THUMS 甲板，与 [Y25] 同源）→
> ② `meshing.stl_to_tet_gmsh` 批量（固定 reduction、单网格 **≤~5 万 tet**）→ ③ 扩展 `joint_loads`
> 输出**各关节反力**（现只有距下）→ ④ 泛化 FE builder 到任意骨 → ⑤ 1–50 m 扫描 → **50×9 骨折矩阵** →
> ⑥ **Logistic**（渐进部位）+ **层次聚类**（Jaccard + 平均联结）→ 对标 §六 **A1–A4**。
> ⚠️ 跖面 BC 已诊断（§3 S3.D）：首轮为与既有可比可仍用 `fixed`，或改用 `spring` 作为"真值"对照。

> **S3.1 步① 侦察 ✅ 完成（2026-10-04）** —— 9 部位 PID 全部实测定位 + 3 骨试点 rc=0：
> tibia_r **38209 tet / 体积 +2.03%**；skull_r_parietal 30714；lumbar_L3 44181。产物
> `results/opensim_fe/{PARTS_MAPPING.md,S31_SCOUT_REPORT.md}`、`temp/opensim_fe/parts/`。
> 🔴 **两个阻塞坑（S3.1 批量前须定案）**：① **gmsh 尺寸控制失效**（`Mesh.MeshSizeMax` 不生效，
> 真骨恒 >5 万 tet）→ 只能走体素 `_delaunay_tet`，质量 `J<0.3`≈17%（差于跟骨 gmsh 6.2%）；
> ② **皮质壳缺失**（中轴骨/颅骨 `*_CORT`/`_external_`/`_shell` 是**壳单元**，`kmesh_io` 只读 solid
> → 只剩松质；长骨/足的 CORT 是**实体 hex**，可用）。其余坑：体素过填（L3 +95%）、薄壳参照体积
> 不可用、整段脊柱多体、骨盆 CORT 9 个厚度变体、单位 mm vs m、FE builder 面名跟骨专用。
>
> **坑①（gmsh 尺寸控制）✅ 已解决（2026-10-04）** —— 根因 = `classifySurfaces(默认
> forReparametrization=False)` 把输入面片钉成强制几何实体；修法 = `classifySurfaces(ang, True, True)`
> + 关点/边界/曲率尺寸源（gmsh t10 recipe）。新 helper `scripts/opensim_fe/part_meshing.py`
> （mm-native `mesh_bone_stl_mm`，**未改 `meshing.py`**）。胫骨 cl=4 → **27670 tet / J<0.3=3.66% /
> volerr 1.11% / 2.4 s**（体素路径 17.4%）；跟骨 cl=3 → 17059 / J2.90%（优于旧 6.2%）；
> 顶骨 cl=3 → 23709 / J3.15%（薄壳增厚仍未解）。报告 `results/opensim_fe/GMSH_SIZE_REPORT.md`。
> **坑② 步①（壳侦察）✅ 完成（2026-10-04）** —— 解析器 `scripts/opensim_fe/shell_io.py` + `shell_scout.py`。
> 实测：`*ELEMENT_SHELL` **452920** 条 / 665 PID；骨皮质厚度在 **`*SECTION_SHELL`**（SECID==PID、ELFORM=16），
> **非** `*ELEMENT_SHELL_THICKNESS`；**壳节点与 SPON 实体 100% 共享**（可直接共节点粘接，**无需 tie**）；
> 骨盆 9 变体是**分区铺满**（面积和 59135 ≈ SPON 自由面 59221，99.85%）⇒ **须全取**（只取 `_CORT_1.5` 仅盖 42%）；
> pyfebio **原生支持壳**（`ShellDomain`/`tri3`/`quad4`，壳厚每域单值 ↔ THUMS 逐 PART 单值）。报告
> `results/opensim_fe/SHELL_SCOUT_REPORT.md`。
>
> **坑② 步②（L3 shell+solid FE）✅ 完成（2026-10-04）** —— pyfebio `ShellDomain` 共节点路径**端到端跑通**
> （`scripts/opensim_fe/l3_shell_feb.py`；`temp/opensim_fe/l3_shell/`）。L3 SPON tet4(E=73.4)+CORT quad4 壳
> (E=15000,t=1.39) **100% 共节点、无 tie**；FEBio 4.13 **rc=0**。**皮质效应：去壳后松质 σ_vm 升 16–27×**
> （max 26.6× / p95 19.2× / mean 15.9×）⇒ **纯松质严重失真、坑② 必须解决**（壳吸收 ~95% 内力）。
> 报告 `L3_SHELL_FEB_REPORT.md`。⚠️ 坑：FEBio 4.13 `.xplt` 省略壳域 `PLT_DOM_NAME` → pyfebio `to_hdf5`
> KeyError，脚本**就地补丁**绕过（未改 pyfebio）；网格是 THUMS 塌缩 hex→tet4（退化、质量 caveat）；
> 报告 §3.5 的 1000 N 段**已作废**。
>
> **坑② 步③（泛化 FE builder）✅ 完成（2026-10-04）** —— `scripts/opensim_fe/bone_feb.py`（+`bone_feb_validate.py`）：
> 任意骨 = 1..n 实体域（**按 PID 读 `*MAT_*` 的 E/ν**）+ 0..n 壳域（各自厚度），自动共节点 / 轴选 / 均布轴压。
> **4 骨 rc=0**：L3(1实体+1壳,tet4) / tibia_r(3实体,**0壳**,tet4+**hex8**) / parietal_r(1+**2壳**,hex8) /
> R_HIPBONE(1+**9变厚壳**)。**关键发现**：`*ELEMENT_SOLID` **混塌缩 tet4 与真 hex8**（须按 distinct 节点数分类；
> 误把 hex8 塌成 tet4 → 负 jacobian 直接失败）。**材料独立抽查已核对**（甲板 card1 原文一致：L3 E=40/0.45、
> cort 13020；tibia cort 18000；hip 15/17300）。L3 回归：无壳与原型差 **<2%**（管线等价）。报告
> `BONE_FEB_BUILDER_REPORT.md`。
>
> **各关节反力 ✅ 完成（2026-10-04）** —— 新模块 `src/climbing/coupling/joint_reactions.py`：自由体法泛化到
> **踝/膝/髋/脊柱**（按关节远端子刚体子树）；每关节 **T2 残差 ~1e-12**（✅）、子树质量单调（✅）。@5 m 峰值：
> subtalar 25.4 / ankle 25.3 / knee 22.4 / **hip 15.5** / **lumbar 22.2 kN**。⚠️ **修正既有 bug**：
> `_realized_states` 从不写 state 时间（`time==0`）→ `PrescribedForce` 恒用 GRF(0) → 无外力的 **lumbar 旧值仅
> 1.5 kN（错，已作废）**；subtalar 旧 26.38 → 修 **25.38 kN**（−4%，不改 ①-④ 结论）。⚠️ **刚性腿伪影**：近端 |F|
> 反低于远端（惯性抵消），非活体趋势。报告 `JOINT_REACTIONS_REPORT.md`。**⇒ 下一步 = 用真实关节载荷驱动 9 骨 → 50×9 骨折矩阵。**
>
> **50×9 骨折矩阵骨架 ✅ 完成（2026-10-04）** —— `scripts/opensim_fe/fracture_matrix.py`：9 代表骨各建 FE（`bone_feb`）
> + **逐高度实测**关节载荷（h=1..50，**未**用 GRF 比例缩放）+ 线性 `σ_vm(h)` 缩放 → 9×50 二值矩阵
> （`results/opensim_fe/fracture_matrix.{json,csv,png}`、`FRACTURE_MATRIX_REPORT.md`）。
> ⚠️ **结果尚不可用**：8/9 部位在 **1 m** 即判骨折（p95 σ_vm **257–1169 MPa** ≫ σ_c）——**单骨点载荷 FE 的 σ 是数值量级**
> （退化塌缩 tet4 + 带状点载荷 → 应力比截面真实值高 ~30×；例：股骨 15 kN 应 ~40 MPa 却 p95 1169 MPa）；
> 且 3 骨（tibia/fibula/T6）@5 m 全载 **FEBio 发散**须降载 5–25× 再线性折回（线性自检 1.84–2.67≠2.0）。
> **⇒ 下一步 = 修 FE 网格质量（退化 tet4 → `part_meshing` gmsh-reparam；壳骨做壳节点投影），再重跑矩阵。**
>
> **FE 网格质量修复 ✅ 完成（2026-10-04）** —— `MESH_FIX_REPORT.md`：4 骨干净网格（gmsh-reparam，min_jac 0.017–0.071、
> 体积误差 0.45–2.47%）。**结论证伪任务前提**：σ 虚高**主因不是退化 tet，而是 BC**（斜置骨沿 bbox 轴加载 → 弯曲；
> 端部 30% 夹持 + 点节点力 → 奇异）。证据：同载荷 tibia 新 p95 **反而更高**（细网格把奇异解更尖）、**median 已物理**
> （0.1–0.6× 名义）、**PCA 轴向后 p95 大降**（tibia 94.6×→26.2×、femur 59.8×→16.9×、跟骨 1.5×）。跟骨达 **32 m**
> （原 1 m）；tibia/femur 仍 1 m（BC 主导）；**fibula 屈曲发散**。**⇒ 下一步 = 修 FE BC（沿解剖骨轴加载 + 端面分布力，替代点/带夹持）。**
>
> **BC-健壮性交叉验证 ✅ 完成（bg_6d584c17）** —— `BC_ROBUST_METRIC_REPORT.md`：**当前矩阵相对排序不可信，必须重做**。三条证据：
> ① 无 FE 的 1D 名义应力排序（Route 1）vs 矩阵 p95 仅 **ρ=0.22**；② **端带剔除后 p95 反而升高**（热点在中段弯曲 → **证伪"端部夹持是 p95 主因"**）；
> ③ **真 bug**：`fracture_matrix.py:570/576` 对 3 根降载骨未乘回 `load_reduction_factor`（tibia 5.05×/fibula 25.27×/T6 2.21×）→ 矩阵混入不同载荷应力、且报告与代码不符
> （已核实：570/576 未用 lrf、`extract_main_domain` 不折算）。**S4 输入候选 = Route 1 的 1D `risk(h)`，但其未建模胫/腓载荷分配（fibula 第一是上界伪影）。**
>
> **降载折算 bug 修正 ✅ 完成（bg_fa8e016f）** —— `MATRIX_SCALING_FIX_REPORT.md`：`fracture_matrix.py::compute_matrix` 补乘 `load_reduction_factor`，
> 新增 `regenerate_fracture_matrix.py` 从缓存重放（**无 FEBio/OpenSim**）：tibia 871.29→**4402.78**、fibula 656.87→**16596.45**、T6 777.24→**1716.77**
> （×5.053/25.266/2.209），其余 6 骨逐字段不变；旧产物已备份 `fracture_matrix_prebugfix_backup.*`。**但 p95 相对排序仍不可信**（8/9 仍在 1 m，不变）
> → 真正要修的仍是 **FE 的载荷方向（沿骨轴）** 与 **统计量选择（p95→BC-健壮）**。
>
> **1D 载荷分配 ✅ 完成（bg_32b364b4）** —— `RISK_1D_LOADSHARE_REPORT.md`：修 Route-1 两个载荷伪影。胫/腓按轴向刚度并联 → **tibia 76.8% / fibula 23.2%**
> （纯并联偏高，另给文献 10%/5% 档）；中轴骨按"上方质量"→ L3 89.5% / T6 52.9% / C5 15.1% / parietal 9.5% × lumbar。**fibula 从 #1 掉到 #2**（纯并联；10% 档 → #6）。
> 修订排序：`calcaneus(#1,8m) > fibula(#2,23m) > tibia(#3,37m) > T6(#4,47m) > L3/femur/C5/R_HIPBONE/parietal(>50m)`。
> ⚠️ **绝对高度偏乐观**：抱石量程(≤4.5 m)内"无一骨折"，与 FE 矩阵(1 m 全折)形成**上下夹逼** —— 1D 低估（忽略弯曲/应力集中）、FE 高估（BC 伪影），真值夹在两者之间。
>
> **论文结论提取 ✅（重解析版 `paper/MinerU_markdown_成人直立位双足坠地…md`）** —— 论文=THUMS AM50 全身 LS-DYNA，1–50 m 双足坠地，9 部位 + von Mises/Logistic/层次聚类。
> 判据（表1，压缩）：足150 / **胫两端70** / 胫骨干200 / **腓两端70** / 腓骨干160 / **股骨颈80** / 股骨干220 / 骨盆180 / 脊柱150 / 颅骨160。
> 结论：①轴向分布 + **双路径**；②**阶跃部位**=足/胫两端/股骨颈/脊柱(Logistic 因完全分离失败)；③**渐进部位**=**腓骨两端 OR1.682 > 颅骨 1.576 > 骨盆 1.236**；④**骨干永不骨折**；⑤**5 群集/3 阶段**(1–6/7–9/10–16/17–34/35–50 m)。
> ⚠️ **决定性错配**：我们 1D 一直用**骨干**强度(胫200/腓160/股220)，论文骨折点是**两端/骨颈**(70/80)—— `bone.py` 已有这些 key，只是用错。
>
> **弯曲感知 1D ✅ 完成（bg_50ec18cf）** —— `RISK_1D_BENDING_REPORT.md`：补 `σ=F_axial/A+M·c/I` 后 4 长骨首折 8m→**1m**，8× 缺口由弯曲解释并闭合；**fibula 升 #1**（`fibula>tibia>calcaneus>femur>…`，ρ=0.80 vs 纯轴向）；**与论文"腓骨最敏感"一致**。真实中段截面 vs 等效圆同量级（0.74–1.27×）。
> **⇒ 下一步 = 按表1 重映射子区域（两端/骨颈）+ 重算 9×50 矩阵 + 对齐清单 + S4。**
>
> **论文对齐重映射 ✅ 完成（bg_40fe69c5）—— 但仍未对齐** —— `PAPER_ALIGNMENT_REPORT.md`：按表1 子区域重映射（胫/腓**两端70**、股骨**颈80**）。
> 1D 结果：足/胫/腓/股 **1 m 即折**、腰椎 6 m / 胸椎 18 m、骨盆/颈/颅 >50 m；**骨干也 1 m 折**（弯曲 1D 是全骨上界）→ 与论文"骨干不折" **✗**；
> Logistic 全为 `all_one`/完全分离 → 论文"渐进三件套 OR 排序" **✗**；Jaccard 聚类仅 **3 类**（矩阵早饱和）→ 5 群集 **部分 ✗**。
> **根因（诚实）**：① 单骨 + 把整条关节反力压在孤骨 → 低高度过折；② **弯曲 1D 是上界**，且我们一直把关节力当"**竖向**"施加在**斜骨**上 → 造出**假弯曲**（真实关节力应**沿骨轴**）；
> ③ FE 被 BC 伪影污染；④ **论文的阶跃/渐进、5 群集是全身动力学"涌现"属性，单骨模型结构性产生不了**。
> **⇒ 下一步候选 = 聚焦尝试：FE 沿骨轴真实关节力 + 分布式接触 BC（去掉"竖向力 vs 斜骨"的假弯曲），看 fibula/胫骨应力是否落到论文量级。
>
> **纯轴向+两端 决断性对照 ✅ 完成（bg_7a2a044f）—— 决定性结果** —— `AXIAL_SUBREGION_REPORT.md`：
> **改用纯轴向 σ + 两端强度(70/80)** → **腓骨首骨折 = 7 m（落进论文 7–9 m ✓）**，足部 8m / 胫骨 8m / 股骨 19m → **论文标志性的"7 m 临界高度"复现**；
> 排序 `腓骨(7m)>足(8m)>胫(8m)>股骨(19m)>胸椎(47m)>腰/颈/骨盆/颅(>50m)`。**弯曲/FE 路线连 fibula 都对不上（1m）→ 纯轴向是目前最接近的。**
> **但仍未全对齐**：❌ 颅骨/骨盆 1–50m 从不骨折（论文=渐进；1D 用整体截面抓不到**枕骨大孔局部/薄壳弯曲**）；❌ 骨干 23–37m 仍折（论文=不折；单骨把整条关节力压在骨干→过应力）；❌ OR 序（`腓>骨盆>颅`，与论文相反）。
> **⇒ 战略结论：单骨管线无法"整篇复现"论文；能对齐的核心 = 7m 阈值 + 双路径 + 相对序，恰好是场景探索所需的。**
>
> **对齐基线 + 场景计划 ✅（C 收口）** —— `docs/论文对齐基线与场景探索计划.md`：冻结**纯轴向+两端**基线（腓7/足8/胫8/股骨19m；7m 阈值 ✓）；
> 局限=颅骨/骨盆渐进 ✗、骨干不折 ✗、5 群集 ~。**建议首个场景 S1 = 垫子/地面刚度（Winkler k）→ 首折高度曲线**（抱石核心问题，复用 `plantar_bc`）。
>
> **场景 S1 垫子刚度 ✅ 完成（bg_61f1ef7b）** —— `PAD_STIFFNESS_REPORT.md`：垫子越软 → 首折越高。FE 实测折减 R_min=**0.7642**(k≈1e3, −23.6% 应力)；
> 最好垫子 vs 刚性：足 8→**11m**、腓 7→**10m**、胫 8→**11m**、股骨 19→**31m**。与既有 FE 跟骨 2.11→5.43m(×2.58) **方向一致、绝对不同源**。给了**有界范围**[8,11]m 与敏感性，未用单一误导数字。**下一步场景候选：S2 落地姿态 / S3 体型 / S4 抱石量程细扫。**
>
> **非垂直扩展 plan-S1（矢量通路）✅ 完成并验真（bg_177ea3ae）** —— `NONVERTICAL_S1_REPORT.md`：GRF 方向化（`ground_normal`/`tilt_deg`）+ 足底三分量力 + `load_transfer` 保留完整 **3D wrench**，
> **全部 opt-in、默认复现旧行为**。**轴向回归硬门槛 PASS**（我独立复跑 `RUN_OPEN_FE=1` → 12 passed，默认 vs 缓存轴向 rel 1e-6；主管线 19/1）。
> 20° 倾角 demo：距下 |F| 25.38→21.49 kN（新增横向 5.54 kN）、|M| 157→274 N·m。⚠️ **纯轴向 1D 只吃纵向分量 → 倾角反而降低 1D 风险**（横向/剪切是 1D 盲区）→ 真正的 3D 影响需 **plan-S4（wrench 上 FE）**。
> **未做**：接触本构/摩擦/CoP/单脚（plan-S2/S4）。**下一步 = plan-S2 单脚 or S3 姿势 or S4 wrench→FE（顺带治跖面奇异 A7）。**
>
> **非垂直扩展 plan-S4（3D wrench 上 FE + 跖面弹性支撑）✅ 完成并验真（bg_313a47bb）** —— `NONVERTICAL_S4_REPORT.md`：新增 opt-in `subtalar_force`/`subtalar_moment`（3D wrench 送 FE）+ `plantar_bc="spring"`。
> **轴向回归 PASS**（我复跑：主管线 19/1、非垂直 S1+S4 **22/1**；默认 gauge_max 逐位复现 186.756369）。
> 三个结论：① **力矩的"压力梯度"路径在本曲面关节面病态**（需 ±193 MPa 带压含拉伸 → 不物理），**必须走幽灵刚体 `RigidMomentLoad`**（且需放开转动 DOF）；② 跖面弹性支撑**缓解轴向 A7 边缘奇异（−18.7%）**，但 **3D wrench 把热点移到关节面棱边 → 不缓解（−7.4%）**；③ 20° 倾斜 gauge ×0.87（方向效应），叠力矩 ×0.68。
> ⚠️ **3D 载荷引入新失效模式（关节面棱边，剪切/力矩驱动）** —— 现有跖面 BC 治不了。**下一步候选 = plan-S2 单脚 / S3 姿势 / S5 FE 接触。**
>
> **非垂直扩展 plan-S2（单脚/不对称）✅ 完成并验真（bg_224c3501）** —— `NONVERTICAL_S2_REPORT.md`：`split` 泛化为任意逐足占比（含单脚）+ `foot_forces` 独立时程；**默认对称逐位复现（0.00e+00）**。
> **单脚 vs 对称 @5m**：承力侧关节反力约 **×2**（subtalar ×2.01 / ankle ×2.02 / hip ×2.45；lumbar ×0.96）；**排序不变**，但**首折大幅前移**（腓 7→1m、跟/胫 8→2m、股 19→4m）；**@5m 越阈骨 0→4**（腓/跟/胫/股）。因子随高度近似恒定。
> **验真：主 19/1、非垂直 S1+S2 23/2。** **下一步候选 = plan-S3 姿势 / S5 FE 接触 / 汇总场景结论。**
>
> **场景结论汇总 ✅（bg_224c3501 之后）** —— `docs/场景探索结论.md`：抱石危险度**单脚 ≫ 3D斜向/力矩 > 直立双足 > 软垫**；单脚载荷 ×2→首折减半（@5m 越阈 0→4）；软垫首折 ×~1.4（有上限）；3D 载荷把热点跖面→关节面（新失效模式，BC 治不了）。
> **临床论文启示 + R1/R2 ✅（bg_5bc28838 / bg_2bd8feb3）** —— 参照 `paper/Boulder_Dissertation_EN.md`(Heck2024/Müller2022, 430 例) + `paper/MinerU_markdown_fspor-7-1609133_*.md`(Beurienne2025, 245 人)，落 `docs/两篇临床论文启示与建议.md`（R1–R5）。
> **R1 踝旋后**（`NONVERTICAL_S5_SUPINATION_REPORT.md`）：opt-in `ground_reaction(roll_deg=)` + `ankle_supination` 模块；默认回归 0；**纯轴向腓骨 7m → 加旋后(15–30°) → 真实高度 2–4.5m 即越阈**。
> **R2 场景矩阵**（`SCENARIO_FIELD_MATRIX_REPORT.md`）：7 运动学→现有旋钮；下肢旗标 100%（vs 临床 61–67%，过度）；踝桶 2/4 场景（腓骨最先）；R5 用姿态分布投影与临床同向。
> **收敛**：**腓骨/外踝 = 现场机制(旋后) + 临床踝伤#1 的交点** → 解释我们此前的"腓骨敏感"。**验真：主 19/1、非垂直 S1+S2+S5 43/3、场景 15。**
> **R4 垫子刚度×旋后权衡 ✅（bg_af224c19）** —— `PAD_SUPINATION_TRADEOFF_REPORT.md`：耦合 S1 R(k)+R1 旋后；新模块 `src/climbing/coupling/pad_supination.py`（opt-in）；**无内部最优 → k*→刚性极限**（L_roll 三档稳健）；**软垫轴向收益(−23.6%) 远小于旋后侧弯(O(10)×)** → **两轴张力**。**验真：主 19/1、非垂直 S1+S2+S5+R4 75/3。**
> **R6 踝韧带扭伤判据 ✅（bg_cd75d981，参考 Tochigi2006 + librarian 失效数据）** —— `docs/踝韧带判据方案.md` + `ANKLE_LIGAMENT_CRITERION_REPORT.md`：新模块 `src/climbing/coupling/ankle_ligament.py`（opt-in）；关节级内翻力矩判据（保守 21 N·m / 高承载 77 N·m 两口径）+ 力/应变（上界）；**β≥10° 两口径下 ATFL 均扭伤**、β=0 安全；**补上临床 #1 踝扭伤 71%**（骨口径仅骨折 27%）。**验真：主 19/1、非垂直 S1+S2+S5+R4+韧带 124/3。**
> **plan-S3 姿势 ✅（bg_3e37bf23）** —— `NONVERTICAL_S3_POSTURE_REPORT.md`：新 opt-in `opensim_fall.LegPose`(锁姿势) + `PassiveStiffness`(解锁+被动刚度)；默认回归 0。**(a) 锁姿势=几何重分布**（远端降/近端随高度升，非吸能）；**(b) 被动刚度不足**——冲击窗口 ~20 ms ≪ 自然周期 → 关节只屈 ~2°、峰值力无净衰减 → **真吸能需主动肌肉**（当前不建模）。**验真：主 19/1、非垂直 S1+S2+S3 35/6。**
> **非垂直分期进度：S1✅ S2✅ S3✅ S4✅**（S5 FE 接触待定）。
> **S5 计划 + 汇总报告 ✅** —— `docs/S5接触方案.md`（FE 脚-垫接触 + 摩擦计划；最小版 = 刚性面 + 摩擦，修 A7 + 摩擦切向力）；`docs/汇总报告_骨与韧带损伤建模.md`（整合 论文对齐 + 临床锚定 + 非垂直 S1–S5 + 场景判据 R1–R6 的总报告）。
> **论文初稿 ✅** —— `docs/论文初稿_抱石落地骨与韧带损伤建模.md`（中文 · 方法+应用论文；摘要/引言/方法/结果/讨论/结论 + 图表清单 + 参考文献节选）。

（以下为历史按序条目，多数已完成）

- ✅ **S3.2 `load_transfer.py`** —— 已实现（关节力 + 跟腱力 → `FeLoadSpec`）。
- ✅ **S3.2b 跟腱力** —— 已实现并检验：**非主因**（生理范围不改变 gauge_max）。
- 🔴 **S3.2c（新，最高优先）—— 距下关节载荷/BC 的边缘奇异** —— 当前 BC 是
  关节面**均匀压力** + 跖面**全固定**，峰值/gauge 恒在关节面棱边、随网格发散
  （voxel 26.3 vs gmsh 48.7 MPa）。这是"偏晚"的**真正嫌疑**。做法：
  (a) 用**分布关节面载荷/接触**替代均匀压力；(b) 判据改用**扩散部位（非棱边）**；
  (c) 校验 gauge 对 BC 的敏感性。**先做这个再继续扫描。**
- **S3.1 批量网格化 9 部位**：复用 `meshing`，固定 `reduction`；**单网格 ≤~5 万 tet**
  （优先 `vtp_to_tet` 体素路径，gmsh 路径会超限）。
- **补下肢被动刚度/肌肉**（S1 锁腿是权宜）——届时跟腱力应由动力学决定，复检 S3.2b。
- **S4** Logistic/聚类对标；**S5** 换软垫。
- **对标口径**：只用**相对模式**（部位排序 / 聚类边界 / OR / 双路径）；
  绝对骨折高度标"不可移植"（§3.5 D5、§六 A5）。
- ✅ **S3.C 方向 C（借 THUMS AM50 跟骨几何+材料）—— 已封档，未复现绝对高度**：
  抽取 + 解剖配准 + 双材料 FE + 38 高度扫描完成（**38/38 rc=0**）；**a=0 首次骨折 ≈2.0-3 m**，
  [Y25] 7-9 m **未复现**（偏高 ~3-4×）。残余差距在 **BC 层面**（跖面全固定奇异 / 150 MPa
  判据 / 跟腱 TractionLoad 边缘 / 多体 F_subt 量级），超出 C 的"不修动力学/BC"边界。
- 🔴 **若要追平 [Y25] 绝对高度 → 必须改 BC**（跖面非全固定 / 接触或分布关节面载荷 /
  重标判据）——**超出方向 C 边界**，需另立方向后再做。
- 🔭 **非垂直落地扩展 → `docs/非垂直落地扩展方案.md`**：把两端 1D 接口升到 3D、
  透传完整 wrench（中段 `subtalar_reaction` 已是 3D）。**S1（斜向矢量）✅ → S2（单脚）✅ →
  S3（姿势）✅ → S4（力矩上 FE）✅ 均已实现并回归；S5（FE 脚-垫接触）计划（`docs/S5接触方案.md`）。** 详见 `results/opensim_fe/NONVERTICAL_S*.md`。
- 🧪 **跖面 BC 敏感性实验（已完成）→ `results/opensim_fe/PLANTAR_BC_REPORT.md`**：把跖面从
  **三向全固定**换成 **三向 Winkler 弹簧**（k≈1e3 N/mm）后，gauge **186.8 → 142.7 MPa（降 24%）**、
  峰值**离开跖面棱边**，首次骨折外推 **2.11 → 5.43 m**（**向 [Y25] 7–9 m 靠近但未到达**）。
  **机制修正**：虚高来自**法向支承刚度**，不是"三向/切向"（真放开切向反而 +14%）。
  默认 BC 仍 `fixed`；详见 `docs/对齐_Y25_假设与结论.md` §2.1。
- 📉 **③④ 次级因素已排除（已完成）**：距下**载荷表示**（均匀/梯度/节点）±10% 内且方向敏感
  （`JOINT_LOAD_REPORT.md`）；**STATIC 充分**（DYNAMIC 瞬态仅 +1.1%，T8 通过）；多体 **77 kg 校准**
  F_subt +1.3% → 首折 2.02→1.99 m（`MASS_CALIB_REPORT.md`）；**跖腱膜**在固定跖面下为 0、
  弹性下 +0.28%@5 kN（`FASCIA_REPORT.md`）。**⇒ 主因仍是跖面法向刚度固定端奇异**。

---

## 8. 已知约束 / 红线（务必遵守）

- **HDF5 上限**：`pyfebio.xplt.to_hdf5` 把网格写成 **attribute**，**>~5 万 tet** 会
  `OSError: object header message is too large`，应力读不回 → **FE 网格必须 ≤~50 k tet**。
- **绝对高度不可移植**：7–9 m 依赖 THUMS 网格/材料/沙漏。
- **落地场景仅覆盖"竖直 + 对称双足 + 腿绷直"**：斜向 / 单脚 / 水平速度 / 力矩**均未建模**
  （见 §3「场景 / 姿态假设」、`docs/非垂直落地扩展方案.md`）——引用结果前先核对场景匹配。
- **跟骨是"结构失效"非"材料强度"**：材料强度×截面 高估 17–34×（`方案更新_v2` P0）→
  应改用 **Voo HIPC 接触力-骨折概率**（Barnes IRCOCB 2019 PMHS 验证）。
- **🔴 不要碰**（并发脚踝 FE 会话拥有）：`scripts/ankle_fe/`、`temp/pyfebio_demo/`、
  `temp/ankle_*`、`.learnings/*` 的 **LRN-010~019**（本项目占用 **020~027**）、
  `FE_PIPELINE_TODO.md`。
- 本项目自留目录：`src/climbing/coupling/`、`scripts/opensim_fe/`、`temp/opensim_fe/`、
  `results/opensim_fe/`、`docs/OpenSim_FE复现方案.md`、`docs/OpenSim_FE交接.md`。

---

## 9. 怎么跑（复制即用）

```powershell
$env:PYTHONPATH="src"
python scripts\opensim_fe\s1_vertical_slice.py --height 5          # 单高度全链路
python scripts\opensim_fe\s2_height_sweep.py                        # 1..50 m 扫描（默认 ℓ=4 mm）
python scripts\opensim_fe\s2_height_sweep.py --heights 5 --char-len 3
python scripts\opensim_fe\s2_plot.py                                # σ_vm(h) 图
python -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q # 回归

# 论文配图（§4b）
python scripts\opensim_fe\make_fracture_figures.py                 # 图 8 + 图 9 全 4 方位

# 坠落仿真动画（§4b）
python scripts\opensim_fe\make_fall_animation.py --height 2.0 --on-pad --fps 25
python scripts\opensim_fe\make_fall_animation.py --height 2.0 --on-pad --passive   # S3b 对比版

# 危险姿势演示（§4b）
python scripts\opensim_fe\make_pad_paradox_figure.py                               # 图 10
python scripts\opensim_fe\make_scenario_animations.py --scenario both              # 动画 A + B（≈93 s）
python scripts\opensim_fe\make_scenario_animations.py --scenario ser --two-foot     # S5 报告口径版
```

---

## 10. 参考

- [Y25] 袁红敏, 高树辉, 魏智彬. 中南大学学报(医学版). 2025;50(11):2051-2061.
- [T22] Li Z, et al. Forensic Sci Res. 2022;7(3):518-527.
- [B25] Beurienne E, et al. Front Sports Act Living. 2025;7:1609133.
- 仓库：`docs/OpenSim_FE复现方案.md`、`docs/OpenSim_FE交接.md`（本文）、`docs/非垂直落地扩展方案.md`（提案）、`docs/方案更新_v2.md`、`阶段总结.md`、`.learnings/`。
