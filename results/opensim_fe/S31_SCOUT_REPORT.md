# S3.1 第一步 · 九部位几何/网格化来源 + 2 骨试点抽取/网格化（Scout Report）

> 目标：把已跑通跟骨的管线（THUMS 抽取 → STL → tet 网格）推广到 [Y25] 的 9 部位，
> 先证明"能推广"，并暴露剩余工程的坑。
> 甲板：`model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`（THUMS AM50 V7.1，ASCII，mm）。
> 新增脚本（**未改任何现有脚本默认行为**）：
> `scripts/opensim_fe/scout_parts.py`（`*PART` 索引）、
> `scripts/opensim_fe/extract_part.py`（任意 PID → STL → tet）。
> 映射表：`results/opensim_fe/PARTS_MAPPING.md`。

---

## 1. 结论（TL;DR）

1. **抽取层：完全推广。** `kmesh_io` 对任意 `*PART` PID 都能抽 solid + 建外表面；
   9 部位全部在甲板里定位到（PID 见映射表），并实测到 solid 单元。
2. **网格化层：能推广，但只能走体素路径（`_delaunay_tet`）。**
   平滑 gmsh 路径（`_gmsh_tet_remesh`）在本环境 **`Mesh.MeshSizeMax` 不生效**（`char_len`
   从 5→12 mm，tet 数几乎不变 ~3.6×10⁵），**任何真骨都压不到 ≤5 万 tet**，无法用。
3. **试点两骨（+1 附赠）在 ≤5 万 tet 约束下全部网格化成功（rc=0）**，但质量是"可用/待改进"：
   - **胫骨（长骨）**：38209 tet，体积误差 **2.03%**，`J<0.3` **17.4%** —— 可接受。
   - **颅骨（取右顶骨 diploë）**：30714 tet，但源面是**薄壳**，体积无可靠参照（见 §3）——
     仅"证明能出网格"，几何保真度存疑。
   - 附赠 **腰椎 L3（左右 SPON 合并）**：44181 tet，但体素路径 **体积过填 +95%**。
4. **最大坑**：中轴骨/颅骨的**皮质是壳单元**（`_CORT`/`_external_`/`_internal_`/`_shell`），
   `kmesh_io` 只读 solid，会**丢皮质**；且薄壳骨用体素法会**欠采样/过填**。

---

## 2. 抽取来源（实测 PID/名称）

完整表见 `results/opensim_fe/PARTS_MAPPING.md`。要点（右侧）：

| 部位 | THUMS PID（右） | 备注 |
|------|-----------------|------|
| 股骨 | 81000000 / 81000001 / 81000100 | 3 段全 solid |
| 胫骨 | 81000600 / 81000601 / 81000700 | 试点；3 段全 solid |
| 腓骨 | 81000800 / 81000801 / 81000900 | 3 段全 solid |
| 足（跟骨） | 81001200 / 81001300 | 已有产物 |
| 骨盆 | 髋 83500200(SPON)+骶 83500100(SPON)；CORT=壳 | — |
| 腰椎 L3 | 89001500+89501500 (L+R SPON) | CORT=壳 |
| 胸椎 T6 | 89000600+89500600 (L+R SPON) | CORT=壳 |
| 颈椎 C5 | 87000500+87500500 (L+R SPON) | CORT=壳 |
| 颅骨 | 右 88000001–51 / 左 88000052–99；diploë solid，皮质板=壳 | 试点取 88000004 |

**左/右**：下肢 `81→82`，躯干 `890xxxxx`(L)↔`895xxxxx`(R)，单椎体需 L+R 两半合并。

**已有可复用几何**：
- `model/opensim/FullBodyModel-4.0/Geometry/`（**VTP，米**，低模）：`r_foot/r_talus/r_tibia/
  r_fibula/r_femur/r_pelvis/sacrum/l_*`；`hat_skull/hat_spine/hat_ribs_scap`（颅/脊柱仅低模）。
- `model/Model/Geometry/`（**STL/OBJ，mm**）：`calcaneus_r/l`、`talus_r/l`、`midfoot_r/l`、
  `metatarsal_*`、`toes_*`、`OK_femur.obj`、`OK_tibia.obj`（无骨盆/脊柱/颅独立骨）。
- `model/FootwithinLowerlimbModel/`：`FootModel_calcn/talus/midfoot/forefoot/digits_*`。

---

## 3. 试点数值（rc / 产物路径 / 可复现命令）

运行约定：仓库根、`$env:PYTHONPATH="src"`、`& .venv\Scripts\python.exe <script>`。
三个试点均 **rc=0**。

### 3.1 胫骨（长骨）——`tibia_r`

```powershell
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name tibia_r --pids 81000600,81000601,81000700 --char-len 4.0
```

| 量 | 值 |
|----|----|
| 抽取节点 / solid 单元 | 7245 / 23550（6651+14482+2417） |
| 外表面 | 4852 tri，边界边=0，**非流形边=17** |
| 散度定理包络体积（参照） | 346130 mm³ = **346.1 cm³** |
| 网格方法 | `delaunay-voxel-repair`（`gmsh-raw` 失败：overlapping facets；`iso-dec-gmsh` 得 409761 tet > 5万被弃） |
| 网格节点 / tet | 7892 / **38209**（**≤5万 ✓**） |
| 网格体积 / 误差 | 353162 mm³ / **+2.03%** |
| `min_jacobian` / `J<0.3` 占比 | 1.02e-3 / **17.41%** |
| 产物 | `temp/opensim_fe/parts/tibia_r/`（`tibia_r_surface.stl`、`tibia_r_volume.vtk`、`tibia_r_mesh.npz`、`tibia_r_elements.npz`、`summary.json`） |

### 3.2 颅骨（不规则/壳多）——`skull_r_parietal`（取右顶骨 diploë）

```powershell
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name skull_r_parietal --pids 88000004 --char-len 3.0
```

| 量 | 值 |
|----|----|
| 抽取节点 / solid 单元 | 2952 / 1430 |
| 外表面 | 6068 tri，**非流形边=174** |
| 包络体积（**不可用参照**） | 3466142 mm³ = 3466 cm³（=颅腔，非骨体积） |
| 网格方法 | `delaunay-voxel-repair`（`gmsh-raw` 访问越界崩溃；`iso-dec-gmsh` 539817 tet > 5万被弃） |
| 网格节点 / tet | 7187 / **30714**（≤5万 ✓） |
| 网格体积 | 119465 mm³（≈119 cm³，**偏大**：顶骨薄板 ~40 cm³，体素把板增厚） |
| `min_jacobian` / `J<0.3` | 1.26e-3 / **17.55%** |
| 产物 | `temp/opensim_fe/parts/skull_r_parietal/` |

> ⚠️ 颅骨是**薄壳穹顶**：`free_surface` 的"包络体积"= 整个颅腔（实测整 5 块颅顶合并 → 2020 cm³）。
> 因此 `volume_error_pct` 对颅骨**无意义**（96.5% 是参照错，不是过填 96.5%）。几何保真只能靠
> 肉眼/包络 bbox 与体素收敛核查。

### 3.3 附赠：腰椎 L3（左右 SPON 合并）——`lumbar_L3`

```powershell
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name lumbar_L3 --pids 89001500,89501500 --char-len 2.0
```

| 量 | 值 |
|----|----|
| 抽取节点 / solid 单元 | 2203 / 8792 |
| 外表面 | 2670 tri，**非流形边=0（干净！）** |
| 包络体积（参照） | 26241 mm³ = 26.2 cm³ |
| 网格方法 | `delaunay-voxel-repair`（`gmsh-raw` 成功但 157556 tet > 5万） |
| 网格节点 / tet | 8600 / **44181**（≤5万 ✓） |
| 网格体积 / 误差 | 51103 mm³ / **+94.7%（过填）** |
| `min_jacobian` / `J<0.3` | 1.96e-3 / **17.16%** |
| 产物 | `temp/opensim_fe/parts/lumbar_L3/` |

### 3.4 抽取可用性普查（`temp/opensim_fe/parts/_verify/summary.json`）

跨 6 部位 22 个代表 PID 一次扫描：**59830 solid 单元**，全部可抽；确认中轴骨/颅骨
**CORT = 0 solid（壳）**。

---

## 4. 剩余工作 · 坑清单（按优先级）

1. **【阻塞级】gmsh 尺寸控制失效**：`_gmsh_tet_remesh` 的 `Mesh.MeshSizeMax` 不生效，
   所有真骨都超 5 万 tet。二选一：(a) 自写 gmsh 重网格（用 `MeshSize` 场 / `CharacteristicLengthMax`
   + `MeshSizeExtendFromBoundary=0`）验证；(b) 全面采用体素路径（现状）。**S3.1 批量前须定案。**
2. **【阻塞级】皮质壳缺失**：中轴骨/颅骨皮质是 `*ELEMENT_SHELL`，`kmesh_io` 不解析 →
   模型只有松质。需要：新增壳解析（`*ELEMENT_SHELL` + `*SECTION_SHELL`）或"壳→实体加厚"。
3. **体素质量/过填**：`J<0.3` 恒 ~17%；瘦小骨（L3 +95%）过填，薄板（顶骨）增厚。
   需**逐骨网格收敛研究**（`char_len` 扫描 + 体积/应力稳定性），不能直接用当前值。
4. **参照体积对壳几何不可用**：`free_surface` 包络体积对薄壳=内腔。判据应改为
   **网格体积收敛** + 骨表面积/内腔一致性，而非 `volume_error_pct`。
5. **整段脊柱合并**：整条腰/胸/颈 = 多椎体 + 椎间盘 + 韧带，多元体非流形（颅顶 5 块实测
   non-manifold=3959）；第一步**只做单代表椎体**（L3/T6/C5），整段另立子任务。
6. **骨盆**：髋骨 `CORT` 有 9 个厚度变体，**只取 1 个**；髋+骶+耻骨联合的合并面需专门处理。
7. **左右/合并规则**：单椎体必须 L+R 两半；选一侧（右）为基准，镜像留给配准步。
8. **单位/入口**：`meshing.stl_to_tet_gmsh` 假定输入**米**（×1000），THUMS STL 是 mm；
   本试点新脚本绕开该入口（用坐标无关的内部函数）。
9. **FE builder 泛化（plan ④）**：`_finalize`/`_write_msh` 的具名面（距下/跖面/跟腱）是
   跟骨专用；其它骨这些集合为空，载荷/BC 需按骨重定义。
10. **HDF5 上限**：单网格 ≤5 万 tet（`pyfebio.xplt.to_hdf5` 用 HDF5 attribute 存网格，
    LRN-027 坑 4）；`extract_part.py --max-tets` 已内建强制。

---

## 5. 诚实结论

- **能推广的部分**：从甲板抽取任意骨（9 部位 PID 全部实测到位）+ 体素四面体化，
  在 ≤5 万 tet 约束下**已对胫骨、颅骨（顶骨）、腰椎 L3 三个新部位跑通（rc=0）**，
  证明跟骨管线不是单骨特例。
- **不能声称的部分**：
  - 平滑 gmsh 网格**当前不可用**（尺寸控制失效，恒 >5 万 tet），故新骨网格质量
    （`J<0.3`≈17%）**低于**跟骨的 gmsh 路径（6.2%）；
  - 中轴骨/颅骨**缺皮质壳**；颅骨试点仅是"能出网格"，几何保真**未被验证**；
  - 体积误差对薄壳不可用，L3 有 +95% 过填。
- **因此**：S3.1"九骨批量网格"可在**单代表骨 + 体素路径 + ≤5 万 tet**下推进并出数量级正确的
  模型，但**在解决坑 #1/#2/#3 之前，不应把新骨网格的绝对应力当作 [Y25] 对标结论** ——
  与项目既有结论（`LEARNINGS.md` LRN-027：绝对高度不可移植，只看相对模式）一致。
