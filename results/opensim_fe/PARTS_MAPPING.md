# [Y25] 九部位 → THUMS 几何/PID 映射表（S3.1 第一步）

> 状态：**实测**。所有 PID / 名称 / 单元数均来自对 THUMS AM50 V7.1 甲板的实际扫描
> （`scripts/opensim_fe/scout_parts.py` + `extract_part.py --report-only`），未臆造。
> 甲板：`model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`
> → `THUMS_model/THUMS_AM50_V71_Occupant_202411.k`（230 MB，ASCII）。
> 实测：1834 个 `*PART`（带 `$HMNAME COMPS` 名称）。
> 机器可读索引：`temp/opensim_fe/parts/parts_index.json`、`temp/opensim_fe/parts/_all_names.txt`。

**9 部位锁定**（与 [Y25]、`docs/OpenSim_FE交接.md` §7 一致）：
足部 / 胫骨 / 腓骨 / 股骨 / 骨盆 / 腰椎 / 胸椎 / 颈椎 / 颅骨。

单位与坐标系：THUMS 全局系，**mm**（`kmesh_io` 不做重定心/缩放）。

---

## 0. 关键实测事实（决定下面所有策略）

1. **THUMS 把"骨"拆成多个 `*PART`**：
   - **长骨**（股/胫/腓）：`end_spon` + `center_spon`（松质，多为退化 hex＝tet）+ `*_CORT`
     （皮质，**实体 hex**）→ **都是 solid**，可合并成一个闭合骨面。
   - **中轴骨 / 颅骨**（骶/髋/椎/颅）：**只有 `SPON` 是 solid**；`*_CORT`、`*_external_*`、
     `*_internal_*`、`*_shell*` 都是 **壳单元（`*ELEMENT_SHELL`，0 个 solid）**，`kmesh_io`
     只读 `*ELEMENT_SOLID`，**会直接丢掉皮质壳**。实测（`_verify`）：`R_SACRUM_CORT`、
     `R_HIPBONE_CORT_1.5`、`L_L3_CORT`、`L_T6_CORT`、`L_C5_CORT` 的 solid 数 **均为 0**。
2. **左右**：下肢/足 `81xxxxxx`=右，`82xxxxxx`=左；躯干 `89xxxxxx`=左，`895xxxxx`=右；
   **单个椎体本身被拆成 `L_`/`R_` 两半**（如 `L_L3_SPON` 89001500 / `R_L3_SPON` 89501500），
   要拼成一整块必须两半都取。
3. **髋骨 `CORT` 有 9 个厚度变体**（`_1.5/_ACETABLUM1.5/_ACETABLUM1.0/_0.75/_2/_3/_2.2/_1.8/_1`）
   ——这是参数化研究的**同一皮质的多个厚度壳**，**绝不能全部合并**（会重叠）；只取 `SPON`
   + 需要时**挑一个**壳厚度。
4. **`meshing` 实测**（本机 gmsh SDK 4.15.2）：
   - `stl_to_tet_gmsh` / `_gmsh_tet_remesh` 的 `Mesh.MeshSizeMax` **在本环境不生效**——
     tibia 从 `char_len` 5→12 mm，tet 数几乎不变（~3.6×10⁵），即 gmsh 路径**无法压到 ≤5 万 tet**。
   - `stl_to_tet_gmsh` 入口会 `×1000`（假定输入是**米**）；THUMS STL 已是 **mm**，直接调用会放大 1000×。
     本次试点新脚本 `extract_part.py` **不复用**该入口，直接调用坐标无关的
     `_gmsh_tet_remesh / _delaunay_tet / _finalize`，全程 mm。
   - 因此按项目既有约定（`docs/OpenSim_FE交接.md` §7、`s3_achilles_sweep.py`）：
     **优先 `vtp_to_tet`/`_delaunay_tet` 体素路径**；gmsh 仅当能 ≤5 万 tet 时才用。

---

## 1. 九部位映射总表

| # | 部位 | 来源（首选 THUMS PID，右/左） | 名称 | 已有独立几何（可直接用） | 合并 / 取代表策略 |
|---|------|------------------------------|------|--------------------------|-------------------|
| 1 | **足部** | 见 §2（跟骨已做） | talus/calcaneus/navicular/cuboid/cuneiform×3/metatarsal×5/phalanx×14 | `model/Model/Geometry/{calcaneus_r, talus_r, midfoot_r, metatarsal_*, toes_*}.stl`；`model/FootwithinLowerlimbModel/FootModel_*.stl` | 跟骨**已单骨**（`temp/opensim_fe/thums_calcaneus/`）；整足=按骨逐个或合并 |
| 2 | **胫骨** | 81000600 / 81000601 / 81000700 ‖ 82000600/01/00 | `right_tibia_end_spon` / `right_tibia_center_spon` / `R_TIBIA_CORT` | `model/opensim/FullBodyModel-4.0/Geometry/r_tibia.vtp`；`model/Model/Geometry/OK_tibia.obj` | **合并 3 PID** 成一根骨（已试点） |
| 3 | **腓骨** | 81000800 / 81000801 / 81000900 ‖ 82000800/01/09 | `right_fibula_end_spon` / `…center_spon` / `R_FIBULA_CORT` | `…/r_fibula.vtp` | 合并 3 PID |
| 4 | **股骨** | 81000000 / 81000001 / 81000100 ‖ 82000000/01/00 | `right_femur_end_spon` / `…center_spon` / `R_FEMUR_CORT` | `…/r_femur.vtp`；`OK_femur.obj` | 合并 3 PID |
| 5 | **骨盆** | 髋骨 83500200（SPON）+ 骶 83500100（SPON）；左 83000200/83000100 | `R_HIPBONE_SPON_15` / `R_SACRUM_SPON` | `…/r_pelvis.vtp`、`…/sacrum.vtp` | **合并髋骨+骶（+耻骨联合 83500300/83000300）**；皮质壳另处理（§3） |
| 6 | **腰椎** | L3 示例：89001500+89501500（SPON）；L1..L5 = 89001300–89001701 ↔ 89501300–89501701 | `L_L3_SPON`/`R_L3_SPON`（CORT 为壳） | 无（仅 `hat_spine.vtp` 低模） | **取单个代表椎体**（建议 L3 或 L1）；整段=L1–L5 全并（含椎间盘） |
| 7 | **胸椎** | T6 示例：89000600+89500600（SPON）；T1..T12 = 89000100–89001201 ↔ 89500100–89501201 | `L_T6_SPON`/`R_T6_SPON` | 无（`hat_spine.vtp`/`hat_ribs_scap.vtp` 低模） | **取单个代表胸椎**（建议 T6/T8）；整段=T1–T12 全并 |
| 8 | **颈椎** | C5 示例：87000500+87500500（SPON）；C1..C7 = 87000100–87000701 ↔ 87500100–87500701 | `L_C5_SPON`/`R_C5_SPON` | 无 | **取单个代表椎体**（建议 C5）；整段=C1–C7 全并 |
| 9 | **颅骨** | 右 88000001–88000051 / 左 88000052–88000099（见 §4） | `frontal_r`/`parietal_r`/`temporal_r`/`occipital_r`/… | 无（`hat_skull.vtp` 极低模） | **取单块代表骨（diploë）**；整颅=多骨合并且**必须处理皮质壳** |

> 左/右镜像 PID 规则：下肢/足仅换前缀 `81→82`；躯干骨 `89 0xxxxx`(L) 与 `895 xxxxx`(R) 成对；
> 颅骨左块 ≈ 右块 PID + 51（如 `frontal_r`=88000001 ↔ `frontal_l`=88000052，实测）。

---

## 2. 足部（部位 1）细目 —— 右足 PID（左足全部 `82xxxxxx` 对应）

| 骨 | SPON PID | CORT PID | solid? |
|----|----------|----------|--------|
| Talus 距骨 | 81001000 | 81001100 | 是 |
| Calcaneus 跟骨 | 81001200 | 81001300 | 是（**已抽取**） |
| Navicular 舟骨 | 81001400 | 81001401 | 是 |
| Cuboid 骰骨 | 81001500 | 81001501 | 是 |
| Medial cuneiform | 81001600 | 81001601 | 是 |
| Intermediate cuneiform | 81001700 | 81001701 | 是 |
| Lateral cuneiform | 81001800 | 81001801 | 是 |
| Metatarsal 1–5 | 81001900–81002300 | 81001901–81002301 | 是 |
| Phalanges | 81002400–81003700 | 81002401–81003701 | 是 |
| 足软组 | `R_FOOT_SKIN` 81201301 / `R_FOOT_FLESH_TISSUES` 81201400 | — | 皮肤/软组织 |

**已有 STL 可用**（`model/Model/Geometry/`，mm）：`calcaneus_r.stl`、`talus_r.stl`、
`midfoot_r.stl`、`metatarsal_*_r2.stl`、`toes_*_r3.stl`；`FootwithinLowerlimbModel/`
另有 `FootModel_calcn_r/talus_r/midfoot_r/forefoot_r/digits_r.stl`。
→ 足部**不需要**从 THUMS 重新抽所有骨，可优先复用这些 STL；与 THUMS 同源的做法是
逐骨按其 PID 抽取（跟骨已示范）。

---

## 3. 骨盆 / 脊柱（部位 5–8）细目与"皮质壳"问题

实测 solid 单元数（`_verify` 扫描）：

| PID | 名称 | solid 单元 |
|-----|------|-----------|
| 83500200 | R_HIPBONE_SPON_15 | 11767 |
| 83500201 | R_HIPBONE_CORT_1.5 | **0（壳）** |
| 83500100 | R_SACRUM_SPON | 4696 |
| 83500101 | R_SACRUM_CORT | **0（壳）** |
| 89001500 / 89501500 | L_L3_SPON / R_L3_SPON | 4396 / 4396 |
| 89001501 / 89501501 | L_L3_CORT / R_L3_CORT | **0（壳）**（上下同） |
| 89000600 / 89500600 | L_T6_SPON / R_T6_SPON | 1593 / 1593 |
| 87000500 / 87500500 | L_C5_SPON / R_C5_SPON | 830 / 830 |

**含义**：从 THUMS 只能拿到**松质实体**。若模型需要皮质骨，必须另走壳单元路径
（`*ELEMENT_SHELL`，`kmesh_io` 目前不解析）或"壳→加厚实体"的转换——列为剩余工作（§报告）。

**椎体取代表 vs 合并该段**：
- **取单代表椎体**（推荐第一步）：`L_k` + `R_k` 两半的 `SPON` 合并即可得一块完整椎体，
  面干净（L3 实测：surface 2670 tri，**non-manifold=0**）。
- **合并该段**（如整条腰椎）：需并联 L1–L5 共 10 个 `SPON` + 对应 10 个 `CORT`(壳)
  + 9 个椎间盘（`L_ANNULUS_OUT(t–t+1)`/`NUCLEUS_PULPOSUS`）。多体、非流形度极高，
  网格需逐椎体或体素融合，**不建议第一步做整段**。

---

## 4. 颅骨（部位 9）细目 —— 复合 + 壳

颅骨在 THUMS 里**不是一块**，而是 ~50 块颅面骨的左右镜像集合（88000001–88000099）。
每块又分**三层**：`*_r`（diploë，实体）/ `*_external_r`、`*_internal_r`（内、外皮质板，**壳**）。

| 右骨 | diploë PID（solid 实测） | 壳 PID |
|------|--------------------------|--------|
| frontal 额 | 88000001（789） | 88000002/03（**0**） |
| parietal 顶 | 88000004（1430） | 88000005/06（**0**） |
| temporal 颞 | 88000007（1598） | 88000008（**0**） |
| occipital 枕 | 88000011（1161） | 88000012（**0**） |
| sphenoid 蝶 | 88000014（1273） | — |
| ethmoid 筛 | 88000017（229） | — |
| mandible 下颌 | 88000033（145） | 88000034/35（壳） |
| 其余（maxilla/zygomatic/palatine/lacrimal/nasal/vomer/teeth…） | 88000021–88000051 段 | 各类 `_external/_internal/_shell` |

**整颅合并的实测陷阱**：把 5 块颅顶 `diploë`（frontal+parietal+temporal+occipital+sphenoid，
6251 solid）合并后，`free_surface_triangles` 给出 25190 tri、**non-manifold=3959**，
且散度定理"包络体积"= **2020 cm³ ≈ 整个颅腔**——因为颅骨是**薄壳穹顶**，其外边界"围住"脑腔，
包络体积不是骨体积。**整颅不能当一块实心骨来网格化**。

→ **试点取单块代表骨**（本次取 `parietal_r` 88000004）。整颅的正确做法（剩余工作）：
逐骨网格 + 皮质壳单独建壳/加厚，或先做"皮质壳 → 闭合实体"的几何缝合。

---

## 5. 网格策略（每部位）

| 部位 | 提取（source） | 主路径 | 兜底 | target reduction | 备注 |
|------|----------------|--------|------|------------------|------|
| 足部（逐骨） | THUMS PID 或既有 STL | 既有 STL→`extract_part.py` | 体素 `_delaunay_tet` | 0.7（iso 路径用） | 跟骨已完成 |
| 胫/腓/股（长骨） | 3 PID 合并 | 体素 `_delaunay_tet`（gmsh 超限） | gmsh（若 ≤5万） | 0.7 | 试点：胫骨 |
| 骨盆 | 髋 SPON+骶 SPON(+耻骨联合) | 逐块体素 | gmsh | 0.7 | 皮质壳另处理 |
| 腰/胸/颈 | 单代表椎体 = L+R SPON | 体素 | gmsh | 0.7 | 皮/壳另处理 |
| 颅骨 | 单块 diploë（如 parietal_r） | 体素 | — | 0.7 | 壳另处理 |

**统一约束：单网格 ≤ 5 万 tet**（`pyfebio.xplt.to_hdf5` 用 HDF5 attribute 存网格，
>~5万 tet 报 `OSError: object header message is too large`，见 `LEARNINGS.md` LRN-027 坑 4）。
`extract_part.py` 已内建 `--max-tets`（默认 50000），超过就换路径/提示加大 `--char-len`。

---

## 6. 复现命令

```powershell
# 仓库根、PYTHONPATH=src
$env:PYTHONPATH="src"

# 1) 甲板 *PART 全量索引（PID + 名称）
& .venv\Scripts\python.exe scripts\opensim_fe\scout_parts.py `
    --out temp\opensim_fe\parts\parts_index.json

# 2) 任意骨抽取（可分部位落盘）
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name tibia_r --pids 81000600,81000601,81000700 --char-len 4.0
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name skull_r_parietal --pids 88000004 --char-len 3.0
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name lumbar_L3 --pids 89001500,89501500 --char-len 2.0

# 3) 只看抽取（不网格化）
& .venv\Scripts\python.exe scripts\opensim_fe\extract_part.py `
    --name apex --pids 81000000,81000001,81000100 --report-only
```

产物：`temp/opensim_fe/parts/<name>/{*_surface.stl, *_volume.vtk, *_elements.npz,
node_ids.npy, nodes_xyz.npy, surface_nodes.npy, *_mesh.npz, summary.json}`。
