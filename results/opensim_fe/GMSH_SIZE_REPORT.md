# gmsh 尺寸控制失效 · 根因与修复（bone tet ≤5 万 + 质量）

> 范围：让 gmsh 路径把真骨压到 **≤5 万 tet**（HDF5 attribute 上限，LRN-027 坑 4）
> 且 `J<0.3` 接近跟骨（目标 ≲8%），替代只能用的体素路径（`J<0.3`≈17%）。
> 新增脚本（**未改** `src/climbing/coupling/meshing.py` 任何默认行为/签名）：
> `scripts/opensim_fe/gmsh_size_probe.py`（诊断/对照扫描）、
> `scripts/opensim_fe/part_meshing.py`（mm-native 修复 helper）。
> 全部数字为本机实测（gmsh SDK 4.15.2，Py 3.12 venv）。

---

## 0. 结论（TL;DR）

1. **根因确定**：`meshing._gmsh_tet_remesh` 调 `gmsh.model.mesh.classifySurfaces(angle)`
   用的是**默认** `forReparametrization=False`。这会让 `createGeometry()` 把输入 STL 的
   顶点/面片当作**必须保留的几何约束**，`generate(2)` 无法真正重网格化表面 ——
   于是 `Mesh.MeshSizeMax` 只封顶、永远够不着，tet 数由**输入面片密度**钉死。
   实测胫骨 `char_len` 3→10 mm，tet 几乎不变（470605→362545）。
2. **修复 = gmsh 官方 `remesh_stl` 的做法**：`classifySurfaces(angle, boundary=True,
   forReparametrization=True)` + 关掉点/边界/曲率尺寸传播（gmsh t10 recipe）。
   修复后 tet 数**随 char_len 单调下降**：胫骨 63115→2547（cl 3→10）。
3. **推荐配置**：`reparam`（forReparam=True + `MeshSizeFromPoints/FromCurvature/
   ExtendFromBoundary=0` + `MeshSizeMin=Max=char_len`）。胫骨 **cl=4 → 27670 tet，
   `J<0.3`=3.66%，体积误差 1.11%，2.5 s**；体素路径同体量下是 17.4%。
4. **不需要自写重网格**：gmsh 内建的 `forReparametrization` 已解决；helper 只是
   把这条路径以 mm-native 方式暴露出来，不动 `meshing.py`。
5. **诚实边界**：非流形/自交的原始骨面（胫骨 17 非流形边、顶骨 174）必须先经
   marching-cubes 修面（`iso-dec`），这会改变几何；薄壳骨（顶骨）仍被"增厚到 ~118 cm³"，
   这是**几何修复**问题，不是尺寸控制问题。皮质壳仍缺失（未纳入本任务）。

---

## 1. 被诊断函数当前设置的全部 gmsh 选项

`meshing._gmsh_tet_remesh(surf, char_len, angle=0.001)`（`meshing.py` L655–722）：

| 选项 | 值 | 说明 |
|------|----|----|
| `General.Terminal` | 0 | 静音 |
| `General.Verbosity` | 0 | 静音 |
| `Geometry.Tolerance` | 1e-6 | 合并容差 |
| `Mesh.Algorithm` | 6 | 表面 Frontal-Delaunay |
| `Mesh.MeshSizeMin` | `char_len*0.4` | 全局尺寸下限 |
| `Mesh.MeshSizeMax` | `char_len` | 全局尺寸上限 |
| `Mesh.Algorithm3D` | 1 | 体 Delaunay |

调用序列：`merge` → `removeDuplicateNodes` → `classifySurfaces(ang)`（**默认**
`boundary=True, forReparametrization=False`）→ `createGeometry` → `synchronize`
→ `addSurfaceLoop/addVolume` → `generate(2)` → `generate(3)`。

**未设置、即用默认值**（关键）：

| 选项 | 默认 | 语义 |
|------|------|------|
| `Mesh.MeshSizeFromPoints` | **1** | 用几何点上的尺寸控制网格 |
| `Mesh.MeshSizeExtendFromBoundary` | **1** | 把边界尺寸向内部延伸 |
| `Mesh.MeshSizeFromCurvature` | 0 | 曲率自适应（已关）|

（本机实测：`Mesh.MeshSizeFromBoundary` **不存在**，`getNumber/setNumber` 均报
`Could not get/set option`；`Mesh.MeshSizeMax` 与 `Mesh.CharacteristicLengthMax`
是同一选项的别名。）

## 2. 根因假设

> **`forReparametrization=False` 让输入 STL 的离散面片/顶点变成强制几何实体，
> 表面无法被粗化；配合 `MeshSizeFromPoints=1` + `MeshSizeExtendFromBoundary=1`，
> 输入点距成为实际尺寸驱动，`MeshSizeMax` 只封顶、永远不生效。**

gmsh 参考手册 t10 明确：网格尺寸优先级里，若 `MeshSizeFromPoints` 置位，则
几何点上的尺寸优先于 `MeshSizeMax`；随后才被约束进 `[Min,Max]` 并乘 `Factor`。
`classifySurfaces(angle)`（`forReparametrization=False`）创建的离散几何把每个输入
顶点保留为角点/强制节点，所以即使 `MeshSizeExtendFromBoundary=0`，表面也被输入点距钉住。

**证据链**（胫骨 iso-dec 面，10632 tri，平均边长 3.853 mm）：

| 配置 | cl=3 | cl=10 | 判读 |
|------|------|-------|------|
| `baseline`（现状） | 470605 | 362545 | **几乎不变** → 尺寸控制失效 |
| `no-boundary`（仅关点/边界源） | 156553 | 102974 | 降下去但**平台化**，表面仍被钉住 |
| `reparam`（forReparam=True） | 63115 | 2547 | **单调下降** → 真正生效 |

只有打开 `forReparametrization=True` 后，`generate(2)` 才能把表面重网格到
`char_len`，tet 才随 `char_len` 下降。

---

## 3. 探针方法

脚本：`scripts/opensim_fe/gmsh_size_probe.py`。它复刻 `_gmsh_tet_remesh` 的管线，
但把 gmsh 选项/分类策略参数化，对 **每个配置 × `char_len∈{3,4,5,6,8,10} mm`** 记录
`tet 数 / 体积误差 / J<0.3 / 用时`，失败/报错**原文**记录。两套输入面：

* `raw`：清理补孔后的原始 STL。
* `iso-dec`：marching-cubes 等值面 → decimate 0.70 → clean（修复自交/非流形，
  与 `extract_part.mesh_stl` 路径 2 一致）。

参考体积取 `temp/opensim_fe/parts/tibia_r/summary.json` 的
`surface.enclosed_volume_mm3 = 346130.49 mm³`（散度定理，可信参照）。

## 4. 配置对比表（胫骨，iso-dec 面，实测）

**tet 数**（`FAIL` = 该面该配置 gmsh 报错；raw 面全部 FAIL，见 §4.3）：

| 配置 | cl=3 | cl=4 | cl=5 | cl=6 | cl=8 | cl=10 |
|------|------|------|------|------|------|-------|
| `baseline`（现状） | 470605 | 409761 | 389244 | 373925 | 362322 | 362545 |
| `reparam-baseline` | 69923 | 30491 | 16982 | 10839 | 6869 | 4189 |
| `no-boundary` | 156553 | 125122 | 113026 | 108122 | 104061 | 102974 |
| `field-constant` | 156553 | 125122 | 113026 | 108122 | 104061 | 102974 |
| `field-matheval` | 156553 | 125122 | 113026 | 108122 | 104061 | 102974 |
| `factor` | 81849 | 82468 | 82321 | 82482 | 82559 | 82477 |
| **`reparam`（推荐）** | **63115** | **27670** | **15439** | **8852** | **4230** | **2547** |

**`J<0.3` 占比（%）**：

| 配置 | cl=3 | cl=4 | cl=5 | cl=6 | cl=8 | cl=10 |
|------|------|------|------|------|------|-------|
| `baseline` | 4.23 | 4.56 | 4.78 | 4.84 | 4.96 | 5.08 |
| `reparam-baseline` | 3.46 | 3.85 | 4.00 | 4.58 | 5.88 | 6.21 |
| `no-boundary` | 7.63 | 9.39 | 10.51 | 11.16 | 11.89 | 12.33 |
| `field-constant` | 7.63 | 9.39 | 10.51 | 11.16 | 11.89 | 12.33 |
| `field-matheval` | 7.63 | 9.39 | 10.51 | 11.16 | 11.89 | 12.33 |
| `factor` | 51.00 | 52.08 | 51.75 | 52.05 | 51.92 | 52.08 |
| **`reparam`** | **3.21** | **3.66** | **4.48** | **4.26** | **5.98** | **8.72** |

**体积误差（%）**（`baseline/no-boundary/field/factor` 恒 2.13；`reparam` 系列见下）：

| 配置 | cl=3 | cl=4 | cl=5 | cl=6 | cl=8 | cl=10 |
|------|------|------|------|------|------|-------|
| `reparam-baseline` | 1.65 | 1.24 | 0.71 | 0.33 | 0.08 | 0.90 |
| `reparam` | 1.54 | 1.11 | 0.93 | 0.11 | 0.87 | 2.18 |

**用时（s，代表性）**：`baseline` 14–18；`no-boundary`≈6–7；`factor`≈10；
`reparam` 0.5–5.5（cl=4 为 2.0 s）。`reparam` 在更少 tet 下**更快**。

### 4.1 各配置判读

* **`baseline`**：尺寸控制失效（平台 36–47 万），不可用。
* **`no-boundary` / `field-constant` / `field-matheval`**：三者**逐位相同** —— 一旦
  `ExtendFromBoundary=0`，Constant/MathEval 背景场就是冗余的。tet 降到 ~10–15 万但
  **平台化**（表面仍被输入点钉住），且 `J<0.3` **恶化到 12%**。⇒ **不是解**。
* **`factor`**：把 `MeshSizeFactor=char_len/平均输入边长`。tet 平台在 ~8.2 万，
  且 `J<0.3`≈**51%**（严重畸变）⇒ 不可用。
* **`reparam-baseline`**：仅把分类改成 `forReparam=True`、尺寸选项仍走 baseline。
  已能压到 ≤5 万（cl≥4），`J<0.3`≈3.5–4%。可用，但比 `reparam` 略差。
* **`reparam`（推荐）**：forReparam=True + t10 关源 + `Min=Max=char_len`。
  **唯一同时**满足"随 char_len 单调缩放 + 质量最好 + 最快"的配置。

### 4.2 最佳 ≤5 万点（胫骨）

| cl | tet | J<0.3 | volerr |
|----|-----|-------|--------|
| 3 | 63115 | 3.21% | 1.54% | ← 超 5 万，弃 |
| **4** | **27670** | **3.66%** | **1.11%** | ← **推荐**（质量最优且 ≤5 万）|
| 5 | 15439 | 4.48% | 0.93% | |
| 6 | 8852 | 4.26% | 0.11% | |

### 4.3 raw 面为何全 FAIL

胫骨 `/tmp`→原始 STL 有 **17 非流形边 / 17 开放边**，gmsh 在
`generate(2)`/`generate(3)` 直接报：

```
Exception: Invalid boundary mesh (overlapping facets) on surface 2 surface 33
```

`reparam`/`reparam-baseline` 在 raw 上被探针**主动跳过**（报错原文）：
`skipped: surface has 17 open edges; forReparametrization=True needs a closed
manifold surface`。原因：`forReparametrization=True` 作用于非流形面会让 gmsh
重参数化**挂死**（本机在 raw+reparam 上实测卡死 >20 min，已加开放边守卫）。
⇒ 非流形骨面**必须**先走 `iso-dec`。

---

## 5. 推荐配置与 helper

**推荐配置 `reparam`**（写进 `part_meshing.py`）：

```python
gmsh.model.mesh.classifySurfaces(angle, True, True)   # forReparametrization=True
# 角度回退：(90, 75, 60, 45, 30, 0.1, 0.001)°，取第一个成功的（实测取到 90°）
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeMin", char_len)
gmsh.option.setNumber("Mesh.MeshSizeMax", char_len)
gmsh.option.setNumber("Mesh.Algorithm", 6)
gmsh.option.setNumber("Mesh.Algorithm3D", 1)
```

**Helper**：`scripts/opensim_fe/part_meshing.py`，公开入口
`mesh_bone_stl_mm(stl_path, char_len, ...)`（mm-native；不走 `stl_to_tet_gmsh` 的
米↔毫米 ×1000 入口）。`meshing.py` **一字未改**。

表面路由：闭合 2-流形（0 开放边 + 0 非流形边）→ 直接用 `raw`；否则
marching-cubes `iso-dec` 修复后再 `reparam`。

## 6. 重网格 before/after（实测）

### 6.1 胫骨 `tibia_r`（cl=4）

| 量 | before（体素路径，`parts/tibia_r/summary.json`） | **after（`gmsh-reparam(iso-dec)`，cl=4）** |
|----|---------------------------------------------------|---------------------------------------------|
| 方法 | `delaunay-voxel-repair` | `gmsh-reparam(iso-dec)` |
| tet | 38209 | **27670** |
| 节点 | 7892 | 6600 |
| 体积误差 | 2.03% | **1.11%** |
| `min_jacobian` | 1.02e-3 | 4.88e-3 |
| **`J<0.3`** | **17.41%** | **3.66%** |
| 用时 | 18.2 s | **2.5 s** |
| surface tri | — | 23904 |

产物体积 349976 mm³（参照 346130 mm³）。

### 6.2 颅骨 `skull_r_parietal`（cl=3）

| 量 | before（体素，cl=3） | **after（gmsh-reparam(iso-dec)，cl=3）** |
|----|----------------------|------------------------------------------|
| tet | 30714 | **23709** |
| 节点 | 7187 | 6373 |
| 网格体积 | 119465 mm³ | 117375 mm³ |
| 体积误差 | 96.55%（**参照错**：3466 cm³ 是颅腔） | 0.72%（对 iso 面自身） |
| `min_jacobian` | 1.26e-3 | 0.1273 |
| **`J<0.3`** | **17.55%** | **3.15%** |
| 用时 | 20.8 s | **2.0 s** |

⚠️ 顶骨的**真实骨体积**仍不可信：`iso-dec` 把薄壳增厚到 ~118 cm³（体素路径同样
~119 cm³）。参照 `3466 cm³` 是整个颅腔的包络，**不是**骨体积。`J<0.3` 的改善是真的，
但**几何保真度仍未验证**（与 S31 结论一致）。顶骨额外扫描：cl=2→73813(<br>超5万)/
J3.15；cl=3→23709/J3.15；cl=4→10168/J6.75；cl=5→5706/J7.92；cl=6→3663/J7.97。

### 6.3 跟骨对照 `calcaneus_r`（干净流形 → `raw` 路径，cl=3）

| 量 | 旧体素（文档） | 旧 gmsh-decimate（文档） | **新 gmsh-reparam(raw)，cl=3** |
|----|----------------|--------------------------|-------------------------------|
| tet | — | — | **17059** |
| 体积误差 | — | 0.00% | 0.49% |
| `J<0.3` | ~17.2% | 6.2% | **2.90%** |
| 用时 | — | — | 0.5 s |

跟骨额外扫描：cl=2→55915/J2.98/0.23%；**cl=3→17059/J2.90/0.49%**；
cl=4→7517/J3.06/0.91%；cl=5→4071/J3.07/1.36%。跟骨 `raw` 干净，无需修面，
证明 `reparam` 在**原始几何**上也稳定，质量优于既有 6.2%。

## 7. 可复现命令

仓库根、`$env:PYTHONPATH="src"`。

```powershell
# 1) 探针：全配置 × char_len 扫描（胫骨，raw + iso-dec）
& .venv\Scripts\python.exe scripts\opensim_fe\gmsh_size_probe.py `
    --surfaces raw,iso-dec `
    --configs baseline,reparam-baseline,no-boundary,field-constant,field-matheval,factor,reparam `
    --chars 3,4,5,6,8,10 --label tibia_r

# 2) 修复 helper：胫骨 cl=4 / 顶骨 cl=3 / 跟骨对照 cl=3
& .venv\Scripts\python.exe scripts\opensim_fe\part_meshing.py `
    --name tibia_r --stl temp\opensim_fe\parts\tibia_r\tibia_r_surface.stl `
    --char-len 4 --ref-volume-mm3 346130.49

& .venv\Scripts\python.exe scripts\opensim_fe\part_meshing.py `
    --name skull_r_parietal `
    --stl temp\opensim_fe\parts\skull_r_parietal\skull_r_parietal_surface.stl `
    --char-len 3

& .venv\Scripts\python.exe scripts\opensim_fe\part_meshing.py `
    --name calcaneus_r --stl temp\opensim_fe\thums_calcaneus\calcaneus_r_anatframe.stl `
    --char-len 3 --ref-volume-mm3 96205.19
```

产物（**新目录**，未覆盖既有产物）：

```
temp/opensim_fe/parts_gmsh/
  gmsh_size_probe_tibia_r.json                # 探针原始记录（含每条报错原文）
  tibia_r/{tibia_r_mesh.npz, tibia_r_summary.json, tibia_r_surface_used.stl}
  skull_r_parietal/{...}
  calcaneus_r/{...}
```

回归测试：

```powershell
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
```

## 8. 诚实边界与遗留

1. **非流形/自交骨面必须修面**：胫骨、顶骨的 raw STL 非流形，`forReparam=True`
   在非流形面上会让 gmsh **挂死**（已加守卫并如实记录）。必须走 `iso-dec`
   （marching cubes），这会改变几何、且对薄壳**增厚**。
2. **薄壳骨几何保真未解决**：顶骨无论体素还是 gmsh 都得到 ~118 cm³（真实板 ~40 cm³）；
   `volume_error_pct` 对薄壳**无意义**。这是"壳→实体"问题，非尺寸控制。
3. **皮质壳缺失**（S31 坑 #2）仍未解决，本任务不含。
4. **选项名**：任务列的 `Mesh.MeshSizeFromBoundary` 在 gmsh 4.15 **不存在**，
   实测报错已记录；等价物是 `Mesh.MeshSizeExtendFromBoundary`。
5. **`factor` 不可用**：`J<0.3`≈51%，明确否决。
6. **`field` 冗余**：背景场在 `ExtendFromBoundary=0` 后与 `no-boundary` 逐位相同，
   且平台化，不能单独当解。
7. **推荐 char_len 依赖骨体积**：胫骨 cl=4、顶骨 cl=3、跟骨 cl=3 是各骨的
   "≤5 万 + 最小 `J<0.3`"点；换骨需按 §7 命令重扫。
8. `meshing.py` **公开函数默认行为/签名未改**（新增 helper 独立模块）。

## 9. 回归测试结果

```
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
...........s........                                                     [100%]
19 passed, 1 skipped, 1 warning in 82.31s (0:01:22)
```

**19 passed / 1 skipped**，与改动前一致（本任务只新增脚本与产物，未触碰被测模块；
`.pytest_cache` 权限 warning 为环境既有，与本任务无关）。