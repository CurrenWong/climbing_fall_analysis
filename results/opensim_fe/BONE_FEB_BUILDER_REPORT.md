# 通用 THUMS 骨 FEBio 构建器（`bone_feb.py`）实测报告

> 状态：**实测**。所有数字来自本机一次完整运行（`scripts/opensim_fe/bone_feb_validate.py`），
> 复现命令见 §7。FEBio = `D:\Program\FEBioStudio\bin\febio4.exe` **4.13.0**，pyfebio **0.3.0**，
> Python 3.12 venv，`PYTHONPATH=src`。源甲板 = THUMS AM50 V7.1
> `model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`。
>
> 本报告**不覆盖**任何既有产物：只新增
> `scripts/opensim_fe/bone_feb.py`、`scripts/opensim_fe/bone_feb_validate.py`、
> `temp/opensim_fe/bone_feb/**` 与本报告。`l3_shell_feb.py` / `thums_feb.py` /
> `kmesh_io.py` / `shell_io.py` / `meshing.py` / `febio_run.py` **均未改动**。

---

## 0. 一句话结论

把"L3 专用 shell+solid 原型"泛化成 **`bone_feb.py`**：任意骨 = N 个独立固体域
（每域按 PID 读 `*MAT_*` 的 E/ν）+ M 个独立壳域（每域独立厚度），自动共节点、自动选轴、
自动均布轴压。**4 块骨全部 rc=0**，覆盖 **0 壳 / 1 壳 / 2 壳 / 9 变厚度壳**、
**纯 tet4 / tet4+hex8 混合** 四种形态：

| 骨 | 固体域 | 壳域 | 单元拓扑 | FEBio rc |
|----|--------|------|----------|----------|
| **L3**（回归） | 1（SPON） | 1（CORT 壳） | tet4 + quad4 | **0** |
| **tibia_r** | 3（SPON_end + SPON_center + CORT 实体） | **0** | tet4 + **hex8** | **0** |
| **parietal_r** | 1（diploë） | 2（external + internal） | **hex8** + quad4 | **0** |
| **R_HIPBONE**（可选） | 1（SPON） | **9**（变厚度 CORT） | tet4 + quad4 | **0** |

---

## 1. 模块 API

`scripts/opensim_fe/bone_feb.py`（可复用模块，非一次性脚本）：

```python
@dataclass
class Material:
    pid: int; E_mpa: float; nu: float
    src_keyword: str          # "*MAT_DAMAGE_2" / "*MAT_VISCOELASTIC" / ...
    src_card: list[float]     # 原始 card1 数字（诊断用）
    note: str                 # 例如 "elastic approximation of ..."

@dataclass
class SolidDomainSpec:
    pids: list[int]           # 该域覆盖的 THUMS PID（可多个）
    mat_name: str             # pyfebio material 名（全 spec 内唯一）
    mesh_name: str            # pyfebio element domain 名（唯一）
    material: Material | None # read_materials 填充
    topology: str = "auto"

@dataclass
class ShellDomainSpec:
    pids: list[int]; mat_name: str; mesh_name: str
    thickness_mm: float
    material: Material | None
    shell_type: str = "elastic-shell"

@dataclass
class BoneSpec:
    name: str
    solid_domains: list[SolidDomainSpec]
    shell_domains: list[ShellDomainSpec] = []

# --- 入口 ---
read_materials(master_k, pids) -> {pid: Material}
extract_mesh(spec, master_k, *, verbose=True) -> dict      # 共节点全局网格
build_bone_feb(spec, out_feb, *, load_n, top_frac=0.30, bot_frac=0.30,
               axis=None, with_shell=True, time_steps=1, mesh=None) -> dict
run_bone_feb(feb_path, *, workdir=None, timeout=900) -> {rc, elapsed_s, ...}
post_process(xplt, hdf5, *, solid_names, shell_specs, label="post") -> dict
```

### 设计要点

- **共节点**：`extract_mesh` 把每个 PID 引用的 **THUMS 原始节点 id 取并集**，
  整体 remap 成 1..N 的局部索引；solid 与 shell 共用同一 `<Nodes>` 块
  （与 `l3_shell_feb.py` 完全一致）。4 块骨实测 `shell_only=0`（壳节点 100% 是固体节点子集），
  **无 tie / 无接触 / 无缝合**。
- **按 PID 读材料**：`read_materials` 流式扫 `*MAT_*` 卡（跨 `*INCLUDE` 树），
  对每 PID 取第一个数字卡的 (E, ν)。支持：

  | LS-DYNA 卡 | 取法 |
  |-----------|------|
  | `*MAT_ELASTIC` | 直接 (E, ν) |
  | `*MAT_PIECEWISE_LINEAR_PLASTICITY` | 取弹性段 (E, ν)，屈服/塑性曲线**忽略**（注明） |
  | `*MAT_DAMAGE_2` / `*MAT_PLASTICITY_WITH_DAMAGE` | 取 (E, ν)，损伤/塑性**忽略**（注明） |
  | `*MAT_VISCOELASTIC` | (K, G0) → 等效 `E=9KG₀/(3K+G₀)`, `ν=(3K−2G₀)/(2(3K+G₀))`（注明） |

  未覆盖的卡（`*MAT_NULL` / `*MAT_LOW_DENSITY_FOAM` 等）**不臆造**：直接不返回该 PID 的条目，
  调用方 `_materialize` 抛错。
- **自动拓扑分类**（见 §5.1）：`*ELEMENT_SOLID` 按去重后节点数分流——
  **4 个 distinct → tet4**（THUMS 塌缩 hex），**>4 个 distinct → hex8**（真实六面体）。
- **通用轴选 + BC**：与 L3 原型同款——最长 bbox 轴（可 `axis=` 覆盖），
  顶/底 `top_frac`/`bot_frac` 带宽节点集，底面 `BCZeroDisplacement` 全固定、
  顶面等分 `NodalForce` 沿 `-axis` 压缩，`STATIC` 单步（`time_steps` 可调）。
- **0 壳**：`with_shell=False` 或 `shell_domains=[]` 都行。
- **多壳**：每个 `ShellDomainSpec` 独立 `ShellDomain(shell_thickness=...)`，互不干扰。

---

## 2. 泛化性证据（实测）

| 维度 | L3 | tibia_r | parietal_r | R_HIPBONE |
|------|----|---------|-----------|-----------|
| 固体域数 | 1 | **3** | 1 | 1 |
| 壳域数 | 1 | **0** | **2** | **9** |
| 壳厚度 (mm) | 1.39 | — | 1.5 / 1.5 | 1.5,1.5,1.0,0.75,2.0,3.0,2.2,1.8,1.0 |
| 固体拓扑 | tet4 | tet4 + **hex8** | **hex8** | tet4 |
| 壳拓扑 | quad4 | — | quad4 | quad4 |
| 全局节点 | 2,203 | 7,245 | 2,952 | 3,175 |
| 壳节点独有 | 0 | 0 | 0 | 0 |
| 材料数 | 2 | 3 | 2 | 2 |
| FEBio rc | 0 | 0 | 0 | 0 |

→ 同一条代码路径同时吃下 **0/1/2/9 壳**、**变厚度**、**tet4/hex8 混合**、
**1/2/3 材料域**，无需分支改造。

---

## 3. 每骨规模 + 材料 + σ_vm（200 N 轴压，实测）

### 3.1 规模与材料（按 PID 读自甲板）

| 骨 / 域 | PID | 材料卡 | E (MPa) | ν | 单元 |
|---------|-----|--------|---------|---|------|
| L3 `L3_spon` | 89001500, 89501500 | `*MAT_DAMAGE_2` | 40.0 | 0.450 | 8,792 tet4 |
| L3 `L3_cort` | 89001501, 89501501 | `*MAT_PLASTICITY_WITH_DAMAGE` | 13,020.0 | 0.300 | 1,340 quad4 (t=1.39) |
| tibia `spon_end` | 81000600 | `*MAT_PIECEWISE_LINEAR_PLASTICITY` | 160.0 | 0.450 | 6,651 tet4 |
| tibia `spon_center` | 81000601 | `*MAT_VISCOELASTIC` (K=2000,G₀=4) | **11.992** | **0.499** | 14,482 tet4 |
| tibia `cort` | 81000700 | `*MAT_PIECEWISE_LINEAR_PLASTICITY` | 18,000.0 | 0.300 | **2,417 hex8** |
| parietal `diploe` | 88000004 | `*MAT_DAMAGE_2` | 1,090.0 | 0.220 | **1,430 hex8** |
| parietal `ext` | 88000005 | `*MAT_PLASTICITY_WITH_DAMAGE` | 14,900.0 | 0.220 | 1,430 quad4 (t=1.5) |
| parietal `int` | 88000006 | `*MAT_PLASTICITY_WITH_DAMAGE` | 14,900.0 | 0.220 | 1,430 quad4 (t=1.5) |
| hip `spon` | 83500200 | `*MAT_DAMAGE_2` | 15.0 | 0.450 | 11,767 tet4 |
| hip `cort`×9 | 83500201..09 | `*MAT_PIECEWISE_LINEAR_PLASTICITY` | 17,300.0 | 0.300 | 2,232 quad4（9 域） |

### 3.2 σ_vm（von Mises，MPa；单元中心值）

**L3（rc=0）**

| 域 | n | max | p99 | p95 | mean | median |
|----|---|-----|-----|-----|------|--------|
| solid `L3_spon` | 8,792 | **0.0507** | 0.0297 | 0.0218 | 0.0088 | 0.0081 |
| shell `L3_cort` | 1,340 | **7.557** | 4.402 | 2.417 | 0.841 | 0.613 |

**tibia_r（rc=0，solid-only）**

| 域 | n | max | p99 | p95 | mean | median |
|----|---|-----|-----|-----|------|--------|
| `spon_end` | 6,651 | **0.0458** | 0.0327 | 0.0269 | 0.0107 | 0.0089 |
| `spon_center` | 14,482 | **0.0343** | 0.0181 | 0.0108 | 0.0034 | 0.0023 |
| `tibia_cort` | 2,417 | **35.76** | 31.86 | 24.46 | 4.691 | 1.202 |

**parietal_r（rc=0，1 solid + 2 壳）**

| 域 | n | max | p99 | p95 | mean | median |
|----|---|-----|-----|-----|------|--------|
| solid `diploe` | 1,430 | **0.533** | 0.374 | 0.229 | 0.082 | 0.068 |
| shell `ext` | 1,430 | **6.567** | 2.480 | 1.586 | 0.790 | 0.838 |
| shell `int` | 1,430 | **3.936** | 3.469 | 2.431 | 1.039 | 1.042 |

**R_HIPBONE（rc=0，1 solid + 9 变厚度壳）**

| 域 | n | max | p99 | p95 | mean | median |
|----|---|-----|-----|-----|------|--------|
| solid `hip_spon` | 11,767 | **0.0182** | 0.0129 | 0.0084 | 0.0024 | 0.0017 |
| shell `1.5` | 938 | **2.736** | 1.982 | 1.220 | 0.301 | 0.209 |
| shell `ac1.5` | 133 | **0.334** | 0.164 | 0.099 | 0.018 | 0.000 |
| shell `ac1.0` | 8 | **0.000** | 0.000 | 0.000 | 0.000 | 0.000 |
| shell `0.75` | 85 | **1.175** | 1.137 | 0.914 | 0.328 | 0.319 |
| shell `2.0` | 557 | **2.638** | 2.339 | 1.613 | 0.571 | 0.408 |
| shell `3.0` | 8 | **0.153** | 0.153 | 0.151 | 0.093 | 0.099 |
| shell `2.2` | 50 | **0.0046** | 0.0033 | 0.0008 | 0.0002 | 0.000 |
| shell `1.8` | 32 | **0.0177** | 0.0136 | 0.0029 | 0.0008 | 0.000 |
| shell `1.0` | 421 | **1.581** | 1.271 | 1.189 | 0.460 | 0.417 |

> 全部 4 骨的 solid 域 σ_vm 都读回成功；**壳域 σ_vm 也全部读回**
> （FEBio 4.13 的 nameless-shell 坑已绕过，见 §5.2），9 个髋壳的
> `n` 与预期逐一相等（938/133/8/85/557/8/50/32/421）。

---

## 4. L3 回归对比（vs 原型 `l3_shell_feb_result.json`，200 N）

原型 `l3_shell_feb.py` **硬编码** trabecular E=73.4/ν=0.45、cortical E=15000/ν=0.30；
新 builder **按 PID 读甲板**得 trabecular E=40/ν=0.45、cortical E=13020/ν=0.30。
几何、轴选、BC、载荷完全相同（同为 L3 L+R、Y 轴、top/bot 30%、200 N）。

### 4.1 有壳 vs 有壳

| 量 | 原型（E=73.4 / 15000） | 新 builder（E=40 / 13020） | 比值 新/原 |
|----|------------------------|----------------------------|-----------|
| trab σ_vm max | 0.0663 | **0.0507** | 0.77 |
| trab σ_vm p95 | 0.0294 | **0.0218** | 0.74 |
| trab σ_vm mean | 0.0120 | **0.0088** | 0.73 |
| cort σ_vm max | 7.381 | **7.557** | 1.02 |
| cort σ_vm p95 | 2.335 | **2.417** | 1.04 |
| cort σ_vm mean | 0.810 | **0.841** | 1.04 |

**量级一致**（全部同阶）。差异**唯一来源 = 材料**：新 builder 用甲板真实 E。
直觉上 trab σ 应 ∝ E 比 40/73.4=0.545，但实测 0.77——因为这是**双材料复合**，
载荷在壳/松质间按**刚度比**分配：原型的 E_cort/E_spon=204，新模型=325。
刚度比升高 → 壳相对更硬 → 松质分到的载荷更少，故其应力降得不如纯 E 比那么多；
壳应力则基本持平。**这正是"壳是载荷主体"的又一佐证。**

### 4.2 无壳（对照，`with_shell=False`）：

| 量 | 原型 无壳 | 新 builder 无壳 | 比值 |
|----|-----------|-----------------|------|
| trab σ_vm max | 1.7649 | **1.7514** | **0.992** |
| trab σ_vm p95 | 0.5633 | **0.5735** | 1.018 |
| trab σ_vm mean | 0.1906 | **0.1936** | 1.016 |

**单材料模型的应力与 E 无关**（线性：σ∝E·ε，ε∝F/EA → σ∝F/A），
所以尽管 E 从 73.4 改成 40，无壳 L3 的 σ_vm 与原型的差异 **< 2%**。
这条**强证据**说明新 builder 的**几何/BC/载荷/求解链路与原型逐字等价**，
有壳时的差异纯粹来自材料取值。

> 原型的"有壳使松质 σ 降 16–27×"结论在本代码路径**同样成立**：
> 无壳 max 1.751 / 有壳 max 0.051 = **34×**（同为 200 N）。

---

## 5. 关键修复 / 实测发现

### 5.1 ⭐ `*ELEMENT_SOLID` 里混着两种单元（最主要发现）

最初照搬 l3 的假设"8 节点卡一律取前 4 个 distinct 当 tet4"，对 L3 成立（L3 SPON 全是
塌缩 hex），但对**长骨皮质和颅骨穹窿**会造成灾难：THUMS 把它们的**真实 hex8**也写在同一
`*ELEMENT_SOLID` 块里。实测去重后节点数分布：

| PID | distinct=4（塌缩→tet4） | distinct=6 | distinct=8（真 hex8） |
|-----|------------------------|-----------|----------------------|
| L3 SPON 89001500 | 4,396 | 0 | 0 |
| hip SPON 83500200 | 11,767 | 0 | 0 |
| tibia spon 81000600/01 | 6,651 / 14,482 | 0 | 0 |
| **tibia CORT 81000700** | 0 | 2 | **2,415** |
| **parietal diploë 88000004** | 0 | 34 | **1,396** |

错误地把 hex8 塌成 tet4 的后果：tibia CORT 出现 **1,307 个负体积 tet + 7 个近零体积
sliver**，FEBio 报 `Negative jacobian detected during mesh initialization`（初始化即失败）；
即使翻正，sliver 仍让刚度阵近奇异，**5 N 就发散**（线性解 displacement≈2×10⁷）。
**修复：按 distinct 节点数分类**（4→tet4，>4→hex8），hex8 原样交给 FEBio。
修复后 tibia CORT 2417 个 hex8 定向全部正确（0 翻转）、200 N 干净收敛。

> 附带：`pyfebio.numpy_to_nodes` 用 `%e`（7 位有效数字）写坐标，会把 sliver tet 的符号
> 翻掉；本模块改为**全精度**（`%.15g`）自建 `<Nodes>`（不改 pyfebio，走公开 `Node/Nodes` API）。

### 5.2 FEBio 4.13 `.xplt` 壳域坑（多壳版）

l3 原型已记录：FEBio 4.13 在 `.xplt` 里**省略壳域的 `PLT_DOM_NAME`**，pyfebio 0.3.0
`to_hdf5` 会 `KeyError('name')`。**单壳**时用 `name=b""` 占位可绕过；但**多壳**时暴露两个新坑：

1. 多个无名壳域都映射到 HDF5 路径 `/meshes/0/domains/` → `ValueError: Unable to
   synchronously create dataset (name already exists)`（实测 parietal_r、R_HIPBONE）。
2. DOMAIN_HDR 的 `id` 字段是**材料 id，不是唯一域号**（实测 parietal 两个皮质壳都
   `id=2`，因共用材料 `bone_cort`）；而 `parse_state` 按**域位置**（`set_id`）查
   `mesh_dict['domains']`，于是 `KeyError(np.int32(3))`（实测 hipbone）。

**修复**（仍**不改 pyfebio 源码**，运行时 swap + `finally` 还原）：
补丁 `_parse_domain_section`，对每个域**按 1-based 位置重编 `id`**、并给无名域取唯一名
`__unnamed_<pos>`。`parse_mesh` 与 `parse_state` 都查 `mesh_dict['domains']`，两处自洽。
读回时：solid 按名字匹配；第 k 个壳域 → 第 k 个 `__unnamed_*`（FEBio 按 MeshDomains
顺序输出，实测 9 个髋壳的 `n` 逐一与预期相等，佐证顺序稳定）。

### 5.3 定向修正（safety net）

`_fix_tet_orientation` / `_fix_hex_orientation` 对负体积单元翻两个节点。**诚实说明**：
在最终 4 块骨的验证中**一次都没触发**（0 翻转）——因为 5.1 的拓扑修复从根上消除了
sliver。它保留为对其它甲板/退化的兜底（修复前它在 tibia CORT 上确实翻了 1,307 个）。

---

## 6. 诚实边界

1. **材料解析覆盖**：只支持 §1 表列的 5 种卡。`*MAT_NULL` / `*MAT_LOW_DENSITY_FOAM` /
   `*MAT_VISCOELASTIC` 之外的粘弹/超弹卡**不支持**（不臆造，缺卡即抛错）。
   塑料/损伤卡只取弹性段，**屈服与损伤被忽略**（每域 `note` 注明）。
2. **`*MAT_VISCOELASTIC` 等效弹性**：tibia `spon_center` 得 E≈11.99 MPa、ν≈0.499（近不可压）。
   这个值**远软于**其它松质，是甲板原样。200 N 下（配合 hex8 + 全精度坐标）能收敛，
   但该域本身是病态的（近不可压 + 位移法 tet4）；若单独用该域、或去掉 cort 支撑，
   低至 ~2 N 就发散——这是**源材料/单元格式的固有问题**，非 builder 缺陷。
3. **tet4 质量**：trabecular 域仍是 THUMS 的塌缩 hex→tet4，形状差、体积可能很小；
   本报告只用其解 σ_vm，未做单元质量统计/重网格。
4. **BC/载荷是通用几何，非解剖生理**：最长 bbox 轴 + 顶/底 30% 带宽均布力。
   目的是隔离变量、验证"构建器"本身；用于临床应力需换解剖级 BC。
5. **壳 σ_vm 读回**：依赖 §5.2 的 pyfebio 运行时补丁 + 按位置映射。若 FEBio 升级后
   改变了域输出顺序，映射会错位（但 `n_elem` 计数校验会暴露）。solid σ_vm 是**直接**
   按名字读，无此风险。
6. **不做**：单元删除/断裂判据；4 积分点应力重构（直接用单元中心值）；网格重划分；
   接触/tie（4 骨均 100% 共节点，`shell_only=0`，无需）。
7. **`with_shell=False` 只支持整骨去壳**，不是逐壳开关（本任务用不到逐壳）。

---

## 7. 复现命令

```powershell
# 仓库根、PowerShell、PYTHONPATH=src
$env:PYTHONPATH = "src"

# 4 块骨一键构建 + FEBio + σ_vm（写 temp/opensim_fe/bone_feb/**）
& .venv\Scripts\python.exe scripts\opensim_fe\bone_feb_validate.py

# 回归套件（应保持 19 passed / 1 skipped）
& .venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q
```

单骨复现（示例：L3 无壳对照，属 §4.2）见运行日志中 `[post:L3_noshell]`。
输出（实测）：

- `temp/opensim_fe/bone_feb/L3/{L3.feb,.xplt,.hdf5,.log,L3_result.json}`
- `temp/opensim_fe/bone_feb/L3_noshell/{L3_noshell.feb,...,L3_noshell_result.json}`
- `temp/opensim_fe/bone_feb/tibia_r/{...}`
- `temp/opensim_fe/bone_feb/parietal_r/{...}`
- `temp/opensim_fe/bone_feb/R_HIPBONE/{...}`
- `temp/opensim_fe/bone_feb/_summary.json`、`_run.log`

---

## 8. 产物清单

| 文件 | 说明 |
|------|------|
| `scripts/opensim_fe/bone_feb.py` | **新模块**（通用 builder） |
| `scripts/opensim_fe/bone_feb_validate.py` | **新**验证驱动（4 骨定义 + 构建 + 求解 + 报告 JSON） |
| `temp/opensim_fe/bone_feb/{L3,tibia_r,parietal_r,R_HIPBONE,L3_noshell}/**` | **新目录**，每骨 `.feb/.xplt/.hdf5/.log/_result.json` |
| `temp/opensim_fe/bone_feb/_summary.json` / `_run.log` | 汇总 / 完整 stdout |
| `results/opensim_fe/BONE_FEB_BUILDER_REPORT.md` | 本报告 |

未改动：`l3_shell_feb.py`、`thums_feb.py`、`kmesh_io.py`、`shell_io.py`、`meshing.py`、
`febio_run.py`（可新增/调用，不动默认行为/签名）；未触碰 `scripts/ankle_fe/`、
`temp/pyfebio_demo/`、`temp/ankle_*`、`FE_PIPELINE_TODO.md`；
未覆盖任何既有产物。
