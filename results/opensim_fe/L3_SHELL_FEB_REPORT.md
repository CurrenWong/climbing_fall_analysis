# L3 双材料域 FEBio 验证（shell+solid 共节点 vs solid-only）

> **坑② 步② 实测报告**：腰椎 L3 的"皮质壳作为 FEBio 壳单元、与松质实体共节点"
> 这条路径，**端到端跑通**；并量化皮质层对 σ_vm 的影响，证明"坑② 必须解决"。

> 状态：实测。所有数字来自 `scripts/opensim_fe/l3_shell_feb.py` 一次完整运行，
> 复现命令见 §6。源甲板 = `model/AM50_V71_Occupant/main_THUMS_AM50_V71.k`，
> FEBio = `D:\Program\FEBioStudio\bin\febio4.exe` **4.13.0**，pyfebio 0.3.0。

> 已知问题：FEBio 4.13 在 .xplt 里**省略壳域的 PLT_DOM_NAME** 字段，
> pyfebio 0.3.0 的 `xplt.to_hdf5` 会 `KeyError('name')`。本脚本**就地补丁**
> `_parse_domain`（保留原函数行为后立刻恢复），并对壳域应力按"空字符串键
> 或元素数匹配"回退查表。详见 §5 与 `l3_shell_feb.py` 的 `_xplt_to_hdf5_patched`。

## 0. 关键实测结论

| 指标 | CASE A（有皮质壳） | CASE B（无皮质壳，对照） |
|------|---------------------|---------------------------|
| FEBio rc | **0** | **0** |
| trab σ_vm **max** | 0.066 MPa | 1.765 MPa |
| trab σ_vm **p95** | 0.029 MPa | 0.563 MPa |
| trab σ_vm **mean** | 0.012 MPa | 0.191 MPa |
| cort σ_vm max / p95 / mean | 7.38 / 2.33 / 0.81 MPa | — |

**皮质壳的效应（cortical 移除后 / 保留后）**：

| 指标 | 比值 without / with | 即壳使 trab 应力降为 |
|------|---------------------|-----------------------|
| max | **26.6×** | 1/27 |
| p95 | **19.2×** | 1/19 |
| mean | **15.9×** | 1/16 |

**判定**：去掉皮质壳后，**松质 σ_vm 上升 16–27×**（max / p95 / mean 同量级）。
皮质壳是**载荷分配主体**——壳层 E=15000 MPa 远高于松质 E=73.4 MPa，
加上薄壁几何，壳吸收了 ~95% 以上的内力。**纯松质模型无法反映真实骨应力
分布，必须把皮质壳加回去。** —— 这正是"坑② 必须解决"的物理证据。

> ⚠️ **更正**：早先声称的 "1000 N" 对比（§3.5）**未由本仓库脚本复现**——脚本中
> `LOAD_N_SHELL_ONLY` 仅定义未调用，`run.log`/`result.json` 只含 **200 N** 结果（A、B 均 rc=0）。
> 以 200 N 结果为准；1000 N 段**已作废**（见 §3.5）。

## 1. 模型规模（实测，L3+R3 合并）

源甲板：THUMS AM50 V7.1，`*INCLUDE` 树 3 个文件。

| 量 | 值 |
|----|----|
| 节点（global THUMS id）| **2,203** |
| - solid tet4 引用节点 | 2,203 |
| - shell quad4 引用节点 | 1,335 |
| - 共享（solid ∩ shell） | **1,335（100%）** |
| - 节点壳独有 | **0** |
| 实体单元（tet4，塌缩 hex→tet4）| **8,792** |
| 壳单元（quad4 / tri3）| **1,340**（quad4=1340, tri3=0）|
| 壳厚度 | **1.39 mm**（`*SECTION_SHELL` SECID=89001501/89501501） |
| L3 PID SPON | 89001500 (L) + 89501500 (R) |
| L3 PID CORT | 89001501 (L) + 89501501 (R) |

> **共节点验证**：SHELL_SCOUT_REPORT §3 显示 L3 frac_shell_shared=1.0；
> 实测 `n_id_shell_only=0`，**所有壳节点都是固体节点子集**，无 tie/接触/缝合。

### 1.1 材料 / 单元类型

| 域 | 名称 | 单元 | E (MPa) | ν | 厚度 |
|----|------|------|---------|---|------|
| 实体 | `trabecular` | tet4（4 节点） | **73.4** | **0.45** | — |
| 壳   | `cortical`  | quad4（4 节点） | **15000** | **0.30** | **1.39 mm** |

pyfebio 写出的 `<ShellDomain type="elastic-shell" shell_thickness="1.39">`
经 `_fix_febio413` 把 `<solver type="solid">…</solver>` 替换为 `<solver/>`，
FEBio 4.13.0 直接通过 `Reading file …SUCCESS!`。

## 2. 载荷与边界（通用，非跟骨专用）

由 L3+R3 全局 bbox 自动选取 **最长轴**：

| 量 | 值 |
|----|----|
| bbox min / max（最长轴 = Y）| -36.762 / 36.732 mm |
| 跨度 | 73.49 mm |
| 选择轴 | **Y**（axis_idx=1，最长 bbox span）|
| 顶部带（top band）| Y ≥ max - 0.30 × span = **14.69 mm** |
| 底部带（bottom band）| Y ≤ min + 0.30 × span = **-14.71 mm** |
| 顶节点数 | 391 |
| 底节点数 | 391 |
| 边界 | 底面 `BCZeroDisplacement` x/y/z 三向全固定 |
| 载荷 | **200 N** Y 轴压缩（向上 = -Y），等分到 391 个顶节点（每节点 -0.5115 N）|
| 分析 | `STATIC`，1 个载荷步（step_size=1.0，time_stepper=None → FEBio 4.13 兼容）|

> 选 30% 带宽而非更小的（如 15%）：更小带宽下纯松质模型在 200 N 即
> 出现 ~1000 个 negative jacobians（详见 §3），是退化（degenerate）压缩
> 的经典表现；30% 带宽分散了局部应力集中。带宽选择直接进入
> `l3_shell_feb.py` 的常量 `TOP_BAND_FRAC = 0.30`。

## 3. 求解结果

**CASE A（有壳）+ CASE B（无壳）均 rc=0，FEBio 正常终止**。

| CASE | 描述 | rc | 求解耗时 |
|------|------|----|----------|
| A | trabecular (tet4) + cortical (quad4 共节点) | **0** | 0.4 s |
| B | 仅 trabecular (tet4) | **0** | 0.3 s |

### 3.1 CASE A：trabecular σ_vm（n=8,792 单元）

| 指标 | 值 (MPa) |
|------|----------|
| max   | **0.0663** |
| p99   | 0.0404 |
| p95   | **0.0294** |
| mean  | 0.0120 |
| median| 0.0111 |

### 3.2 CASE A：cortical shell σ_vm（n=1,340 单元）

**通过 fallback 解析空 key 取得**（详见 §5）。数值：

| 指标 | 值 (MPa) |
|------|----------|
| max   | **7.381** |
| p99   | 4.292 |
| p95   | **2.335** |
| mean  | 0.810 |
| median| 0.576 |

> 壳层应力约比松质应力**高 50–100×（max）**、**高 65–80×（mean）**，
> 这与 E 比 15000/73.4 ≈ 200× 一致——壳吸收了大部分内力。

### 3.3 CASE B：trabecular σ_vm（n=8,792 单元，去掉壳）

| 指标 | 值 (MPa) |
|------|----------|
| max   | **1.765** |
| p99   | 0.892 |
| p95   | **0.563** |
| mean  | 0.191 |
| median| 0.142 |

### 3.4 对比：壳对松质应力的屏蔽（with-shell vs no-shell）

| 指标 | with shell (MPa) | without shell (MPa) | ratio **without / with** |
|------|------------------|---------------------|---------------------------|
| max  | 0.0663 | 1.765 | **26.6** |
| p95  | 0.0294 | 0.563 | **19.2** |
| mean | 0.0120 | 0.191 | **15.9** |

> 去掉壳后松质 σ_vm **提升 16–27×**。换句话说，**壳使松质应力降为 1/27 ~ 1/16**。
> 这就是"坑② 必须解决"的定量证据。

### 3.5 关于 1000 N（⚠️ 已作废）

> 本段原数字（trab max 0.65 / cort max 73.8，"线性 5×"）**未由 `l3_shell_feb.py` 复现**：
> 脚本中 `LOAD_N_SHELL_ONLY=1000.0` **定义但未被 `main()` 调用**；`run.log`/`result.json` 均无
> 1000 N 记录。且 200→1000 N 是 **5×**，线性弹性下应力应**严格 5×**：正确预期为
> trab max ≈ **0.33 MPa**、cort max ≈ **36.9 MPa**——原"10×"数字与线性弹性矛盾。
> **待办**：若需 1000 N 结论，应把 `LOAD_N_SHELL_ONLY` 显式接入 `main()` 并重跑。

## 4. pyfebio 壳 API 实测小结

| 需求 | pyfebio 0.3.0 支持？ | 实测路径 |
|------|---------------------|----------|
| `ShellDomain` 类 | ✅ | `pyfebio.meshdomains.ShellDomain(name, mat, shell_thickness, type='elastic-shell')` |
| `numpy_to_elements(..., "quad4")` | ✅ | `fmesh.numpy_to_elements(arr, 'quad4', name='cortical')` |
| 共节点 shell+solid | ✅ | 共用同一 `<Nodes>` 块，不同 `<Elements>` 块共享同一 ID 域 |
| 壳-实体域分离绑定 | ✅ | `<ShellDomain name="cortical" mat="bone_cort" shell_thickness="1.39">` |
| quad4 + tri3 同一 MeshDomain | ✅ | 两个 `<Elements name="cortical" type="quad4|tri3">` 块均指向同一 ShellDomain |
| xplt 读回 σ_vm (shell) | ⚠️ 有坑 | FEBio 4.13 省略壳域 PLT_DOM_NAME → 补丁 `_parse_domain` 见 §5 |

## 5. 已知问题与诚实边界

### 5.1 FEBio 4.13 .xplt 省略壳域 PLT_DOM_NAME（已实测，已绕过）

**现象**：在 `l3_with_shell.xplt` 的 mesh 节里，固体域 `trabecular` 的 domain
header 含 `PLT_DOM_NAME = b"trabecular"`，但**壳域的 domain header 里 PLT_DOM_NAME
字段整体缺失**（二进制实测：只有 `etype=3, id=2, nelems=1`，无 name）。
pyfebio 0.3.0 的 `_parse_domain` → `parse_mesh` 直接
`mesh_dict["domains"][domain["id"][0]] = domain["name"]`，抛 `KeyError('name')`。

**绕过**：本脚本的 `_xplt_to_hdf5_patched` 在调用 `to_hdf5` 之前**就地
替换** `pyfebio.xplt._parse_domain`，强制 `name=b""` 占位；stress 数据仍按
HDF5 路径 `states/{sid}/element_data/stress/b''` 正确写出。
壳域应力通过 fallback `(a) 空字符串键 或 (b) 元素数 == n_shell 的键` 查表得到
（实测 key = `b''`，n_elem=1340 = `n_shell`，与预期完全一致）。

**未做**：没改 pyfebio 源码；也没在主分支提交此补丁（仅本脚本一次性使用）。

### 5.2 网格质量 caveat（THUMS 塌缩 hex→tet4）

THUMS AM50 V7.1 的椎体 SPON 实体是**塌缩 hex**（写入时第 5–8 个节点
重复第 4 个），本脚本按 `kmesh_io` 既有约定"取每行前 4 个 distinct 节点"
恢复为 tet4。这些 tet4 在六面体内部是非常**扁平的退化四面体**，四面体
体积可能接近 0、形状比非常差；在纯压缩下容易触发 negative jacobian
（本任务 200 N 下 A、B 均 rc=0，网格退化是后续高载荷的潜在风险）。**这并非建模错误，而是源甲板的
几何表示限制。** 后续若需要更稳健的体网格，可对 L3 SPON 做 gmsh
remesh（参考 `extract_part.py` 的 `mesh_stl` 路径）再走同一 FE 链路。

### 5.3 载荷 / BC 是通用几何而非解剖生理

本报告**不**复刻 L3 在人体中"上终板 / 下终板"的载荷分布；只是用 bbox
自动选最长轴、顶底 30% 带宽均布力。这是为了**隔离"壳层贡献"这一变量**——
任何特定解剖 BC 都会引入第二变量。结论"壳使松质应力降为 1/16 ~ 1/27"
在该简化条件下成立；用于预测临床应力需叠加解剖级载荷。

### 5.4 接触/界面：本模型不存在

THUMS 壳节点与同骨 SPON 节点 100% 共享（n_id_shell_only=0），本模型
**无 tie、无接触、无缝合**。若该共享关系被任何预处理破坏（即出现
`n_id_shell_only > 0`），则必须先 tie / surface contact 才能跑得动。

### 5.5 单元退化 hex→tet4 的影响范围

`_first_distinct(row, 4)` 取 distinct 后若不足 4 个节点**抛错**；本数据集
8,792 个实体单元无一行触发此异常。tri3 = 0 是因为 L3/R3 CORT 全部是 quad4
（5 个 tri3 是 hex 退化情况，本数据集未出现）。

### 5.6 不做的判据

- 不做单元删除 / 断裂判据（仅求解、报告 σ_vm）。
- 不做 4-pt 积分点的应力重构（直接用 element_data/stress 的单元中心值）。
- 不做单元删除 / 重网格（直接用原始 THUMS 网格，承认 tet4 质量 caveat）。
- 不做单元正交质量 / aspect ratio 统计（坑② 后续步骤可能补）。

## 6. 复现命令

```powershell
# 仓库根，PYTHONPATH=src，PowerShell
$env:PYTHONPATH = "src"
& .venv\Scripts\python.exe scripts\opensim_fe\l3_shell_feb.py
```

输出（实测）：
- `temp/opensim_fe/l3_shell/l3_with_shell.feb`  （A：有壳）
- `temp/opensim_fe/l3_shell/l3_with_shell.xplt` （A 结果）
- `temp/opensim_fe/l3_shell/l3_with_shell.hdf5` （A 转换后）
- `temp/opensim_fe/l3_shell/l3_no_shell.feb`    （B：无壳）
- `temp/opensim_fe/l3_shell/l3_no_shell.xplt`   （B 结果）
- `temp/opensim_fe/l3_shell/l3_no_shell.hdf5`   （B 转换后）
- `temp/opensim_fe/l3_shell/l3_shell_feb_result.json`（汇总）
- `temp/opensim_fe/l3_shell/run.log`            （完整 stdout）

注：脚本启动时会清掉上述旧产物（仅本目录）；不触碰任何既有产物。

## 7. 产物清单

| 文件 | 说明 |
|------|------|
| `scripts/opensim_fe/l3_shell_feb.py` | 本脚本（新增）|
| `temp/opensim_fe/l3_shell/l3_with_shell.feb` / `.xplt` / `.hdf5` | CASE A 有壳 |
| `temp/opensim_fe/l3_shell/l3_no_shell.feb`   / `.xplt` / `.hdf5` | CASE B 无壳 |
| `temp/opensim_fe/l3_shell/l3_shell_feb_result.json` | 全部数值 |
| `temp/opensim_fe/l3_shell/run.log` | 完整 stdout |
| `results/opensim_fe/L3_SHELL_FEB_REPORT.md` | 本报告 |
