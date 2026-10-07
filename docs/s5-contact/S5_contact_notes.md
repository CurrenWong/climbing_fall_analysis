# S5-contact 配方提取笔记（§3.4 / T4 — 只读研究）

> 写于 2026-10-05。**只读**：未触碰 `scripts/ankle_fe/**`、`temp/pyfebio_demo/**`。
> 唯一新增文件即本文。
> 目标：下一位 agent 拿到此文件即可写 `src/climbing/coupling/pad_contact.py`，**不需要再读 ankle_fe**。
>
> 来源标注约定：
> - **proven in repo** = `scripts/ankle_fe/build_feb.py` 在 FEBio 4.13 上跑通过（脚踝 FE 冒烟已 rc=0）；
> - **demo-only / unverified** = 仅出现在 `temp/pyfebio_demo/**` 的探针脚本里、尚未在生产 .feb 中实际产出；
> - **NOT FOUND IN REPO** = 文档里看到、但本仓库任何位置都没出现过的字段。
>
> 单位 mm–N–MPa–s；FEBio 路径 `D:\Program\FEBioStudio\bin\febio4.exe`。

---

## 0. 一次性快照

| 项 | 状态 | 主出处 |
|---|---|---|
| 1. working sliding-elastic block | **proven in repo** | `scripts/ankle_fe/build_feb.py:674-685` |
| 2. tied-elastic block | **proven in repo** | `scripts/ankle_fe/build_feb.py:662-673` |
| 3. section ordering | **proven in repo**（踩坑注释） | `scripts/ankle_fe/build_feb.py:600,620-621` |
| 4. surface = real facet | **proven in repo**（踩坑注释） | `scripts/ankle_fe/build_feb.py:300-306,348-352,586-588,116` |
| 5. plotfile 变量（无 contact traction） | **proven in repo** | `scripts/ankle_fe/build_feb.py:687-693` |
| 6. pyfebio API + contact 类清单 | **proven in repo** + demo 探针 | `temp/pyfebio_demo/gen_load_xml.py:8-12,45-50` + `dir(pyfebio.contact)` 探针 |
| 7. Material 块（hyper + raw XML） | 双向（demo + raw proven） | `temp/pyfebio_demo/demo.py:56-67` + `scripts/ankle_fe/build_feb.py:561-569,74-88` |
| 8. minimal contact deck | 组合 §1–§7 | 同上 |

---

## 1. Working sliding-elastic contact block（**proven in repo**）

出处：`scripts/ankle_fe/build_feb.py:674-685`，由 `_pn`（= `SurfacePair` 名）驱动循环写出。
紧接其后是 `<Output>`（`build_feb.py:687`）和 `<Control>`（`build_feb.py:696-703`），
对应 FEBio 4.13 单步简写。

```xml
<contact name="tibiotalar" surface_pair="tibiotalar_pair" type="sliding-elastic">
  <laugon>PENALTY</laugon>
  <penalty>1</penalty>
  <auto_penalty>1</auto_penalty>
  <two_pass>0</two_pass>
  <node_reloc>0</node_reloc>
  <symmetric_stiffness>0</symmetric_stiffness>
  <tolerance>0.02</tolerance>
</contact>
```

> **逐字引用**：`scripts/ankle_fe/build_feb.py:674-685`（`build_feb.py:676` 注释写明：
> "contact 的 surface_pair 引用的是 **SurfacePair 的名字**（User Manual 3.14）"）。

### 1.1 子标签逐项表

| 子标签 | 类型 / 域 | 取值（**proven**） | 备注（来源） |
|---|---|---|---|
| `name` | 属性 | "tibiotalar"（原 `_pn[:-5]` 去掉 `_pair`） | `build_feb.py:676` |
| `surface_pair` | 属性 | **引用 `<SurfacePair name=…>` 的名字**（不是 primary/secondary 名） | `build_feb.py:675-676` 注释 + User Manual 3.14 |
| `type` | 属性 | `"sliding-elastic"`（字面量） | `build_feb.py:676-677` |
| `<laugon>` | 子元素 | `"PENALTY"`（环境变量 `LAUGON`，默认 `PENALTY`） | `build_feb.py:34,678` |
| `<penalty>` | 子元素 | 浮点（默认 `1`） | `build_feb.py:37,679` |
| `<auto_penalty>` | 子元素 | `"1"` 当 `laugon=="PENALTY"`；否则 `"0"` | `build_feb.py:680`（**proven 联动规则**） |
| `<two_pass>` | 子元素 | `"0"` 或 `"1"`（env `TWO_PASS`，默认 `0`） | `build_feb.py:36,681` |
| `<node_reloc>` | 子元素 | `"0"` 或 `"1"`（env `NODE_RELOC`，默认 `0`） | `build_feb.py:35,682` |
| `<symmetric_stiffness>` | 子元素 | `"0"`（写死） | `build_feb.py:683` |
| `<tolerance>` | 子元素 | `0.02`（写死） | `build_feb.py:684` |

### 1.2 没在生产 XML 里出现的字段（**demo-only / unverified**）

`dir(pyfebio.contact)` 探针给出 `SlidingElastic.model_fields` 还含：

- `gaptol`, `minaug`, `maxaug`, `knmult`, `seg_up`, `update_penalty`, `smooth_aug`,
  `tension`, `fric_coeff`, `flip_primary`, `flip_secondary`,
  `shell_bottom_primary`, `shell_bottom_secondary`, `offset`。

**这些都没有在 `build_feb.py` 的生产 XML 里被写出**——也就是说
"FEBio 4.13 上 friction-less sliding-elastic 跑通"是 **proven**；
而任何含 `<fric_coeff>` 的滑动接触配方在本仓库内是 **demo-only / unverified**。

> 引申：写 `pad_contact.py` 的摩擦（μ）参数时，**不要从 pyfebio 的 kwargs 套默认值**；
> 必须先在最小冒烟 .feb 里确认 FEBio 4.13 接受并能写反应力。

---

## 2. Tied-elastic block + 何时用（**proven in repo**）

出处：`scripts/ankle_fe/build_feb.py:662-673`，**手册 §3.14.6** 引用见 `build_feb.py:654`。

```xml
<contact name="tie_tibiotalar_tibia" surface_pair="tie_tibiotalar_tibia" type="tied-elastic">
  <penalty>1.0</penalty>
  <auto_penalty>1</auto_penalty>
  <two_pass>0</two_pass>
  <laugon>PENALTY</laugon>
  <tolerance>1.0</tolerance>
  <symmetric_stiffness>0</symmetric_stiffness>
  <search_tol>0.01</search_tol>
  <search_radius>1.0</search_radius>
</contact>
```

> **逐字引用**：`scripts/ankle_fe/build_feb.py:662-673`；注释行 654–656 明确写
> "参数取手册默认值：penalty=1.0 / auto_penalty=0 / two_pass=0 / tolerance=1.0 /
> gaptol=0(off) / search_tol=0.01 / search_radius=1.0"。

### 2.1 tied vs sliding 选择

仓库原文（`scripts/ankle_fe/build_feb.py:335-339`）：

> 手册 §3.14.5/§3.14.6：tied 用来「连接两个非共形网格」，约束
> primary 的**节点**连到 secondary 的**面**上；tied-elastic 在此基础上
> 强制界面两侧**位移连续**（solid-solid）。这里 primary = 脱开的软骨底面，
> secondary = 粗骨关节面。目的：软骨可独立加密，骨保持已验证的粗网格。

判定准则：

| 场景 | 用 |
|---|---|
| 两片网格**几何共形**（同节点编号或同坐标） | **不需要**任何 contact |
| 两片网格**不共形**，但**不允许相对滑动**（骨-软骨界面、内/外关节囊） | **tied-elastic** |
| 两片网格**不共形**，**允许相对滑动**（足底-垫、关节软骨对软骨） | **sliding-elastic**（可加 μ） |

仓库里 S5-contact 的目标场景（足底-垫）→ **sliding-elastic**（§1）。
tied-elastic 仅在同模型内部、骨与软骨分块加密时才需要。

### 2.2 tied 与 sliding 子标签差异（**proven**）

`Tied*` 多了 `search_tol` / `search_radius`，但**没有** `node_reloc`、`fric_coeff`、
`smooth_aug`、`tension`、`flip_primary/secondary`、`shell_bottom_*`、
`offset`——这些是 sliding 系列专属。

---

## 3. .feb 段序规则（**proven in repo**，踩坑注释）

### 3.1 顶层段序：Loads → Boundary → Contact

仓库原文（`scripts/ankle_fe/build_feb.py:620-621`）：

> ★ 段序必须为 Loads -> Boundary -> Contact（pyfebio 生成的 XML 即此顺序）；
> 顺序错会报 "unrecognized tag" 这种误导性错误。

FEBio 4.13 实际接受的顶层顺序（来自 `build_feb.py:557-704`，按出现顺序）：

```
<febio_spec version="4.0">
  <Module type="solid"/>
  <Globals><Constants>…</Constants></Globals>
  <Material>…</Material>
  <Mesh>…</Mesh>
  <MeshDomains>…</MeshDomains>
  <LoadData>…</LoadData>
  <Loads>…</Loads>
  <Boundary>…</Boundary>
  <Contact>…</Contact>     ← 必须紧跟 <Boundary>
  <Output>…</Output>
  <Control>…</Control>     ← 顶层，不是 <Step> 内
</febio_spec>
```

### 3.2 Mesh 内部：`<SurfacePair>` 必须跟在 `<Surface>` 之后

仓库原文（`scripts/ankle_fe/build_feb.py:600`）：

> ★ SurfacePair 必须定义在 Surface **之后**（User Manual §3.6.8）

`build_feb.py` 写顺序（节选）：

```python
# 595-599：先写所有 <Surface>
for nm, sel in all_surfs.items():
    w(f'    <Surface name="{nm}">') ...
# 600-610：再写所有 <SurfacePair>
for _tn, _tp, _ts in tie_pairs: ...
for _pn, _pa, _pb in PAIRS: ...
```

### 3.3 其他（仅注释，无 .feb 验证）

`build_feb.py:582`：`<NodeSet>` 的文本**直接就是节点列表**，不能嵌套 `<node>`（User Manual §3.6.3）。
`build_feb.py:616-618`：`<SolidDomain mat=…>` 填**材料名**而不是数字 id，否则
`invalid value for attribute "mat"`。两者均属 §3 衍生，但**不是 contact 专属**——只在骨架里复用。

---

## 4. Surface 构造：必须是真实单元 facet（**proven in repo**，3 处踩坑注释）

### 4.1 规则

`<Surface>` 内嵌的**子标签必须按单元类型命名**（`<tri3>` / `<quad4>` / `<tri6>` / `<quad8>`），
且**面节点必须真实属于该单元**（即网格里的 tet/hex/penta facet，不是 npz 里的"虚拟四边形"）。

### 4.2 仓库三处同源注释

| 位置 | 原文摘录 |
|---|---|
| `scripts/ankle_fe/build_feb.py:300-306`（`cart_outer` docstring） | "必须从 **tet 单元的真实外表面** 取，不能从 npz 的原始四边形取！……接触面若仍用原始 offset 四边形（其三角化方式与 tet 分解不一致），面就**与网格实际的面不重合**，FEBio 的接触搜索关联不上 ⇒ 接触完全不检测" |
| `scripts/ankle_fe/build_feb.py:348-352`（tied-elastic 注释） | "面必须是**单元的真实 facet**，否则 FEBio 报 'N invalid facets' 并静默丢弃（User Manual §3.6.5 的老坑，cart_outer 踩过一次）。软骨已全部拆成 tet4 ⇒ facet 只能是**三角形**；npz 里底面却是 quad/tri 混合，直接照抄会丢掉全部 quad（实测 93 面里 61 个 invalid = 正好是 hex8 的个数）" |
| `scripts/ankle_fe/build_feb.py:586-588`（`<Surface>` 写块） | "★ Surface：子标签**按单元类型命名** `<quad4>`/`<tri3>`（User Manual §3.6.5）。既不是 `<elem>`，也不是 `<Elements type="quad4">` —— 用错会报 'invalid value for attribute surface'" |
| 契约 §5 红线 | "接触面必须是**真实单元面**（tet facet），不能是 npz 里的四边形。" |

### 4.3 接触面来源算法（伪代码）

依据 `cart_outer`（`build_feb.py:299-328`）总结的"通过 elem_blocks 反查 facet"模式：

```python
# 给一个 domain 名，遍历该域所有 tet4 单元，挑"全节点属于指定集合"的面 → 朝外定向
def outer_facets_from_real_tets(domain: str, must_have_nodes: set[int]) -> list[list[int]]:
    out, seen = [], set()
    for _, ids in elem_blocks[(domain, "tet4")]:        # 真实装配好的 tet 单元
        for f in TET_FACES:                             # tet4 的 4 个三角面索引
            fv = [ids[z] for z in f]
            if not all(x in must_have_nodes for x in fv):
                continue
            k = tuple(sorted(fv))
            if k in seen or not facet_ok(fv):           # 3 distinct nodes, area > MIN_FACET_AREA
                continue
            seen.add(k)
            # 朝外定向（法向背离对顶点）
            opp = ids[[z for z in range(4) if z not in f][0]]
            pts = XYZ[[i - 1 for i in fv]]
            n = sum(np.cross(pts[i], pts[(i + 1) % 3]) for i in range(3))
            if np.dot(n, pts.mean(0) - XYZ[opp - 1]) < 0:
                fv = fv[::-1]
            out.append(fv)
    return out
```

⚠ 这段**仅作为算法意图复刻**，**不是**直接 copy 自 `cart_outer`——里面
`elem_blocks` / `XYZ` / `MIN_FACET_AREA` 都是 `build_feb.py` 的局部变量。
新模块应当按同思路自己建索引（详见 §8）。

---

## 5. plotfile 变量（**proven in repo**）

### 5.1 实际产出的 XML

`scripts/ankle_fe/build_feb.py:687-693`：

```python
w('  <Output>')
# ★ 不要输出 "contact traction"：它会写出 PLT_FACE_DATA(surface_data)，
#   而 pyfebio.xplt.to_hdf5 的 parse_state 在查 mesh_dict["surfaces"][set_id] 时
#   会 KeyError（它不认我们这些自定义 Surface 集合）→ 整份 .xplt 都读不出来。
w('    <plotfile type="febio"><var type="displacement"/>'
  '<var type="stress"/><var type="reaction forces"/></plotfile>')
w('  </Output>')
```

即，**唯一可读回的三项**：

| var type | HDF5 路径 | 用途 |
|---|---|---|
| `displacement` | `states/{last}/node_data/displacement/all` | 节点位移 (N,3) mm |
| `stress` | `states/{last}/element_data/stress/...` | 单元 Cauchy 应力 |
| `reaction forces` | `states/{last}/node_data/reaction forces/all` | 节点反力 (N,3) N |

### 5.2 `contact traction` 不可用（**proven in repo**，被显式禁用）

`build_feb.py:688-690` 注释明确：

> 不要输出 "contact traction"：它会写出 PLT_FACE_DATA(surface_data)，
> 而 pyfebio.xplt.to_hdf5 的 parse_state 在查 mesh_dict["surfaces"][set_id] 时
> 会 KeyError（它不认我们这些自定义 Surface 集合）→ 整份 .xplt 都读不出来。

**后果**：要在后处理拿"接触压力分布 / 接触合力"，必须用以下任一方法
（**均未在本仓库实现**，属于 S5 的 T3 工作）：

1. 节点反力求和（脚底 BC 节点的反力即地面支承反力的镜像）—— 简单但不区分分布；
2. **关闭** `contact traction` 之后自己从 `stress`（element_data）反算（接触面相邻的
   单元 → 应力的法向分量做面积分）—— 复杂但保留分布；
3. 改用 `febio4.exe` 之外的 xplt 解析器（`mesh_dict["surfaces"]` 需手动补 set_id）。

⚠ 在未解决解析器问题之前，**禁止加 `<var type="contact traction"/>`**，
否则 `.xplt` 整文件报废。

---

## 6. pyfebio API 形态 + contact 类清单

### 6.1 contact 类清单（**proven via runtime probe**）

运行：

```
$env:PYTHONPATH="src"
& "D:\Project\climbing_fall_analysis\.venv\Scripts\python.exe" -c "
import pyfebio.contact as C
for n in sorted(dir(C)):
    o = getattr(C, n)
    if isinstance(o, type) and n[0].isupper() and not n.startswith('_'):
        print(n, '->', list(getattr(o, 'model_fields', {}).keys()))"
```

→ 输出（**proven** 探针即 .venv 当前版本）：

```
BaseXmlModel           → []
Contact                → ['all_contact_interfaces']                              # 容器类
ContactPotential       → ['type','name','surface_pair','kc','p','check_intersections',
                           'R_in','R_out','R0_min','w_tol','integration_rule']
Sliding2               → ['name','surface_pair','laugon','two_pass','penalty',
                           'auto_penalty','update_penalty','tolerance','gaptol',
                           'minaug','maxaug','search_tol','search_radius','knmult',
                           'seg_up','node_reloc','type','ptol','pressure_penalty',
                           'symmetric_stiffness','smooth_aug','dual_proj']
SlidingBase            → ['name','surface_pair','laugon','two_pass','penalty',
                           'auto_penalty','update_penalty','tolerance','gaptol',
                           'minaug','maxaug','search_tol','search_radius','knmult',
                           'seg_up','node_reloc']
SlidingBiphasic        → […,'type','ptol','pressure_penalty','symmetric_stiffness',
                           'fric_coeff','contact_frac','smooth_aug','smooth_fls',
                           'flip_primary','flip_secondary',
                           'shell_bottom_primary','shell_bottom_secondary']
SlidingElastic         → […,'type','symmetric_stiffness','smooth_aug','tension',
                           'fric_coeff','flip_primary','flip_secondary',
                           'shell_bottom_primary','shell_bottom_secondary','offset']
SlidingFacetOnFacet    → […,'type','smooth_aug']
SlidingNodeOnFacet     → […,'type','fric_coeff','fric_penalty','ktmult']
TiedBase               → ['name','surface_pair','laugon','tolerance','penalty',
                           'minaug','maxaug']
TiedBiphasic           → […,'type','gaptol','ptol','auto_penalty','update_penalty',
                           'two_pass','search_tol','search_radius','knmult',
                           'pressure_penalty','symmetric_stiffness']
TiedElastic            → […,'type','auto_penalty','update_penalty','two_pass',
                           'knmult','search_tol','search_radius','gaptol',
                           'symmetric_stiffness','flip_primary','flip_secondary']
TiedFacetOnFacet       → […,'type','search_tolerance','gap_offset']
TiedNodeOnFacet        → […,'type','search_tolerance','offset_shells',
                           'max_distance','special','node_reloc']
```

> `xml_tag` 在探针里都返回 `None`（pyfebio 用 pydantic 字段名做 tag，不挂元数据）。
> 实际生成时由 `m.save()` 按字段名输出。

### 6.2 API 调用形态（**proven in repo**）

`temp/pyfebio_demo/gen_load_xml.py:45-50`：

```python
m.contact_.all_contact_interfaces.append(
    C.SlidingElastic(name="tib", surface_pair="topsurf,topsurf",
                     laugon="PENALTY", penalty=1.0, tolerance=0.1))
```

⚠ **demo 与生产的 divergence（务必看清）**：

| 维度 | demo (`gen_load_xml.py:48`) | 生产 (`build_feb.py:606-610,676`) |
|---|---|---|
| `surface_pair` 取值 | **字符串** `"topsurf,topsurf"`（逗号拼接 two surface 名） | **字符串** `"tibiotalar_pair"`（= `<SurfacePair name>`） |
| `<SurfacePair>` 元素 | demo 没写 | 生产先写 `<SurfacePair><primary>…</primary><secondary>…</secondary></SurfacePair>`，再用其名 |

**结论**：生产写法（引用 `<SurfacePair>` 的 `name`）是 **proven in repo**；
demo 的逗号拼接写法是 **demo-only / unverified**——别抄 demo。

### 6.3 怎么选 contact 类

仓库实际只用到 `SlidingElastic` 与 `TiedElastic`（`build_feb.py:663, 677`）。
其它类（`SlidingBiphasic` / `ContactPotential` / `SlidingFacetOnFacet` 等）**本仓库未用过**，
S5 接触（足底-垫）走 **SlidingElastic** 即可，**TiedElastic** 仅在骨-软骨加密网格拼接处用。

---

## 7. Material 块语法（hyper + raw XML）

### 7.1 pyfebio 高层：`CoupledMooneyRivlin`（**demo-only / unverified**）

`temp/pyfebio_demo/demo.py:56-67`：

```python
# 材料: CoupledMooneyRivlin 软骨
my_material = pyfebio.material.CoupledMooneyRivlin(
    id=1,
    name="cartilage",
    c1=pyfebio.material.MaterialParameter(text=10.0),
    c2=pyfebio.material.MaterialParameter(text=1.0),
    k=pyfebio.material.MaterialParameter(text=1000.0),
)
my_model.meshdomains_.add_solid_domain(
    pyfebio.meshdomains.SolidDomain(name="box", mat="cartilage")
)
my_model.material_.add_material(my_material)
```

> 注：`demo.py` 只生成 `.feb` 不解算（除非带 `--run`），所以 `CoupledMooneyRivlin`
> 是否在 FEBio 4.13 上**真的收敛**是 **unverified**。生产模型（`build_feb.py:74-88`）
> 用的是 `isotropic elastic`，**未用过任何 hyperelastic 材料**——S5 的"垫"是
> hyperfoam / ogden / 自写 XML，**没有任何 proven 先例**。

### 7.2 原生 raw XML `<Material>` 形状（**proven in repo**）

`scripts/ankle_fe/build_feb.py:561-569`（紧接 `<Globals>` 之后）：

```python
w('  <Material>')
mid = {}
for i, (name, (mt, pr)) in enumerate(MATERIALS.items(), 1):
    mid[name] = i
    w(f'    <material id="{i}" name="{name}" type="{mt}">')
    for k, v in pr.items():
        w(f'      <{k}>{v}</{k}>')
    w('    </material>')
w('  </Material>')
```

→ 实际写出的 XML（以 `mat_cort_tibia` 为例，`build_feb.py:74-78`）：

```xml
<Material>
  <material id="1" name="mat_cort_tibia" type="isotropic elastic">
    <density>2e-09</density>
    <E>18000.0</E>
    <v>0.3</v>
  </material>
  <material id="2" name="mat_spon_long" type="isotropic elastic">
    <density>8.615e-10</density>
    <E>160.0</E>
    <v>0.45</v>
  </material>
  …
</Material>
```

### 7.3 给 S5 的 hyper 材料怎么落

S5 计划里"垫"是 σ(ξ) = s_pl·(1−e^(−ξ/ξ_pl)) + s_d·(ξ/ξ_d)^n + E_bot·max(ξ−ξ_d, 0)
（`docs/S5_contact_contract.md:36-46`），含"平台+压实"两段。仓库 **未尝试过**
任何 hyperelastic 拟合——S5 的 T1 `pad_foam.py` 必须**先做最小冒烟**：
用最简单的 `neo-Hookean`（pyfebio 的 `NeoHookean`）跑通 single-element 受压，
再换 `hyperfoam` / 手写 XML 拟合曲线。

---

## 8. Minimal contact-deck skeleton（组合 §1–§7，**未在本仓库运行**）

> **状态**：本骨架是从 §1–§7 的 proven 片段**拼接**出的"参考实现"，**未在仓库跑通**。
> 写 `pad_contact.py` 时请把它当模板，**首次跑前必须先做最小冒烟**。

### 8.1 目标

1 个足部（cartilage-on-bone 的 hex 块）+ 1 个垫（hyperfoam 的 hex 块）+ 1 对
sliding-elastic contact（FEBio 4.13 单步求解）。

### 8.2 脚手架（伪 Python → raw XML，可直接用 §6 的 pyfebio 包装替换）

```python
# src/climbing/coupling/pad_contact.py  （示意，未在本仓库运行）
from pathlib import Path
import os, math, numpy as np

# ---- 0. 单位：mm - N - MPa - s ----
WORK = Path(r"D:\Project\climbing_fall_analysis\temp\opensim_fe\pad_contact\pad_contact.feb")
WORK.parent.mkdir(parents=True, exist_ok=True)

# ---- 1. 网格（占位：用 meshing.vtp_to_tet / pad_mesh.py 输出） ----
# 期望 schema：{"nodes": (Nn,3), "tets": (Ne,4), "surfaces": {"foot_bot":[(3,)[…]],
#                                          "pad_top":[(3,)[…]]}, "node_sets":{...},
#                "volume_mm3":..., "min_jacobian":..., "jacobian_frac_lt_0_3":...}
foot = ...   # §3.2 pad_mesh.py 产物（占位）
pad  = ...

# ---- 2. 节点全局编号拼装 ----
nodes_all = np.vstack([foot["nodes"], pad["nodes"]])
N_foot = len(foot["nodes"])
N_pad  = len(pad["nodes"])

# ---- 3. 真实 facet 取法（§4.3 算法） ----
def real_tri_facets(domain_elem_ids: np.ndarray, must_in: set[int]) -> list[list[int]]:
    """domain_elem_ids: (Ne, 4) tet4 全局节点号；must_in: 该面的所有节点必须 ∈ 此集"""
    TET_F = [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]
    out, seen = [], set()
    for ids in domain_elem_ids:
        for f in TET_F:
            fv = [int(ids[z]) for z in f]
            if not all(x in must_in for x in fv): continue
            k = tuple(sorted(fv))
            if k in seen: continue
            seen.add(k)
            out.append(fv)
    return out

foot_bot = real_tri_facets(foot["tets"], set(foot["node_sets"]["bot_nodes"]))
pad_top  = real_tri_facets(pad["tets"],  set(pad["node_sets"]["top_nodes"]))

# ---- 4. 写 .feb ----
w = (WORK).write_text  # 一次性写完，下面 demo 用 list 拼接
L = []
a = L.append
a('<?xml version="1.0" encoding="ISO-8859-1"?>')
a('<febio_spec version="4.0">')
a('  <Module type="solid"/>')
a('  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>')

# Material：脚部 = iso-elastic 占位，垫 = ★hyperfoam 待 T1 落定
a('  <Material>')
a('    <material id="1" name="mat_foot" type="isotropic elastic">'
  '<density>1e-09</density><E>1000.0</E><v>0.3</v></material>')
a('    <material id="2" name="mat_pad"  type="hyperfoam">'
  '<density>5e-10</density><mu1>0.05</mu1><alpha1>10</alpha1>'
  '<mu2>0.0</mu2><alpha2>1</alpha2><k>0.1</k></material>')  # 占位参数，**未 proven**
a('  </Material>')

# Mesh
a('  <Mesh>')
a('    <Nodes name="all">')
for i, p in enumerate(nodes_all, 1):
    a(f'      <node id="{i}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
a('    </Nodes>')

# 元素块（**块名要唯一**：`build_feb.py:576-578` 踩坑）
a(f'    <Elements type="tet4" name="foot__tet4">')
for e, ids in enumerate(foot["tets"], 1):
    a(f'      <elem id="{e}">{",".join(str(int(x)) for x in ids)}</elem>')
a('    </Elements>')
a(f'    <Elements type="tet4" name="pad__tet4">')
for e, ids in enumerate(pad["tets"], 1):
    eid = e + len(foot["tets"])
    # 注意：pad 的节点要偏移 N_foot
    shifted = [int(x) + N_foot for x in ids]
    a(f'      <elem id="{eid}">{",".join(str(x) for x in shifted)}</elem>')
a('    </Elements>')

# NodeSet（**文本就是节点列表**，不要套 <node>：`build_feb.py:582`）
ns_bot_foot = foot["node_sets"]["bot_nodes"]
ns_bot_pad  = [x + N_foot for x in pad["node_sets"]["bot_nodes"]]
a(f'    <NodeSet name="foot_bot_ns">'
  f'{",".join(str(x) for x in sorted(ns_bot_foot))}</NodeSet>')
a(f'    <NodeSet name="pad_bot_ns">'
  f'{",".join(str(x) for x in sorted(ns_bot_pad))}</NodeSet>')

# Surface：子标签按单元类型命名 **tri3**（`build_feb.py:586-588`）
a('    <Surface name="foot_bot">')
for n, fv in enumerate(foot_bot, 1):
    a(f'      <tri3 id="{n}">{",".join(str(x) for x in fv)}</tri3>')
a('    </Surface>')
a('    <Surface name="pad_top">')
for n, fv in enumerate(pad_top, 1):
    # 面节点也要整体偏移 N_foot
    fv = [x + N_foot for x in fv]
    a(f'      <tri3 id="{n}">{",".join(str(x) for x in fv)}</tri3>')
a('    </Surface>')

# SurfacePair：**必须在 <Surface> 之后**（`build_feb.py:600, §3.2`）
a('    <SurfacePair name="foot_pad_pair">')
a('      <primary>foot_bot</primary>')
a('      <secondary>pad_top</secondary>')
a('    </SurfacePair>')
a('  </Mesh>')

# MeshDomains：**mat 填材料名而非数字 id**（`build_feb.py:616-618`）
a('  <MeshDomains>')
a('    <SolidDomain name="foot__tet4" mat="mat_foot"/>')
a('    <SolidDomain name="pad__tet4"  mat="mat_pad"/>')
a('  </MeshDomains>')

# LoadData + Loads + Boundary + Contact — **段序固定**（`build_feb.py:620-621, §3.1`）
a('  <LoadData>')
a('    <load_controller id="1" type="loadcurve">'
  '<interpolate>LINEAR</interpolate><extend>CONSTANT</extend>'
  '<points><pt>0.0,0.0</pt><pt>1.0,1.0</pt></points></load_controller>')
a('  </LoadData>')
a('  <Loads>')
# 顶面（足部）压缩位移/力：略，按 build_feb.py:628-636 的 <surface_load type="traction"> 写
a('    <surface_load type="traction" surface="pad_bot">'
  '<traction>0,0,-0.1</traction></surface_load>')   # 占位方向
a('  </Loads>')
a('  <Boundary>')
# 垫底面固定：参照 build_feb.py:638-641
a('    <bc name="fix_pad" type="zero displacement" node_set="pad_bot_ns">'
  '<x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof></bc>')
a('  </Boundary>')
a('  <Contact>')
# sliding-elastic：**逐字抄 §1**，contact 的 surface_pair 引用的是 SurfacePair 名
a('    <contact name="foot_pad" surface_pair="foot_pad_pair" type="sliding-elastic">')
a('      <laugon>AUGLAG</laugon>')               # ★ 与脚踝默认 PENALTY 不同；接触对足底-垫
                                                 #   标准做法是 AUGLAG（`build_feb.py:32-33` 注释）
a('      <penalty>1</penalty>')
a('      <auto_penalty>0</auto_penalty>')         # AUGLAG 时 auto_penalty 必须 0（联动规则见 §1.1）
a('      <two_pass>0</two_pass>')
a('      <node_reloc>1</node_reloc>')             # 垫被压实时把节点投回从面（`build_feb.py:33` 注释）
a('      <symmetric_stiffness>0</symmetric_stiffness>')
a('      <tolerance>0.02</tolerance>')
# ★ 摩擦 μ = 0.6 等参数在本仓库**未 proven**：见 §1.2
a('    </contact>')
a('  </Contact>')

# Output：**禁** <var type="contact traction"/> —— `build_feb.py:688-693`
a('  <Output>')
a('    <plotfile type="febio"><var type="displacement"/>'
  '<var type="stress"/><var type="reaction forces"/></plotfile>')
a('  </Output>')

# Control：**顶层**单步简写（`build_feb.py:694-703`）
a('  <Control>')
a('    <analysis>STATIC</analysis><time_steps>50</time_steps>'
  '<step_size>0.02</step_size>')
a('    <time_stepper type="default"><max_retries>10</max_retries>'
  '<opt_iter>15</opt_iter><cutback>0.25</cutback></time_stepper>')
a('    <solver type="solid"><linear_solver type="pardiso"/>'
  '<max_refs>25</max_refs><reform_each_time_step>1</reform_each_time_step></solver>')
a('  </Control>')
a('</febio_spec>')

WORK.write_text("\n".join(L), encoding="latin-1")
print("wrote", WORK)
```

### 8.3 骨架与 §1 的两处关键差异（务必注意）

1. `laugon="AUGLAG"` 而非仓库默认的 `PENALTY`：足底-垫是"标准接触对 + 摩擦"，
   仓库注释 `build_feb.py:32-33` 明确写 "PENALTY 收敛差, AUGLAG 是标准解法"。
   仓库默认 `PENALTY` 仅因为 **脚踝**是关节对、约束了刚性模态。脚底-垫是**真正
   的接触**，用 `AUGLAG`。**auto_penalty 联动规则**保持（`laugon=="PENALTY"` 才 1，
   否则 0）。
2. `node_reloc=1`：仓库默认 `0`，但 `build_feb.py:32-33` 注释指出"垫被压实时
   关键"——S5 必须打开。

---

## 9. 反向映射：从 §1–§7 到仓库源码

| § | 文件 | 行 | 摘要 |
|---|---|---|---|
| §1 sliding | `scripts/ankle_fe/build_feb.py` | 674–685 | `<contact type="sliding-elastic">` 写块 |
| §1 LAUGON/penalty env | `scripts/ankle_fe/build_feb.py` | 32–37 | 接触参数环境变量 |
| §2 tied | `scripts/ankle_fe/build_feb.py` | 662–673 | `<contact type="tied-elastic">` 写块 |
| §2 tied 取舍 | `scripts/ankle_fe/build_feb.py` | 335–339, 654 | 注释 |
| §3 段序 | `scripts/ankle_fe/build_feb.py` | 600, 620–621 | 注释 + 顺序 |
| §3 Control 顶层 | `scripts/ankle_fe/build_feb.py` | 694–703 | 注释 + 写块 |
| §4 真实 facet (1) | `scripts/ankle_fe/build_feb.py` | 300–306 | `cart_outer` docstring |
| §4 真实 facet (2) | `scripts/ankle_fe/build_feb.py` | 348–352 | tied-elastic 段 |
| §4 Surface 子标签 | `scripts/ankle_fe/build_feb.py` | 586–588 | 注释 |
| §4 facet_ok | `scripts/ankle_fe/build_feb.py` | 100–108 | 函数（3 distinct + area>MIN） |
| §5 plotfile | `scripts/ankle_fe/build_feb.py` | 687–693 | `<Output>` + face data 注释 |
| §6 pyfebio API | `temp/pyfebio_demo/gen_load_xml.py` | 4, 8–12, 45–50 | import + 探针 |
| §6 contact 类清单 | runtime probe | — | `dir(pyfebio.contact)` 见 §6.1 |
| §7 hyper demo | `temp/pyfebio_demo/demo.py` | 56–67 | `CoupledMooneyRivlin` 构造 |
| §7 raw XML | `scripts/ankle_fe/build_feb.py` | 561–569 | `<Material>` 写块 |
| §7 iso elastic 例 | `scripts/ankle_fe/build_feb.py` | 74–88 | `MATERIALS` 表 |

---

## 10. 已知 gap（**NOT FOUND IN REPO**）

下列项**在本仓库任何位置都未出现过**——S5 接触实现时必须自己定，且要先做最小冒烟：

1. **`<fric_coeff>` 在生产 .feb 中的可接受写法**——`build_feb.py` 里的
   `sliding-elastic` 没有任何 friction。仓库内**无 friction 的 proven 先例**。
2. **hyperelastic / hyperfoam / ogden 任何一种**在 FEBio 4.13 + 本仓库网格上的
   proven 收敛曲线。`CoupledMooneyRivlin` 只见 demo，未实际求解。
3. **接触压力 / 接触合力后处理**——`<var type="contact traction"/>` 已被 §5
   显式禁用；其它读法（element_data stress 反算、parse_state 补丁）**均未实现**。
4. **rigid wall**（垫 = 解析刚性面，§3.1 最小版路径）——`build_feb.py` 全是
   solid-solid，**无 rigid wall / rigid body** 配方。`pyfebio.material.RigidBody`
   存在于目录列表里（§6.1 探针），但**本仓库无任何引用**。
5. **pad mesh generator** `pad_mesh.py`（T2 输出 schema）：见契约 §3.2，
   本仓库**仅有契约签名，无实现**。
6. **多步 `<Step>`**（时间步分段加 contact）：`build_feb.py:694` 注释明确
   "`Control` 是顶层元素（单步简写），不包在 `<Step>` 里"。**多步+接触的
   段序**是 **NOT FOUND IN REPO**。

---

> **读毕即用规则**：写 `pad_contact.py` 之前先对照 §9 行号把 §1–§7 的
> **逐字片段**贴到 IDE，做一次最小冒烟（单 tet + 单 sliding-elastic +
> 顶面 traction / 底面 fixed），**不接触 hyper 也不接触 μ**。冒烟通过后再
> 按 §8 扩到两域 + 垫。
