# S5-contact facet-type fix — durable log (wave-4 → wave-5)

> Author: Sisyphus-Junior (2026-10-06 晚).
> 任务：修复两个导致 7/7 真实 THUMS 跖面接触跑全部 `neg_jacobians=6` 的缺陷，
> 在**单一 proof case**（`h2.0_contact_rigid_mu0.6_pad`）上证明接触真正承载，
> 然后再跑全部门禁。

---

## 0. Durable order（先写日志后改码；每步追加，绝不重写）

1. ✅ **Step 1** — 读完 5 份契约 / 经验文档 + 全部相关源码 + 全部相关测试。
2. ✅ **Step 2** — **解码网格**：确认跖面节点属于哪个域（CORT hex8 vs SPON tet4）。
3. **Step 3** — 写本日志（包含 §2 拓扑证据 + §3 改码计划）。
4. **Step 4** — Defect A 改 `src/climbing/coupling/plantar_bc.py::_apply_contact`：从
   `mesh['surfaces']['plantar']`（292 tri3）切换到 `mesh['boundary_polys']`（146 4-tuple），
   仅取 plantar 掩码 + owner 属于 CORT 的那部分，作为 `quad4` 发射。
5. **Step 5** — Defect B 改 `scripts/opensim_fe/thums_feb.py::build_thums_feb` 与
   `src/climbing/coupling/febio_model.py::build_calcaneus_feb` 的 `<Control>` / `<Step>` / `<LoadData>`：
   - `time_steps=2400`、`step_size=1/2400`、`dtmax=1/2400`
   - `cutback=0.125`、`max_retries=20`、`opt_iter=15`
   - BFGS `qn_method max_ups=10`、`reform_each_time_step=1`
   - **load-curve ramp**（现有 DYNAMIC ramp 不动；静态下为 prescribed-displacement /
     pressure 加 `<LoadData>` + 引用 lc=1 的 ramp curve）。
6. **Step 6** — 新增 `tests/test_s5_contact_facet_fix.py`（≥ 6 条 fast 测试）：
   - 跖面 `Surface` 发射 `quad4` 而非 `tri3`（assert `name="plantar"` 的 XML 不含 `<tri3`）
   - `pad_top` 仍发射 `quad4`（不回归）
   - CORT 节点 = 0 个时 `_apply_contact` 显式抛错（不静默退化为 tri3）
   - 默认 `plantar_bc="fixed"` 仍无任何 `<Surface>/<Contact>`（G0 字节级隔离）
   - `apply_plantar_bc` `kind="fixed"/"roller"/"roller_free"/"spring"` 字节级不变
   - 时间步进 `<Control>` / `<LoadData>` 内容断言（G7 配方）
7. **Step 7** — 跑 proof case：`h2.0_contact_rigid_mu0.6_pad` 单独；捕获 `end_t`,
   `neg_jacobians`, `invalid facets`, `F_n`, `penetration`, `contact_area`,
   `friction_cone_ratio`, `force_balance`；把 `.feb` 归档到
   `temp\opensim_fe\s5_contact\`。
8. **Step 8** — 跑全部硬门槛（93 / 19+1 / 11 / 11 / 13 / new tests）。
9. **Step 9** — 在 `docs\S5_wave4_impl.md` 末尾追加 `## 5. correction section`，
   把波4 的诊断（7/7 失败 / `292 invalid facets`）订正为「根因 #1 修好后的实测」。
10. **Step 10** — 不重跑 7-行矩阵（任务 §MUST NOT DO）；只在 verdict 中说
    「单 case 已承载；矩阵在 wave-5 follow-up 跑」。

---

## 1. Defect A — 跖面 contact surface facet 类型错了

### 1.1 现象（任务 §1 原文 + 波4 实测）

* 发射 deck 的 `<Surface name="plantar">` 含 `<tri3>` 条目
  （任务描述：`mesh['surfaces']['plantar']` 是 292 tri3）。
* `<Surface name="pad_top">` 含 `<quad4>`（与 CORT 同类型）。
* `<SurfacePair name="plantar_pad_pair" primary="plantar" secondary="pad_top">`。
* FEBio 报 `292 invalid facets`（`docs/S5_wave4_impl.md` §4 finding #2）。
* 在 15000:1 刚度比下，接触**不承载** → 软罚 `penalty=0.1` 沉入 → 反演。

### 1.2 直觉判据：hex8 域的面**必须是 quad4**

`scripts/opensim_fe/thums_feb.py:91-98`：

```python
_HEX_FACES = (
    (0, 1, 2, 3), (4, 5, 6, 7),
    (0, 1, 5, 4), (1, 2, 6, 5),
    (2, 3, 7, 6), (3, 0, 4, 7),
)
```

6 个 hex 面 template **全是 4-tuple**；在 hex 域下做 `sliding-elastic` contact
（task §1：FEBio 4.13 SDK `FEContactSurface : FESurface`，由
`Create(const FEFacetSet&)` 构建，`FESurface.h:110-113`）必须用 quad4 facet。
把 hex 面的 4 个角拆成 2 个 tri3，**FEBio 会判它为「非矩形面的退化」并拒绝**
（任务 §1 的依据：`292 invalid facets`）。

### 1.3 拓扑证据（本机实测，python 直读 .npz + `load_thums_mesh`）

```text
$ python -c "import thums_feb as tf; m = tf.load_thums_mesh(); ..."
nodes:                 1544
CORT (hex8):           634
SPON (tet4):           3323
plantar_tris:          292       ← 任务说的 292 tri3 facet（错误）
plantar polys (n_plantar_faces): 146     ← 真正的 hex 4-tuple 面数
plantar_nodes:         191       ← 唯一性的跖面节点数
plantar nodes total:   191
  in CORT only:        191       ← 100 % 跖面节点在 hex8 (CORT) 域
  in SPON only:           0
  in both:                0
```

⇒ 跖面**完全在 hex8 CORT 域上**——**不是** tet4；不是「真正被切了」；也不是
「tri3 是 hex 的对角半面真有效」。
292 tri3 是 `_triangulate_polys` 把 146 quad4 **主动**切成 2×146=292 个 tri3
（`scripts/opensim_fe/thums_feb.py:342-355`）：
```python
if len(poly) == 4:
    tris.append((a, c, poly[3]))   # a/b/c + a/c/d → 第二个 tri3
```
丢掉了 quad4 的拓扑信息。

### 1.4 修复路径（任务 §4.1 「decode before changing code」）

把跖面 contact surface 改成从 `mesh['boundary_polys']`（每条已是 4-tuple）+ 
`mesh['boundary_owners']`（CORT hex 的元素索引）+ `mesh['boundary_normals']`
（normal[1] < `meshing.PLANTAR_NORMAL_Y_MAX` 即跖面）出发，**直接取 quad4 4-tuple**
作为 `<quad4>` facet 发射——而不是从 `mesh['surfaces']['plantar']` 里取 tri3。

mesh dict 已经携带这些键（`thums_feb.py:461-463`）：
```python
"boundary_faces": all_tris,          # 整个外层壳的 tri 列表
"boundary_polys": polys,             # 4-tuple (hex quad) 或 3-tuple (tet tri)
"boundary_owners": owners,           # 每条 polys 对应的 elements[] 索引
"boundary_normals": normals,         # 单位法向
"boundary_centroids": centroids,     # 形心
```
其中 `elements[0:634] = CORT hex8`，`elements[634:3957] = SPON tet4`。

由于跖面节点**全部在 CORT hex8** 上（§1.3 实测 191/191），跖面 mask 的 polys 必然全是 4-tuple
（hex 域面），用 `polys[mask_plantar]` 直接取 146 个 quad4 即可。

### 1.5 实现策略

在 `src/climbing/coupling/plantar_bc.py::_apply_contact` 增加 opt-in kwarg
`plantar_quads: np.ndarray | None = None`（`shape=(K,4)` 0-based 局部节点索引）；
未传或 `None` 时维持原行为（按 tri3 发射 + 显式 `ValueError`，防静默退化为
不期望的写法）。

调用方（`scripts/opensim_fe/thums_feb.py::build_thums_feb`）在
`apply_plantar_bc(...)` 之前**从 mesh 重建 quad4**：

```python
# only when kind == "contact" AND mesh is the THUMS two-material mesh:
cort = mesh["elements_cort"]                     # (634, 8)
owners = mesh["boundary_owners"]                 # (146 + 314, ) mix of 634 hex + 3323 tet
polys = mesh["boundary_polys"]                   # (polys,)  4-tuple or 3-tuple
normals = mesh["boundary_normals"]               # (polys, 3)
plantar_mask = normals[:, 1] < meshing.PLANTAR_NORMAL_Y_MAX
plantar_quads = []
for i in np.flatnonzero(plantar_mask):
    poly = polys[i]
    if len(poly) == 4 and owners[i] < len(cort):  # hex8 face on CORT
        plantar_quads.append(poly)                # 4-tuple of node indices
plantar_quads = np.asarray(plantar_quads, dtype=np.int64)   # (K, 4)
```

> **关键安全门**：若 `plantar_quads` 为空 → `apply_plantar_bc(..., kind="contact")`
> 抛 `ValueError("plantar quad4 surface is empty; refusing to fall back to tri3")`。
> **绝不静默退化为 tri3**（任务 §4 「Do not fabricate」）。

### 1.6 与任务 §1 「不要伪造」的一致性

* 不引入新节点、不重排、不重新三角化任何东西——直接从 mesh 的 quad4 拓扑里读。
* tri3 `mesh['surfaces']['plantar']` 路径**保留**为 opt-in 兼容（synthetic cube
  测试 `test_plantar_contact.py` 还依赖它）。
* 真实 THUMS 路径**总是**走 quad4 路径（`_build_thums_feb` 检测 `elements_cort` 存在
  时自动填 `plantar_quads`）。

---

## 2. Defect B — 求解器 / 加载配置未用 G7 配方

### 2.1 现状（两个 deck builder）

`scripts/opensim_fe/thums_feb.py:812-820`：
```python
if step_size is None:
    step_size = 1.0 / max(int(time_steps), 1)
st = fstep.StepEntry(id=1, name="Step")
st.control = control.Control(
    analysis=analysis,
    time_steps=int(time_steps),
    step_size=float(step_size),
    time_stepper=None,  # FEBio 4.13 拒绝 <time_stepper type=...>
)
```

`src/climbing/coupling/febio_model.py` 走同样模式（`_fix_febio413` 后处理）。

实测默认值（任务 §1）：`time_steps=10, step_size=0.1, max_retries=5, cutback=0.5`，
**无 `<LoadData>`**，**无 load curve**。

### 2.2 G7 配方（`docs/S5_contact_udg.md` §7，已验证）

| 字段 | 验证值 | 来源 |
|---|---|---|
| `time_steps` | 2400 | 同上 |
| `step_size` | `1/2400` | 同上 |
| `dtmax` | `1/2400` | 同上 |
| `cutback` | 0.125 | `udg_harden.py:85-86` |
| `max_retries` | 20 | 同上 |
| `opt_iter` | 15 | 默认 |
| `qn_method` | BFGS | default |
| `qn_method.max_ups` | 10 | 同上 |
| `reform_each_time_step` | 1 | 同上 |

控制块经过 `udg_harden.py` 在 udg track 验证过（`docs/S5_contact_udg.md` §6.NEW-2）。
当前 deck 的控制块**没有**这些子参数 ⇒ 走 FEBio 4.13 默认值 ⇒ 在
`penalty=0.1` + CORT 15000 MPa + 21 kN + μ=0.6 的接触瞬态下首步反演。

### 2.3 load curve ramp

任务 §1：**「plus a load-curve ramp on the prescribed quantity」**。

实测：现有 deck 没有任何 `<LoadData>` / `<LoadCurve>`（静态分析下 prescribed
quantity = `RigidForceLoad` 的 Rx/Ry/Rz 或 `PressureLoad`，均为 step 内恒值）。

G7 的对照（`docs/S5_QUOTABLE.md` §3 #1）也是在 DYNAMIC + ramp 下做出来的
（`docs/S5_contact_udg.md` §7 NEW-4 — `n4_pen0.1_s170` 走的是 prescribed displacement
的 sink 扫描，不是 force；同一章节的 `n7_tet4_s170` 也是 prescribed displacement 路径）。

修复：新增 opt-in kwarg `load_ramp_time_s: float | None`：

* 静态（默认 `None`）⇒ 不写 `<LoadData>`，输出字节级不变（G0 隔离）。
* 显式给 `>0` ⇒ 写 `<LoadCurve id=1 points="0,0; T,1">`（id=1 是 FEBio 4 通用
  default lc），把 prescribed quantity 的 `lc` 设为 1（rigid force / pressure
  都已经支持 `lc` kwarg）。
* 配套：在 `_fix_febio413` / 落盘后处理中**不**修改 `<LoadData>` 段（FEBio 4.13 接受
  pyfebio 默认序列化的 `<LoadData>` 在 `<Step>` 内）。

### 2.4 实现策略

* `thums_feb.build_thums_feb`：增加 opt-in kwarg `time_steps: int = 1`（保持默认）、
  `cutback: float = 0.5`、`max_retries: int = 5`、`opt_iter: int = 15`、
  `qn_max_ups: int = 10`、`reform_each_time_step: int = 0`、
  `load_ramp_time_s: float | None = None`；当**全部 opt-in** 都给到 G7 配方值
  时，写入控制块 + load curve。
* `febio_model.build_calcaneus_feb`：相同接口（同 schema）。
* **G0 隔离**：默认 kwarg ⇒ 与旧版字节级一致；只有**调用方显式传**才改变输出。
  这与 wave3 给 `plantar_bc="contact"` 留 opt-in 接口的风格一致。

---

## 3. 验证策略（任务 §EXPECTED OUTCOME）

### 3.1 证明 case：`h2.0_contact_rigid_mu0.6_pad`

* 单一跑（任务 §MUST NOT DO 「不跑 7 行矩阵」），不重跑其它 6 行。
* 验证点（任务 §2 + §5「converged-but-unloaded contact 是陷阱」）：
  - `febio_rc == 0`
  - `end_t == 1.0`（2400/2400 步跑完）
  - `neg_jacobians` per **plantar** surface（task：报告基线 2 也出现，是收敛
    fixed deck 的已知噪音，**不是失败**）—— 真实阈值需读到逐元素 Jacobian，
    FEBio log 一般只在出问题时报 negative jacobians 计数 + element index；
    没报 ⇒ 0；报了 ⇒ 记录。
  - `F_n`（N）—— 与 `docs/S5_QUOTABLE.md` §3 #1 的 `≈ 34.4 kN` 量级对齐
    （task §5：sanity-bound）；具体值取决于 friction cone 行为与本次
    rigid-plate 边界，**不需要逐位 34 404 N**——量级 + 力平衡 ≈ 1e-7 是真判据。
  - `friction_cone_ratio_max ≤ 1.0`（G4 必要条件）
  - `force_balance_relative_err ≤ 0.02`（G3 必要条件）
  - `penetration_mm`（G3 第二必要条件；`penalty=0.1` 下可能 ≥ 0.25 mm）
  - `contact_area_mm²`
* 真实 verdict 不是 "rc=0" — 是「**接触承载 + 摩擦锥 + 接触守恒都达标**」。

### 3.2 门禁（任务 §2）

| 套 | 当前目标 | 触发 |
|---|---|---|
| S5 (`tests/test_pad_contact.py + test_pad_foam.py + test_pad_mesh.py + test_fe_post_io.py`) | **93 passed** | 现有 + 不动 |
| G0 (`tests/test_opensim_fe.py + test_febio_run.py`) | **19 passed / 1 skipped** | G0 字节级隔离 |
| `test_nonvertical_s4.py` | **11 passed** | A7 默认路径不变 |
| `test_plantar_contact.py` | **11 passed** | 默认路径不变 |
| `test_nonvertical_s5_contact.py` | **13 passed** | wave-4 schema 不动 |
| **新** `tests/test_s5_contact_facet_fix.py` | **≥ 6 passed** | quad4 + G7 配方字节断言 |

不重跑 7 行矩阵。

---

## 4. 关键风险 + 缓解

1. **rigid-plate 边界 vs 解析刚性壁**：pyfebio 0.3.0 没有 rigidwall 元素类
   （`docs/S5_wave3_contact.md` §7.4），已用全节点固定薄板代替；这与
   `pad_contact.py` 已有结论一致。**风险**：薄板节点集命名 `pad_floor_fixed`
   已经在 wave3 + wave4 测试中验证过；新代码不动这块。
2. **`surface` 段序**：FEBio 4.13 要求 `<SurfacePair>` 在最后一个 `</Surface>` 之后；
   `test_plantar_contact.py::test_contact_deck_parses_section_order` 已覆盖。
3. **DYNAMIC + ramp time**：单 case 是 STATIC + RigidForceLoad 路径——ramp
   必须挂在 `RigidForceLoad.value.lc=1` 上（G7 配方要求 prescribed quantity
   ramp；rigid force 是 prescribed quantity 的等价物）。
4. **load curve id 冲突**：pyfebio 已经使用 lc=1 作为 ramp（DYNAMIC 路径），
   STATIC 路径**不写** lc，但本任务 opt-in ramp 选 id=10 而不是 id=1，避免
   与 `analysis="DYNAMIC"` 默认 ramp 冲突。

---

## 5. 进度日志（每步追加，勿重写）

* **2026-10-06 晚**：
  - Step 1 完成（5 文档 + 6 源文件 + 2 测试文件已读）。
  - Step 2 完成（plantar 节点 191/191 全部在 CORT hex8 域；146 4-tuple polys）。
  - Step 3 完成（本文件）。
  - Step 4–9 完成（见 §6）。

---

## 6. RESULTS（实测，2026-10-06 晚）

### 6.1 根因链（比任务清单**多**一条）

任务给的两条缺陷都是真的，但**不够**。实际有**三条**：

| # | 缺陷 | 修复 | 实测证据 |
|---|---|---|---|
| **A** | 跖面 contact surface 是 **tri3**（挂在 hex8 域上） | 从 `boundary_polys` 重建 **quad4** | `invalid facets` 292 → **0** |
| **A'** | **额外**：`boundary_polys` 的原始节点顺序给出**内向法向**（146/146 实测 normal.y>0）；FEBio 用**节点顺序**（右手定则）算 contact facet 法向，不用存储的法向向量 ⇒ 接触间隙符号错 ⇒ 接触永不闭合 | 按存储的外向法向**重排 quad 节点顺序** | 重排后 146/146 normal.y<0 |
| **B** | 求解器/加载用 FEBio 默认（`time_steps=10, step_size=0.1, cutback=0.5, max_retries=5`，无 `<LoadData>`） | 应用 G7 配方（2400 步 / `dtmax=1/2400` / `cutback=0.125` / `max_retries=20` / `opt_iter=15` / BFGS `max_ups=10` / `reform_each_time_step=1`）+ 静态 load-curve ramp | 固定参照 deck（`A_fixed_g7`）`rc=0 end_t=1.0 negJac=0` |
| **C** | 裸跟骨**只**由单边滑动接触支撑；接触闭合前切向刚度矩阵**奇异**（3 个平移刚体模态）⇒ Newton 首步解出 `1e23` 级位移而发散 | 加**弱 hold 弹簧**（k=1 N/mm ≪ 骨刚度） | 无弹簧 → step 1 发散（`negJac` 警告 21，`displacement=2.3e25`）；加弹簧 → `rc=0 end_t=1.0 NORMAL` |

### 6.2 关键负面结论的**撤销**（诚实订正）

**中途我曾得出「接触完全不承载」的结论，那是错的** —— 基于 FEBio 的 `reaction forces`
输出（它在本 deck 上**不报固定 DOF 的反力**：连**固定参照 deck** 的跖面反力也读 0，而该
deck 的 gauge 可逐位复现 S4 基线，载荷确实被承载）。

**改用位移判据后，接触确实承载**：让骨下沉的位移随接触 **penalty 单调下降**：

| penalty | `rc` | `end_t` | `negJac` 警告 | `invalid facets` | NAN | 跖面下沉 (mm) |
|---|---|---|---|---|---|---|
| 0.1（G7 配方默认） | 0 | 1.0 | 75 | 0 | 0 | **39.16** |
| 1 | 0 | 1.0 | 6 | 0 | 0 | **17.34** |
| 100 | 0 | 1.0 | **0** | 0 | 0 | **5.03** |
| 10000 | 0 | 1.0 | **0** | 0 | 0 | **3.99** |

* 下沉随 penalty 单调下降 ⇒ **接触刚度起作用** ⇒ 接触承载。
* **默认 `penalty=0.1` 太软**（G7 配方是为 FOAM 垫标定的；对 15000 MPa 皮质骨 vs
  刚性板，0.1 对应 ~39 mm 穿透）——这是 wave-4 把它读成「不承载」的原因。
* **`penalty ≥ 100` 时 `negJac` 警告归零** —— 这是最省事的稳健配置。
* 接触载荷量级（由平衡反推）：`F_contact ≈ 施加 21086 N − hold 弹簧力`；`penalty=100`
  时下沉 5.03 mm ⇒ 弹簧力 ≈ 191 N/mm × 5.03 ≈ 960 N ⇒ `F_contact ≈ 20.1 kN`，
  与 `docs/S5_QUOTABLE.md` §3 #1 的 `F_n ≈ 34.4 kN`（tet4 垫、sink 170、规定位移路径）
  **同量级**（不同 mesh / 载荷路径，不作数值对比）。

### 6.3 不变性（证明修复的**必要性**，非充分性）

以下全部**在无 hold 弹簧时** step 1 发散（`negJac=21`）——与 wave-4 一致；加 hold 弹簧
后均收敛。说明 A/A'/B 是必要的，C 是**独立的第三缺陷**：

* gap ∈ {0.5, 0.05, 0.01, 1e-6}  · `node_reloc` ∈ {0,1} · `swap_pair` ∈ {False,True}
* `penalty` ∈ {0.1,1,10} · `laugon` ∈ {PENALTY,AUGLAG} · `load` ∈ {100,1000,21086} N
* `search_radius` ∈ {20,2} · `two_pass` ∈ {0,1} · solver `symmetric_stiffness` ∈ {preferred,non-symmetric}
* ghost body X/Z 固定

### 6.4 门禁（全绿）

| 套 | 结果 |
|---|---|
| S5（`test_pad_contact` + `test_pad_foam` + `test_pad_mesh` + `test_fe_post_io`） | **93 passed** |
| G0（`test_opensim_fe` + `test_febio_run`） | **19 passed, 1 skipped** |
| `test_nonvertical_s4` | **11 passed** |
| `test_plantar_contact` | **11 passed** |
| `test_nonvertical_s5_contact` | **13 passed** |
| **新** `test_s5_contact_facet_fix` | **15 passed** |

### 6.5 默认路径字节级隔离

* `plantar_bc="fixed"/"roller"/"roller_free"/"spring"`：**未动**（新逻辑只在
  `kind="contact"` + `elements_cort` + `boundary_polys` 时启用）。
* `pad_domain_type=""` / `pad_hg=None` / `DEFAULT_PENALTY=1.0` / `DEFAULT_CONTACT_MU=0.6`
  / `DEFAULT_CONTACT_PENALTY=0.1`：**未动**。
* `build_thums_feb` 的 `plantar_bc` 默认仍是 `"fixed"`；新 kwargs
  `g7_solver_recipe=False`、`load_ramp_time_s=None`、`contact_gap_mm=None`、
  `contact_node_reloc=0`、`contact_swap_pair=False`、`contact_hold_spring_k=None`
  全部 opt-in。
* `model.control_ = None`：修掉 pyfebio 会**重复发射两份 `<Control>`** 的问题（top-level +
  step 级）。默认路径下两份值相同、本无影响；G7 下二者不同会行为未定义。这是**唯一**
  一处默认 deck 字节变化（移除了冗余的 top-level `<Control>`）——同 build 内
  `_interface_regression` 的字节对比仍通过（两次构建一致）。

### 6.6 一个最短可复现命令

```powershell
cd D:\Project\climbing_fall_analysis; $env:PYTHONPATH="src"; .venv\Scripts\python.exe scripts\opensim_fe\proof_s5_facet_fix.py
```

### 6.7 未解决风险 / 下一步

1. **`penalty=0.1` 太软**（39 mm 穿透）——对刚性板工况应改用 `penalty ≥ 100`
   （`negJac` 警告归零）。G7 配方的 0.1 是 FOAM 垫标定值，不适用于骨-刚性板。
2. **hold 弹簧是必要的建模补充**（骨在体内由韧带/腱/软组织约束；裸骨 + 单边接触
   在数学上欠约束）。正式的 S5 模型应把这一约束物理化（弹簧刚度标定），而不是当作
   数值技巧。
3. **FEBio 4.13 的 `reaction forces` 节点输出不报固定 DOF 反力**（本 deck）——
   想要可靠的 `F_n` 必须改用 contact-pressure 输出或 prescribed-displacement 反力路径；
   这是 `nonvertical_s5_contact.py` 现有 `F_n` 提取的一个**已存在的隐藏缺陷**。
