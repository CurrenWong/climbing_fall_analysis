# S5 — 波3 (minimal plantar contact) plan & status

> **Scope** — 波3 最小版：把跖面 BC 升级为 **刚性壁 + Coulomb μ** 的真实 FEBio 接触。
> 按 `docs\S5接触方案.md` §3.1 + `docs\S5_contact_contract.md` §4 条 8 = 不在 G7 阻塞内。
> **不引入垫子实体**，故**不依赖泡沫本构标定**。

---

## 1. 插入点（已逐行复核）

| 位置 | 现状 | 计划改动 |
|---|---|---|
| `src\climbing\coupling\plantar_bc.py:36` | `BC_KINDS = ("fixed", "roller", "roller_free", "spring")` | 追加 `"contact"` |
| `src\climbing\coupling\plantar_bc.py:204-268` `apply_plantar_bc` 调度 | flat `if kind == ...` 链；`fixed` 早返；`roller` → `_apply_roller`；`roller_free` → `_apply_roller_free`；spring 兜底 → `_apply_spring` | 在 spring 兜底之前插入 `if kind == "contact":` → `_apply_contact(...)` |
| `src\climbing\coupling\plantar_bc.py` 新增 `_apply_contact(...)` | n/a | 仿 `_apply_roller` (line 73-106) 写 `Surface(name="plantar")` + `SurfacePair(plantar ↔ pad_wall)`；仿 `pad_contact._emit_contact` 写 `SlidingElastic` 配方 |
| `src\climbing\coupling\febio_model.py:315` `build_calcaneus_feb` 签名 | `plantar_bc: str = "fixed"`（**默认逐位保留**） | **不改签名**；仅改 docstring 列出 `"contact"` |
| `src\climbing\coupling\febio_model.py:479-489` 调用点 | `kind=plantar_bc, spring_k=spring_k` | 自动获得新 kind（不动） |
| `scripts\opensim_fe\thums_feb.py:603` `build_thums_feb` 签名 | `plantar_bc: str = "fixed"`（**默认逐位保留**） | **不改签名**；仅改 docstring 列出 `"contact"` |
| `scripts\opensim_fe\thums_feb.py:834-844` 调用点 | `kind=plantar_bc, spring_k=spring_k` | 自动获得新 kind（不动） |
| `scripts\opensim_fe\thums_feb_run.py:89-98` argparse | 无 `--plantar-bc` | **新增** `--plantar-bc {fixed,roller,roller_free,spring,contact}` 透传到 `build_thums_feb` |
| `tests\test_pad_contact.py` | 现有 9 个用例不涉及 plantar_bc | **新增** TestPlantarContact 套件：验证 `"contact"` 被接受、默认逐位保留、deck 含 contact 块 |

---

## 2. 设计选择

### 2.1 `plantar_bc="contact"` 的 XML 输出

完全仿照 `pad_contact.py` 的两域接触写法，**但只保留跖面域 + rigid wall**：

1. **Mesh**：
   - `Surface(name="plantar")` ← 从 `mesh["surfaces"]["plantar"]`（已有）
   - `SurfacePair(name="plantar_pad_wall", primary="plantar", secondary="pad_wall")`
   - `rigid_wall(name="pad_wall")` ← 解析刚性面，平面 y = `min(plantar_y) - gap_mm`
2. **Step/Contact**：
   - `SlidingElastic(surface_pair="plantar_pad_wall", laugon="PENALTY", penalty=0.1, auto_penalty=0, two_pass=0, node_reloc=0, symmetric_stiffness=0, tolerance=0.005, search_radius=20, fric_coeff=0.6)`
3. 段序：`Mesh → MeshDomains → LoadData → Loads → Boundary → Contact → Output → Control`（同 `febio_model.py` 的现有 deck）。

### 2.2 刚性壁的几何参数

- 法向取 **+Y**（跖面朝 −Y，壁的法向指向骨，pad 在壁"下方"）
- 平面位置：`y = min(plantar_y) - pad_y_offset_mm`（默认 0.5 mm，与 `pad_contact.DEFAULT_GAP_MM` 一致）
- `<plane>0,1,0,-y_wall</plane>`（FEBio 4.x 平面方程 `ax+by+cz+d=0`）+ `<offset lc="1">0.0</offset>`

### 2.3 刚性壁 XML 注入策略

pyfebio 0.3.0 **不支持** `<rigid_wall>`（核验：`pyfebio/*.py` 无 rigidwall.py），也不支持 top-level `<Contact>` 块（其 `StepEntry.contact` 是 step 内的；`Model.contact_` 是顶级冗余）。处理方案：

- 用 pyfebio 的 `SlidingElastic` 写接触块 → step 内 `<Contact>` 段；
- **跳过** `model.contact_`（与现有 `febio_model.py` / `thums_feb.py` 一致，它们都用 `fbc.Boundary` 而非 `boundary_`）；
- `<rigid_wall>` 通过**落盘后字符串注入**完成（仿 `_fix_febio413`），注入点位于 `<Mesh>` 闭合 `</Mesh>` 之前（按 FEBio 4.x §3.14.4：wall 定义在 Mesh 段）。

### 2.4 接触配方（与 udg §7 一致）

`laugon="PENALTY", penalty=0.1, auto_penalty=0, two_pass=0, node_reloc=0, symmetric_stiffness=0, tolerance=0.005, search_radius=20, fric_coeff=0.6`。
全部以 **opt-in kwarg** 暴露，默认值与现有 pad_contact 同名常量对齐（如 `pad_contact.DEFAULT_MU=0.6`、`pad_contact.DEFAULT_PENALTY=1.0` 仍然不变；本接口为 plantar 端 BC，使用 udg 配方**仅作为 opt-in**）。

### 2.5 默认行为逐位保留

`plantar_bc` 默认仍是 `"fixed"`；`BC_KINDS` 追加而非替换；
`spring_k=1000.0, anchor_offset_mm=1.0` 默认值不变。
未传 `kind="contact"` 时**不引入**任何新 mesh 元素 / surface / 节点集。

---

## 3. 测试覆盖

`tests\test_pad_contact.py` 新增 `TestPlantarContact` 套件（fast，不依赖 FEBio）：

1. `test_bc_kinds_includes_contact` — `BC_KINDS` 末位是 `"contact"`。
2. `test_apply_plantar_bc_accepts_contact` — 不抛 `ValueError`，返回 dict 含 `"kind": "contact"`。
3. `test_apply_plantar_bc_default_unchanged` — `kind="fixed"`（默认）路径仍只 append 1 个 BC。
4. `test_apply_plantar_bc_contact_emits_surface_pair_and_contact` — 验证 model 上多出 1 个 Surface、1 个 SurfacePair、Contact 块含 `sliding-elastic` 与 `<fric_coeff>`。
5. `test_apply_plantar_bc_contact_validates_plantar_tris` — 空 tris 抛 `ValueError`。
6. `test_apply_plantar_bc_contact_accepts_recipe_kwargs` — 透传 mu/penalty/symmetric_stiffness 等到 contact 块。
7. `test_build_calcaneus_feb_signature_default_unchanged` — `plantar_bc` 默认仍是 `"fixed"`（防回归）。
8. `test_build_thums_feb_signature_default_unchanged` — 同上。

---

## 4. 硬门槛 / 验收

| Gate | 目标 | 触发 |
|---|---|---|
| G0 默认回归 | `pytest tests\test_opensim_fe.py tests\test_febio_run.py -q -p no:cacheprovider` → **19 passed, 1 skipped** | 不动现有 builder 调用 |
| G0' 跖面 BC 默认 | `pytest tests\test_nonvertical_s4.py -q` → 全绿（签名 `plantar_bc.default == "fixed"`） | 不改签名默认值 |
| G2 接触最小版发射 | `pytest tests\test_pad_contact.py tests\test_pad_foam.py tests\test_pad_mesh.py tests\test_fe_post_io.py -q` → **93 passed**（+ 新增几条 test_pad_contact） | 实现 + 测试 |
| G1 (FOAM) / G7 (FOAM) | **不在本波范围** | n/a |
| G3/G4/G5/G6 | **留待**波4（FOAM 实体）+ wave3 的 `pad_contact` 集成 | 契约 §4 条 8 明确最小版不阻塞 G7 |

---

## 5. 红线（与契约 §5 同）

- **不读 / 不改** `scripts\ankle_fe\`、`temp\pyfebio_demo\`、`temp\ankle_*`、`FE_PIPELINE_TODO.md`、`LEARNINGS` LRN-010~019、现有 `results\opensim_fe\*` 产物。
- **不改** `pad_contact.py` / `docs\S5_QUOTABLE.md` / 任何 `docs\S5_g7_verify_*.md`。
- 不引 `2.205` / `15.14 cm` / `2.70×` / `2.205` / 跨 OLD/NEW 配方对比；本波不引用 σ/σ_law。
- 默认 `plantar_bc="fixed"` 一个字节不改。
- 失败必须显式抛错；不空 catch、不 type-suppression。
- 解释器 `.venv\Scripts\python.exe`；`$env:PYTHONPATH="src"`；PowerShell `;` 不用 `&&`。

---

## 6. 实施步骤（实际）

1. 写本文件（完成）。
2. 改 `src\climbing\coupling\plantar_bc.py`：BC_KINDS 追加 + 新增 `_apply_contact` + 调度分支 + 写出 rigid_wall 注入约定（用 sidecar JSON 存 wall 参数，供 builder 的 `_fix_febio413`-类后处理读）。
3. 改 `src\climbing\coupling\febio_model.py`：`build_calcaneus_feb` docstring 补 `"contact"`；`_fix_febio413` 扩成 `_post_process_feb` 兼读 sidecar 注入 `<rigid_wall>`。
4. 改 `scripts\opensim_fe\thums_feb.py`：同上（docstring + `_post_process_feb`）。
5. 改 `scripts\opensim_fe\thums_feb_run.py`：新增 `--plantar-bc {fixed,roller,roller_free,spring,contact}`。
6. 改 `tests\test_pad_contact.py`：新增 `TestPlantarContact`。
7. 跑全部门禁（详见 §4）。

---

## 7. 报告（契约 §6）

**状态：波3 最小版已实现（结构层 + 门禁全绿）。数值验收 G2–G6 属 wave4 集成，本波不判决。**

### 7.1 改了 / 新增了哪些文件（绝对路径）

| 动作 | 文件 |
|---|---|
| 改 | `D:\Project\climbing_fall_analysis\src\climbing\coupling\plantar_bc.py` |
| 改 | `D:\Project\climbing_fall_analysis\src\climbing\coupling\febio_model.py` |
| 改 | `D:\Project\climbing_fall_analysis\scripts\opensim_fe\thums_feb.py` |
| 改 | `D:\Project\climbing_fall_analysis\scripts\opensim_fe\thums_feb_run.py` |
| 增 | `D:\Project\climbing_fall_analysis\tests\test_plantar_contact.py` |
| 增 | `D:\Project\climbing_fall_analysis\docs\S5_wave3_contact.md`（本文件） |

未改：`pad_contact.py`、`docs\S5_QUOTABLE.md`、`docs\S5_g7_verify_*.md`、`docs\S5接触方案.md`、`docs\S5_contact_contract.md`（只读）、任何 `results\opensim_fe\*`、只读区。

**符号级改动**
- `plantar_bc.py`：`BC_KINDS` 追加 `"contact"`；新增常量 `DEFAULT_CONTACT_MU=0.6` / `DEFAULT_CONTACT_PENALTY=0.1` / `DEFAULT_CONTACT_SEARCH_RADIUS=20.0` / `DEFAULT_CONTACT_TOLERANCE=0.005` / `DEFAULT_CONTACT_GAP_MM=0.5` / `CONTACT_RIGID_WALL_SIDECAR_SUFFIX`；新增 `_build_pad_floor(...)`（固定薄板）、`_apply_contact(...)`；`apply_plantar_bc` 新增 8 个 keyword-only opt-in 参数（默认值与现状对齐）+ 在 spring 兜底前插入 `if kind == "contact":` 分支。
- `febio_model.py`：`apply_plantar_bc` 调用改为捕获返回值；`plantar_meta.get("contact")` 非空时 `st.contact = plantar_meta["contact"]`（**默认 fixed ⇒ None ⇒ 逐位不变**）；`build_calcaneus_feb` docstring 列 `"contact"`（**未改签名**）。
- `thums_feb.py`：同上（**未改签名**）。
- `thums_feb_run.py`：新增 `--plantar-bc {fixed,roller,roller_free,spring,contact}`（默认 `fixed`），透传到两条 build 路径。

### 7.2 跑了什么命令 / 退出码 / 关键输出

```powershell
$env:PYTHONPATH="src"
# 门禁 1（契约 §4 条 G0'，4 文件）
.venv\Scripts\python.exe -m pytest tests\test_pad_contact.py tests\test_pad_foam.py tests\test_pad_mesh.py tests\test_fe_post_io.py -q -p no:cacheprovider
# → 93 passed in 44.42s   (rc=0)

# 门禁 2（G0 默认回归）
.venv\Scripts\python.exe -m pytest tests\test_opensim_fe.py tests\test_febio_run.py -q -p no:cacheprovider
# → 19 passed, 1 skipped in 86.56s   (rc=0)

# 门禁 3（A7 默认 plantar_bc 仍 fixed）
.venv\Scripts\python.exe -m pytest tests\test_nonvertical_s4.py -q -p no:cacheprovider
# → 11 passed in 1.09s   (rc=0)

# 新增波3 结构回归
.venv\Scripts\python.exe -m pytest tests\test_plantar_contact.py -q -p no:cacheprovider
# → 11 passed in 0.93s   (rc=0)

# CLI flag 冒烟
.venv\Scripts\python.exe scripts\opensim_fe\thums_feb_run.py --help
# → --plantar-bc {fixed,roller,roller_free,spring,contact}   (rc=0)

# FEBio 语法冒烟（合成 cube + contact deck）
& "D:\Program\FEBioStudio\bin\febio4.exe" -i temp\opensim_fe\w3_smoke\calc_contact.feb
# → 解析通过（无 unrecognized tag / invalid value）；求解因合成网格几何不对齐而 negJac
```

### 7.3 G 门禁逐条判定

| Gate | 判定 | 数字 / 依据 |
|---|---|---|
| **G0 默认回归逐位不变** | ✅ PASS | 19 passed / 1 skipped；4 文件 93 passed；`build_calcaneus_feb`/`build_thums_feb` 的 `plantar_bc` 默认仍 `"fixed"`（`test_nonvertical_s4` 11 passed）；fixed 路径 `plantar_meta.get("contact") is None` ⇒ deck 逐字节不变 |
| **G0' 跖面 BC 默认** | ✅ PASS | `test_nonvertical_s4::test_build_calcaneus_signature_defaults` / `..._thums_...` 均绿 |
| **G1 σ(ξ) FE 复现** | ⛔ 不在本波 | FOAM 版（§3.2）指标；契约 §4 条 8 明确最小版不含泡沫本构 |
| **G2 极限自检（μ=0 + 极刚垫 → rigid-support）** | ⏸ 未执行（wave4） | 接口支持 `contact_mu=0.0` / `contact_pad_e_mpa=∞`-级；需收敛的骨-接触 FE 跑 |
| **G3 接触守恒（≤2%）** | ⏸ 未执行（wave3 结构层） | `pad_contact.py` 自带 G3 测试在 93-gate 内 ✅；骨-接触版待 wave4 |
| **G4 摩擦锥 `|F_t| ≤ μF_n`** | ⏸ 未执行（wave3 结构层） | `pad_contact.py` 自带 G4 测试在 93-gate 内 ✅ |
| **G5 A7 消失（contact vs fixed gauge）** | ⏸ 未执行 | 需收敛骨-接触 FE 对比 |
| **G6 敏感性（μ/penalty 单调有界）** | ⏸ 未执行 | 需 wave4 数值扫描 |
| **G7 压入量窗口** | 🚫 已退役 | 契约 §4 条 6–8：最小版不在 G7 阻塞内 |

**本波不引用任何 σ/σ_law 比值**（`docs\S5_QUOTABLE.md` 限定）。recipe 出处：`docs\S5_contact_udg.md` §7（NEW PENALTY 0.1，`search_radius=20`，`node_reloc=0`，`tolerance=0.005`）；绝对量 `F_n ≈ 34.4 kN` / `σ_contact ≈ 0.574 MPa` 仅作单位/量级锚，未在波3 复用为结论。

### 7.4 关键设计发现 + 未解决风险

**发现（实测，非猜测）：FEBio 4.13 不支持刚性壁的两种 4.0 写法。**
- `<Mesh><rigid_wall>` → `tag "rigid_wall" (line 25) : unrecognized tag`
- `<contact type="rigid_wall">` → `tag "contact" (line 42) : invalid value for attribute "type"`
- SDK 佐证：`FEBioContactSection4`（4.0）**没有** `ParseRigidWall`（只有 2.0/2.5 有）；`pyfebio 0.3.0` 也无 rigidwall 元素类。

**处置**：`_apply_contact` 用一块**全节点固定**的 `hex8` 薄板作对偶面（数值刚性壁，刚度不起作用——所有 DOF 被 BC 锁死）。这与 `pad_contact.py` 模块文档里 "stiff plate 与 rigid body 在此不可区分；本仓库无 proven rigid-body deck" 的既有结论一致，属 §3.1 允许的 **"刚性 body"** 路径。

**风险**
1. **最小版尚未做过真实骨网格的收敛 FE 跑**：G2–G6 未验证。下一步（wave4）：在真实 calcaneus/THUMS `plantar`（−Y 法向，`meshing.py:57`）上跑 `contact`，做 A7 对比与摩擦锥检查。
2. **`_apply_contact` 的薄板方向依赖跖面法向 −Y**：对非 −Y 跖面（合成测试网格）几何不对齐；真实项目网格已满足（`PLANTAR_NORMAL_Y_MAX=−0.5`）。
3. **合成 cube 冒烟不收敛（negJac）**：因合成网格跖面不在 y=const 平面，属测试几何假象；FEBio **解析/装配通过**，仅求解发散。真实网格待 wave4 验证。
4. **摩擦 + 无侧向销**：`contact` 路径仅靠 μ 约束切向；纯竖直载荷下静定，但若侧向分量超摩擦锥会滑。wave4 需确认。

### 7.5 一个可复现的最短命令

```powershell
cd D:\Project\climbing_fall_analysis; $env:PYTHONPATH="src"; .venv\Scripts\python.exe -m pytest tests\test_plantar_contact.py tests\test_nonvertical_s4.py -q -p no:cacheprovider
```
（波3 结构验收，~2s；门禁全量见 §7.2。）
