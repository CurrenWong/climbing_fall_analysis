# 附页 · 跖面接触不收敛 —— **穷举排除与两个具体缺陷**（2026-10-16 续）

> **性质**：`V1_UNDERCONSTRAINT_METHODS_SURVEY.md` 的续篇（用户令「继续尝试」）。
> **本件新增**：① **手册罚因子的量级公式与我们差 5 万倍**；② **对偶板被完全固定（`pad_e_mpa` 形同无效）**；
> ③ **穷举排除**：失败对 8 个接触/求解参数**全部免疫**。

---

## §0 · 三条新结论

| # | 发现 | 强度 |
|---|---|---|
| **1** | **接触罚因子比手册公式小约 5 万倍**（`ε=0.1` vs `ε=E·A/V ≈ 5000`），且 `auto_penalty=1` 时 `penalty` **只是缩放因子** | 硬（手册逐字公式） |
| **2** | **对偶薄板的上、下两排节点**全部**被固定** ⇒ 板是**刚体块**，`pad_e_mpa`（含 0.185 MPa 跟垫）**形同无效** | 硬（代码 + 实测逐位相同） |
| **3** | **失败对 penalty / gap / 载荷 / 板模量 / 板固定 / 分析类型 / auto_penalty 全部免疫** ⇒ **不是接触参数问题** | 硬（穷举） |

**⇒ 三条合起来的推论**：**接触在本配置下没有提供有效约束**，而**唯一能改变结果的是加 `hold_spring_k`**。

---

## §1 · ★ 罚因子的量级（手册 §3.14 comment 2，逐字）

> 「an algorithm has been implemented in FEBio that attempts to calculate a good initial
> value for the penalty factor ε for a particular node/integration point *i*:
> **ε_i = E·A / V**.
> Here, **A** is the area of the element the integration point belongs to, **V** is the element
> volume and **E** is a measure of the elasticity modulus…
> When the **auto-penalty flag is on, the value of the penalty parameter serves as a
> **scale factor** for the automatically-calculated penalty factor.**」

**量级核对**：皮质 `E = 15000 MPa`；单元尺度 ~3 mm ⇒ `A/V ≈ 1/h` ⇒ **`ε ≈ E/h ≈ 5000`**。

| | 值 | 与正解之比 |
|---|---:|---:|
| 手册正解（auto 计算） | **≈ 5000** | 1× |
| **本项目历史配方** | **0.1** | **1/50000** |
| `auto_penalty=1` + `penalty=0.1` | **≈ 500**（0.1 是缩放因子） | 1/10 |

**⇒ `penalty = 0.1`（`auto_penalty=0`）比正解小约 5 万倍** —— 按手册，罚因子直接乘穿透量给出接触力
（comment 1：「the contact traction is determined by the gap (i.e. penetration distance)
**multiplied by the user-defined penalty factor**」）⇒ **接触力 ≈ 0 ⇒ 表面自由穿透 ⇒ 等价无支撑。**

**⚠️ 但本轮实测：把 penalty 提到 5000（= 正解量级）仍失败，签名一字不变 ⇒ 罚因子不是（唯一的）原因**（§3）。

---

## §2 · ★ 对偶板被「完全固定」—— `pad_e_mpa` 形同无效

`src/climbing/coupling/plantar_bc.py::_build_pad_floor`（**改前**）：

```python
node_ids = sorted(top_ids.values()) + sorted(bot_ids.values())
model.mesh_.add_node_set(NodeSet(name="pad_floor_fixed", ...))
```

⇒ **上表面（`wall_y`）与下表面（`bot_y`）两排节点全部进入固定集**（deck 实测
`pad_floor_fixed` = **全部 50 个节点**）⇒ **板是刚体块，不能变形**
⇒ **`pad_e_mpa` 无论取 200000（数值刚性）还是 0.185（跟垫）都不改变任何东西。**

**实测佐证**：`pad_e = 0.185` 与 `200000` 给出**逐位相同**的失败签名
（`neg_first=8`, 6 次重试, 末 t=0.1667, ~0.6 s）。

**⇒ 处置（已做，opt-in）**：新增 `pad_fix_bottom_only`（默认 `False` ⇒ **原路径逐位不变**）——
只固定**底面** ⇒ 板成为「**可变形层 + 刚性地基**」= 真正的跟垫层。
plumb 到 `apply_plantar_bc(contact_pad_fix_bottom_only=)` 与
`build_thums_feb(contact_pad_fix_bottom_only=)`。**实测 deck 固定节点数 50 → 25 ✓**
（**但失败签名仍一字不变** —— 见 §3）。

---

## §3 · ★★ 穷举排除：失败对**接触参数**全部免疫

h=5 m，anat 网格，`hold_spring_k=None`（无旋钮），STATIC。

| 变量 | 扫过的范围 | 失败签名 |
|---|---|---|
| **penalty** | 0.1 / 1 / 100 / 1000 / **5000** | **全部 `neg_first=8`，步 6，末 t=0.1667** |
| **auto_penalty** | 0 / 1（× penalty 0.1/1/10） | 同上 |
| **gap** | −0.05 / −0.005 / 0.0 / +0.001 / +0.01 / +0.05 | **同上** |
| **外加载荷** | **×1e-6（0.026 N！）/ ×1e-3 / ×0.1 / ×1** | **同上**（×1e-6 首现 631 个，仍失败） |
| **板模量** | 0.185 / 5 / 200000 MPa | **同上** |
| **板固定** | 全 50 节点 / 只底面 25 节点 | **同上** |
| **分析类型** | STATIC / DYNAMIC(+`<solver>` 阻尼) | STATIC 秒退；DYNAMIC 变成**刚体飞走** |
| **hold 弹簧** | `None` vs **100** | **★ 唯一改变结果的变量**（`contact_hold100` 全高度收敛） |

**⇒ 无歧义：这不是接触参数的调参问题。**

**⇒ 而且「载荷缩到 0.026 N 仍失败」说明：失败与外加载荷无关** ⇒ 不是力学问题。

### 参考构型是干净的（排除「网格本身坏」）

直接读 deck 的参考坐标算 `det J`：

| 域 | 单元 | `det J` | 倒置 |
|---|---|---:|---:|
| `calcaneus`（骨皮质） | 634 hex8 | min **1.539** / 中位 3.467 | **0/634** |
| `pad_floor`（对偶板） | 16 hex8 | min/中位 **1279** | **0/16** |

**⇒ 网格参考构型无倒置。** 而日志显示 `Total number of equilibrium iterations = 0`
⇒ 错误发生在**第一次 reform 之后、迭代之前**。

### 已排除的其它嫌疑

- **接触法向**：`pad_top` 面法向 `(0, +1, 0)`（朝骨 ✓，正向比例 **100%**）；
  `plantar` 法向 Y 中位 **−0.796**（朝板 ✓）⇒ **不是 B5「接触法向」坑**。
- **`subtalar_joint` 238 invalid facets**：该警告在**能收敛的 `fixed` deck 里同样出现**
  ⇒ **既有/无害**（不是本缺陷的来源）。
- **`min_residual = 1e-20`**（默认）；`rtol = 0`、`dtol = 0.001`、`etol = 0.01` —— 未见异常。

---

## §4 · 现在的状态（诚实）

**「`hold_spring_k=None` 时跖面接触不收敛」是一个已被穷举到接触参数之外的问题。**
按手册 §8.12.4 的诊断，它的性质仍是 **under-constrained**（force-driven contact，
接触界面本应定义良约束却未能提供）。

**⇒ 可选路径只剩两条**

| 路径 | 说明 |
|---|---|
| **(i) 接受 `hold_spring_k` 是必要约束** | 但它同时是自由旋钮（`σ` 随 `k` 单调趋向 `fixed`）⇒ **只能当诊断口径**；若要交付，必须给 `k` 一个**物理依据**（否则违反「不调参凑靶」） |
| **(ii) 换界面类型：`tied` / `sticky` 首步** | 手册 §3.14.9 `sticky`：「similar to a tied interface except that it allows **initial separation** of the tied surfaces and **breaking of the tie** after a user-defined normal traction is exceeded」⇒ **t=0 强 tie（消刚体模态）+ 超限自动脱开** = **用物理界面替代旋钮** | 
| **(iii) 给骨一个非旋钮的最小约束** | 手册 §8.3.1.5「eliminate it by **applying additional constraints**」；FEBio 的 `<rigid>` BC 或把骨的**转动**锁住（平移仍由接触平衡） |

**本件未做 (ii)/(iii) —— 它们是下一轮的第一顺位。**

---

## §5 · 诚实边界

1. 全部 h=5 m / a=0 / anat / `PressureLoad` / STATIC（或 DYNAMIC+阻尼）口径的**诊断跑**，
   **无一条构成可交付口径**。
2. **§1 的「5 万倍」是按 `A/V ≈ 1/h` 的量级估算**，未逐积分点算 `E·A/V`；**量级结论成立**。
3. **对偶板 `bottom_only` 的改动未产生行为差异**（§2 末）⇒ 它修的是**表里不一的实现**
   （`pad_e_mpa` 本该有效），**不是**收敛性的解药。
4. 未尝试：`contact_type` 换 `tied`/`sticky`、`rigid` BC、`<Control>` 级 `symmetric_stiffness` 的
   pyfebio 正规接口、`two_pass=1`、`node_reloc=1`。
5. **不动任何判定**：V1 仍为 ❌【口径不可复现】。

---

## §6 · 溯源

| 项 | 位置 |
|---|---|
| 罚因子公式 | FEBio 手册 **§3.14 comment 2**（`temp/thums/febio_manual.txt` L5063-5081） |
| 罚方法/AL 的取舍 | 手册 §8.5.1 / §8.5.2；§3.14 comment 4 |
| sticky 界面 | 手册 §3.14.9 |
| 板固定集（缺陷 2） | `src/climbing/coupling/plantar_bc.py::_build_pad_floor`（`node_ids = sorted(top)+sorted(bot)`） |
| 本轮改动（opt-in） | `plantar_bc.py`（`pad_fix_bottom_only`、`contact_penalty`、`contact_pad_thickness_mm`、`contact_auto_penalty`）＋ `thums_feb.py` 同名透传 |
| 穷举数据 | `results/opensim_fe/v1_penalty_scan.json`、`v1_load_scale_scan.json`、`v1_gap_window_scan.json`、`v1_pad_fix_bottom_only_scan.json`、`v1_pad_material_scan.json` |
| 参考构型 det J | `temp/v1_padjac.py` |
| 接触法向 | `temp/v1_normal.py` |
| 前篇 | `results/opensim_fe/V1_UNDERCONSTRAINT_METHODS_SURVEY.md` |
| 零回归 | `pytest -k "plantar or contact or thums or bone_feb"` ⇒ **72 passed** |
