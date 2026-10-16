# 🔴 **重大口径发现**：V1 的「跖面接触」是**惰性的** —— 它逐位等于 `spring` BC

> ## ⚠️ **2026-10-16 更正（后件 `V1_PLANTAR_CONTACT_WINDING_FIX.md`）**
>
> **本件「接触本性惰性」的结论**不准确**。真机制已定位 = **两个接触面的面片绕向（法向）
> 不满足 FEBio 的 `sliding-elastic` 约定**，使 FEBio **恒算出正值（= 分离）的间隙
> ⇒ 接触永不啮合**。**把 `plantar` 与 `pad_top` 的绕向都反转后，接触生效**
> （`max` 152.49 ≠ 151.90；且所需 `hold_k` 从 **≥100 降到 ~1**）。
>
> **⇒ 本件的正确读法**：「**在当前（as-built）绕向配置下**，那个 `<Contact>` 块贡献恰好 0」。
> §1 的全部**观测**仍成立（逐位相同、参数全免疫），**只是原因换了**（从「机构惰性」改为「未啮合」）。
> §2 的**扩散影响**相应地改为：**凡用 `plantar_bc="contact"` 得到的数，都是「未啮合的空接触」，
> 实际等于对应的 `spring` 口径。**

> **发现日期**：2026-10-16（用户令「继续」）
> **性质**：**口径发现**（不是模型改进）。**不动任何判定**。
> **强度**：**逐位相同（bit-for-bit identical）**，非近似。

---

## §0 · 一句话

**`plantar_bc="contact"` + `contact_hold_spring_k = k` 的 FE 结果，
与 `plantar_bc="spring"` + `spring_k = k` 的 FE 结果**逐位相同** ⇒
那个 `<Contact>` 块**贡献恰好 0** ⇒ **项目从未有过真正生效的跖面接触**；
所有 `contact_hold*` 口径的数**实际都是 `spring` 口径的数**。**

---

## §1 · 证据

### 1.1 逐位对照（h=5 m，F_sub = 26383.1 N，anat 网格）

| 配置 | deck 里有 `<Contact>`？ | `max` (MPa) | 中位 (MPa) |
|---|---|---:|---:|
| `plantar_bc="contact"` + `contact_hold_spring_k=100` | **是** | **151.896519** | **65.185651** |
| `plantar_bc="spring"` + `spring_k=100` | **否** | **151.896519** | **65.185651** |

**`np.allclose(a, b, rtol=0, atol=0)` ⇒ `True`（逐位相同）**

**⇒ 同一 k 下，有没有那个 `<Contact>` 块，结果一模一样。**

### 1.2 接触参数**全部无效**（若接触真在起作用，它们必须有效）

h=5 m，`hold_spring_k=None`，STATIC：

| 变量 | 范围 | 结果 |
|---|---|---|
| `penalty` | 0.1 → **1e8** | **失败签名一字不变**（`neg_first=8`，6 步，末 t=0.1667） |
| `gap` | **−0.5**（深过盈）→ +0.05 | **一字不变** |
| `auto_penalty` | 0 / 1 | 一字不变 |
| 板模量 `pad_e_mpa` | 0.185 → 200000 | **逐位相同** |
| 板固定（全 50 节点 / 只底面 25） | — | **逐位相同** |
| **外加载荷** | **×1e-6（0.026 N）** → ×1 | **一字不变** |

**最硬的一条：深过盈 0.5 mm + `penalty = 1e8`（巨硬接触）仍与 `penalty=0.1` 同款失败**
⇒ **接触界面根本没有参与计算。**

### 1.3 其它已排除的嫌疑

- **参考构型无倒置**（读 deck 参考坐标算 `det J`）：`calcaneus` 634 hex8 倒置 **0**（min 1.539）、
  `pad_floor` 16 hex8 倒置 **0**（min 1279）
- **接触法向正确**：`pad_top` 法向 `(0,+1,0)`（朝骨，正向比例 **100%**）；
  `plantar` 法向 Y 中位 **−0.796**（朝板）⇒ 不是项目已知的 B5 坑
- **`subtalar_joint` 238 invalid facets**：该警告在**能收敛的 `fixed` deck 里同样出现** ⇒ 既有/无害
- **FEBio 认出了接口**：日志有 `CONTACT INTERFACE DATA / contact interface 1 - Type: sliding-elastic`
  及全部参数 ⇒ **认了但没用**
- **`laugon` 写法有效**：日志回显 `laugon : PENALTY (0)` ⇒ 解析为「罚方法」

---

## §2 · 对下游的**扩散影响**（必须一并更正）

| 受影响处 | 现写法 | 应改为 |
|---|---|---|
| `V1_FIXABILITY_ASSESSMENT.md` §1.1 数据列「`contact+hold100`」 | 当作**接触口径**的一列 | **它是 `spring k=100` 的数**（逐位相同） |
| `V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` §1「第 1 层：不再焊死 ✅ 已实现」 | 认为 `contact`/`spring` 两条路都实现了承载 | **`contact` 路不存在**；承载**只**由 `spring`（或 `hold` 弹簧）提供 |
| 同上 §1「`contact + hold_k` 的 `hold_k` 是自由旋钮」 | 把它当作接触路径的参数 | **它就是 `spring` 的 `k`** —— 所谓「接触路径」=「弹簧路径 + 一个空块」 |
| `PLANTAR_CONTACT_FIX_20261015.md` 全篇 | 分析「跖面接触为何不收敛」 | 分析的是**`spring` BC 在 k 小/无支撑时不收敛**，以及**一个惰性接触块** |
| 本轮前件 `V1_CONTACT_UNCONVERGED_ADDENDUM.md` §1（罚因子差 5 万倍） | 归因到接触罚因子 | 该归因**不成立** —— 接触惰性 ⇒ **罚因子怎么写都无所谓** |
| 本轮前件同 §3「穷举排除」 | 结论「不是接触参数问题」 | **✅ 该结论被本件坐实并给出原因**（接触惰性） |

**⇒ 「跖面接触不收敛」这个议题的**正确表述**是**：

> **`spring` BC（含 `hold` 弹簧）在 k 偏小/无支撑时不收敛（欠约束）；
> 而 `contact` BC 从未真正生效，因此从未被检验过。**

---

## §3 · 仍未回答的问题（诚实）

**为什么那个 `<Contact>` 块惰性？** 本件只证明「它惰性」，**未定位机制**。已知的排除项见 §1.3。

**候选（未验证，按怀疑度排序）**
1. **`<SurfacePair>` / `<Contact>` 的放置或引用方式**（`<SurfacePair>` 在 `<Mesh>` 内、
   `<Contact>` 在 `<Step>` 内 —— 均符合常规，但未与一个**最小可复现算例**对照过）
2. **对偶面 `pad_top` 整体被固定**（虽然试过 `bottom_only` 放开了上排，仍惰性 ⇒ 证据偏弱）
3. **FEBio 4.13 对 `sliding-elastic` 的某个前置条件未满足**（未定位）
4. **`<Contact>` 块被写在 `<Loads>` 之后**——顺序是否影响未测

**⇒ 定位手段（下一轮第一顺位）**：
做**最小可复现算例**（两个 cube 的 sliding-elastic 接触），
再逐项加回（`SurfacePair` 位置 / 固定面 / 板厚 / 分格），**找出哪一项让接触失效**。

---

## §4 · 本件做了什么

| 动作 | 结果 |
|---|---|
| 双接触（距下帽 + 跖面板）、双初始接触、零旋钮 | ❌ 全部失败（同款签名） |
| `contact` vs `spring` 逐位对照 | ✅ **逐位相同 ⇒ 接触惰性** |
| 深过盈 + 硬罚探针（0.5 mm / 1e8） | ❌ 无影响 ⇒ 再证惰性 |
| 已知能收敛配置的**控制实验**（`hold100`） | ✅ **仍收敛**（151.8965）⇒ **本轮代码改动零破坏** |

**代码改动（全部 opt-in，默认路径逐位不变；零回归 `pytest -k "plantar or contact or thums or bone_feb"`
⇒ 72 passed）**：`plantar_bc.py` 的 `pad_fix_bottom_only`、`contact_penalty`、
`contact_pad_thickness_mm`、`contact_auto_penalty`；`thums_feb.py` 同名透传（含 `contact_pad_e_mpa`）。

---

## §5 · 诚实边界

1. 「逐位相同」是**本机同一 FEBio 4.13、同一 deck 生成器**下的结果（`rtol=0, atol=0`）；
   跨机器可能有末位差异（手册 §8.12.5 提示多核可给略不同结果）。
2. 本件**没有**证明 `<Contact>` 块在所有情形下都惰性 —— 只说**在本模型/本配置下**它贡献为 0。
3. 本件的**归因**（惰性）是**否证式**的：接触参数全无效 ⇒ 接触未参与。
   **机制未定位**（§3）。
4. 全部 h=5 m / a=0 / anat / `PressureLoad` 口径的**诊断跑**，无一条构成可交付口径。
5. **不动任何判定**：V1 仍为 ❌【口径不可复现】。

---

## §6 · 溯源

| 项 | 位置 |
|---|---|
| 逐位对照 | `temp/v1_cmp.py`（`contact_hold100` vs `spring_100`） |
| 控制实验 | `temp/v1_ctrl.py`（`hold100` 三变体 ⇒ 151.8965 相同） |
| 深过盈探针 | `temp/v1_jam.py`（gap −0.5 / penalty 1e8） |
| 穷举 | `results/opensim_fe/v1_penalty_scan.json`、`v1_load_scale_scan.json`、`v1_gap_window_scan.json`、`v1_pad_fix_bottom_only_scan.json`、`v1_pad_material_scan.json` |
| 双接触 | `results/opensim_fe/v1_dual_contact_scan.json` |
| 参考构型 det J | `temp/v1_padjac.py` |
| 接触法向 | `temp/v1_normal.py` |
| 前件 | `V1_CONTACT_UNCONVERGED_ADDENDUM.md`、`V1_UNDERCONSTRAINT_METHODS_SURVEY.md`、`V1_PLANTAR_LOAD_BEARING_FEASIBILITY.md` |
| 已知收敛基线 | `temp/opensim_fe/thums_sweep_anat_contact_hold100/` |
| 铁律 | `[MEMORY gotcha #26]`；`docs/notes/FE_PRESOLVE_GATES.md` **闸门 3**「你传的参数 ≠ deck 里写的 ≠ FEBio 读到的」 |
