# P2 · V1 结果 —— 修法 1 单变量重跑 T6 能力门：**0/4 → 3/4 PASS**

> **日期**：2026-10-15 ｜ **性质**：**单变量**重跑（只改压板 z 注册基准，其余逐字不动）
> **脚本**：`scripts/opensim_fe/p2_t6_regfix.py`（加法式；monkeypatch 注入注册修正）
> **产物**：`results/opensim_fe/P2_V1_GATE_VERDICT.json`、`P2_V1_REGFIX.json`、`P2_V1_REGFIX_post.json`
> **前置**：`P2_T6_ROOTCAUSE.md`（根因定位）、`P2_BC_SCOPE_INVENTORY.md` §5ter/§6（原结案）
> **判定**：**逐字调用官方 `p2_head_contact.judge_case`**（BEAR∈[0.80,1.20]、`R2_MIN=0.20`、`bend=(r2≥0.20) OR (szz_max>0)`）

---

## 0. 结论

**本件直接推翻 §5ter 的「能力门 0/4」**：修掉**一处**几何注册基准后，**4 档中 3 档通过官方能力门**，接触承载偏差 ~2e-06。

| h | solve load | 求解 | `bear_frac` | `r2` | **官方 judge_case** |
|---:|---:|---|---:|---:|---|
| **10** | 19 675.9 N | **rc=0** | **1.00000185** | 0.099 | **PASS** |
| **20** | 29 794.1 N | **rc=0** | **1.00000248** | 0.067 | **PASS** |
| **30** | 38 257.3 N | **rc=0** | **1.00000319** | 0.034 | **PASS** |
| 50 | 52 618.6 N | **rc≠0** | — | — | FAIL（无 h5） |

**对照实录**（`P2_HEAD_CONTACT.json`）：**4 档全 FAIL**，且**最低诊断档 1220.8 N** 即 1888–2515 负 Jacobian、首步 t=0.1 全单元倒置。

> **⇒ 「与载荷无关的全体倒置」已消失。** 该现象正是 `P2_T6_ROOTCAUSE.md` 判定的**纯刚体模态**签名 ——
> 注释所依赖的前提（「t=0 必须由接触提供 z 锚位」）原先从未成立，现已成立。

## 1. 单变量声明（可审计）

| 改了什么 | 改前 | 改后 |
|---|---|---|
| 压板 z 注册基准 | `z_wall = geom["a_max"] + GAP_MM`（`a_max` 来自**顶带节点集**） | `z_wall = max(node z over **保留 cap 面片**) + GAP_MM` |

**T6 上二者差 1.2342 mm（实测下移量）。**

**其余逐字不动**：`GAP_MM=−0.3`、`PLATE_DISP_MM=0.5`、`TOP_FRAC=0.30`、`LOAD_FRAC=0.20`、
`CAP_AXDIV_MIN=0.5`、`penalty = 50×F_ref/57941.3`、`fric=0`、`lateral_support=[(0.05,0.95)]`、
`N_ANCHOR=3`、`TIME_STEPS=10`、`bottom_fix_axial_only=False`、`with_shell=True`。

**⚠️ 对股骨（正对照）是恒等变换**：实测 `a_max ∈ cap` ⇒ 下移量 **0.0000 mm** ⇒ 板几何数值不变
⇒ 生成的 `.feb` 不变 ⇒ **P5.0 的既有结论不受影响**（该论证是**几何解析**的，未重跑股骨 FE）。

## 2. 两条必须记账的保留

### 2.1 `h=50` 仍不收敛（诚实的剩余缺口）

- h=50 = **52 618.6 N** 时 FEBio 非零退出，无 h5。
- ⚠️ **不能据此说「T6 不如股骨」**：股骨 P5.0 在 **57 941.3 N** 通过；两者**载荷不同、几何不同**，
  本件未做同载荷对照 ⇒ **差异未归因**。
- 但可确定：**能力门不是 4/4**，引用时必须写 **3/4**。

### 2.2 🔴 σ 幅值**不自洽**（能力门不检查这一项）

| h | σ median | σ p95 | σ max | σ_zz ∈ | σ_c（皮质域） |
|---:|---:|---:|---:|---|---:|
| 10 | 209.6 | 897.6 | **1421.4** | [−686.7, +29.2] | 150 |
| 20 | 319.8 | 1410.7 | **2276.9** | [−1013.3, +52.7] | 150 |
| 30 | 416.3 | 1919.9 | **3298.9** | [−1268.7, +89.7] | 150 |

**⇒ h=10 的 `|σ_zz|max = 686.7 MPa` 已是皮质骨强度（150 MPa）的 4.6×；σ_max 1421 MPa 达 9.5×。**

- **这不在能力门判据内**（判据只查 收敛 / 接触承载 / 弯曲路径存在），**故 3/4 PASS 是真的**，
  但**不等于「模型物理可信」**。
- **性质与 V8b 的既有发现同类**：`SEGMENTS_6` 载荷段纯骨柱 ⇒ σ 超骨强度 3.5–5.6×（**模型不自洽**），
  真根因在**载荷路径缺串联柔顺**（G6 方向对、形式不可锚定；G7 待立项）。
- **本件不做归因**。数值已列全，供后续判断。**不得**引用 T6 的 σ 幅值作结论。

## 3. 🐞 途中发现的**第二个 harness bug**（与物理无关，但会制造假阴性）

- `p2_head_contact.official_domain("T6")` 返回 **`shell_t6_cort`**（.feb 里的域名）；
- 但 **h5 里该壳域被 `post_process` 命名为 `__unnamed_2`**（实测 h5 域 = `ghost / plate / t6_spon / __unnamed_2`）；
- ⇒ `domain_stat(h5, "shell_t6_cort")` **直接 KeyError** ⇒ `bear_frac` / `r2` / `szz_max` 全为 `None`
  ⇒ `judge_case` 判 **FAIL**，**即使求解完美收敛**。
- **本次 V1 首轮就撞上它**，产生过一次假 FAIL。本件按**元素数 712** 把官方域名解析到 h5 键后复判。
- ⚠️ **含义**：即便 §5ter 当时收敛了，`domain_stat` 也会失败 ⇒ **T6 的皮质域判定在现 harness 下从未真正跑通过**。
  建议修 `official_domain`/`domain_stat` 的域解析（或给 .feb 的壳域补 `<name>`）。

## 4. 自错记录（如实）

| # | 自错 | 修正 |
|---|---|---|
| 1 | V1 首版把 `R2_MIN` **臆造成 0.90**（官方为 **0.20**，且闸门是 **OR**：`r2≥0.20 或 szz_max>0`） | 改为**直接引用**官方常量 + 调官方 `judge_case`，不再二次发明判据 |
| 2 | `P2_T6_ROOTCAUSE.md` §3 的**「修法 2（恢复底端 z 约束）」已作废** | 底端一固定 ⇒ 底端反力分掉载荷 ⇒ `C = F − R ≠ F` ⇒ **`bear_frac` 判据当场失效**。P5 设计**必须**让骨轴向自由、由接触供 z 锚位 ⇒ **注册修正是唯一正确修法** |
| 3 | `temp/` 下脚本 `ROOT` 又写成 `parents[1]`（应为 `parents[2]`），两次 | 已修（MEMORY 已有此坑记录，仍复发） |

## 5. 对既有判定的影响

| 项 | 状态 |
|---|---|
| **A2** | **不变**：仍 **4/5 PASS（混 BC）**，胸椎仍唯一遗留 `常值1`（**本件未重跑 A2 矩阵**） |
| **统计维度** | **不变**：仍 **2/4** |
| **A1** | **不变**：仍 FAIL（`A1_BOUNDARY_WINDOW_REPORT.md`） |
| **P2 结案（§6）** | ⚠️ **直接受冲击**：结案理由是「T6 能力门 0/4、根因未定位」。本件证 **0/4 源于一个注册 bug**，且根因已定位、单变量修复后 **3/4 PASS** ⇒ **§6 的「当前几何下不可交付」不再成立**，应改记为**「修法后 3/4 可交付，h=50 待归因」**。**该改记须你拍板**（本件只提供证据，不自行改判）。 |
| **冻结输入** | **未改**（本件全部产出为新件；`p2_head_contact.py` 一行未动，修法经 monkeypatch 注入） |

## 6. 下一步（建议，待拍板）

1. **改记 P2 结案**（§6）：0/4 → 3/4，根因与修法入账。
2. **h=50 归因**：单变量（如降 penalty / 加步数 / 查 FEBio 日志第几步发散），**不与股骨比**（载荷不同）。
3. **σ 幅值问题单独立项**：与 V8b/G6/G7 同类（载荷路径缺串联柔顺）⇒ 很可能**同一把钥匙**。
   若成立，则 T6 的 σ 与 V8b 一并解决 —— **这是当前最高价值的一条线索**。
4. **修 harness 域解析**（§3），否则后续 T6 判定仍会假 FAIL。

## 7. 溯源

| 项 | 位置 |
|---|---|
| 根因件 | `results/opensim_fe/P2_T6_ROOTCAUSE.md` |
| 单变量脚本 | `scripts/opensim_fe/p2_t6_regfix.py` |
| 前置闸门 | `scripts/opensim_fe/contact_registration_audit.py` |
| 能力门判定 | `scripts/opensim_fe/p2_head_contact.py:79-80`（阈值）、`:312-328`（`judge_case`） |
| 修法注入点 | `p5_head_contact.plate_geom_from_mesh:114`（原 `z_wall = a_max + GAP_MM`） |
| 原结案 | `results/opensim_fe/P2_BC_SCOPE_INVENTORY.md` §5ter / §6 |
| V1 产物 | `results/opensim_fe/P2_V1_GATE_VERDICT.json`、`P2_V1_REGFIX.json`、`P2_V1_REGFIX_post.json` |
| h5 与 .feb 域名不符 | 实测 h5 域 = `ghost / plate / t6_spon / __unnamed_2`（壳 712 单元） |

**复现**：

```powershell
cd D:\Project\climbing_fall_analysis
& .venv\Scripts\python.exe scripts\opensim_fe\contact_registration_audit.py   # 零 FE 前置闸门
& .venv\Scripts\python.exe scripts\opensim_fe\p2_t6_regfix.py                 # V1：修法1 单变量
#   诊断单档：--loads 2000   （⇒ solve load 1220.8 N，与实录同档）
```
