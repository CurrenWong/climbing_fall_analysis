# 审计 · **「标定前载荷链」残留**（2026-10-15）

> **触发**：`V1_PLANTAR_BC_SENSITIVITY.md` 的 7.6% 未归因差异 ⇒ 查出**归档 sweep 基线
> 是标定前载荷链的产物**。
> **性质**：**审计件**（只报，不改活脚本 —— 是否重跑属**拍板**事项）。

---

## §0 · 结论

**`results/opensim_fe/s2_thums_height_sweep.json` 与 `s2_thums_height_sweep_anat.json`
是 2026-10-08 刚性载荷链标定「之前」的产物**；**5 个活脚本仍把它们当输入** ⇒
**它们的结论建立在已被取代的载荷链上**。

---

## §1 · 证据链（四步，每步可复现）

| # | 步骤 | 工具 | 结果 |
|---|---|---|---|
| 1 | 探针路径 vs sweep 路径的 deck **sha256** 对比 | `consistency_probe_vs_sweep.py` | **逐位相同**（只差 plotfile 文件名）⇒ 差异不在 deck |
| 2 | **同一份 deck** 过两套求解+后处理 | `..._step2.py` | 同一 deck 给 **362.561133**（＝探针值）⇒ **归档的 336.92 才是异常** |
| 3 | 用**现在的链**重跑 `fixed` 基线（两网格） | `s2_thums_sweep.py --tag _refixed_now` | `now/存档` = **1.0431 … 1.1202** |
| 4 | 比对**载荷列** `subtalar_kN` | `chk_load` | **载荷比 ≈ gauge 比（差 <0.3%）** ⇒ σ 的变化 **100% 由载荷解释**（线性响应） |

### 1.1 归因表（`now / 存档`）

| h | 载荷 `subtalar_kN` 比 | `gauge_max` 比 | 差 |
|---:|---:|---:|---:|
| 1 | 1.120159 | 1.1182 | 0.2% |
| 2 | 1.077650 | 1.0761 | 0.15% |
| 3 | 1.101389 | 1.0994 | 0.2% |
| 4 | 1.105446 | 1.1030 | 0.25% |
| 5 | 1.081017 | 1.0790 | 0.2% |
| 6 | 1.044283 | 1.0431 | 0.1% |
| 8 | 0.974597 | 0.9753 | −0.07% |

**⇒ 载荷与 σ 同比例变化** ⇒ 差异**不是**后处理、**不是**网格、**不是**数值 —— **就是载荷链**。

### 1.2 仓库内的书面证据（两处独立）

- `V1_CALIB_TODO.md:18`：**「`_anat` 扫描用**标定前载荷链**（10-04 vs 10-08 标定）」**
- `ALIGN_A3_RERUN_TODO.md:14`：**「L2 刚性标定 `RIGID_STANDING_CHAIN`（2026-10-08 采纳）
  改了默认刚性链路，缓存 `s1_h5_opensim.json` 是**标定前快照**」**

---

## §2 · **受影响的活脚本（5 个，均未改动）**

| 脚本 | 读什么 | 影响 |
|---|---|---|
| `scripts/opensim_fe/nonvertical_s5_contact.py:51` | `s2_thums_height_sweep.json` 的**逐高度载荷** | 线 2 / S5 全部结果按旧链 |
| `scripts/opensim_fe/joint_load_probe.py:64` | `s2_thums_height_sweep_anat.json`（h=5, act=0 的**载荷**） | 探针载荷旧链 |
| `scripts/opensim_fe/proof_s5_facet_fix.py:54` | **硬编码** `SUBTALAR_LOAD_N = 21086.2`（旧 h=2；**现链 22723.6**） | 单点证明用旧载荷 |
| `scripts/opensim_fe/a1_gauge_check.py:16` | `s2_thums_height_sweep_anat.json` | A1 侧数字旧链 |
| `scripts/opensim_fe/plantar_bc_probe.py:57` | `s2_thums_height_sweep_anat.json` | 跖面 BC 探针旧链 |

**另**：`docs/reports/总体进度.md:196`、`docs/paper/论文v2_出处索引.md:120`、
`V1_FOOT_TODO.md`、`results/opensim_fe/A1_SOLUTION_PATH.md`、`MASS_CALIB_REPORT.md`、
`S2_THUMS_REPORT{,_ANAT}.md`、`NONVERTICAL_S5_CONTACT_REPORT.md` 等**引用了这些旧数字**。

**⚠️ 注**：`论文v2_出处索引.md:120` 与 `总体进度.md` 的 **V1 行**已**显式标注**〔旧链；当前链 2 m …〕
⇒ 这两处**无污染**（已如实标注）。**其余未标注者**需逐个判定。

---

## §3 · 建议（**待拍板**，本件不执行）

1. **给所有 `*_sweep*.json` 加「载荷链指纹」字段**（如 `subtalar_kN` 的输入哈希或
   标定版本号），并在**使用处断言指纹匹配** ⇒ 否则「拿旧基线比新结果」还会再犯。
2. **逐个判定 §2 的 5 个脚本**：仍要用的 ⇒ 同链重跑；已归档的 ⇒ 加**显式旧链标注**。
3. **文本更正**：§2 末列出的未标注引用处，按**替换式**补「〔旧链〕」或换成新链数字。

---

## §4 · 方法论教训（**建议升格为闸门 10**）

> **比较两个数之前，先证它们出自同一条链。**
> 「基线」若来自旧产物，必须先**同链重跑**，再比。

这是「**检查器自身要有已知正对照**」（纪律 4）在**基线**上的对应物 —— 本件就是
**先信了旧基线**（v1 用 0.842–0.919），才得出「符号相反」这一**错误**结论。

**⇒ 建议闸门 10 · 「基线同链」**：任何「A 与 B 相比 ⇒ 差异 x%」的结论，
必须给出**同一链内**的双侧数值；跨链比较一律作废重跑。

---

## §5 · 溯源

| 项 | 位置 |
|---|---|
| 触发件 | `results/opensim_fe/V1_PLANTAR_BC_SENSITIVITY.md`（v2 §3） |
| 审计脚本 | `scripts/opensim_fe/consistency_probe_vs_sweep{,_step2}.py` |
| **同链基线**（应用） | `s2_thums_height_sweep_refixed_now.*` · `..._anat_refixed_now.*` |
| **旧链产物**（标注、勿直接比） | `s2_thums_height_sweep.json` · `s2_thums_height_sweep_anat.json` |
| 书面证据 | `V1_CALIB_TODO.md:18` · `ALIGN_A3_RERUN_TODO.md:14` |


---

## §6 · **工具已落地**：载荷链指纹（闸门 10 的日常工具）

### 6.1 新增件

| 件 | 作用 |
|---|---|
| `src/climbing/coupling/loadchain.py` | **纯函数**：`fingerprint()`（sha256 前 16 hex，6 位有效数字归一 ⇒ 吸收浮点末位噪声）· `table_from_json()`（读 sweep JSON 的载荷列，`_kN` 自动归一到 N）· `compare()`（**先在共有高度上**比指纹 + 逐高度比值跨度） |
| `scripts/opensim_fe/loadchain_audit.py` | CLI：扫 `results/opensim_fe/*sweep*.json`，逐份判 ✅同链 / ⚠️常数缩放 / ⛔形状不同。**带磁盘缓存** `LOADCHAIN_CURRENT.json`（`_load_at` 每高度跑一次 OpenSim ⇒ 必须缓存，且**增量落盘**可中断续跑） |
| `tests/test_loadchain.py` | 8 个单测，**含正/负对照**（纪律 4）与**真实数据回归锚** |

### 6.2 **审计结果（`loadchain_audit.py`，当前链指纹 `06cee88a51412a3c`）**

| 文件 | 载荷指纹 | 比值跨度 | 判定 |
|---|---|---|---|
| `s2_thums_height_sweep.json` | `538abf07260d11c3` | 17.49% | ⛔ **异链·形状不同**（0.9400 → 1.1202） |
| `s2_thums_height_sweep_anat.json` | `538abf07260d11c3` | 17.49% | ⛔ **异链·形状不同**（同上） |
| **`s2_thums_height_sweep_plantar_bc_base.json`** | `c7349f65f3afcaf8` | 0.00% | ⚠️ **异链·常数缩放 ×1.0810** ← **本件新发现** |
| `s2_thums_height_sweep_anat_contact_hold100.json` | `b5fdb95d1026d678` | 0.00% | ✅ 同链 |
| `s2_thums_height_sweep_anat_refixed_now.json` | `b5fdb95d1026d678` | 0.00% | ✅ 同链 |
| `s2_thums_height_sweep_contact_hold100.json` | `7fa3e02c69cc06e0` | 0.00% | ✅ 同链 |
| `s2_thums_height_sweep_refixed_now.json` | `7fa3e02c69cc06e0` | 0.00% | ✅ 同链 |
| `s2_thums_height_sweep_v1rerun20261015.json` | `3d796dc1a40b898f` | 0.00% | ✅ 同链 |

**读法**

1. **⛔ 形状不同** ⇒ 不能只乘一个系数修正；**必须同链重跑**（§2 的那两份）。
2. **⚠️ 常数缩放 ×1.081** ⇒ `..._plantar_bc_base.json` **也是标定前的**，但它只是**整体 ×1.081**
   （＝标定因子 +8.1%）⇒ 若某处只关心**比值**（同一文件内的相对比较），它可能仍可用；
   若关心**绝对值**，须 ×1.081 或重跑。
3. **✅ 同链** ⇒ 可安全用作基线（本件新落的 `*_refixed_now` 与 `*_contact_hold100` 都在此列
   ⇒ **验证了上一轮补救的正确性**）。

**⚠️ 指纹可比性 caveat**：指纹建在**该文件与当前链的共有高度集**上 ⇒
**只有高度覆盖相同的两份文件才共享指纹**（如 `anat_*` 两份都是 h=1..6 ⇒ 同 `b5fdb95d…`）。
**判定看 ✅/⚠️/⛔，不要把不同覆盖的指纹直接等值比较。**

### 6.3 载荷链变更的**形状**（`loadcurve_shape_audit.py`）

`current / archived` 比值**不是常数**：**1.1202（h=1）→ 0.9400（h=10）→ 1.0117（h=50）**，
**两次穿过 1.0**；归一化曲线在中段明显更平（h=8：归档 1.7149 vs 当前 1.5509）。

**⇒ 这不是「一个标定系数」，是**另一张载荷-高度表** ⇒ 佐证闸门 10 必须查**形状**而不只是均值。

### 6.4 建议（仍未做，属拍板）

1. 用 `loadchain_audit.py` 把 §2 的 5 个活脚本**逐个判定**（它们读的正是 ⛔ 那两份）。
2. 在需要长期依赖的产物上加**指纹字段**（而非只靠外部审计）。
