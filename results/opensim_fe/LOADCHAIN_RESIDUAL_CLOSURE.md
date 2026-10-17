# 载荷链残留：**暴露面收口**（2026-10-16）

> **触发**：优先 2「载荷链一致性（O）」。
> **性质**：**工具审计 + 全库检索 + 影响量化**。**不改任何判定。**
> **结论**：暴露面**远小于** MEMORY 记录（「5 个活脚本」→ **实际 1 个探索性探针**）；
> **当前权威 V1 扫描全部同链** ⇒ **无已发布数字受影响**；仍需给那 1 个探针加标注。

---

## §1 · 审计（`scripts/opensim_fe/loadchain_audit.py`，闸门 10 的日常工具）

扫描 13 份 `*sweep*.json`（5 份**无载荷列** ⇒ skip），当前链指纹 = `06cee88a51412a3c`：

| 文件 | 载荷指纹 | 比值跨度 | 判定 |
|---|---|---|---|
| **`s2_thums_height_sweep.json`** | `538abf07260d11c3` | **17.4921%** | ⛔ **异链·形状不同**（0.9400→1.1202） |
| **`s2_thums_height_sweep_anat.json`** | `538abf07260d11c3` | **17.4921%** | ⛔ **异链·形状不同** |
| `s2_thums_height_sweep_anat_contact_hold100.json` | `b5fdb95d1026d678` | 0% | ✅ 同链 |
| `s2_thums_height_sweep_anat_refixed_now.json` | `b5fdb95d1026d678` | 0% | ✅ 同链 |
| `s2_thums_height_sweep_contact_hold100.json` | `7fa3e02c69cc06e0` | 0% | ✅ 同链 |
| `s2_thums_height_sweep_refixed_now.json` | `7fa3e02c69cc06e0` | 0% | ✅ 同链 |
| `s2_thums_height_sweep_v1rerun20261015.json` | `3d796dc1a40b898f` | 0% | ✅ 同链 |
| **`s2_thums_height_sweep_plantar_bc_base.json`** | `c7349f65f3afcaf8` | 0% | ⚠️ 异链·**常数缩放 ×1.0810** |

**⇒ 只有 2 份是异链，正是 MEMORY 标的「标定前产物」，且**比值形状不同** ⇒
**不是常数重标定可得**。

## §2 · 暴露面（全库检索，**更正 MEMORY**）

MEMORY 记「**5 个活脚本**仍读它们」。本次逐库检索（`temp/`、`scripts/opensim_fe/`、
`src/`、`tools/`、`tests/`）：

| 命中 | 性质 |
|---|---|
| **`temp/fe_out/a1_gauge_check.py`** | **唯一**真正**读取**异链表的脚本；且是**探索性探针**（自述「零 FE」、what-if：能否把 A1 的 #2 推到 16–18 m） |
| `temp/a3_report_regen/NONVERTICAL_S5_CONTACT_REPORT.md` | 仅**正文文字提及** |
| `temp/fe_out/{smoke,thums_smoke,v1regress,v1rerun}/S2_THUMS_REPORT*.md` | 仅**正文文字提及** |
| `s2_thums_height_sweep_plantar_bc_base.json` 的读者 | **无** |

**⇒ 「5 个活脚本」应更正为「1 个探索性探针（+4 处文字提及）」。**

## §3 · 影响量化（同链重算该探针的结论）

把 `a1_gauge_check.py` 的算法（gauge 口径，阈值 150 MPa）分别跑在

- 异链 `s2_thums_height_sweep_anat.json`（38 行）
- 同链 `s2_thums_height_sweep_anat_refixed_now.json`（6 行）
- 同链 `s2_thums_height_sweep_anat_contact_hold100.json`（6 行）

| h (m) | 异链 gauge_max | 同链 gauge_max | Δ% |
|---|---|---|---|
| 1 | 246.550 | 128.482 | **−47.89** |
| 2 | 276.430 | 161.299 | **−41.65** |
| 3 | 282.246 | 172.167 | −39.00 |
| 4 | 292.678 | 186.053 | −36.43 |
| 5 | 308.506 | 201.710 | −34.62 |
| 6 | 327.404 | 217.656 | −33.52 |

| 口径 | 首折 |
|---|---|
| 异链（标定前） | **1.000 m** |
| 同链 `anat_refixed_now` | **1.603 m**（**+60.3%**） |
| 同链 `anat_contact_hold100` | **3.865 m** |

### ⚠️ 归因边界（**不得过度归因**）

**−34~−48% 的 σ 差 **不能**归因于载荷链单项** —— 该文件对是
`anat`（标定前）vs `anat_refixed_now`，**同时**含：① 载荷链 ② **BC re-fix**
③ 可能的配准变化。三者未分离。

**但可确证的两点**：

1. **载荷比形状差 0.94–1.12（跨度 17.5%）** 经 `首折 ∝ σ^(−1/p)`（`1/p ≈ 3.16`）
   映射 ⇒ **σ 差 ≈ ±14%、首折差 ≈ ×1.5–1.6** ⇒ **与观测的 ×1.60 同量级**
   ⇒ **载荷链残留对首折的影响正好落在「不可忽略」的量级**
   （与本轮「需要 13% 级精度」的语境正好冲突 ⇒ 必须标注）。
2. **该探针的结论必须重跑**：其数字（1.000 m）建立在异链上。

## §4 · 裁定与动作

1. **无已发布数字受影响** —— 全部**当前权威 V1 扫描**（`anat_contact_hold100`、
   `anat_refixed_now`、`v1rerun20261015`、两个 `*_refixed_now`）**均 ✅ 同链**。
2. **暴露面 = 1 个探索性探针** ⇒ 给 `a1_gauge_check.py` 加**异链警示**（本件即其档），
   并把「必须改用同链表重跑」写入其 docstring 视野内。
3. **`s2_thums_height_sweep{,_anat}.json` 与 `..._plantar_bc_base.json` 标记为
   【标定前·禁引】**（本件 + 审计 JSON `temp/loadchain_audit.json` 为凭）。
4. **MEMORY 更正**：「5 个活脚本」→「**1 个探索性探针（+4 处文字提及）**」。

---

## §5 · 溯源

| 项 | 位置 |
|---|---|
| 审计工具 | `scripts/opensim_fe/loadchain_audit.py`；链定义 `src/climbing/coupling/loadchain.py` |
| 审计输出 | `temp/loadchain_audit.json` |
| 影响量化 | `temp/loadchain_impact.py` → `temp/loadchain_impact.json` |
| 唯一读者 | `temp/fe_out/a1_gauge_check.py` |
| 原始记录 | `STALE_LOADCHAIN_AUDIT.md`；`[MEMORY gotcha #1/#10 表]` |
