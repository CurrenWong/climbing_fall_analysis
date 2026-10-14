# A2 判定修订（入账件）—— **4/5 PASS**（混 BC 条件下）

> **日期**：2026-10-14 ｜ **性质**：**判定入账**（含条件声明；原统一-BC 判定**保留**）
> **依据**：`results/opensim_fe/B1_A2_STAT_CALIBER_RULING.md` §4.4（**官方管线**复核）
> **触发**：用户拍板「接受混 BC 入账」+「全链统一 BC 列为后续」

---

## 0. 修订后的判定

| 口径 | 判定 | 命中 | 说明 |
|---|---|---|---|
| **统一 BC（BC2，原判定）** | **FAIL** | **3/5** | **保留不覆盖**；仍是不做任何 BC 改动时的官方数字 |
| **混 BC（股骨换 P5 接触 BC）** | **PASS** | **4/5** | **本次入账**；逐部位口径 5/6 |

**逐目标（混 BC）**：足部 ✅｜胫骨 ✅｜**股骨 ✅（常值1 → 阶跃）**｜腰椎+胸椎 ❌｜颈椎 ✅

**唯一遗留**：**胸椎 `常值1`**（独立模型缺口，与本轮 BC 改动**无关**）。

## 1. 依据链（含方法正对照）

1. **口径钉死**：官方 A2 = `s4_alignment_bc2.py` Route 乙 + **median**
   （`fracture_matrix.json` 一代与 `paper_alignment.py` 用 p95，为参考/敏感性）。
2. **A2 股骨规则**（`s4_alignment_bc2.py:176`）：`(ends ≥ 80) | (shaft ≥ 220)` ——
   弱子区域是**股骨颈 80**，非骨干 220。
3. **官方管线复核（`B1_A2_STAT_CALIBER_RULING.md` §4.4）**：`s4_alignment_bc2.py` Route 乙
   **原封执行**（仅外层 monkeypatch `BC2_DIR`/`BC2_JSON`，femur h5 换 `p5_pen_h50`、
   `fe_solve_load_n` 15227.65 → 57941.3 同步）。
   - **正对照（未 patch）**：`n_hit=3/5`、`n_loose=4`、`FAIL`、股骨 `常值1`
     ⇒ **逐位复现仓库记录** ✅ ⇒ 驱动可信。
   - **Patch 后**：股骨 `ends_stat=103.76`、`first_ends=33.0`、behavior **阶跃**；
     **`n_hit` 3/5 → 4/5、`verdict` FAIL → PASS**；逐部位 4/6 → 5/6。
4. **支撑复核**：线性度 σ/load 离散 ends **2.76%** / shaft **0.53%**（单点缩放合法）；
   解剖三区 + 弱区比例敏感性（f∈{0.20…0.35} 全给阶跃；唯一反例 f=0.15 距阈仅 6.7%）；
   **弱区排除头球后 weak_ref 反升**（112.5–120.6 vs 82.5–91.3）⇒ 非载荷伪热点。

## 2. 本入账的**条件性**（必须与判定并列引用）

1. **混 BC**：仅**股骨**采用 P5 接触 BC，其余 8 部位仍 BC2。**本判定不是单一一致模型的结果。**
2. **4/5 恰等于阈值 ≥4/5（无余量）** —— 余量为零，任一处波动即翻回 FAIL。
3. 判定基于**单一 h50 解 + 线性缩放**（线性度已验 2.76%/0.53%，但非逐高度重解）。
4. **胸椎仍是真实漏判**；A2 的「完整复现」未达成（论文 6 目标 vs 我方 6 目标中 1 漏）。

## 3. 对统计维度的影响

| 代号 | 状态 |
|---|---|
| A1 | ❌ FAIL（真差异，`A1_BOUNDARY_CONVENTION_REPORT.md`） |
| **A2** | **FAIL（统一 BC）→ PASS（混 BC，本次入账）** |
| A3 | ✅ 3/3（G9 收口） |
| A4 | 【架构不可判定】 |

⇒ **统计维度：1/4 → 2/4**（A3 ✅ + A2 ✅；A1 ❌、A4 不可判定）。

## 4. 后续（用户指定：(B) 全链统一 BC 列为后续）

- **(B) P5 主工程**：9 部位**统一** BC 重解 ⇒ 把本件的「混 BC 上界」升级为「单一模型结果」。
  **成本大**，且会同时动另外 8 部位（风险远大于收益）⇒ 列为**后续**，不阻塞本入账。
- **胸椎模型缺口**（独立于本轮）：`常值1`，需新的机制假设 + 预注册。
- **解剖颈的骨性标志定义**（替代比例口径）—— 本件用比例，已在 §1.4 做敏感性。
- **`MEMORY.md` 归档**：已 ~21.6 KB，超 10 KB 硬限一倍余。

## 5. 溯源

| 项 | 位置 |
|---|---|
| 修订依据（含官方管线复核） | `results/opensim_fe/B1_A2_STAT_CALIBER_RULING.md`（§3.4/§3.5/§4.1–4.4/§5/§6） |
| 官方管线驱动脚本 | `temp/b1_official_rerun.py` |
| 独立复算脚本 | `temp/b1_caliber_recompute.py`、`b1_three_zone.py`、`b1_three_zone2.py`、`b1_full_verdict.py`、`b1_three_zone_verdict.py`、`b1_linearity.py` |
| 原判定（保留） | `results/opensim_fe/s4_alignment_bc2.json`（`A2.n_hit=3`、`verdict=FAIL`）、`S4_ALIGNMENT_BC2_REPORT.md`、`A2_MISS_DIAGNOSIS.md` |
| 口径裁定 | `A2_CONVENTION_RULING.md`（严格/逐部位）、`A2_AGGREGATION_RULING.md`（维持 median） |
| P5 源数据 | `temp/opensim_fe/p5_head_contact/femur_r/p5_pen_h{10,20,30,50}/femur_r.hdf5` |
