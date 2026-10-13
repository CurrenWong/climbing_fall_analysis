# A2 口径决议 —— strict vs loose（由 [Y25] 编码判定）

> **日期**：2026-10-13 ｜ **性质**：口径裁定（**独立于我们的结果**，由 [Y25] 的既有编码决定）
> **触发**：`results/opensim_fe/A2_MISS_DIAGNOSIS.md` §3 的口径警示
> （strict 3/5 vs loose 4/5，差异全来自「腰椎+胸椎」一个 target）。

---

## 0. 裁定

> **采用【严格 / 逐部位】口径。宽松计数废弃。**
> **A2 = FAIL 确认**（胸椎是**独立的、真实的漏判**，宽松计数把它掩盖了）。

---

## 1. 依据（**代码内已有两个编码，取更贴近 [Y25] 的那个**）

### 1.1 编码 A：`PAPER_STEP`／`paper_alignment.py:74`（**逐部位**，与 [Y25] 同名注释）

```python
#: 论文 [Y25] 的「阶跃响应」部位（Logistic 因完全分离而拟合失败）。
PAPER_STEP    = {"足部", "胫骨", "股骨", "腰椎", "胸椎", "颈椎"}
PAPER_GRADUAL = {"腓骨", "骨盆", "颅骨"}
```

⇒ **[Y25] 的阶跃部位是 6 个独立部位，腰椎与胸椎各自独立**。

### 1.2 编码 A 的判定逻辑（`paper_alignment.py:475-483`）—— **逐部位比对**

```python
for cn, bone, _, subs in REGIONS:
    b = _behav(logi1d[cn])
    target = "阶跃" if cn in PAPER_STEP else "渐进"
    if b == target: verdict = "✓"
    elif target == "渐进": verdict = "✗"    # 论文渐变、我方却恒定 ⇒ 未复现
    else: ...
```

⇒ **参考实现本身就是逐部位判**；没有「组内任一即命中」的逻辑。

### 1.3 编码 B：`A2_TARGETS_STRICT`／`s4_alignment_bc2.py`（**分组**，5 个 target）

把 **「腰椎+胸椎」合成一个 target**，再配 strict（组内全中）/ loose（组内任一）两种计数。

**⇒ 该分组在 [Y25] 编码里没有依据**（`PAPER_STEP` 明列两个独立部位）。
**⇒ 宽松计数的「过」全部来自这个无依据的分组**：腰椎=阶跃使命中，掩盖了**胸椎=常值1 的独立漏判**。

---

## 2. 逐部位判定（采用编码 A 的口径）

| [Y25] 阶跃部位 | 我方行为 | 命中 |
|---|---|:--:|
| 足部 | 阶跃 | ✅ |
| 胫骨 | 阶跃 | ✅ |
| 股骨 | **常值1** | ❌ |
| **腰椎** | 阶跃 | ✅ |
| **胸椎** | **常值1** | ❌ |
| 颈椎 | 阶跃 | ✅ |

**⇒ 4/6 命中**（不是「3/5」，也不是「4/5」——那两个数都是**分组口径**的产物）。

> **✅ 已实测验证**（逐部位跑 `paper_alignment.PAPER_STEP` 比对，读冻结的 `logistic_yi`）：
> ```
> 足部  complete_separation  阶跃  OK      腓骨  all_one  渐进  MISS
> 胫骨  complete_separation  阶跃  OK      股骨  all_one  阶跃  MISS
> 腰椎  complete_separation  阶跃  OK      胸椎  all_one  阶跃  MISS
> 颈椎  complete_separation  阶跃  OK      骨盆  all_zero 渐进  MISS
> 颅骨  complete_separation  渐进  MISS
> PAPER_STEP 命中 = 4 / 6
> ```

**⇒ 无论用哪个计数法，`胸椎` 这个独立漏判都在** ⇒ **A2 判定不变：FAIL。**

---

## 3. 顺带更正：`A2_TARGETS_STRICT` 与 `PAPER_STEP` 的内部不一致

| 编码 | 阶跃部位的表示 | 部位数 |
|---|---|---|
| `PAPER_STEP`（`paper_alignment.py`） | `{足部, 胫骨, 股骨, 腰椎, 胸椎, 颈椎}` | **6** |
| `A2_TARGETS_STRICT`（`s4_alignment_bc2.py`） | `[足部], [胫骨], [股骨], [腰椎,胸椎], [颈椎]` | **5** |

**同一仓库内两条编码不一致** —— 这是「strict 还是 loose」这个歧义的**真正来源**。
⇒ **建议**：以 `PAPER_STEP` 为准（它带 [Y25] 出处注释、且早于 S4），
把 `A2_TARGETS_STRICT` 的「腰椎+胸椎」拆开，并把阈值从「≥4/5」相应重述为**逐部位**门限（如 **≥4/6**）。

> ⚠️ **本件不擅自改代码/阈值**——`A2_TARGETS_STRICT` 与「≥4/5」是**已冻结的判据**，
> 改动须走预注册。本件只**裁定口径**并**指出内部不一致**。

---

## 4. 影响

| 项 | 结论 |
|---|---|
| **A2 判定** | **不变：FAIL**（3/5 或 4/6，胸椎独立漏判在两种口径下都算漏） |
| **宽松计数** | **废弃**——它用一个无依据的分组掩盖真实漏判 |
| **修复路径** | 不变：**解剖端分布载荷（改 FE BC，大件）**，见 `A2_MISS_DIAGNOSIS.md` §4 |
| **待办** | 是否把 `A2_TARGETS_STRICT` 对齐 `PAPER_STEP`（拆开脊柱）—— 须预注册，**待用户拍板** |

---

## 5. 溯源

| 项 | 来源 |
|---|---|
| [Y25] 阶跃部位编码 | `scripts/opensim_fe/paper_alignment.py:73-76`（`PAPER_STEP` / `PAPER_GRADUAL`） |
| 参考判定逻辑（逐部位） | `scripts/opensim_fe/paper_alignment.py:475-483` |
| 分组编码 | `scripts/opensim_fe/s4_alignment_bc2.py`（`A2_TARGETS_STRICT` / `a2_verdict`） |
| 实测行为 | `results/opensim_fe/s4_alignment_bc2.json`（`logistic_yi`） |
| 诊断 | `results/opensim_fe/A2_MISS_DIAGNOSIS.md` |
