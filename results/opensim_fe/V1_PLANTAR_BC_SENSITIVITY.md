# V1 跖面 BC 敏感性 —— **仅 BC 一项就把首折从 3 m 推到 4 m**（2026-10-15）

> **性质**：**敏感性/诊断结果**（**不**改 V1 判定；V1 仍为 ❌【口径不可复现】）。
> **动机**：`PLANTAR_CONTACT_FIX_20261015.md` §9.4 —— 两步法未收口，改走「把 `hold_spring_k`
> 限定为**诊断口径**」这条路，用**同一条 `s2_thums_sweep` 流水线**量化跖面 BC 对 V1 的影响。
> **脚本**：`scripts/opensim_fe/s2_thums_sweep.py`

---

## §0 · 结论（一句话）

**在同一网格、同一流水线、同一载荷链、同一判据下，只把跖面 BC 从 `fixed`（三向全固定）换成
`contact`（真实接触 + 弱 hold 弹簧 k=100 N/mm），解剖配准网格的首折高度就从 3 m 变成 4 m**
（σ 系统性下降 **−8% ~ −16%**）。

**⇒ 跖面 BC 是 V1 的一个「口径级」自由度，其影响量级（×~2 in 高度）与 V1 的目标带
（6–9 m）同阶。** 这与 `V1V2_CALIBER_RULING.md` 的 ❌【口径不可复现】是一致的、且更硬。

---

## §1 · 数据（阈值 150 MPa；`gauge_max`，ℓ = 4 mm，域 `calcaneus`/CORT，act = 0）

### 1.1 解剖配准网格（**V1 现行口径所用的网格**）

| h (m) | `fixed` | `contact + hold k=100` | 比 | 越过 150？ |
|---:|---:|---:|---:|---|
| 1 | 114.8806 | 105.6015 | 0.919 | — / — |
| 2 | 149.8816 | 131.7919 | 0.879 | — / — |
| **3** | **156.6101** | 140.4150 | 0.897 | **✔ / —** |
| **4** | 168.6184 | **151.3479** | 0.898 | ✔ / **✔** |
| 5 | 186.7564 | 163.5426 | 0.876 | ✔ / ✔ |
| 6 | 208.5384 | 175.6669 | 0.842 | ✔ / ✔ |
| **首折** | **3 m** | **4 m** | | |
| 内插首折 | **2.02 m** | **≈3.96 m** | **×1.96** | |

### 1.2 `calcnframe` 网格（对照；**不是** V1 的口径）

| h (m) | `fixed` | `contact + hold k=100` | 比 |
|---:|---:|---:|---:|
| 1 | 258.4155 | 285.6689 | **1.106** |
| 2 | 336.9152 | 357.9889 | 1.063 |
| 3 | 352.0211 | 382.0972 | 1.085 |
| 4 | 379.0165 | 412.7406 | 1.089 |
| 5 | 419.5998 | 446.8582 | 1.065 |
| 6 | 467.8925 | 481.6212 | 1.029 |
| 8 | 570.4416 | 548.8225 | **0.962** |

**⇒ 两条网格上比值**符号相反**（anat 恒 <1、calcnframe 恒 >1，除 h=8）**
⇒ **跖面 BC 的效应与网格/坐标系强耦合** ⇒ 又一个「口径」自由度。

---

## §2 · 判据与边界（**必须与数字并列引用**）

**✅ 可主张**：

1. **同流水线内** `fixed → contact(+hold k=100)` 使解剖网格首折 **3 m → 4 m**、σ **−8%~−16%**。
2. `hold_spring_k = 100 N/mm` 是**消刚体模态的数值化装置**，**非物理接触刚度**；
   已在 k ∈ [10, 1e5] 上验过 σ 对 k 的**敏感性**（`PLANTAR_CONTACT_FIX_20261015.md` §3.2：
   h=2 下 0.987→1.000 单调趋向 `fixed`）⇒ **k 越大越像 `fixed`**。
3. 故本件给的是**区间**而非点：**首折 ∈ [3 m, 4 m]**（两端分别是 k→0 的纯接触倾向与
   k→∞ 的 `fixed` 极限的**上/下界方向**，**不做外推**）。

**❌ 不可主张**：

1. ❌ 不能把它当「V1 的新实测值」—— V1 已是 ❌【口径不可复现】（`V1V2_CALIBER_RULING.md`），
   本件**只增加一条口径敏感性的证据**。
2. ❌ 不能说「真实跖面 BC 使骨折风险降低 8–16%」—— `contact + hold` 是**近似**，
   且 h=1 下 calcnframe 网格的比值是 **+10.6%**（符号相反）。
3. ❌ 不能引「跖面接触 BC 下 V1 = 4 m」而不带「hold k=100、anat 网格、gauge_max、ℓ=4mm」
   这四个口径标签（＝闸门 9）。

---

## §3 · 为什么这**加强**了 V1 的 ❌【口径不可复现】（而不是削弱）

`V1V2_CALIBER_RULING.md` 的论据是「**统计量口径**（`median`/`max`/…）跨度太大」。
本件补上**第二个自由度：边界条件**：

| 自由度 | 跨度 | 来源 |
|---|---|---|
| **统计量口径** | `gauge_max` 2 m ↔ 矩阵 `median` 13 m | 裁定件 §1.2 |
| **跖面 BC**（本件） | 首折 **3 m ↔ 4 m**（anat） | 本件 §1.1 |
| 网格/配准 | `anat` 3 m ↔ `calcnframe` 1 m | 本件 §1.1/§1.2 |

**⇒ 三个自由度各自都能把 V1 的答案动 1 m 量级甚至更多** ⇒ 「口径不可复现」这一判定
**由本件进一步支持**，且**给出量级**。

---

## §4 · 溯源与复现

| 项 | 位置 |
|---|---|
| 脚本 | `scripts/opensim_fe/s2_thums_sweep.py` |
| 本件数据（anat） | `results/opensim_fe/s2_thums_height_sweep_anat_contact_hold100.{json,csv}` |
| 本件数据（calcnframe） | `results/opensim_fe/s2_thums_height_sweep_contact_hold100.{json,csv}` |
| 基线（anat / fixed） | `results/opensim_fe/s2_thums_height_sweep_anat.json` |
| 基线（calcnframe / fixed） | `results/opensim_fe/s2_thums_height_sweep.json` |
| 自动报告 | `results/opensim_fe/S2_THUMS_REPORT_ANAT_CONTACT_HOLD100.md`、`S2_THUMS_REPORT_CONTACT_HOLD100.md` |
| 接触/`hold` 口径来源 | `PLANTAR_CONTACT_FIX_20261015.md` §3（探针）、§9（bisection） |
| V1 口径裁定 | `results/opensim_fe/V1V2_CALIBER_RULING.md` |

**一键复现**：

```
.venv\Scripts\python.exe scripts\opensim_fe\s2_thums_sweep.py ^
  --mesh temp\opensim_fe\thums_calcaneus\calcaneus_r_anatframe.npz ^
  --heights 1 2 3 4 5 6 --acts 0 ^
  --plantar-bc contact --contact-hold-spring-k 100 ^
  --tag _anat_contact_hold100
```

**⚠️ 已知未解**：探针（`plantar_contact_fix_probe.py`）同一网格同一高度下
`ref_fixed` 给 **362.56 MPa**（h=2、calcnframe），而本流水线给 **336.92 MPa** —— **相差 7.6% 未归因**。
本件的**所有比较都在本流水线内**，故该差异不影响 §1 的**同流水线比值**，但**两条路不能互引**。
