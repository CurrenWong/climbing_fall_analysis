# 足部 FE 模型盘点 —— 「现在做了完整脚底吗？」

> **触发**：用户问「现在 FE 做了不光根骨，还有完整脚底吗」（2026-10-16）。
> **性质**：**取证件**（查实际文件/摘要，不跑求解）。**不改任何判定**。
> **本件新增检查**：`temp/v1_wholefoot.py`、`temp/v1_footpart.py`、`temp/parts_sum.py`（只读）。

---

## §0 · 三句答案

1. **V1 现行：没有。** 就是**单骨跟骨**（1544 节点）。
2. **踝子线：也没有前足。** 是「踝 + 后足」（tibia / fibula / **talus** / **calcaneus** / **navicular** + 软骨 + 韧带），**无跖骨/趾骨**。
3. **但「整只脚」的骨块其实都在手上，且都已提取/网格化** —— 只是**从未组装、从未求解**。
   另外**旧管线曾用过"整足"网格**（`r_foot.vtp`，**首折 26.5 m**），Option C 时被换成单骨跟骨。

---

## §1 · 现在实际在求解的 FE（三个）

| # | 模型 | 组成 | 状态 |
|---|---|---|---|
| **1** | **V1 现行**（`calcaneus_r_anatframe.npz`） | **单骨跟骨**：1544 节点，CORT hex8 **634** + SPON tet4 **3323**；面：距下关节面 2495.33 mm²、跖面 2953.78 mm²、achilles | ✅ 在用（V1/F 线/矩阵） |
| **2** | **踝子线**（`thums_lowerlimb.npz`） | `tibia(+spon+marrow) / fibula(+spon) / talus(+spon) / calcaneus(+spon) / navicular(+cort+spon)` **+ 软骨层**（`tibiotalar_*`、`subtalar_*`）**+ 韧带**（`ligaments.npz`） | ⚠️ 在用但**量程 ≤519 N**（真实 25.3 kN，差 49×）；**无前足** |
| **3** | `joint_bc="contact"` 变体 | 上述 (1) **+ `E=200 GPa` 刚性帽**（数值刚性板，**不是**距骨） | ⚠️ opt-in；补跑未复现（`V1_LOADSHARE_VERIFICATION.md` §3.2） |

---

## §2 · 手上「已经有」但**没组装**的足部骨块（THUMS AM50，同一 deck）

`temp/opensim_fe/parts/` —— 全部来自 **同一个** `model\AM50_V71_Occupant\main_THUMS_AM50_V71.k`，
帧 = **`THUMS global, mm`**：

| part | 节点 | 单元 | 三角片 | 含哪些骨（PID） |
|---|---:|---:|---:|---|
| **`foot_r`** | **2662** | **1314** | 4152 | **48 PID = 24 块骨**：舟骨 / 骰骨 / 内·中·外楔骨 / **5 跖骨** / **5 近节 + 4 中节 + 5 远节趾骨** |
| **`talus_r`** | 916 | 2160 | 796 | 距骨（spon + cort，2 PID） |
| `tibia_r` | 7245 | 23550 | 4852 | 胫骨 |
| `lumbar_L3` | 2203 | 8792 | 2670 | L3 |
| `skull_r_parietal` | 2952 | 1430 | 6068 | 顶骨 |
| `_verify_parts` | 18537 | 59830 | 19164 | 验证合并体 |

**⇒ 拼一只「完整脚」的全部骨块都在**：
`calcaneus(81001200/1300)` ＋ `talus(81001000/1100)` ＋ **`foot_r`（81001400 → 81003701，中足+前足）**
—— **同一 THUMS deck、同一全局坐标系** ⇒ **几何上可直接拼**。

**⚠️ 但（三条必须并报）**
1. **`foot_r` 从未被任何求解脚本使用** —— 全仓引用它的只有**画图/可视化**脚本：
   `scripts/opensim_fe/make_fracture_figures.py:33`（+ `temp/` 下的临时可视化草稿）。
   **在 `thums_feb.py` 里搜 `foot_r` 零命中**（只有一处注释提到 `calcaneus<->metatarsal`
   跖腱膜两端连接做不到，`:244`）⇒ **确认：从未组装、从未求解。**
2. **`foot_r` 只有 `spon` 有单元**：摘要里 **`*_CORT` 的 `n_elements` 全为 0**
   （24 块骨的 cortical PID 均 0 单元）⇒ **皮质壳不在这个 part 里**（THUMS 里皮质走 `*ELEMENT_SHELL`，
   而本 part 的提取只取了 solid）⇒ **拼起来会缺皮质壳**（与中轴骨 `T6` 同款的坑）。
3. **`foot_r` 摘要报的 `enclosed_volume_cm3 = 1245.1` 明显偏大**（中足+前足的**骨**体积应远小于此；
   同目录 `lumbar_L3` 报 26.2 cm³ 属正常）⇒ 该值是**表面包围体积**（含骨间空隙/外壳），
   **不是骨体积**，**引用时须注明**。

---

## §3 · 历史：**曾经的「整足」FE**（这就是 V1 口径混乱的来源之一）

`scripts/opensim_fe/s2_thums_sweep.py:421-424` 逐字：

> **目标几何替换**：旧 S1–S3 的 `"calcaneus"` 域**实际是 `r_foot.vtp` 的整足网格（121.7 cm³）**，
> 不是单一跟骨；本模型的 `"calcaneus"` 是真实 AM50 跟骨。两者数值**不可直接比较**，只能看趋势/方向。

`scripts/opensim_fe/thums_feb.py:3-4` 逐字：

> Option C …: **swap the crude whole-foot deformable geometry (`r_foot.vtp`)** for the
> registered THUMS AM50 right calcaneus.

**⇒ 所以：**

| | 旧管线 S1–S3 | 现行（Option C 之后） |
|---|---|---|
| 几何 | **`r_foot.vtp` 整足**（OpenSim Rajagopal2015 表面，**1000 点**，包围体积 **121.7 cm³**） | **单骨跟骨**（THUMS，双域） |
| 质地 | **"crude" 粗网格**（可视化面三角化） | 真实 AM50 跟骨几何 |
| 首折 | **26.5 m**（`s2_thums_sweep.OLD_FIRST_FRACTURE_M`） | **2 m**（固定口径）/ 4 m / 5.43 m |
| 现状 | **退役** | ✅ 在用 |

> ⚠️ `r_foot.vtp` 的**命名有歧义**：`registration.py:12` 说它是
> 「the **`calcn_r` body** geometry in the Rajagopal model」——
> 即**文件名是 r_foot，但它挂在 `calcn_r` body 上、覆盖整只脚**。两处说法本质一致。

**⇒ 教训**：V1 的「多口径 1/2/3/13/26.5 m」里，**26.5 m 就是「整足粗网格」那条口径**
（`V1V2_CALIBER_RULING.md`）。**这也是 V1 被判 ❌【口径不可复现】的一部分原因。**

---

## §4 · 还有一套**完整足部几何**（只有表面，没网格）

`temp/ankle_geom/` —— **IITD 足踝模型**（`IITD_Ankle_foot_model_NLSRclass.OSIM`）的导出：

| 文件 | 说明 |
|---|---|
| `calcaneus_r.stl` / `calcaneus_l.stl` | 跟骨 |
| `talus_r.stl` / `talus_l.stl` | 距骨 |
| **`midfoot_r.stl`**（5.9 MB） | **中足** |
| **`metatarsal_{index,middle,ring,thumb,tiny}_r2.stl`** | **5 块跖骨** |
| **`toes_{index,middle,ring,thumb,tiny}_r3.stl`** | **5 组趾骨** |
| `talus_cartilage_r.stl`、`tibial_talus_cartilage_r.stl` 等 | 软骨 |
| `femur_cartilage_*`、`*_meniscus_r`、`OK_tibia.obj` | 膝（同批导出） |

**⇒ 这是**第二套完整的足部几何**（含骨 + 软骨），但**只有 STL 表面、没有 FE 网格、从未求解**。

---

## §5 · 所以「完整脚底」的可达性（诚实评估）

| 要做「完整脚底 FE」需要 | 现状 |
|---|---|
| 足部各骨的**几何** | ✅ **有两套**：THUMS（跟骨/距骨/`foot_r`，同坐标系）＋ IITD（全套含软骨） |
| 中足/前足的**体网格** | ⚠️ THUMS `foot_r` **已有 hex8**（1314），但**只有松质、无皮质壳**；IITD **无网格** |
| 各骨**之间的关节**（距下、跟骰、楔骰、跖跗、跖趾…） | ❌ **完全没有**（现有最多做到距下 + 胫距） |
| **跟垫 / 足底软组织** | ❌ **没有**（G5 的 1D 跟垫模型未接进 FE） |
| **足底接地 BC** | ⚠️ 只有「单骨跖面焊死」这类；**没有真正的足底接触** |
| **收敛能力** | ⚠️ 项目记录：无纯实体重网格器、无网格收敛证据（`R2 风险`），且 `foot_r` 皮质壳缺失会重现 `T6` 那类问题 |

**⇒ 结论：材料（几何）基本齐，**但「组装 + 关节 + 软组织 + 收敛」四件事都没做**。**
现有能做的最多到 **`(1) 跟骨 + (2) 距骨 + `foot_r`（中足+前足）`** 三块拼合的**骨性**脚 ——
**没有关节连接、没有跟垫、没有皮质壳**。

---

## §6 · 溯源

| 项 | 位置 |
|---|---|
| V1 单骨模型 | `temp/opensim_fe/thums_calcaneus/calcaneus_r_anatframe.npz`（1544 节点） |
| 载荷路径与跖面 BC | `results/opensim_fe/V1_LOAD_PATH_FACTS.md` |
| 刚性帽（≠距骨） | `src/climbing/coupling/joint_contact.py:47`（`DEFAULT_CAP_E_MPA = 200_000.0`） |
| 多骨踝+后足 | `scripts/ankle_fe/ankle_lig_feb.py:13-17`；`temp/thums/thums_lowerlimb.npz` |
| **THUMS `foot_r`（中足+前足，hex8）** | `temp/opensim_fe/parts/foot_r/{summary.json, foot_r_elements.npz, foot_r_volume.vtk}` |
| `foot_r` 仅用于画图 | `scripts/opensim_fe/make_fracture_figures.py:33` |
| **旧「整足」管线与 26.5 m** | `scripts/opensim_fe/s2_thums_sweep.py:421-424`（`OLD_FIRST_FRACTURE_M = 26.5`）；`thums_feb.py:3-4` |
| `r_foot.vtp` 身份（calcn_r body） | `scripts/opensim_fe/registration.py:12,133,536` |
| IITD 完整足几何（仅 STL） | `temp/ankle_geom/` |
| V1 口径裁定 | `results/opensim_fe/V1V2_CALIBER_RULING.md` |
