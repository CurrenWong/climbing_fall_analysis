# PENALTY_FIX_RESEARCH — 攻踝 tied/关节接触罚因子上限：可行方案调研

**日期**: 2026-10-09　**状态**: 调研完成 + **E1 已实跑（结论已改写，见 §0b/§7）**　**触发**: 用户拍板走路线 b
**关联**: `docs/notes/TIE_DIVERGE_TODO.md`（P1d2b 链）· `[LRN-20261005-056]` · `ANKLE_ROADMAP_TODO.md` A3

---

## 0b. E1 实测结果（2026-10-09，取代本文档初版对 H1/H2 的判断）

**驱动**: `scripts/ankle_fe/penalty_auto0_probe.py`（单变量：deck 里 `auto_penalty` 1→0，行级断言保证只动该位）
＋`temp/pyfebio_demo/k_joint.py`（加 opt-in `AUTO0=1`，沿用**原 δ 口径**保证苹果对苹果）。

### (1) 收敛行为：2×2 因子 @ p=2
| 单元 | auto_penalty | dt(=1/TIME_STEPS) | 失败于 t | 发散次数 | 耗时 |
|---|---|---|---|---|---|
| A | 1 | 0.05 | ~0.002 | 11 | 0:30 |
| C | 1 | 0.005 | **0.00045** | 11 | 0:29 |
| B | **0** | 0.05 | **0.132** | 17 | 8:34 |
| D | **0** | 0.005 | **0.133** | 82 | 12:32 |

- p=1：两臂**都收敛**（auto=1 rc=0 17.6s；auto=0 rc=0）⇒ 关 auto 不破坏已知好点。
- ⇒ **可收敛上限在 p∈(1, 2]**，与 auto 开关、dt 都无关；**dt 杠杆被证伪**（缩小 10× 无效，只多烧一倍时间）。

### (2) 刚度：同一组 penalty、只翻 auto_penalty
| PENALTY | k（auto=1，历史表，本会话逐位复现） | k（auto=0，项目铁律口径） | F(auto=1) | F(auto=0) |
|---|---|---|---|---|
| 0.001 | 49.1 | **0.2** | 6.99 N | 0.03 N |
| 0.01 | 150.9 | **1.7** | 50.65 N | 0.25 N |
| 0.1 | 451.3 | **18.1** | 220.08 N | 2.36 N |
| 1.0 | **880.0** | **70.5** | 419.56 N | 17.54 N |

### (3) 三条改写原判断的结论
1. **k 表的"物理含义"作废**：auto=1 下 `<penalty>` 是 scale factor，49→880 不是 k 的收敛序列，**"k 单调趋近 1059、p≥2 挡住验收"这条叙事不成立**。
2. **但 auto=1 不是病灶，而是"通往同一刚度的更好条件化路径"**：auto=1 @ p=1 的等效 k≈**880 N/mm**（距 Anderson 1059 仅 −17%）**且收敛**；auto=0 要靠线性外推 **p≈15** 才够到同一刚度，而 **p=2 就已发散**。
   ⇒ 真正的堵点不是"罚值非法"，而是：**要用原始罚参数表达物理接触刚度，就必须给到 ~15–30× 于周围网格的刚度比，而 penalty 接触恰在此区间崩溃。**
3. **生产链同样带 auto=1**：`temp/thums/febio_lig/ankle_lig_WK.feb` 三个接触（tibiotalar sliding-elastic + tie_tibiotalar_tibia/talus 两个 tied-elastic）**全部 `auto_penalty=1`** ⇒ 436.8 N 是在此语义下得到的；`ankle_lig_feb.py:1828` 只对 AUGLAG 分支置 auto=0。

### (4) 靶标口径差被量化（E6 的硬数据）
脚本自算列 `E·A/t` 用**实测闭合面积**（p=1: A=344mm²）+ t=3.0 ⇒ **1377 N/mm**；auto=0 档更 soft 时该列升到 1980–2327。
而文档靶标 1059 用的是 **Anderson 的 300mm² / t=3.4** ⇒ 同一次实现在两套口径下差 **−36% ~ −55%**。

### (5) 诊断开关
`<output_negative_jacobians>` 在 **FEBio 4.13 是 unrecognized tag**（1.9KB 日志秒退，表现为"无负 J"的**假阴性**）⇒ 想定位翻转单元需另找标签或用 `-g2` 存失败状态。

## 0d. E8 实测：AUGLAG 增广整定（2026-10-09，取代 §0c 中"原始罚参数路线不可能成功"的悲观收口）

**驱动**: `temp/pyfebio_demo/k_joint_auglag.py`（k_joint 的变体，`LAUGON=AUGLAG` + 可配 `PENS`/`FEB_TIMEOUT`）。
AUGLAG 分支下 `build_feb.py:696` 自动写 `auto_penalty=0` ⇒ 天然符合项目铁律，无需 patch。

| aug_penalty | 状态 | F (N) | δmax (mm) | k=F/δ (N/mm) | 闭合面积 (mm²) | vs Anderson 1059 |
|---|---|---|---|---|---|---|
| 2 | 收敛 | 339.59 | 0.488 | 695 | 343 | −34% |
| 8 | 收敛 | 405.85 | 0.455 | 893 | 311 | −16% |
| **16** | **收敛** | **426.37** | 0.451 | **944.9** | 304 | **−11%** |
| 32 | **超时 >30 min** | — | — | — | — | 上限在此 |

### 结论
1. **AUGLAG + `auto_penalty=0` + `aug_penalty≈16` 给出 k≈945 N/mm，收敛**，落在 Anderson 1059 的 **−11%**；
   本模型几何自算 `E·A/t`（A=304mm², t=3.0）= 1215 ⇒ 差 −22%。
2. **接触刚度上限 ≈950–1000 N/mm，与表述方式无关**：PENALTY+auto=1@p=1 → 880；AUGLAG@16 → 945；
   再往上（AUGLAG@32 / PENALTY p≥2）不是发散就是超时 ⇒ **残余 −11% 是表述上限，不是标定误差**。
3. 🔑 **力值对表述稳健（retro-validation）**：同一 smoke deck 上 `PENALTY+auto=1@p=1` = **419.56 N**
   vs `AUGLAG+auto=0@16` = **426.37 N**，差 **1.6%** ⇒ **生产链的 436.8 N 不是"罚参数语义"的产物**；
   作废的只是「`<penalty>` 就是 k」这个**解读**，不是力值本身。
4. **新基线候选**：`AUGLAG / auto_penalty=0 / aug_penalty=16`（合规、可解释、k≈945、收敛、快）。
   已排除项：**PENALTY 原始罚参数路线 p≥2 不可用**（auto=0 时 p=2 崩于 t=0.132；auto=1 时 p=2 秒崩）。
   ⚠️ 注意：subtalar 历史上因"AUGLAG 更差（0.46→0.65）"选了 PENALTY（`ANKLE_ROADMAP_TODO.md:329`），
   该结论产生于**未整定增广参数**时 ⇒ 需按 aug_penalty 扫描**复验**后才可翻案。

## 0e. R13 实测：subtalar「AUGLAG 更差」——维持原判，但归因已精确化（2026-10-09）

**驱动**: `temp/pyfebio_demo/r13_subtalar_auglag.py`（只改 deck 里 **subtalar 那一块**，tibiotalar 保持生产原样）
＋`r13b_stall_probe.py`（步长/reform 两因子）＋`r13c_auglag_tol.py`（增广容差/node_reloc）。
配置 = A1+A2 生产档：`SUBTALAR=1 / TALUS_SUPPORT=spring / TALUS_K=500 / MAX_UPS=0 / CART_DIR=cart_m1 / DISP_MM=+1.0`。

| case | rc/NORMAL | F_top_z (N) | fail_t | 负J | 残余 max gap (mm) |
|---|---|---|---|---|---|
| **base：PENALTY + auto=1（生产现状）** | 0 / ✅ | **449.6312** | — | — | — |
| AUGLAG auto=0 pen=1 | 1 | 423.16 | — | — | — |
| AUGLAG pen=2 / 8 / 16 | ✗ | 390.5 / 377.6 / 372.7 | 0.829 / 0.773 / 0.758 | **0** | — |
| + `tolerance` 0.02→**0.2** | ✗ | 396.33 | **0.841** ↑ | 0 | 0.110 |
| + `tolerance`→1.0 | ✗ | 387.05 | 0.831 | 0 | 0.135 |
| + `time_steps` 20→40 / `max_refs` 50→100 / 两者 | ✗ | 373~378 | **0.773（逐位不变）** | 0 | — |
| + `node_reloc=1` | ✗ | 396.3334 | 0.841288（与 tol0.2 **逐位相同**） | 0 | 0.110 |

**基线逐位复现** roadmap 的 `+1.0 → 449.63 N` ✅（同日同链可复现性确认）。

### 结论
1. **原决策维持**（subtalar 保持 `PENALTY + auto=1`），但**理由必须改写**：
   ~~"AUGLAG 表述更差"~~ ⇒ 真因是**两条叠加**且都不是罚参数能解决的：
   - **(i) 增广容差过硬**：deck 硬写 `<tolerance>0.02</tolerance>`（`ankle_lig_feb.py:1830`）＝日志里的
     `D multiplier REQUIRED 2.000000e-02`；放松到 0.2 把停滞点从 t=0.773 推到 **0.841**，再放松到 1.0 无收益。
   - **(ii) 界面残余间隙合不上**：`maximum gap` 0.043 → 0.110 → 0.135 mm（**随容差放松反而变大**），
     且 `node_reloc=1` **实测零影响**（与 W1 逐位相同）⇒ 非节点重定位可解，属几何/咬合问题
     （与 `ankle_lig_feb.py:109` 记的"t=0 未咬合（中位间隙 0.97mm）⇒ 接触病态"同族）。
2. **罚参数不是关键**：p=2/8/16 只把停滞点 0.829→0.773→0.758，**罚越硬停得越早**（stiffness-contrast）。
3. **逐项排除清单**（subtalar）：步长细分 ✗｜reform 预算 ✗｜node_reloc ✗｜罚值提高 ✗｜增广容差部分有效但不足 ✗。
4. ⚠️ **tag 名教训（第二次）**：`<max_aug_iterations>` 在 FEBio 4.13 **不是合法 tag**（秒退）；
   与 `<output_negative_jacobians>` 同类 ⇒ **FEBio tag 必须查手册/实测，不能猜**。

## 0c. 由此改写的下一步（E8 已完成见 §0d；E1' 已被 §0d 取代；R13 已完成见 §0e）

| 轮 | 内容 | 理由 |
|---|---|---|
| **E1'** ⭐ | **auto=1 下改变 auto 估计的"依据"**：auto_k = f(primary surface 刚度, 单元尺寸) ⇒ 试 cartilage E↑ 或接触区网格减半，看 k 能否从 880 走向 1059–1377 且保持收敛 | 这是**唯一能把 k 顶到物理值而不用原始罚参数爆掉**的杠杆（项目从未用过） |
| **E7** | 生产链 `ankle_lig_feb.py` 全接触 auto=0 对照（436.8 N 会变成多少） | 决定"绝对力值是否需要重标定"这条上游问题 |
| **E8** | AUGLAG 增广参数整定（`aug_penalty`/`max_aug_iterations`/`update_tolerance`）作为"近刚性接触"的标准解法 | p 路线走到 ~15 必崩，换表述才有机会 |
| **E6** | 口径修复：k 改 `F·t/(A·δ̄)`、靶标按实测 A/t 重算；文档替换式更正 | 无论后续结果都该做 |

---

## 0. 一句话结论

**「罚因子上限」这个堵点的定义本身可能是错的**：整张 `k=49/151/451/880 N/mm` 表是在 **`auto_penalty=1`** 下扫出来的，
而 FEBio 官方语义下 `auto_penalty=1` 时 `<penalty>` 只是**乘性 scale factor**（基值由 primary surface 刚度 + 网格尺度自动估计），
官方并明确警告该模式可能生成**过大的罚而翻转单元**。⇒ 先用 `auto_penalty=0` 重扫，再谈"抬上限"；
同时发现 **dt 从未被与罚联合扫描**（失败恰好发生在第一个满步 t=0.05）。

---

## 1. 现状取证（R1，均为文件实测）

| 事实 | 证据 |
|---|---|
| FEBio **4.13.0**，求解器 STATIC，`dtmin=1e-4 / dtmax=0.05 / max_retries=10 / opt_iter=10` | `temp/thums/febio/_fast/t_p8.log:20`；`scripts/ankle_fe/ankle_lig_feb.py:1863-1889` |
| 罚扫描驱动 = `pen_converge.py`/`k_joint.py` → 调 **`build_feb.py`**，配置 `TIE=0, LAUGON=PENALTY, TIME_STEPS=20, DRIVER=disp, DISP_MM=0.1472` | `temp/pyfebio_demo/pen_converge.py:14-15`、`k_joint.py:79-81` |
| `TIME_STEPS=20` ⇒ `step_size=1/20=0.05` = `dtmax`，**第一步即满步 50 μm** | 同上 + `ankle_lig_feb.py:1866,1887` |
| 产出 deck 实测含 **`<auto_penalty>1</auto_penalty>`**（`laugon=PENALTY` 分支自动开） | `temp/thums/febio/ankle_tibiotalar.feb:35629`；生成逻辑 `scripts/ankle_fe/build_feb.py:696` |
| 项目**自订铁律**：auto_penalty=1 会让 `<penalty>` 失效 ⇒ 罚因子扫描必须关掉 | `scripts/ankle_fe/ankle_lig_feb.py:127-128`（**扫描走的 build_feb.py 分支未执行此铁律**） |
| p≥2 失败形态 = **求解发散的下游症状**，非网格初始翻转 | `_fast/t_p8.log:215-275`：t=0.05、steps completed=0、位移范数 1.0e2→1.07e5→3.09e2→1.89e3（要求 9.6e-2）、两次 diverging/reform、末了才 "Negative jacobian detected." |
| 失败算例 ≠ 生产算例 | 扫描用 `build_feb.py`（仅 tibiotalar、TIE=0、位移驱动）；生产 436.8 N 来自 `ankle_lig_feb.py`（韧带+subtalar+tie+刚体距骨） |
| 另一条链（S5）已独立判定"负 J 是派生物" | `docs/s5-contact/S5_contact_udg.md:172,191,209`：`ls_check_jacobians=no(0)` 时无负 J 记录 |

## 2. 根因候选（排序，每条配可判别测试）

| # | 候选机制 | 判据 | 成本 |
|---|---|---|---|
| **H1**（高） | `auto_penalty=1` 使罚由 auto 基值×scale 决定；auto 基值按 primary（较硬侧）刚度+网格尺度估计，天然过大 ⇒ 高 scale 时翻转/发散 | `auto_penalty=0` 重扫 p∈{1e-3…16}：rc 与 k(p) 是否变样 | ~10 min |
| **H2**（高） | **dt 不足**：满步 50 μm + 高罚 = 接触冲击，Newton 首步即崩；`max_retries=10` 打满 | 同一 p 下 `TIME_STEPS` 20→200/500（dt 0.005/0.002）能否再抬 1–2 档 | ~15 min |
| **H3**（中） | primary surface 选择不当（较硬侧做 primary ⇒ auto 罚更大；手动罚时影响 search/投影） | 交换 surface_pair 主从 + `node_reloc=1`（本项目注释称其专治接触病态） | ~15 min |
| **H4**（中） | 无 continuation：高罚点若能"从 p=1 的解续起"则可跨过发散区 | `MULTISTEP` 分段罚：step1 p=1 → step2 p=2 → step3 p=4（FEBio step 内参数仅该 step 生效，见 `ankle_lig_feb.py:1412`） | 中 |
| **H5**（低） | AUGLAG 增广参数未调（初值扫描 2→315 N 说明增广本身没收敛） | 调 `aug_penalty` / `max_aug_iterations` / `update_tolerance` | 中 |
| **H6**（诊断） | 不知哪些单元翻转 | `-g2`（负 J 时存状态）+ `output negative jacobians=1`（打印全部翻转单元/高斯点），官方 webinar 推荐 | ~5 min |

**官方口径**（FEBio Webinar Contact Modeling Part 1/2 + McGill 文档）：
- "负 J ... failed to converge ... Max. nr of retries reached" = **跑到一半失败/变形过大**，官方处方是降载荷/降罚 ⇒ 对我们方向相反，必须换杠杆（dt / continuation / 关 auto / 换 primary）；
- auto penalty 以 primary surface（较硬侧）为准，"very large penalty ... going to invert elements"；
- auto 开启时用户给的 penalty 变成 scale factor。

## 3. ⚠️ 靶标/口径自洽性（R4：比罚上限更根本）

1. **k 的定义不对齐**：代码 `k = F_total / δmax`（`k_joint.py:100,117`），F 是全接触面反力积分、δmax 是单点最大闭合；
   而靶标 `E·A/t` 是**均匀压缩**刚度（应 `F·t/(A·δ̄)`）。用 δmax 作分母 ⇒ **880 系统性偏高**（若 δ̄/δmax≈0.6–0.8）。
2. **靶标不同源**：文档写 `12×300/3.4≈1059`（A=300、t=3.4 取自 **Anderson**），脚本自己算的是 `12·A1/3.0`（`k_joint.py:118`，
   A1=**本模型实测闭合面积**；加载后实测 **344 mm²**）。按实测几何：`12×344/3.0≈1376`、`12×344/3.4≈1214` ⇒ 880 的差距从 −17% 变 −36%/-27%。
3. ⇒ **"k 必须收敛到 1059"未必是本模型应收敛到的值**；正确判据 = 按本模型 A、t 重算的 k_phys + 统一 δ 口径。
4. 方向性提示：`TIE_DIVERGE` 记"真值可能高 1.5–2× ⇒ 436.8N 偏低"。若 k 定义偏高（见 1）与 F 偏低同向叠加，**该倍数可能高估** ⇒ 需实测而非推断。

## 4. 方案排序（建议执行序）

| 轮 | 内容 | 期望 | 成本 | 若失败的价值 |
|---|---|---|---|---|
| **E1** | `auto_penalty=0` 重扫 p∈{1e-3,1e-2,1e-1,1,2,4,8}（同 deck） | 拿到**真** k(p) 与真实可收敛上限；k 表或成立或作废 | ~10 min | k 表可信度结论（替换 🟠 定量基础） |
| **E2** | 在 E1 的最高可收敛 p 上，`TIME_STEPS` 20→200→500 | 抬上限 1–2 档 | ~15 min | dt-罚交互曲线（项目空白） |
| **E3** | `node_reloc=1` + 交换 primary/slave | 再抬档 / 降穿透 | ~15 min | 主从敏感性数据 |
| **E4** | MULTISTEP 罚 ramp（1→2→4 续解） | 跨过发散区 | 中 | continuation 可行性结论 |
| **E5** | AUGLAG 增广参数整定 | 换算法路径 | 中 | 备选算法可行性 |
| **E6** | 口径修复：k 改 `F·t/(A·δ̄)`、靶标按实测 A/t 重算；文档替换式更正 | 判据正确 | 零算例 | **无论 E1–E5 结果如何都该做** |

**正对照（判据必备）**：新增"接触刚度补丁测试"——均匀块压刚性板，E/ν 已知 ⇒ 解析 `k=E·A/t`；
扫 penalty 看 k 是否收敛到解析值 + 记录可收敛的最高 p（与 B1 的 BC patch test 同款纪律：`B12_PATCH_TEST_REPORT.md`）。
**无此正对照，"罚收没收敛"永远说不清**（SOUL：检查器本身要有已知正对照）。

## 5. 风险与红线

- **k 表在 E1 出结果前应视为「不可引用」**（deck 事实：`auto_penalty=1`）；但**暂不改写** 🟠 行，等 E1 实测后**替换式**更正（不叠加相反论断）。
- 扫描算例（简化 smoke deck）≠ 生产链（韧带+subtalar+tie）；E1–E3 的上限结论须在 `ankle_lig_feb.py` 生产链上**复验**才可写入结论。
- 全程沿用 `.venv\Scripts\python.exe`；FEBio 抢核 ⇒ 串行（历史 FEBio 慢 3×/超 600s）。
- 不动 pad 条件位恒等锁；不改 `436.8 N` 相关已发布数字（待 E1 判定后再联动）。

## 6. 未做/待办
- [ ] E1–E5 实跑（用户放行后）
- [ ] 生产链复验（E1 结论落 `ankle_lig_feb.py`）
- [ ] 接触刚度正对照脚本
- [ ] 文档替换式更正：`TIE_DIVERGE_TODO.md` 下一步 / 总体进度 🟠 行 / `[LRN-20261005-056]` 增补