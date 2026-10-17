# 冻结 BC2 数据保全 —— 清单与恢复路径

> **日期**：2026-10-13 ｜ **机器可读**：`FROZEN_BC2_MANIFEST.json`
> **本地路径（不入库）**：`temp/fe_out/fracture_matrix_bc2/`
> **体量**：83 文件 / **51.59 MB**

---

## 1. 为什么它不在 git 里（**政策 vs 实测的冲突**）

> `temp/fe_out/fracture_matrix_bc2/` 是 4 个沿用部位判定的唯一数据副本。
> 但 `.gitignore:72` 明确把 `*.feb`/`*.xplt`/`*.hdf5`/`*.npz` 与 `temp/` 归为
> 「**Regenerable simulation artefacts (safe to delete — reproduced by scripts)**」，
> 理由是仓库**对外发布**、只收 code + docs + **可再生** results。

> **⇒ 不可入库（会违反该发布政策）；但「可再生」这个前提已被实测证伪**（同条件重跑差 **6–22%**）。
> **⇒ 本件的作用：把「静默丢失」变成「可检测丢失 + 可恢复」。**

## 2. 冻结后确有一次代码漂移（**候选，非已证**）

冻结登记日 **2026-10-09**；此后以下模块被改动 ⇒ 用当前代码重跑得到的是**另一套 FE**：

| 模块 | 冻结期最后改动 |
|---|---|
| `scripts/opensim_fe/bone_feb.py` | 37f782bd45d8 (2026-10-09) |
| `scripts/opensim_fe/fracture_matrix.py` | 392d5f0e4789 (2026-10-06) |
| `scripts/opensim_fe/rerun_matrix_bc2.py` | 37f782bd45d8 (2026-10-09) |
| `scripts/opensim_fe/paper_alignment.py` | 392d5f0e4789 (2026-10-06) |
| `scripts/opensim_fe/s4_alignment_bc2.py` | 37f782bd45d8 (2026-10-09) |

**冻结后被改动的提交**（恢复时须回退这些）：

- `scripts/opensim_fe/bone_feb.py` `30d67a2fa0dc` (2026-10-14) — F 路线 W4-a: BC 稳定性门槛达成(50/50 高度全域可解, 此前 13/50) —— 但翻出可识别性问题: A2-股骨 不可由修 BC 解决 —— (改动, 加法式默认关断) bone_feb.py 新增 bottom_fix_axial_only 参数(底部夹持只约束载荷轴 DOF, 释放横向) + 复用既有 lateral_support(默认 None); (稳定性扫描, 最难高度 h=50/57.9kN) 6 变体: A0_baseline(现状全3DOF+无侧向)=FAIL, A1_axial_only(仅轴向+无侧向)=FAIL, A2_axial_lat90=OK(ends 13.99/shaft 81.74), A3_full_lat90=OK(ends 13.99/shaft 81.74 与A2逐位相同!), A4_axial_lat80=OK(17.14/82.75), A5_axial_latfull=OK(13.07/81.50) => A2 与 A3 逐位相同 => 发散根因是侧向屈曲而非固定端弯矩, 只要侧向支撑覆盖>=80%自由长则底部夹持DOF无关(与 §9.14 腓骨'导向须覆盖自由长>=80%'一致); (全50高度, 侧向90%, 干净重跑——首跑混入上轮无侧向缓存 h=1-7/10/12/13/16, 已把缓存移开为 femur_r_unbraced 后重跑): 50/50 全收敛, shaft 13.90(9.78kN)->81.74(57.94kN), ends 2.34->13.99, sigma(h) 光滑单调无跳变, 判定=全不骨折(all_zero, shaft 始终<220) => 门槛达成; ★核心发现: 三种 BC 给三种分类 — 现状无侧向=all_one; 无侧向+逐高度=屈曲分支(h=10 达 781MPa, h=16 达 531MPa, h=43/44 仅 137-160MPa)分支跳变无定义; 侧向90%+逐高度=all_zero => 股骨在 20-58kN 处屈曲载荷附近, 静力解在分支间跳变, sigma(h) 在现 BC 族下不是良定义的量, 分类由侧向支撑强度这一一阶参数控制 => A2-股骨 不可由'修 BC 修数值'解决, 需先给'坠落中股骨受多少侧向支撑'物理锚定; 诚实边界: 只测股骨/覆盖率是扫描非推导/未与独立参照核对/A2 现行判定 4/6 不变; 建议下一步: 给侧向支撑找物理锚定(可核 g9_fibula_fe_support.py 的 >=80% 论证据否移植), 或锚定后仍不可定义则记【判据核心输入未定】(同 V8b/A4 待遇); 新增 results/opensim_fe/FE_BC_STAB.md + FE_BC_STAB.json + scripts/opensim_fe/fe_bc_stab.py
- `scripts/opensim_fe/bone_feb.py` `03858010c879` (2026-10-11) — feat(线2 OpenSim_FE): 结构件 G1/G3/G4/G8/G9 + [Y25] OR 口径重审 + 文档同步
- `scripts/opensim_fe/fracture_matrix.py` `03858010c879` (2026-10-11) — feat(线2 OpenSim_FE): 结构件 G1/G3/G4/G8/G9 + [Y25] OR 口径重审 + 文档同步

**⚠️ 6–22% 偏差的成因仍未定案 —— 不要把「代码漂移」当结论：**
`03858010c879`（2026-10-11）早前已逐字核查为**对 `bone_feb.py` 与 `fracture_matrix.py` 的纯加法**
（见 `FE_BC_PROBE.md` 排障记录）⇒ **纯加法本身不应改变数值**。
**⇒ 剩余候选（均未验）**：网格数据变化 / 材料参数来源 / FEBio 求解器版本 / 求解随机性。
**⇒ 因此第 3 节第 1 条的回退只是「已知的必要条件」，不保证复现出冻结值。**

## 3. 恢复路径（**任一条即可**）

1. **回到冻结期代码**：上表提交哈希 ⇒ `git checkout <commit> -- <modules>` ⇒ 重跑
2. **逐字节校验副本**：本清单含**逐文件 sha256** ⇒ 任何外部副本都能验证一致
3. **降级使用**：若已丢失 ⇒ A2 那 4 个部位**退回「不可复现数据」**状态，须登记

## 4. 建议（**待拍板**）

- 把 `temp/fe_out/fracture_matrix_bc2/` 另存到**仓库外**（或私有仓库），本清单入 git 作校验基准
- 恢复后，**重跑并与本清单逐文件比对**，差异本身是有价值的证据（可定位是哪一段代码漂移造成）

## 5. 溯源

| 项 | 来源 |
|---|---|
| 不可复现实测 | `FE_BC_PROBE.md`（`V0_control` 未复现冻结值，差 6–22%） |
| 代码漂移核查 | `git log` on 冻结期相关模块 |
| 「沿用 4 骨」的判据 | `A2_RECOMPUTE_SCOPE.md` |
| A2 新计数 | `A2_RECOMPUTE_RESULT.md` |

