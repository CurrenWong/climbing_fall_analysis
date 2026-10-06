# S5 接触收敛性：可靠区间与 G7 工况的不可解（旧 AUGLAG 配方下）（2026-10-06）

> 本文记录接触稳健化的**决定性与可复现**结论。它**取代**此前两条中间结论：
> ①「不存在既能收敛又传力的设置」；②「2.70× 是大应变平冲头压入物理，可作机制/量级陈述」。
> 新的结论是：**接触在中等压入区间稳健可信；G7 工况（ξ≈0.85）解不在收敛包线内且不可复现。**

> ### 🔄 2026-10-06 晚间补注 —— **本文标题/TL;DR 的「不可解」应读作「旧配方下不可解」**
>
> 本文结论**仅对当时使用的 OLD AUGLAG 配方**成立（该配方本身是真实、可复现的实验记录）。
> 此后两件事分别变了：
>
> | 事项 | 现状 |
> |---|---|
> | G7 **可复现性** | ✅ **已解除** —— `laugon="PENALTY"` + `penalty=0.1`（**保留 `search_radius=20`，大而非小**）⇒ **复现 ×4**、`rc=0 NORMAL / end_t=1 / negJac=0`、力平衡 ≈1e-7、`F_n ≈ 34.4 kN` 跨 3 网格漂移 **<0.6 %**。依据：`S5_contact_udg.md` §7、`S5_g7_verify_mesh.md` §3 |
> | G7 的 **σ/σ_law 收敛性** | ❌ **仍失败** —— 4 项独立 FAIL：自由 ξ **+30.62 %**、固定 ξ **+43.13 %**、小应变 −0.75 %…−5.47 %、跨配方 −5.85 %…−9.38 % |
>
> **比值发散的真凶也已改判**：不是平头边缘奇异，是 **`penalty=0.1` 软罚压力外溢**（压力溢出刚性足迹约
> 一个单元宽 ⇒ 穿透随网格变 ⇒「固定行程」与「固定压入量」无法同时收敛）。依据：`S5_g7_verify_metric.md` §5/§6。
>
> ⇒ **本文的 `15.14 cm / 2.70×` 禁用结论、ξ ≤ 0.25 可靠包线、OLD 配方记录，全部继续有效。**
> ⇒ 但**现行可引用数字一律以 `docs/S5_QUOTABLE.md` 为准**；
> 波3/波4 的启动口径见 **`docs/S5_contact_contract.md` §4 条 6–8**（已改为「**可启动但限定可信口径**」）。

## 1. 结论（TL;DR）

| 项目 | 结论 |
|---|---|
| **可靠区间** | 压头**下沉 ≤ 50 mm**（名义应变 **ξ ≤ 0.25**）：`rc=0`、**穿透 ≈ 0**（<0.25 mm）、FE/law 平滑 |
| **可信的平冲头增强** | **FE/law ≈ 1.56 → 1.64**（随压入单调缓增） |
| **G7 工况** | 2 m 坠落需 ξ≈0.85（压到 15% 厚度）→ **远在可靠区间之外**，tet 泡沫单元反演，**3 次相同配置给出 3 个不同结果**（全部失败） |
| **因此** | **15.14 cm / 2.70× 不可引用**（既非标定值、也非可复现值） |

## 2. 收敛配方（已固化进库）

FEBio 4.13 的 **AUGLAG 默认子参数**（`gaptol=0` 关闭、`minaug=0`、`maxaug=10`、`smooth_aug=0`）在本 Ogden 软垫 deck 上**发散**。收敛配方：

```python
build_pad_contact_feb(..., laugon="AUGLAG", penalty=1.0,
                      tolerance=0.005,
                      aug_controls=AUGLAG_CONVERGENT_CONTROLS)  # gaptol=0.001, minaug=3, maxaug=200, smooth_aug=1
```

- FEBio 4.13 **没有** `<augtol>` 标签 —— 增广容差就是 `<tolerance>` 本身。
- 库已暴露：常量 `pad_contact.AUGLAG_CONVERGENT_CONTROLS`，参数 `tolerance=` / `aug_controls=`（`build_pad_contact_feb` 与 `run_pad_contact` 都有）。**默认行为不变**（不传 = 逐字同旧输出）。
- 验证：`temp/opensim_fe/s5_t9/verify_lib.py`（纯库路径，无裸 XML 注入）。

## 3. 证据

### 3.1 可靠区间扫描（AUGLAG pen=1 + 配方，n=1200）

> **⚠️ 引用限定（2026-10-06）**：下表的 1.56 / 1.59 / 1.62 / 1.64 **来自 OLD 配方**
> (`laugon="AUGLAG"`, `penalty=1.0`, `tolerance=0.005`, `aug_controls=AUGLAG_CONVERGENT_CONTROLS`,
> `node_reloc=1`, `search_radius=20`, `n=1200`, default `<Control>`)；
> **(i) 配方**：OLD AUGLAG（不是 NEW PENALTY 0.1）；
> **(ii) 网格收敛性**：在 tet4 `(3,3,2)` 单网格上重现历史值至 ≤0.1%，但 NEW 配方同区间锚点
> 在 3 网格 (3,3,2)→(6,6,4)→(12,12,8) 上漂移 **−0.75 %…−5.47 %**，refinement step 非单调递增 ⇒
> **未证明 OLD 配方在该网格族上做了 mesh-convergence 扫描**，本节只是「同 OLD 配方下数值重现」，
> 不是「与同龄段该配方已网格收敛」；
> **(iii) 不可与 NEW 配方数字并列**：NEW vs OLD、其 σ/σ_law 在 anchor sinks 上单调偏移
> **−5.85 %…−9.38 %**（`S5_g7_verify_anchors.md` §4.1），约为 ≤0.1% 重现噪声的 60–100× ⇒ 两条曲线
> 不在同一根 y 轴上；OLD 配方装上 NEW stepper 会因 closure `negJac=55` 立即失败 ⇒ 两配方**不可互换**；
> **(iv) 权威裁决**：`docs\S5_QUOTABLE.md`。

| sink (mm) | ξ | F_n (N) | **FE/law** | 穿透 (mm) | 结果 |
|---|---|---|---|---|---|
| 20 | 0.098 | 1589 | 1.56 | −0.164 | ✅ |
| 30 | 0.149 | 2163 | 1.59 | −0.228 | ✅ |
| 40 | 0.199 | 2700 | 1.62 | −0.235 | ✅ |
| 50 | 0.248 | 3232 | 1.64 | −0.086 | ✅ |
| 60 | — | — | — | — | ❌ negJac=53 |
| 70 | — | — | — | — | ❌ negJac=20 |

穿透为负 = 接触面有**亚毫米微间隙**（健康）。FE/law 从 1.56 平滑升到 1.64 ⇒ 这是**真实的平冲头几何效应**，但**量级只有 ~1.6，不是 2.7**。

### 3.2 G7 工况不可复现（决定性）

同一配置（`(3,3,2)`、sink=170 mm、n=2400、AUGLAG pen=1 + 配方）**连跑 3 次**：

| 重复 | 结果 |
|---|---|
| repro_0 | rc≠0，**negJac=27** |
| repro_1 | rc≠0，**negJac=3** |
| repro_2 | rc≠0，**negJac=15** |

⇒ 不仅失败，**失败方式本身都是非确定的**。此前一次 rc=0 的记录是**抽中了好签**。原因：FEBio/PARDISO 多线程 + 自适应时步，在单元濒临反演时结果随机。

### 3.3 为何 G7 不可解（物理/离散原因）

- 2 m/77 kg 落 0.20 m 垫 ⇒ 需 ~580 kPa ⇒ 1D 律给出 ξ≈0.85（压掉 85% 厚度，接近触底）。
- 该应变下：① 泡沫 Ogden 极陡；② 垫用 **tet4**（`pad_mesh` 由 hex 裂解），大应变易体积锁定/畸变；③ `(3,3,2)`=**100 mm 单元**，冲头边缘是应力奇异。
- 细化到 (6,6,4)/(8,8,6)/(10,10,8) **仍未解决**（同样反演）⇒ 不是单纯粗网格，而是**该材料+该应变超出 tet 泡沫的稳健域**。
- 升 penalty（3/5/10）更早反演；载荷驱动（pressure 0.58 MPa）在第一步即反演。

## 4. 对 S5 / 论文的含义

1. **S5 的定量结论仍不可用**，但原因**已彻底查清**（超稳域 + tet 泡沫 + 极端压入），不再是"谜"。
2. **可写的正面结论**：在 ξ ≤ 0.25 的可靠区间，3D 接触 FE 相对 1D 单轴应力律的增强为 **~1.6**，且穿透≈0 —— 这是"平冲头 + 摩擦 + 周围材料约束"的真实几何效应。
3. **验收口径仍须"同边界 FE vs FE"**；**G7 的 15.14 cm / 2.70× 一律不得引用**。
4. 若要 S5 覆盖 G7 工况，需换数值策略：**六面体单元**（而非 hex→tet 裂解）、近不可压大应变泡沫专用公式、或**自适应网格/重划分** —— 属于后续工作。

## 5. 复现命令

```powershell
$env:PYTHONPATH="src"
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\verify_lib.py     # 配方生效 + 可靠区间
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\envelope.py      # 区间边界
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\repro.py         # G7 不可复现（3 次）
```

产物在 `temp/opensim_fe/s5_t9/`（`finish/`、`mesh_refine/`、`repro/`、`envelope/`、`verify_lib/`、`h2_auglag/`）。

## 6. 与早前结论的关系（避免误引）

> **⚠️ 引用限定（2026-10-06）**：本节谈到的所有**OLD 配方数字**（1.56–1.64）**只在该配方族内可引用**——
> **(i) 配方**：OLD AUGLAG pen=1 + AUGLAG_CONVERGENT_CONTROLS，n=1200，default `<Control>`；
> **(ii) 网格收敛性**：在 tet4 `(3,3,2)` 单网格上重现历史值 ≤0.1%，但未做 mesh-convergence 扫描
> （见上节 §3.1 引用限定）；
> **(iii) 不可与 NEW PENALTY 配方数字并列**：NEW vs OLD σ/σ_law 在 anchor sinks 单调偏移
> **−5.85 %…−9.38 %**（`S5_g7_verify_anchors.md` §4.1），两条曲线不可叠加；
> **(iv) 权威裁决**：`docs\S5_QUOTABLE.md`。

- `pad_contact.py` docstring 中 T8 的**材料验证**（均质自由侧面 FE ≡ 1D 律 ≤0.5%）**仍然成立**。
- 但其中**接触**部分在高 sink 下得到的 2.6–3.5× / 3–6.4× 上界，取自**收敛包线之外**的解，**不应再引用**；可信值以本文 §3.1 的 ~1.6 为准。

## 7. 六面体（hex8）路线：已试，**失败**（2026-10-06）

假设：G7 的失败是 `tet4` 在大应变下反演 → 换 `hex8` 应当更稳。

**实现**：`pad_mesh.pad_block_hex(nx,ny,nz,size_mm)`（新；与 `pad_block_tet` 同 schema，
`hexes` (M,8) + quad4 面；`pad_contact.py` 的发射层已泛化为同时支持 `tet4`/`hex8` 与
`tri3`/`quad4`）。默认路径**逐字不变**。

**结果：hex8 比 tet4 更差。**

| 测试 | tet4 | hex8 |
|---|---|---|
| 带接触 sink=50 mm | ✅ (3,3,2) | ❌ 6,6,4 在 ~15 mm 反演；3,3,2 在 **4.3 mm** 反演 |
| **纯单轴压缩（无接触、侧向自由）** | ✅ 撑到 **ξ≈0.80**（`pad_foam` 已验证 2.54%） | ❌ **ξ=0.50 即反演**（negJac 24），ξ=0.85（negJac 8），粗网格 3×3×3 同样 |

**判定**：单轴测试**不含接触**，hex8 仍在 ξ≈0.5 反演 ⇒ 限制来自 **hex8 单元公式本身**（全积分
三线性 hex8 对该软泡沫的**体积/剪切锁定**），**与接触算法无关**。换单元类型**不能**解决 G7。

**含义**：
- "G7 不可解"的根因更精确地表述为：**该 Ogden 泡沫在 ξ≳0.5 超出全积分 tet4/hex8 单元族的稳健域**。
- 覆盖 G7 需要**换公式**而非换网格：混合 u/p（近不可压专用）单元、专用泡沫大应变本构、
  或自适应重划分（remeshing）。
- `pad_block_hex` 保留为工具（有回归测试），但**不用于 S5 主线**。

**复现**：
```powershell
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\hex_probe.py     # 带接触：hex 早早反演
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\hex_diag.py      # 小压入量扫描 + tet 对照
.venv\Scripts\python.exe temp\opensim_fe\s5_t9\hex_uniaxial.py  # 无接触单轴：hex ξ=0.5 反演
```

> 踩坑记录：手写最小 FEBio deck 必须用
> `<?xml version="1.0" encoding="ISO-8859-1"?>` + `<febio_spec version="4.0">` + `<Module type="solid"/>`
> 作根；写了 `<febio>` 会得到 **"FATAL ERROR: Failed opening input file"**（不是力学报错，极易误判）。
