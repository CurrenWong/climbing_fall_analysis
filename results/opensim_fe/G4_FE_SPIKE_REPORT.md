# G4 ② FE 深水区 — P0 可行性 spike：FEBio 原生损伤材料（fibula_r）

> **筛选级 spike，不是标定模型。** 目的只有一个：确认 FEBio 能否解析/求解`elastic damage` + `DC Drucker shear stress` + `CDF Simo` 物料，以及它对一根代表骨（`fibula_r`）产生什么损伤场。所有数字本机实测；未对任何参数调参。

_生成时间：2026-10-11 21:27:48；FEBio `D:\Program\FEBioStudio\bin\febio4.exe` (probe=4.13.0)。_

## 0. 一句话结论

- **解析测试：ACCEPT**（FEBio 4.13 接受该物料类型；已收敛并产生损伤场）
- 受损模型在 **load = 125 N** 收敛；单次求解 wall-clock **1.0 s**；该载荷重复求解 **4/5 收敛**（非确定性，见 §4.3）。
- 产生损伤场：`Dmax` = 0.9807，D>0 单元占比 86.3%。

- 主承载域 `fib_cort`（PID 81000900）：E=18500.0 MPa，ν=0.3，σ_y=37.0 MPa，E_t=16700.0 MPa → Simo `a` = 21.362（c=0.0）。

## 1. 最小解析测试（先跑，决定性）

- 文件：`D:\Project\climbing_fall_analysis\temp\opensim_fe\g4_fe_spike\_parse_test.feb`（1 个 hex8，8 节点，单物料 `elastic damage`，底固定 + 顶节点力 -9.5 N/节点，10 载荷步）。
- 结论：**ACCEPT**。物料摘要 出现于 FEBio 日志（解析器接受 neo-Hookean + CDF Simo + DC Drucker shear stress）

FEBio 日志中的物料摘要（节选，证明解析成功）：

```
  1 - cube (type: elastic damage)
	density ........................................... : 1
	  density ......................................... : 1e-06
	  E ............................................... : 18500
	  v ............................................... : 0.3
	damage: (type: CDF Simo)
	  Dmax ............................................ : 1
	  a ............................................... : 21.362
	  b ............................................... : 0
	criterion: (type: DC Drucker shear stress)
	  c ............................................... : 0
```

完整日志：`D:\Project\climbing_fall_analysis\temp\opensim_fe\g4_fe_spike\_parse_test.log`。

## 2. 注入的 `<material>` XML（原文）

`pyfebio` 没有 damage/reactive 物料类，故用**文本替换**注入（`<material name="fib_cort" ...>…</material>` → 下块；`name=`/`id=` 保持不变，`<MeshDomains>` 仍解析到同一物料）。同时在 `<plotfile>` 加 `<var type="damage"/>`。

```xml
    <material name="fib_cort" id="3" type="elastic damage">
      <elastic type="neo-Hookean">
        <density>1e-06</density>
        <E>18500.0</E>
        <v>0.3</v>
      </elastic>
      <damage type="CDF Simo">
        <a>21.361959960016154</a>
        <b>0.0</b>
        <Dmax>1.0</Dmax>
      </damage>
      <criterion type="DC Drucker shear stress">
        <c>0.0</c>
      </criterion>
    </material>
```

**映射依据**（FEBio User Manual 4.13）：§4.6.1 `elastic damage`（未约束弹性，对 ν=0.3 用 plain 而非 `uncoupled`）；§4.6.2.1 `CDF Simo`（`a`=α，与 Ξ 同量纲；`b`=β，最大损伤 1−β；`Dmax`≤1）；§4.6.3.7 `DC Drucker shear stress` （`c` ∈ [−27/8, 9/4]，`c=0` ⇒ von Mises，k=σ_y/√3）。

## 3. THUMS 参数（逐 PID，实测卡）

| domain | PID | src_keyword | E (MPa) | ν | σ_y (MPa) | E_t (MPa) | E_t/E |
|---|---|---|---|---|---|---|---|
| `fib_spon_end` | 81000800 | *MAT_PIECEWISE_LINEAR_PLASTICITY | 145.6 | 0.45 | 18.6 | 29.0 | 0.1992 |
| `fib_spon_center` | 81000801 | *MAT_VISCOELASTIC | 11.992 | 0.499001 | — | — | — |
| `fib_cort` | 81000900 | *MAT_PIECEWISE_LINEAR_PLASTICITY | 18500 | 0.3 | 37.0 | 16700.0 | 0.9027 |

主承载域 = **`fib_cort`**。Simo α = σ_y/√3·(1−(4/27)c)^(1/6) = **21.362** MPa（c=0.0）。E_t 仅记录，不用于本次 spike（损伤演化只由 `a` 控制速率；`b=0` ⇒ Dmax≤1）。

## 4. 实测数字（rc / wall-clock / σ_vm / 损伤场）

| 运行 | load (N) | rc | wall-clock (s) | 收敛步 | 状态 |
|---|---|---|---|---|---|
| elastic twin | 125 | 0 | 0.5 | 1 | OK |
| damaged | 125 | 0 | 1.0 | 1 | OK |

受损模型降载搜索（同 load 下 elastic twin 亦按需重跑）：

| load (N) | config | ok | rc | wall (s) | FEBio error (verbatim) |
|---|---|---|---|---|---|
| 1000 | default | no | 1 | 0.9 | 98 negative jacobians detected. | failed to converge at time : 1 |
| 1000 | time_stepper+10steps | no | 1 | 0.0 | 98 negative jacobians detected. | failed to converge at time : 1 |
| 500 | default | no | 1 | 1.4 | 324 negative jacobians detected. | failed to converge at time : 1 |
| 500 | time_stepper+10steps | no | 1 | 0.1 | 324 negative jacobians detected. | failed to converge at time : 1 |
| 250 | default | no | 1 | 2.7 | 265 negative jacobians detected. | failed to converge at time : 1 |
| 250 | time_stepper+10steps | no | 1 | 0.0 | 265 negative jacobians detected. | failed to converge at time : 1 |
| 125 | default | yes | 0 | 1.0 | (no ERROR block) |

### 4.1 σ_vm（MPa）逐域：elastic vs damaged

| domain | elastic max | elastic p95 | elastic mean | damaged max | damaged p95 | damaged mean |
|---|---|---|---|---|---|---|
| `fib_spon_end` | 0.0310 | 0.0176 | 0.0055 | 0.8766 | 0.5268 | 0.1883 |
| `fib_spon_center` | 0.1556 | 0.0442 | 0.0149 | 17.9304 | 0.2919 | 0.1799 |
| `fib_cort` | 111.6490 | 65.9914 | 21.2641 | 36.7131 | 34.0448 | 16.7192 |

主承载域 σ_vm 偏差（(damaged−elastic)/elastic）：

- max:  -67.1%（111.649 → 36.713 MPa）
- p95:  -48.4%（65.991 → 34.045 MPa）
- mean: -21.4%（21.264 → 16.719 MPa）

### 4.2 损伤场 D（末状态，主承载域）

- 单元数：842
- `Dmax` reached：**0.9807**
- mean D：0.3520；median D：0.3594
- D>0 单元占比：**86.3%**；D>0.5 占比：30.3%

### 4.3 收敛稳健性（对**同一** .feb 重复求解）

在 load = 125 N 对同一 `fibula_r_damage_c0.feb` 重复求解 **5** 次：**4/5 收敛**。

| # | ok | rc | wall (s) | d_max | d_mean | D>0 | FEBio error |
|---|---|---|---|---|---|---|---|
| 0 | yes | 0 | 0.8 | 0.9865 | 0.2245 | 86.3% |  |
| 1 | yes | 0 | 0.8 | 0.9865 | 0.2252 | 86.3% |  |
| 2 | yes | 0 | 0.8 | 0.9865 | 0.2261 | 86.3% |  |
| 3 | yes | 0 | 3.2 | 0.9954 | 0.3582 | 86.3% |  |
| 4 | no | 1 | 4.3 | — | — | — | 260 negative jacobians detected. | 78 negative jacobians detected. | 3 negative jacobians detected. |

**关键发现（诚实边界）**：软化型损伤在**力控制**下处于稳定性边界，FEBio 的非线性迭代路径对（并行）线性求解器的浮点舍入敏感 —— **同一输入文件在不同次运行会收敛/发散不同**，收敛状态下的损伤场也有散布。故下方「单次求解 ≤ 5 min」的判定须理解为「**在收敛的那些运行里**」的耗时，而非「每次都收敛」；这是本 P0 最重要的可行性风险，后续若要在力控制下稳定求解，需考虑位移控制 / 粘性正则 / 弧长法（P1）。

## 5. P0 gate 判定

| gate | 阈值 | 实测 | 判定 |
|---|---|---|---|
| 单次求解 ≤ 5 min | ≤ 300 s | 1.0 s | **PASS** |
| 解析/收敛（能跑） | rc=0 且收敛 | 是 | **PASS** |
| 收敛稳健性（新增观察） | ≥ 1 次收敛（原计划未设） | 4/5 | **PASS（弱）** |
| 网格客观性 | 粗网格重跑 σ/首发点偏差 | **未测** | **未能廉价执行** |

**网格客观性说明（诚实边界）**：本 spike 直接复用 THUMS AM50 固定网格，仓库内没有为「纯实体 fibula」接线的重网格器；在 P0 时间内无法廉价生成一套更粗的同胚网格（不同单元阶/拓扑会改变约束与点载荷带），故此项 **未测**，留待后续（P1）。不作任何「已满足网格客观性」的声称。

## 6. 阻塞与复现

受损模型在记录载荷下收敛（rc=0）；**唯一阻塞/风险** = §4.3 的收敛非确定性（软化本构在力控制下处于稳定性边界）。本载荷 4/5 次收敛。

```powershell
$env:PYTHONPATH="src"
& .venv\Scripts\python.exe scripts\opensim_fe\g4_fe_spike.py --bone fibula_r --c 0
```

产物：

- `results/opensim_fe/G4_FE_SPIKE_REPORT.md`（本文件）
- `temp/opensim_fe/g4_fe_spike/`（`_parse_test.feb/.log`、`*_elastic.feb`、`*_damage_*.feb`、`injected_material.xml`、`.xplt/.hdf5/.log`、`g4_fe_spike.json`）

**未修改任何既有模块**：本 spike 只新增 `scripts/opensim_fe/g4_fe_spike.py`，其余均以 import + 调用复用。

