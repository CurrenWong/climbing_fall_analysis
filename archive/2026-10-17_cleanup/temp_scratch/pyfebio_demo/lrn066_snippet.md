
---

## [LRN-20261005-066] 弱弹簧**已实现并验证**：FEBio `<Discrete>` 语法 + 3 个实测坑；**但力控制未被弹簧解决**（真因=关节初始张开）

- **承接**：`[LRN-20261005-065]`（弱弹簧 vs 位移控制的本质差别）、`[LRN-20261005-064]`（TIE 移植成功、600 N 达标）

### ✅ FEBio 弹簧语法（**实测权威，别再猜**）

`FEBio_User_Manual.pdf` §3.6.7 + §3.16：

```xml
<!-- Mesh 段内 -->
<DiscreteSet name="tal_spring_x">
  <delem>3186,4030</delem>          <!-- 2 节点弹簧：i,g（g = 接地哑节点） -->
</DiscreteSet>
<!-- ★ 顶层，且必须在 <MeshDomains> 之后 -->
<Discrete>
  <discrete_material id="1" type="linear spring">
    <E>1.548</E>                    <!-- 弹簧刚度 -->
  </discrete_material>
  <discrete dmat="1" discrete_set="tal_spring_x"/>
</Discrete>
```

**最小算例校验**（1 个 hex8 立方体 + 1 根接地弹簧）：
`Δk=99、δ=0.1` ⇒ ΔF 实测 **104.70 N** = 理论 `Σ_t k·δ(t) = 10×10.5 = 105 N` ✓

### 🐛 三个实测语法坑（都不是靠读手册能避免的）

1. **`<NodeSet>` 的值是「逗号列表」，不是 `<node id>` 子元素**
   `<NodeSet name="bottom">1,2,3,4</NodeSet>` ✓ ；写成多个 `<node id="1"/>` ⇒
   `tag "NodeSet" : invalid value`
2. **两种 BC 的 dof 写法不同**（手册 §6.2 示例）
   - `zero displacement` ⇒ `<x_dof>1</x_dof>`（**每分量一个标签**）
   - `prescribed displacement` ⇒ `<dof>x</dof><value lc="1">-0.1</value>`
     （**必须给 `lc` 指向 load curve ⇒ 需要 `<LoadData>` 段**）
   - 混用 ⇒ `tag "dofs"/"z_dof" : unrecognized tag`
3. **`<Discrete>` 必须写在 `<MeshDomains>` 之后**
   - 放在 `</Mesh>` 与 `<MeshDomains>` 之间 ⇒ `tag "discrete" : invalid value for attribute`
   - （我按"逻辑顺序"放，踩了这个；最小算例里恰好放在后面才没暴露）

**⇒ 通用教训**：小算例验证要**照抄结构顺序**，只能证明"该顺序可行"，**不能证明顺序无关**。

### ✅ 实现（`ankle_lig_feb.py`，env 开关，默认关 ⇒ 零回归）

- 开关：`TALUS_SUPPORT=fix`（默认）| `spring`；`TALUS_K`（默认 500 N/mm/方向）；
  `TALUS_SPRING_L`（哑节点臂长，定向用，默认 5 mm）
- `spring` 时：对 `talus_low_ns` 的每个节点、每个横向方向建 1 根接地弹簧
  （哑节点 = 原节点 + L·e_d，固定）；单根 `k = TALUS_K / n_nodes`
- 距骨横向硬约束**跳过**，改固定哑节点
- 实测：**323 节点 × 2 方向 = 646 根弹簧，单根 k=1.548 N/mm，收敛 `converged at time : 1`**
- **零回归**：旧格式 `209.5589 N`、新格式+TIE `592.6068 N` 逐字一致

### 📊 硬结果：弹簧 vs 硬约束（轴向，`DRIVER=dispx`）

| DISP (mm) | 硬约束 Fz | 弱弹簧 Fz | 差 |
|---|---|---|---|
| 0.25 | 124.04 | 111.26 | **−10.3%** |
| 0.50 | 230.90 | 203.78 | −11.7% |
| 1.00 | 418.10 | 369.26 | −11.7% |
| 1.60 | 592.61 | 505.27 | −14.7%（且该点不稳定） |

**⇒ 两者差了 10~12%，不是 bug、不是误差 —— 是 `[LRN-20261005-065]` 说过的：**
**边界条件真的不同**（横向释放后关节在轴向受压时能横向让开 ⇒ 同位移下力更小）。
**⇒ 任何"弱弹簧 vs 硬约束"的数字对比，必须先声明用哪种边界条件。**

### ⛔ 关键负面结果：**弹簧并没有解决力控制**

`DRIVER=force FORCE_N=600 TALUS_SUPPORT=spring`，扫 `TALUS_K = 50 / 500 / 5000`：
**全部 `ERROR / 11 negative jacobians`**（与弹簧刚度无关）。

**⇒ 弹簧确实解掉了「刚体模态奇异」（`[LRN-20261005-057]` 的判断），但力控制仍失败。**
**真因是 `[LRN-20261005-058]` 早已定位的那条**：**关节初始张开**
（t=0 中位间隙 ~1 mm、接触面积仅 ~31 mm²）⇒ 力加载时 Newton 过冲 ⇒ 单元翻转。

**⇒ 正解 = Anderson 的「先落位、后加载」**，而手册 §6.1.2 明确支持：
> "This makes it possible, for example, to **switch between load and displacement
> controlled boundary conditions in multi-step analyses**."

⇒ **多步：(step 1) 位移控制把关节压合到位 → (step 2) 换成力控制加载到 600 N。**
（`<Step><step id="1">…</step><step id="2">…</step></Step>`，
空 `<Control/>` 表示沿用；`<Boundary>`/`<Loads>` 可放在 step 内只对该步生效。）

### ⏭️ 下一步（明确）

1. **多步落位**：step1 位移落位（消除 ~1 mm 初始间隙）→ step2 力控制到 600 N
2. 之后才是「用 OpenSim 的力/力矩驱动 ⇒ 姿态作为输出」的姿势判据
3. 仍待核实：**Anderson §2.2 的弹簧布置口径**（挂点/总刚度 vs 单点/释放哪些方向）
   —— 本实现是「每节点、每横向方向各一根、方向合力 = TALUS_K」的合理猜测，**未对标原文**
