
---

## [LRN-20261005-068] FEBio `<Rigid>` 机制**语法已实测验证**（路线 A 可行）—— 刚体 + 参考点 + 接地弹簧 + 转轴解除

- **承接**：`[LRN-20261005-067]`（Anderson §2.2 口径 = 刚性骨 + reference node；可变形模型做不到转轴解除）
- **来源**：`FEBio_User_Manual.pdf` §3.11.16 / §3.12 / §4.18 / §3.7.1
- **验证方式**：最小算例（1 个 rigid hex8 + 1 个单节点"地面"刚体 + 1 根 rigid spring）

### ✅ 实测通过的完整语法

```xml
<!-- 1) 刚体材料（§4.18）—— E/ν 对结果无影响，但 contact 的 auto-penalty 会用到 ⇒ 给合理值 -->
<material id="1" type="rigid body" name="bodyA">
  <density>1e-9</density>
  <E>12000</E><v>0.3</v>
  <!-- center_of_mass 可省（省则自动算，但必须给 density；两者都省 = 报错） -->
</material>

<!-- 2) 域声明为刚体（§3.7.1：type="rigid-solid" 是 rigid 材料的默认） -->
<SolidDomain name="partA" type="rigid-solid" mat="bodyA"/>

<!-- 3) 把节点集挂到刚体（§3.11.16）—— 这些节点**不再有自由度** -->
<bc type="rigid" node_set="nsA"><rb>bodyA</rb></bc>     <!-- rb = 材料 name（也可用 id） -->

<!-- 4) 刚体段的约束 / 载荷 / 连接（§3.12） -->
<Rigid>
  <rigid_bc type="rigid_fixed">           <!-- 锁平动+转动 DOF：Rx,Ry,Rz,Ru,Rv,Rw -->
    <rb>bodyG</rb>
    <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>1</Rz_dof>
    <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>
  </rigid_bc>
  <!-- ★★ 转轴解除/固定 = 把对应 R*_dof 写成 0 —— 这就是 Anderson 的
       "fixed against rotation about the inversion-eversion / internal-external axes,
        these restraints were removed" -->
  <rigid_connector type="rigid spring">   <!-- §3.12.4.7：连**两个刚体** -->
    <body_a>1</body_a><body_b>2</body_b>
    <insertion_a>0.5,0.5,0.5</insertion_a>   <!-- 任意点，不必是节点 -->
    <insertion_b>0,0,-1</insertion_b>
    <k>500</k>                                <!-- free_length 省略 ⇒ 取初始长度 -->
  </rigid_connector>
  <rigid_load type="rigid_force">         <!-- ★ 用 <dof>Rz</dof>+<value>，**不是** <force> -->
    <rb>bodyA</rb><dof>Rz</dof><value lc="1">100</value>
  </rigid_load>
</Rigid>
```

### 🔑 实测证据（数值验到底）

| 检验 | 结果 |
|---|---|
| 解析 + 初始化 | ✅ 全段被接受 |
| 自由度计数 | `Nr of equations = 1`（只有 Rz 自由）✅ 约束精确生效 |
| **弹簧刚度** | `Δz = 0.235768` 精确满足非线性平衡 `k(L−L₀)(1.5+Δz)/L = 100`（积=**0.2000**）✅ |
| **转轴解除** | `Ru_dof=0, Rw_dof=0` ⇒ 两个角点位移完全不同（node1 `z=1.602` vs node5 `z=0.024`）⇒ **刚体真的转动** ✅ |

### 🐛 三个实测坑

1. **`rigid_force` 的参数是 `<dof>Rz</dof>` + `<value>`，不是 `<force>`** ⇒ 写 `<force>` 报
   `tag "force" : unrecognized tag`（手册 §3.12.3.1 的示例里 `</dofs>` 是**别字**，应为 `</dof>`）
2. **刚体的"接地"必须造一个哑刚体**：`rigid spring` 连的是**两个刚体**；
   "connected to the ground" ⇒ 用一个（可以是单节点的）**全约束哑刚体**当 ground
3. **`<plotfile>` 里没有 `rigid displacement`/`rigid rotation` 这两个输出变量** ⇒ 写了会
   `FATAL ERROR: Output variable ... is not defined`。刚体位移要从**节点位移**（`<node_data>`）读。
   另：1 个自由度的 rigid 问题 `pyfebio.xplt.to_hdf5` 读不了（`Invalid FEBio file`）

### ⚠️ FEBio 4 的关键行为

**手册 §4.18 明写：rigid body 的自由度「initially unconstrained」（可自由平动 + 绕任意轴转动），
必须自己加约束**（FEBio 3 时代是初始全约束 ⇒ 这是一个版本行为差异，照着旧经验写会错）。

### 🧠 方法论

- **"我以为的期望值"要按精确物理算一遍**：我先把期望写成 `F/k = 0.2`，实测 0.2358，
  差点误判成"没生效"。**代入精确非线性平衡后 4 位吻合** ⇒ 弹簧是对的，是我的期望是错的。
  **⇒ 判据必须先自证，再用来判别人。**
- 手册示例可能有别字（`</dofs>`），**以 FEBio 实际报错为准**。
