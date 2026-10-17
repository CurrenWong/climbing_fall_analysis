
---

## [LRN-20261005-069] 路线 A **已落地**：距骨刚体 + 接地弹簧跑通（位移控制）；力控制仍卡在**关节初始张开 ⇒ 必须落位**

- **承接**：`[LRN-20261005-067]`（Anderson 口径）、`[LRN-20261005-068]`（`<Rigid>` 语法实测）

### ✅ 已实现（`ankle_lig_feb.py`，env 开关，默认关 ⇒ 零回归）

`TALUS_RIGID=1`：
- 距骨 `talus_cort`/`talus_spon`/`talus_marrow` 域 → `<SolidDomain type="rigid-solid" mat="RIGIDtalus">`
- 916 个距骨节点经 `<bc type="rigid" node_set="talus_rigid_ns"><rb>RIGIDtalus</rb></bc>` 挂到刚体
- 地面 = 2 个哑节点（`ground_rigid_ns`）挂 `RIGIDground`，全 6 DOF 固定
- `<rigid_connector type="rigid spring">` k=`TALUS_SPRING_K`（默认 500 N/mm）
- 转轴锁定开关 `TAL_ROTLOCK`（1=落位期固定 / 0=解除）
- `<Rigid>` 段**必须写在 `</Boundary>` 之后**（`_RIGID_LINES` 缓冲后 flush）

### 📊 位移控制下跑通

| 配置 | 状态 | Fz |
|---|---|---|
| `TAL_ROTLOCK=1`（转轴锁） | **NORMAL** | **202.29 N** |
| `TAL_ROTLOCK=0`（解除） | **NORMAL** | **187.90 N** |

⇒ **转轴解除使结果变化 −7.1%** ⇒ **`rigid_rotation` 的"固定/解除"机制真的在起作用。**

### 🐛 两个真 bug（各花了很久）

1. **`<rigid_connector>` 的 `body_a`/`body_b` 是「材料 ID」，不是名字**
   我按最小算例写死 `1/2`，但真模型里 1 号材料是 `mat_tibia_cort`（可变形）⇒ 报
   **`Nonlinear constraint 1 (<unnamed>) failed to initialize`**。
   正解：刚体材料追加在末尾 ⇒ `body_a = len(MATERIALS)+1`。
   **⇒ 最小算例"恰好成立"的常量是陷阱：写死常量前先问"它在新上下文里还是这个值吗"。**
2. **`<Rigid>` 段的位置**：放在 `</MeshDomains>` 之后（手册 §3.1 推荐序）会初始化失败；
   放在 **`</Boundary>` 之后**（与跑通的最小算例同序）才行 ⇒ 现用 `_RIGID_LINES` 缓冲后 flush。

### ⛔ 力控制仍未打通（本轮到此）

`DRIVER=force`：负 Jacobian 已消失（刚体化有效），但**变成"不收敛"**：
| 配置 | 结果 |
|---|---|
| `TAL_ROTLOCK=0`, k=500 | t=**0.227** 失败（229.6 s） |
| `TAL_ROTLOCK=1`, k=500 | t=**0.0** 失败 |
| `TAL_ROTLOCK=0`, k=5000 | t=0.235 失败 |
| `TAL_ROTLOCK=0`, k=50000 | t=0.172 失败 |

**症状**：施加 **5 N** 却在胫骨顶读出 **177.5 N**、距骨约束反力 **5004 N** ⇒ 远未平衡。

**⇒ 根因（与 `[LRN-20261005-058]` 一致）：关节初始张开**（t=0 中位间隙 ~1 mm、接触仅 ~31 mm²）。
加载时胫骨与距骨仅靠 **500 N/mm 弹簧 + 韧带**相连 ⇒ 近似机构 ⇒ 发散。

**⇒ 正解 = Anderson 的落位协议（`[LRN-20261005-067]` 已核实原文）**：
> step1 把胫骨带到与弱约束距骨的**最小接触** → step2 **施加 200 N 预载**（解除转轴约束之后）
> ⇒ 关节自行沉降到解剖对位 → 然后才是 duty cycle

**⇒ 下一步 = 多步分析**（`<Step>`：`<step id="1">` 落位 + `<step id="2">` 力加载到 600 N）。
手册 §6.1.2 明确支持跨步在位移控制 ↔ 载荷控制间切换。

### 🧠 方法论

- **"最小算例通过"只能证明语法，不能证明把它搬到真模型里也成立**：本轮两个 bug
  （材料 ID、段序）都是"最小算例恰好没暴露"的类型。**⇒ 搬到大模型后必须重跑并逐项核对。**
- **失败态会给出"看起来很合理"的数值**（如恒定 239.19 N、距骨反力 5004 N）。
  **⇒ 读数前先确认 `t=1.0` 且 `NORMAL`；否则数值无意义。**
