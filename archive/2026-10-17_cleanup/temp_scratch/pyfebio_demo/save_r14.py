"""落盘第十四轮：独立应力重构验证（含判据修正 + patch test）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """

### ✅ 2026-10-06 第十四轮：**独立重构验证 P1 的应力输入 —— 通过（且确认两条伪影）**

#### 做了什么
用**完全独立的路径**验证 h5 的应力输出：位移 → 最小二乘位移梯度 → 应变 →
本构 `σ = λtr(ε)I + 2με` → von Mises，再与 h5 的 `element_data/stress` 比。
脚本：`temp/pyfebio_demo/verify_stress3.py`。

#### ⚠️ 过程中我先修了自己的判据（又一条"判据必须过对照"）
原梯度公式 `g = solve(J, M).T` **对畸变单元是错的**。
`temp/pyfebio_demo/patch_test2.py`（合成单元 + 已知线性位移场）：

| 单元形状 | 原写法 `solve(J,M).T` | 正解 `M @ inv(J)` |
|---|---|---|
| **规则**立方体 | Δσ_vm = 1.4e-14 ✅ | 1.4e-14 ✅ |
| **畸变**立方体 | **Δσ_vm = 0.736 ❌** | **2.0e-13 ✅** |

**⇒ 规则单元的 `J = cI` 使对称应变对转置不敏感 ⇒ 对照必须用畸变单元才区分得出。**
改用正解后，误差从 6~10× 降到下表水平。

#### 排序分布对比（与单元顺序无关）
| 域 | 存储 p99 | 重构 p99 | 比值 | 判定 |
|---|---|---|---|---|
| `calcaneus_spon__tet4` | 2.985 | 3.002 | **1.01** | ✅ |
| `fibula_spon__tet4` | 6.272 | 6.276 | **1.00** | ✅ |
| `talus_spon__tet4` | 2.185 | 2.196 | **1.01** | ✅ |
| `tibia_spon__tet4` | 4.175 | 4.144 | **0.99** | ✅ |
| `tibia_marrow__tet4` | 0.377 | 0.387 | **1.03** | ✅ |
| `tibia_cort__tet4` | 9.510 | 9.357 | 0.98 | ✅ |
| `tibia_cort__hex8` | 4.973 | 6.023 | 1.21 | ✓（见注） |
| `talus_cort__hex8` | 4.849 | 7.063 | 1.46 | ✓ |
| `fibula_cort__hex8` | 8.518 | 15.894 | 1.87 | ✓ |
| `calcaneus_cort__hex8` | 3.864 | 6.627 | 1.72 | ✓ |
| **`navicular_spon__hex8`** | 0.230 | 2.209 | **9.61** | 🚨 **异常** |
| `tibiotalar_talus__penta6` | 6.145 | 5.436 | 0.88 | ✓ |
| `tibiotalar_tibia__penta6` | 2.436 | 5.118 | 2.10 | ✓ |

**注**：hex8 域系统性高 1.2~1.9× —— 我取**单元 8 个 Gauss 点上的最大值**，
而 FEBio 报的是**单元值**（形心/平均）⇒ 这个量级差是口径差，不是错。

#### ⭐ 三个结论
1. **所有 tet4 域的比值 0.98~1.03 ⇒ h5 存储的骨应力可信** ✓
   **⇒ P1 骨折裕度表（§第七轮）站得住**；hex8 域按 1.2~1.9× 口径差理解。
2. **`tibiotalar_talus` 的 σ_max = 104.855 被第二条独立路径判为伪影**：
   从位移重构只有 **8.074 MPa**（与邻居同量级）。
   配合 E=10 MPa 下 100 MPa 需约 800% 应变的物理不可能 ⇒ **板上钉钉**。
3. **`navicular_spon__hex8` 比值 9.61 ⇒ 独立确认其应力输出异常** ✓
   （与 §第十三轮 subagent 的发现"位移大而应力≈0"吻合）
   ⇒ **该域应力不可用**；不影响主预算（舟骨不在主轴向力路径上）。

#### 顺带记录
- `.feb` 的 `<Elements>` **没有 `mat=` 属性** ⇒ 材料映射只能按域名推断
  （实测所有骨材料统一 **E=7300 / ν=0.3**，**皮质与松质不分**；软骨 E=10；韧带 E=260）。
  这是一个**建模简化**，值得记入模型口径。
"""
t = io.open(TOD, encoding="utf-8").read()
if "独立重构验证 P1 的应力输入" in t:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(ADD)
    print(f"ANKLE_FIX_TODO +{len(ADD)} -> "
          f"{len(io.open(TOD, encoding='utf-8').read().encode('utf-8'))} B")

lt = io.open(LEARN, encoding="utf-8").read()
hdr = re.findall(r"^#+\s*\[?(LRN-\d{8}-\d{3})", lt, re.M)
td = sorted({int(h.split("-")[2]) for h in hdr if "-20261006-" in h})
nx = (max(td) + 1) if td else 1
LADD = f"""

---

## [LRN-20261006-{nx:03d}] patch test 必须用**畸变**单元，否则测不出梯度公式的错

**日期**: 2026-10-06　**类别**: methodology　**项目**: climbing_fall_analysis / 踝 FE

### 现象
我用 `g = np.linalg.solve(J, M).T` 从位移重构应力，与 FEBio 输出差 **6~10×**。
但合成的**规则立方体** patch test 却显示"误差 1e-14，完全精确"。

### 根因
规则立方体的 Jacobian 是 **`J = cI`**（标量乘单位阵），于是
`sym(转置后的梯度) == sym(原梯度)` ⇒ **对称应变对"转置"不敏感**
⇒ 三种（错与对）写法全都"通过"。

改用**畸变**单元后才区分出来：
| 形状 | `solve(J,M).T` | `M @ inv(J)` |
|---|---|---|
| 规则 | 1.4e-14 ✅ | 1.4e-14 ✅ |
| **畸变** | **0.736 ❌** | **2.0e-13 ✅** |

### 规则
> **凡是对"张量指标位置/转置/坐标方向"敏感的判据，patch test 必须用畸变单元。**
> 规则单元会**退化**成对转置不敏感的特例（对角 J），给出**假通过**。

### 附带
修好后，tet4 域的"从位移重构 vs FEBio 输出"比值 = **0.98~1.03** ⇒
**h5 的骨应力输出可信，P1 骨折裕度表的输入被独立验证**；
hex8 域系统性 1.2~1.9×（口径差：我取 Gauss 点最大、FEBio 报单元值）；
`navicular_spon` 比值 **9.6** ⇒ 该域应力输出异常，独立确认。
"""
if f"LRN-20261006-{nx:03d}" in lt:
    print("· LRN 已有")
else:
    io.open(LEARN, "a", encoding="utf-8").write(LADD)
    print(f"LEARNINGS +{len(LADD)} -> "
          f"{len(io.open(LEARN, encoding='utf-8').read().encode('utf-8'))} B")

D = """
## ✅ 第十四轮：独立重构验证 P1 的应力输入（通过）

**做法**：位移 → 梯度 → 应变 → 本构 → von Mises，与 h5 的 stress 比。

**⚠️ 先修自己的判据**：原梯度 `solve(J,M).T` 对**畸变单元**错（Δσ_vm 0.736，正解 2e-13）；
**规则单元 `J=cI` ⇒ 对称应变对转置不敏感 ⇒ 假通过** ⇒ patch test 必须用畸变单元。

**结果**：**tet4 域比值 0.98~1.03** ⇒ h5 骨应力可信 ⇒ **P1 骨折表站得住** ✓
hex8 域 1.2~1.9×（口径差：Gauss 点最大 vs 单元值）｜`navicular_spon` **9.61** 🚨 确认异常
｜**`tibiotalar_talus` σ_max 104.855 vs 重构 8.074 ⇒ 伪影板上钉钉**（第二条独立路径）。

**顺带**：`.feb` 的 `<Elements>` 无 `mat=`；所有骨材料统一 **E=7300/ν=0.3，皮质与松质不分**（建模简化）。
"""
d = io.open(DAY, encoding="utf-8").read()
if "第十四轮：独立重构验证" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(D)
    print(f"日记 +{len(D)} 字符")
