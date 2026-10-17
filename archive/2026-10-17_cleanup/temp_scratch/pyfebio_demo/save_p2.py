"""落盘 §第八轮（P2 收口：三假设推翻 + 柔度在韧带）。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
OVR = os.path.join(R, "docs", "总体进度.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

ADD = """

### ⚠️ 2026-10-06 第八轮：**P2 收口 —— 三个假设全被推翻，靶子找错了**

#### 结论：**当前软骨网格是健康的，"修退化软骨"这个任务本身是错的**

| 假设 | 实测（脚本 `p2_detj.py` / `p2_thick.py` / `p2_end2.py`） | 判定 |
|---|---|---|
| 软骨网格退化 | `tibiotalar_talus__penta6` det(J) **min 1.2245**、`tibiotalar_tibia__penta6` **0.7746**；全模型 det(J)≤0 的单元 **0 个**、<1e-2 的 **0 个** | ❌ **推翻** |
| 超薄单元 | 软骨厚度 min **0.6828 mm**、中位 1.298、max 1.491；<0.1 mm 的 **0 个** | ❌ **推翻** |
| 软骨被压穿 | +1.6mm 末态厚度只减 **0.0142 mm**（中位），det(J) 末态 min 0.6116/1.2039 | ❌ **推翻** |
| σ_max=104.85 是网格奇点 | **E=12 MPa 下 100 MPa 需要约 800% 应变 ⇒ 物理不可能** | ✅ **改判 = 罚因子接触尖峰** |

**⇒ 之前携带的「软骨近退化单元 det(J)=2.39e-3」是旧网格的数字**，对当前
`cart_patch m=1` 网格不成立 ⇒ **状态表问题 #2/#3 的归因作废。**

#### ⭐ 本轮真正的发现：**柔度全在韧带，不在软骨**

逐域 dz 跨度（`p2_by_domain2.py`；弯曲度量，dz = 轴向）：
| 类别 | 域 | dz 跨度 | z 跨度 |
|---|---|---|---|
| 韧带 | **lig_MD2_tibiocalcaneal** | **1.6480 mm** | 36.5 mm |
| 韧带 | **lig_CFL_calcaneofibular** | **1.5290 mm** | 28.7 mm |
| 韧带 | **lig_SPRING_calcaneonavicular** | **1.2661 mm** | 5.0 mm |
| 韧带 | MD1 / PTCL / MD4 / ATFL / PTFL | 0.85~1.05 mm | 17~28 mm |
| 骨 | fibula_cort | 0.0984 mm | 70.7 mm |
| 骨 | calcaneus_cort | 0.0388 mm | 67.1 mm |
| 骨 | tibia_cort | 0.0206 mm | 51.4 mm |
| 骨 | talus_cort | 0.0100 mm | 43.2 mm |

⇒ 骨只变形 **0.01~0.10 mm**，关节只闭 **0.014 mm**，**韧带吸收 0.85~1.65 mm**
⇒ **300 倍差距** ⇒ 与 `[LRN-080]`「韧带承担 77% 载荷」**独立吻合** ✓

**⇒ 状态表 #4「接触刚度 ~150 N/mm」重新解释**：那是**韧带主导的系统刚度**，
不是接触刚度；**关节在本模型里近似刚性**。

**⇒ 状态表 #3「压缩上限 +1.6~1.7 mm」也重新解释**：那是**韧带**给出的极限
（不是软骨被压穿）。

#### 对计划的影响（重要）
- **P2「修退化软骨」删除** —— 网格健康，没有可修的。改记为：
  「**软骨应力口径用 p99，不用 max**（max 是罚因子接触尖峰）」。
- **新增 P2'：刚度/压缩上限是韧带问题** —— 任何刚度标定都要针对韧带
  （`LIG_MODE=shell` 用固定 `E=260`、不用 `a_cross`；这直接决定系统刚度）。
  这也是「tet4 vs penta6 差 2.2×」（问题#6）的**新嫌疑方向**：单元类型改变骨的
  分担比例 ⇒ 韧带/关节的载荷分配变化 ⇒ 系统刚度变。
- **P1 骨折裕度表不受影响** ✓ —— 骨域 max/p99 = 1.34~1.35（无尖峰），
  与软骨的 17× 尖峰完全不同。
- **P3/P4/P5 不变**；**P2' 现在比 P3 更该先做**（它同时管刚度可信度与压缩上限）。
"""
t = io.open(TOD, encoding="utf-8").read()
if "P2 收口 —— 三个假设全被推翻" in t:
    print("· TODO 已有")
else:
    io.open(TOD, "a", encoding="utf-8").write(ADD)
    print(f"ANKLE_FIX_TODO +{len(ADD)} -> "
          f"{len(io.open(TOD, encoding='utf-8').read().encode('utf-8'))} B")

O = io.open(OVR, encoding="utf-8").read()
OLD = "| **骨折判据** | ❌ **未建** | **用户目标的真正缺口** ⇒ `ANKLE_FIX_TODO.md` 计划 P1 |"
NEW = ("| **骨折判据** | ✅ **已建（首版）** | 见 §第七轮：592.6N 下腓骨 11.40 MPa / 裕度 0.38；"
       "内翻 ⇒ 腓骨载荷份额 +57% |\n"
       "| **软骨网格健康度** | ✅ **健康** | det(J) min 0.77、无 <0.1mm 单元、末态只压 0.014mm "
       "⇒ §第八轮：原「退化单元」归因**作废** |")
if OLD in O:
    O = O.replace(OLD, NEW)
    io.open(OVR, "w", encoding="utf-8").write(O)
    print(f"总体进度 已更新 -> {len(O.encode('utf-8'))} B")
else:
    print("⚠ 总体进度未命中")

D = """
## ⚠️ P2 收口（第八轮）：三个假设全被推翻，柔度在韧带

**三推翻**：①软骨网格**健康**（det(J) min 0.77、≤0 的 0 个）②无超薄单元（最薄 0.683mm）
③软骨**没被压穿**（只压 0.014mm）。**改判**：σ_max=104.85 是**罚因子接触尖峰**
（E=12MPa 下 100MPa 需 800% 应变，物理不可能）。

**⭐ 真发现：柔度全在韧带** —— 逐域 dz 跨度：韧带 **0.85~1.65 mm** vs 骨 **0.01~0.10 mm**、
关节 **0.014 mm** ⇒ **300 倍** ⇒ 与 LRN-080（韧带 77% 载荷）独立吻合。
⇒ 状态表 #4「接触刚度 150 N/mm」其实**是韧带主导的系统刚度**；#3「压缩上限」也是**韧带极限**。

**计划变化**：P2「修退化软骨」**删除**；新增 **P2' 刚度/压缩上限=韧带问题**（优先于 P3）；
P1 骨折表不受影响（骨域 max/p99=1.34，无尖峰）。
"""
d = io.open(DAY, encoding="utf-8").read()
if "P2 收口（第八轮）" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(D)
    print(f"日记 +{len(D)} 字符")
