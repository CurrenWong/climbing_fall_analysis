"""追加 §七：收敛研究（查清"不是什么"+ 定位"是什么"）。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"

ADD = """

---

## 七、收敛研究：2.2 倍刚度差 —— 查清了「不是什么」（2026-10-04）

> 详见 **`CONV_TODO.md`**；learnings LRN-033/034/035

### ✅ 结论 1：**不是单元格式问题**（单元级测试证明）

单楔形单向压缩 vs 解析解 `F = E·A·ε`（E=10, **ν=0** ⇒ 去掉体积耦合）：

| 单元 | F (N) | F/F_exact |
|---|---|---|
| **penta6** | 9.702 | **0.970** |
| **tet4×3**（WEDGE_T 拆法） | 9.702 | **0.970** |
| penta15 | 0.791 | 0.079 ← 我的实现坏了 |

**薄/歪斜扫描**（h 0.2~5 mm、顶面平移 0~1.5 mm）：
**5 种情形 penta6/tet4×3 比值全部 = 1.000**（含 h=0.2、skew=(1.5,0.5)）。

⇒ **两者等价，「哪个更准」是伪命题。**
**此前把 2.2 倍归因于「tet4 体积锁死」是错的** —— ν 扫描也支持：
**ν=0 时比值已是 1.99**，锁死只贡献 1.99→2.38。

### ✅ 结论 2：最可能来自**软骨网格退化 + 差别剔除**

| 单元 | 软骨总数 | `DET_MIN` 剔除 | 剔除率 | 接触面 |
|---|---|---|---|---|
| tet4 | 1,052 | 22 | **2%** | 148/191 |
| penta6 | 315 | 43 | **14%** | 138/171 |

- **同一个 `DET_MIN` 对不同单元不等效**：3 个 tet 各有自己的 det(J)、
  只剔最差那个；楔形是**整块一个** det ⇒ 一脚退化丢整块。
- 剔除位置在**偏移面最薄处** ⇒ 打穿**接触承载路径**。
- 但阈值 0.5→0.02 时剔除 43→34、**载荷不变**（117 N）
  ⇒ 被剔的是**真退化**（det≤0.02），闸门不冤。
- 关韧带后比值 2.22 → **3.10**（212.9 vs 68.8 N）⇒ 韧带也在掩盖差异。

### 🔜 下一步：根因在**软骨生成**，不在单元类型

`scripts/ankle_fe/cartilage.py` 用**骨面顶点法向**偏移；曲率大/尖锐处
偏移方向与面不垂直 ⇒ 楔形扭曲退化。
1. 诊断偏移质量（每楔形 det(J) + 偏移方向与面法向夹角）
2. 改用**面法向**偏移，或对退化处平滑/重投影
3. 目标：`DET_MIN=0.5` 下软骨剔除率降到 tet4 同级（~2%）
4. 再重跑对比 —— 若一致，**600 N 档数值才可信**

### 📌 本轮确定/修好的 FEBio 语法与自 bug

- `penta15` 节点序 = **VTK 序**（角 0-5 + 边中点 6=01,7=12,8=20,9=34,10=45,11=53,12=03,13=14,14=25）
- prescribed displacement：`<dof>x</dof><value>V</value><relative>0</relative>`，
  **每个自由度一条 BC**；`<x_dof>1</x_dof>` 只属 `zero displacement`
- 探针 .feb 必须带 `<solver type="solid">`（否则报 solver 缺失，易误判）
- `<value lc="1">` 无 `<LoadData>` 时报 `Invalid load curve ID`
- **pyfebio 补丁抽成 `temp/pyfebio_demo/fepatch.py`** 复用（忘了打补丁会
  `KeyError: np.int32(4)`）

### ⚠️ 遗留
- **penta15 实现有 bug**（单元测试 0.079），暂不能用它做 p 收敛
- 2.2/3.1 倍的**定量归属**未定（剔除 vs 接触面形状）

### 📁 本轮产物
`CONV_TODO.md`（含全部结论表）· `temp/pyfebio_demo/fepatch.py`（复用补丁）
· `elem_test2.py` + `thin_test.py`（单元级验证）· `bisect_contact2.py`
· `temp/pyfebio_demo/elem_test2/` 与 `thin_test/`（探针 .feb 存档）
· `scripts/ankle_fe/ankle_lig_feb.py` +`CART_E`/`CART_NU`、+`prism15`（有 bug）
"""

t = M.read_text(encoding="utf-8")
M.write_text(t.rstrip() + ADD, encoding="utf-8")
print(f"{M.stat().st_size:,} bytes")
