"""追加 §六：600 N 打通 + 双靶验证 + 两个自查 bug。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"

ADD = """

---

## 六、✅ 600 N 打通 + 双靶验证（2026-10-04，承接 §五 的"下一步 1"）

> 详见 **`ANKLE_600N_TODO.md`**（含遗留问题）；结果并入
> `docs/吴恺论文_重解析关键信息.md` §八

### 🎯 成果：模型从「0.378 mm 就翻」到「645.7 N」

```
DRIVER=disp TOPFIX=lateral TRUNC_MM=60 FIX_TALUS=all
CART_ELEM=prism TIME_STEPS=30 DISP_MM=1.20
→ rc=0 / 7.4 s / 支反力 645.7 N
```

| 靶子 | 吴恺 | 本模型 | 偏差 |
|---|---|---|---|
| 最大接触应力 | 3.74 MPa | 软骨 von Mises max **2.96 MPa** | **−21%** |
| 接触面积 | 335.5 mm² | **287.4 mm²**（间隙<0.20 mm） | **−14%** |

**两个靶子都在同一量级、偏差 ~20% 以内** —— 比只比软骨应力 p99 硬得多。

### 🔧 根因与修法

**根因**：薄软骨（0.703 mm）的 6 节点楔形被**拆成 3 个 tet4** ⇒ 细长 sliver
⇒ 末态 det(F) 跌到 **0.523**（tibiotalar_tibia）⇒ 再压就翻负。
**修法**：`CART_ELEM=prism` —— 直接产出 **penta6**，不拆。

排查路径（重要）：
1. FEBio **不报单元号** ⇒ 自己读末态位移**逐单元算 det(F)** 定位
2. 试 `hex8`（四边形底）：歪斜 det(J)≈0.5，剔 42/61 ⇒ 软骨破洞、
   接触面 150/198→**51/64**，弃用
3. 改「四边形也拆 2 个 penta6」⇒ 接触面回到 **150/198** ✓

### 📌 本轮修好**两个我自己的 bug**（教训比结果更值钱）

**① hex8 的 det(F) 基向量必须从 node0 出发**（`(0,1),(0,3),(0,4)`）
我原先用 `(0,1),(2,3),(4,7)`，两对近乎平行 ⇒ det(G)≈0
⇒ **纯刚体平移都算成 det = −1.34**，**误报腓骨单元翻转**，白追一轮。
★ **判据**：喂一个纯刚体平移的假位移，**det(F) 必须恰好 = 1**。

**② 楔形的 Jacobian 必须采样**（沿 r,s,t 扫）
底面形心 det>0 但三角内部/沿 t 退化 ⇒ FEBio 在积分点报
`Negative jacobian detected during mesh initialization`。
且 `dt` 列必须是**上下三角形心差** `(P3+P4+P5-P0-P1-P2)/3`，写成 `P3-P0` 会算错。

### 📌 FEBio 4.13 的 `penta6` 节点序（实测，别再猜）

**底三角 `0,1,2` + 顶三角 `3,4,5`（同向对应）**
6 种排列实测，只有 `[1,2,3,4,5,6]` 通过。
（探针 .feb 必须自带 `<solver type="solid">`，否则报
`Component "" needs to have property "solver" defined`，会误判成元素序问题。）

### ⚠️ 遗留：楔形 vs 四面体刚度差 **2.2 倍**（未解决）

同位移 0.1472 mm：

| 软骨单元 | 支反力 | 软骨 p99 |
|---|---|---|
| tet4 | **259.3 N** | 0.765 MPa |
| prism | **116.9 N** | 0.273 MPa |

同几何/材料/边界 ⇒ 差异只能来自**离散格式**（线性 tet4 在 ν=0.4 下体积锁死偏刚；
线性 penta6 偏柔/可能沙漏）。**哪个更接近真值未定 ⇒ 需网格收敛研究**。
在此之前 600 N 档数值只能当量级参考。

（已排除：`DET_MIN` 0.5→0.02 几乎无变化，117 N 不变 ⇒ **元素剔除不是主因**。）

### 📁 本轮新增/修改

| 文件 | 变化 |
|---|---|
| `scripts/ankle_fe/ankle_lig_feb.py` | +`CART_ELEM`（tet4/prism）、+`DET_MIN` 可调、修 `detJ_penta6`（采样）、修 hex8 det(F) 基向量 |
| `temp/pyfebio_demo/diag_detF3.py` | 逐单元 det(F) 诊断（定位翻转单元） |
| `temp/pyfebio_demo/final_600.py` | 600 N 档 + 接触面积 + 应力 一站式 |
| `temp/pyfebio_demo/bench.py` | 配置扫描基准（纯 Python，**不用 .bat**） |
| `temp/pyfebio_demo/probe_penta2.py` | penta6 节点序实测探针 |
| `ANKLE_600N_TODO.md` | 本轮计划/结果/遗留 |
| `.learnings` | LRN-030/031/032 |

### 📌 工程教训：**别用含中文的 .bat**

cmd 按 GBK 读批处理，UTF-8 的中文会把行解析错乱（实测两次：
`'G_MODE' is not recognized`）。⇒ **改用 Python 脚本设 `os.environ` +
`subprocess`** 驱动，彻底绕开。
"""

t = M.read_text(encoding="utf-8")
M.write_text(t.rstrip() + ADD, encoding="utf-8")
print(f"{M.stat().st_size:,} bytes")
