"""记录第二轮的新教训到 .learnings/。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
E = ROOT / ".learnings/ERRORS.md"
L = ROOT / ".learnings/LEARNINGS.md"

et = E.read_text(encoding="utf-8")
adde = """

---

## ERR-20261004-004 — 求解器配置问题被误当成物理问题（"每步 97 次迭代磨不完"）

**严重度**: 中（浪费大量时间）
**现象**: 加了韧带的踝关节模型在力控制下：
- 600 N 时 **首个载荷步（1/220 载荷）** 就报 22 个 negative jacobian；
- 降到 6 N 能跑 3 步，但报 **"Max nr of iterations / Max nr of reformations reached"**；
- 每步 **~97 次牛顿迭代**，3 步花 **407 s**。

**误判过程**（我的错误链条）：
1. 以为是**刚体模态**（LRN-016/017 的教训）⇒ 去固定舟骨、去写转动约束 BC —— 都没用
2. 以为**韧带壳与骨共用节点**引入自由转动自由度 ⇒ 改成只用骨节点的条带/三角网 —— 也没用
3. 以为 **`elastic-truss`** 是正解 ⇒ 真实模型读文件就崩（见 ERR-20261004-001）

**真因**（读 FEBio 迭代日志才看出来）：
```
convergence norms :     INITIAL         CURRENT         REQUIRED
   residual            9.840343e+06    3.239042e+01    0.000000e+00   ← 要求降为 0，永不可能
   displacement        1.870091e-01    4.064515e-02    1.257888e-07   ← 容差过严，差 5 个数量级
```
- **残差判据 REQUIRED = 0** ⇒ 该判据形同废除（**同一个坑 MEMORY 里已记过一次**）
- 实际只能靠位移判据，而它 **1.26e-07 过严** ⇒ 死磨
- 缺 `<qn_method type="BFGS">`（FEBio 的标准拟牛顿加速）

**修法**：按手册 §3.3 写全 `<solver>`：`max_refs/diverge_reform/reform_each_time_step/
dtol/etol/rtol/lstol/min_residual/rhoi/qn_method(BFGS)` + `<time_stepper>` 的
`dtmin/dtmax/max_retries/opt_iter`。

**教训**：
1. **"失败在 1/220 载荷"+"每步上百次迭代" ⇒ 先查求解器配置，别先改物理/几何**
2. **必须读 FEBio 的迭代日志**（`convergence norms` 那几行），它直接告诉你哪条判据没满足；
   只看"ERROR TERMINATION"会盲猜
3. **`REQUIRED = 0` 是个红旗**：任何"要求降到 0"的收敛判据都不可能满足

---

## ERR-20261004-005 — pyfebio 读 shell/beam 域的 .xplt：**两处** KeyError

**严重度**: 低（有补丁）
**现象**: `pyfebio.xplt.to_hdf5()` 对含 shell/beam 域的模型失败：
1. `parse_mesh` → `KeyError: 'name'`
   （`domain["name"]` —— 注意是**双引号形式**，只替换 f-string 里的 `{domain['name']}`
   **不会命中**，本轮因此白试一轮）
2. `parse_state` → `KeyError: np.int32(4)`
   （`mesh_dict[set_lut[key]][set_id]` —— **`set_id` 与 `set_id + 1` 两处都要兜底**）

**修法**（`temp/pyfebio_demo/read_lig.py` 已实现）：用 `inspect.getsource` 取源码、
字符串替换加 `.get(...)` 兜底、`exec` 回写模块。

**附带**：`<Nodes>` 的坐标是**结构化数组** `(id, x, y, z)`，不是 `(n,3)`；
要 `np.column_stack([c["x"], c["y"], c["z"]])`。
"""
E.write_text(et + adde, encoding="utf-8")

lt = L.read_text(encoding="utf-8")
addl = """

---

## LRN-20261004-027 — 提速的正解是「把模型做小且做对」，不是「调大迭代上限」

**类别**: best_practice
**场景**: 加了韧带的踝关节模型每轮求解要 **900 s 都跑不完**，无法迭代。

**我一开始想的**：调大 `<max_refs>`（25→100）、放宽容差、加 BFGS —— 有效但**仍是分钟级**。

**真正的杠杆**：**模型里有 25,716 个单元是"整根胫骨"**（松质 6651 + 骨髓 14482）。
而**踝关节模型根本不需要整根胫骨** —— 吴恺的模型就是"足踝部"。

**截断到踝关节近端 60 mm 后**：**32,488 → 9,109 单元**（3,588 节点），
求解 **900+ s → 5.4 s**（**>160 倍**），且载荷加在截面上正对应吴恺的
"胫腓骨**上截面** 600 N"。

**教训**：
1. **先问"这个模型真的需要这些单元吗"**，再去调数值参数。
   领域直觉（踝关节模型不需要整根胫骨）比求解器调参有用 100 倍。
2. **迭代太慢时，优先做"建模范围"的裁剪**，而不是买更大的算力或调容差。
3. 顺带：**截断必须在节点编号之前做**，否则被删单元会留下孤立节点。

---

## LRN-20261004-028 — 收敛失败时**必须读迭代日志的 `convergence norms`**

**类别**: best_practice
**场景**: FEBio 报 `ERROR TERMINATION`，只看到"failed to converge"。

**只看结论**会盲猜（我猜了刚体模态、猜了壳的转动自由度、猜了 truss 语法，
全都是错的）。**读日志里这几行**立刻定位：
```
convergence norms :     INITIAL         CURRENT         REQUIRED
   residual            9.840343e+06    3.239042e+01    0.000000e+00
   energy              1.636477e+02    6.190982e-02    1.636477e+00
   displacement        1.870091e-01    4.064515e-02    1.257888e-07
```
- 哪一行 `CURRENT > REQUIRED` ⇒ 就是那条判据没过
- **`REQUIRED = 0` ⇒ 该判据不可能满足**（FEBio 里 `rtol=0` 时残差判据形同废除）

**教训**：**求解器的日志比错误信息信息量大得多**。
错误信息只说"失败了"，日志说"哪条判据、差多少"。
（附带：本轮还踩了 `"N O R M A L   T E R M I N A T I O N"` 字母间带空格、
用 `==` 判断导致把**4 个成功算例全判成失败** —— 判断前先 `replace(" ", "")`。）

---

## LRN-20261004-029 — 复现文献数字前，先确认**指标定义**与**运动学方向**

**类别**: knowledge_gap
**场景**: 复现吴恺 2012 的"韧带位移（拉伸长度）：胫跟 5.24 > 胫舟 5.04 >
跟腓 2.05 > 下胫腓前 1.92 mm"。

**本模型（严格照他所述的边界：600 N 压缩 + 固定跟骨 + 约束距骨）得到**：
ΔL = **0.001 ~ 0.08 mm**，且 **`胫跟` 缩短最多（最后一名）** ——
与他的"胫跟拉伸最大（第一名）"**方向相反**。

**物理解释**：纯**轴向压缩**下，胫骨朝被固定的足骨**靠近** ⇒
所有"胫骨→足骨"的韧带**都在缩短**。要让 `胫跟` **拉伸** 5 mm，
胫骨必须**远离**跟骨（关节分离/内侧张开）—— **纯压缩不可能做到**。

**教训**：
1. **差 2 个数量级 + 名次颠倒 ⇒ 先怀疑"指标定义不同"，而不是"我的模型不准"**。
   "位移（拉伸长度）"至少有三种读法：①长度变化 ②节点位移模 ③两点相对位移向量模。
2. **论文正文的数字可能来自图**，纯文本抽取拿不到（本次 `图 4` 的数值就不在文本里）
   ⇒ **复现前先确认能拿到原始数值**，否则没有锚点。
3. **照抄文献的边界条件时，要问"这个边界下运动学会往哪个方向走"**：
   若算出的方向与文献相反，先检查是否缺了一个不在文字里的约束/加载细节。
"""
L.write_text(lt + addl, encoding="utf-8")
print("已追加 ERR-004/005 与 LRN-027/028/029")
print(f"ERRORS {len(E.read_text(encoding='utf-8').encode('utf-8')):,} | "
      f"LEARNINGS {len(L.read_text(encoding='utf-8').encode('utf-8')):,}")
