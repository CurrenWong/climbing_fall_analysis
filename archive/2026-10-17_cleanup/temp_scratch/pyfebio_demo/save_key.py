"""压缩前收口：把「关键问题」三层判断写进 CONV_TODO.md，并接好指针。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

# ── ① CONV_TODO.md 顶部插入「关键问题」节 ──
T = ROOT / "CONV_TODO.md"
s = T.read_text(encoding="utf-8")
KEY = """# CONV — 收敛研究：tet4 vs penta6 的 2.2 倍刚度差

---

## 🎯 关键问题（一句话 + 三层，**先读这一节**）

**模型的「刚度」有约 2 倍不确定度、且查不出原因 ⇒ 无法判断 600 N 那组结果可不可信。**

### 第一层（表层）
`tet4 259.3 N` vs `penta6 116.9 N`，同几何/同材料/同边界，仅软骨单元类型不同。
五个假设已排除（单元格式 / 体积锁死 / 几何 / 剔除 / 接触面积）
⇒ **至少有一个是错的，而我不知道是哪个。**

### 第二层（影响面）
下游结论全部挂在它上面：

| 想回答的问题 | 依赖 |
|---|---|
| 这套 FE 管线能不能当工具用？ | 载荷-位移关系可信吗 |
| 韧带排序差异是"机理"还是"量级不够"？ | 位移尺度可信吗 |
| 双靶验证（2.96 MPa / 287 mm²）算通过吗？ | 应力近似自归一，但**接触状态**受刚度影响 |
| 要不要去验 §3 的内旋机理？ | 在错误刚度上验机理没意义 |

### 第三层（最根本）
**没有独立的地面真值。**
- 吴恺韧带位移（1.92~5.24 mm）复现不了（差 30~40 倍），已判定不可用
- 应力/接触面积只能说"量级对上、偏差 ~20%"
⇒ **即使解决 2.2 倍，验证仍然是薄的。**

---

## 🎯 最尖锐的、还没查的候选：**接触面的「定向」**

**矛盾特征**（这是整个排查里唯一的硬矛盾）：

| | tet4 | penta6 |
|---|---|---|
| 接触面积 @0.20 mm | 238.3 mm² | 228.6 mm²（差 4~9%） |
| 同位移软骨应力 p99 | 0.765 MPa | 0.273 MPa（**差 2.8 倍**） |

**面积对、行为不同** ⇒ 典型症状是**法向定向 / 主从面配错**
（FEBio `sliding-elastic` 接触对 primary/secondary 定向敏感）。

**线索**：`cart_outer` 里两条分支算朝向用的参考点**不一样**
```python
tet4 路径：rest = 对面那个 tet 节点（单点）
prism 路径：rest = 底面三角形质心
```
两者算出的法向**理论上应一致，但从没验证过**。

---

## 🎯 方法论问题：我一直在用「拆开排查」

单元素 patch test / 切出体积 / 单独算面积 / 扫 ν —— 每个各排除一个因素，
但**如果真因是某种交互**（接触 × 单元 × 装配），**拆开就永远找不到**。

⇒ 下一步必须是**不拆的对比**：把 tet4 导出的接触面强行喂给 prism 模型
（加 `CONTACT_FROM` 开关）。这是唯一能一刀切开
「几何生成」还是「接触行为」的实验。

---

## 🔜 建议顺序

1. **交换接触面**（最便宜、最可能一击解决）
2. 若还不是 ⇒ 逐单元刚度贡献对比（同位置单元在原模型里单独加载）
3. 若仍查不出 ⇒ **承认第三层**，去找真正的地面真值
   （吴恺模型文件 / 其他足踝 FE 基准数据），
   而不是继续在内部自洽性上打转
"""
s = s.replace("""# CONV — 收敛研究：tet4 vs penta6 的 2.2 倍刚度差

**状态：五个假设全部逐一排除，仍未归因。**（2026-10-04）
**这是负面结果，如实记录；下面每个"排除"都有实测支撑。**
""", KEY + """
---

**状态：五个假设全部逐一排除，仍未归因。**（2026-10-04）
**这是负面结果，如实记录；下面每个"排除"都有实测支撑。**
""")
T.write_text(s, encoding="utf-8")
print(f"CONV_TODO.md -> {T.stat().st_size:,} bytes")

# ── ② FE_PIPELINE_TODO.md：接手入口补一行指针 ──
FP = ROOT / "FE_PIPELINE_TODO.md"
f = FP.read_text(encoding="utf-8")
if "CONV_TODO" not in f:
    L = f.split("\n")
    hi = next(i for i, l in enumerate(L) if "接手入口" in l or "\U0001f9ed" in l)
    L.insert(hi + 1, "")
    L.insert(hi + 2,
             "> ⚠️ **刚度有 2 倍量级不确定度（未归因）**：同几何/材料/边界下"
             "tet4 259.3 N vs penta6 116.9 N。五个假设已排除，最尖锐候选是"
             "**接触面定向**。600 N 档的**应力/接触面积**双靶可用，"
             "**载荷-位移关系不可信** → `CONV_TODO.md`（顶部「关键问题」节）")
    FP.write_text("\n".join(L), encoding="utf-8")
    print(f"FE_PIPELINE_TODO.md -> {FP.stat().st_size:,} bytes（已插指针）")
else:
    print("FE_PIPELINE_TODO.md 已有 CONV 指针")

# ── ③ MEMORY.md：看余量，必要时等量替换 ──
M = WS / "MEMORY.md"
raw = M.read_bytes()
print(f"\nMEMORY.md 磁盘 {len(raw):,} / 10240  (余量 {10240 - len(raw):+,})")
print("  含 CONV 指针:", "CONV_TODO" in M.read_text(encoding="utf-8"))
