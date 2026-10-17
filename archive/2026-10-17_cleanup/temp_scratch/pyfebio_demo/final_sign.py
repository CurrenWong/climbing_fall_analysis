"""收尾：TODO 追加 + 日记 + 全量校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
TOD = os.path.join(R, "ANKLE_FIX_TODO.md")
DAY = os.path.join(W, "memory", "2026-10-06.md")

TADD = """

### 🚨 2026-10-06 第三轮：**符号真相推翻上轮**（`[ERR-20261006-003]`）+ 真压缩曲线

**`DRIVER=dispx` 下 `DISP_MM > 0` 才是压实**。上轮把「正解」定成 `DISP_MM<0`，**方向判反了**。

**三条独立证据**
1. `.feb` 实体：写入的 z 值 = **`−DISP_MM`**
2. 几何：顶面 z̄=−653.5、距骨 z̄=−709.4（距骨在**下方**）⇒ `DISP=+1.0` 顶面 dz=**−0.927**（朝距骨）
3. 接触力（距骨全固定·无韧带）：`+1.0` ⇒ **149.59 N**；`−1.0` ⇒ **0.0000 N**
4. 闭环：历史 `DISP_MM=+1.6` 的 **592.6020 N 逐位复现** ⇒ 历史约定本来是对的

**⇒ 真·压缩曲线**：+0.25/+0.5/+1.0/**+1.6**/+2.5 mm → 124.0/230.7/418.1/**592.6**/615.4 N
⇒ **600 N 在 `DISP_MM≈+1.89 mm` 达成** —— **「600 N 到不了」作废**；`[LRN-072]` 的 −132…−1017 N
是**拉伸**值（由韧带撑）。

**载荷路径**（`[LRN-080]`）：压缩下**韧带 77% / 接触 28%**；`FIX_TALUS=transverse` 时几乎全走韧带
⇒ 路径对照必须 `FIX_TALUS=all`（否则无韧带时近奇异）。`shell` 模式用固定 `E=260`、**不用 `a_cross`**
⇒ 本轮的韧带几何修复（`[ERR-20261006-004]`）**只影响 `truss` 模式**。

**下一步**
1. 用**正确符号**重跑 C2/C3 姿势矩阵（现有姿势数据是在**拉伸**方向测的，需重测）
2. 查**接触为何只有 ~150 N/mm**（预期 1000+）—— 这是让接触成为主路径的前提
3. 再接 OpenSim 姿势输入
"""
t = io.open(TOD, encoding="utf-8").read()
if "符号真相推翻上轮" in t:
    print("· TODO 已存在")
else:
    io.open(TOD, "a", encoding="utf-8").write(TADD)
    print(f"ANKLE_FIX_TODO +{len(TADD)} 字符")

DADD = """
## 🚨 第三轮：符号真相 —— 我上一轮的"修正"是把方向搞反了

**做 ②（修韧带插入点几何）时意外挖出更根本的错误。**

### 起点：`DISP=+1.6` 在修复韧带前后都给 592.6020（逐位相同）
追下去发现三件事：

**1. 韧带几何修复本身有效，但只影响 `truss` 模式**
- `ligaments.py` 两处真 bug：① 纤维轴用「点云主成分」⇒ 对宽短韧带与骨-骨方向差
  **PTCL 89.3°、SYND 84.0°**（改成骨-骨方向后，PTCL 刚度修正 **2.81×**、SYND **1.37×**；
  CFL/ATFL/MD* 只差 0.7~4.8° ⇒ 1.00×）② 配对退化成**扇形**（ATFL 4 根杆全汇到同一腓骨节点）
  ⇒ 改为 1:1 + 就近补足，10 条韧带足迹覆盖 **100%**
- **但 `ankle_lig_feb.py` 里 `a_cross` 只在 `LIG_MODE=truss` 用**；默认 `shell` 用固定 `E=260`
  ⇒ **FE 结果不变**（592.6020 逐位复现的原因）

**2. 顺手纠正自己的错判**：CFL 其实一直是对的（轴夹角 0.7°、刚度 1.00×）——
我先前引用的"CFL 34.33 vs 解剖 28.5 ⇒ 稀释 20%"是把**杆长**（弯曲膜的直线弦）当成了
**膜沿纤维的延伸**（面积公式的正确分母）⇒ 该 caveat 撤回。

**3. 真正的大发现：`DISP_MM > 0` 才是压实**
| 证据 | 结果 |
|---|---|
| `.feb` 写入值 | `z = −DISP_MM` |
| 顶面 dz（距骨在下方）| `+1.0` ⇒ **−0.927**（朝距骨=靠近）；`−1.0` ⇒ **+1.000**（远离）|
| 接触力（距骨全固定·无韧带）| `+1.0` ⇒ **149.59 N**；`−1.0` ⇒ **0.0000 N** |
| 历史闭环 | `DISP=+1.6` ⇒ **592.6020 N 逐位复现** |

**⇒ 上轮 `[ERR-20261005-009]` 的"正解 `DISP_MM<0`"方向判反、已作废；`LOAD_SIGN` 对 `dispx` 无效。**

### 真·压缩曲线（全 NORMAL）
| DISP_MM | +0.05 | +0.10 | +0.25 | +0.50 | +1.00 | **+1.60** | +2.50 |
|---|---|---|---|---|---|---|---|
| 力 (N) | 36.08 | 58.25 | 124.04 | 230.71 | 418.10 | **592.60** | 615.35 |

**⇒ 600 N 在 `DISP_MM≈+1.89 mm`** ⇒ **「600 N 到不了」的旧结论作废**；
`[LRN-20261005-072]` 的 −132…−1017 N 是**拉伸**值（韧带撑），**已作废**。

### 载荷路径分解（距骨全固定·压缩 +1.0mm）
| 配置 | 力 |
|---|---|
| 韧带开·接触开 | 541.61 N |
| **韧带关** | 149.59 N ⇒ 接触贡献 |
| **接触关** | 414.58 N ⇒ 韧带贡献 |

**⇒ 压缩下韧带占 77%、接触 28%**；`FIX_TALUS=transverse` 时几乎全走韧带（去韧带 ⇒ 0.0042 N，
且小 `DISP` 时 `|u|max` 冲到 10⁶ ⇒ **近奇异，该对照只能用 `FIX_TALUS=all`**）。

### 落盘
`[ERR-20261006-003]`（符号）·`[ERR-20261006-004]`（ligaments 两 bug + 影响范围）·
`[LRN-20261006-079]`（真压缩曲线）·`[LRN-20261006-080]`（载荷路径）·
MEMORY 改两处 + 压到 **10231/10240 ✅** · 使用说明 §2 换成真压缩曲线。

### 下一步
1. **用正确符号重跑 C2/C3 姿势矩阵**（现有姿势数据是在**拉伸**方向测的，需重测）
2. **查接触为何只有 ~150 N/mm**（预期 1000+）—— 让接触成为主路径的前提
3. 再接 OpenSim 姿势输入
"""
d = io.open(DAY, encoding="utf-8").read()
if "符号真相" in d:
    print("· 日记已存在")
else:
    io.open(DAY, "a", encoding="utf-8").write(DADD)
    print(f"日记 +{len(DADD)} 字符")

# ── 全量校验 ──
print("\n" + "=" * 74)
print("最终校验（解码字节；⚠️ 与 getsize 的差 = 行数 CR）")
print("=" * 74)
ok = True
import ast
CHK = [(os.path.join(W, "MEMORY.md"), 10240),
       (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
       (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
       (os.path.join(R, ".learnings", "ERRORS.md"), None),
       (os.path.join(R, "scripts", "ankle_fe", "ligaments.py"), None),
       (os.path.join(R, "scripts", "ankle_fe", "ankle_lig_feb.py"), None),
       (os.path.join(R, "scripts", "ankle_fe", "check_c2.py"), None),
       (os.path.join(R, "scripts", "ankle_fe", "lig_strain.py"), None)]
for p, cap in CHK:
    if not os.path.exists(p):
        print(f"  ✗ 缺失 {os.path.basename(p)}")
        ok = False
        continue
    n = len(io.open(p, encoding="utf-8").read().encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    print(f"  {os.path.basename(p):<34} {n:>8} B{f}")
for s in ("ligaments.py", "ankle_lig_feb.py", "check_c2.py", "lig_strain.py",
          "acceptance.py"):
    ast.parse(io.open(os.path.join(R, "scripts", "ankle_fe", s),
                      encoding="utf-8").read())
print("  AST ✓ ×5")
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    t2 = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", t2, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup
u = io.open(os.path.join(R, "docs", "踝FE可用模型使用说明.md"),
            encoding="utf-8").read()
need = ["DISP_MM > 0", "真·**压缩**曲线", "592.60", "+1.89 mm", "载荷路径",
        "149.59", "414.58"]
miss = [k for k in need if k not in u]
print(f"  使用说明关键内容缺失：{miss or '无 ✓'}")
ok &= not miss
mm = io.open(os.path.join(W, "MEMORY.md"), encoding="utf-8").read()
need2 = ["符号真相", "DISP_MM>0", "真·压缩曲线", "+1.9 mm", "韧带 77%"]
miss2 = [k for k in need2 if k not in mm]
print(f"  MEMORY 关键内容缺失：{miss2 or '无 ✓'}")
ok &= not miss2
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
