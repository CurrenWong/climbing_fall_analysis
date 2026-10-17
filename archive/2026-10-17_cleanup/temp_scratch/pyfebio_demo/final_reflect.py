"""反思落盘 + 最终校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
DAY = os.path.join(W, "memory", "2026-10-06.md")

DADD = """
## 🪞 反思：用户问"这条路是不是之前走过" —— 是，而且是我自己的老毛病

### 事实：方向符号翻了两次，最后**回到了最初约定**
| 时间 | 我下的结论 | 结果 |
|---|---|---|
| 10-04 | 「力控制病态」 | 后被推翻 |
| 10-05 | **`DISP_MM<0` 才是压实**（`[ERR-20261005-009]`）+ 拉伸标定表 + 「600 N 到不了」+ 整套姿势/判据 | **错** |
| 10-06 | **`DISP_MM>0` 才是压实**（`[ERR-20261006-003]`） | **回到最初约定** |

**⇒ 我离开了本来正确的约定，在一个"反过来的世界"里做了一到两轮完整工作，然后回到原点。**

### 根因
1. **用"约定相关量"推"方向"**：10-05 我的"决定性证据"是反力正负
   （`−132.56 N` 标成"压缩"、`+187.77 N` 标成"拉伸"）—— **标签是从符号读出来的**。
   当天本可做三条**与约定无关**的检查（各 1~2 次运行）：`.feb` 写了什么／材料点往哪动／
   固定距骨时哪个方向有接触力。
2. **不验证自己的验证器**：今天我的检查脚本自己把 `DISP=+1.0` 判成"拉开"（脚本内符号写反）；
   `compress_curve.py` 的 NORMAL/ERROR 是假阴性。**没跑已知正对照 ⇒ 不算证据。**

### 三级放大（比错本身更糟）
- 结论进了 **7 个文件的"务必先读这里"横幅** + 使用说明 §1 把 `LOAD_SIGN=+1` 标成「**必开**」⇒ **主动误导**
- 我"修"它时只**前置**新横幅、**没删旧的** ⇒ 3 个文件出现**两条相反横幅叠在一起**

### ⭐ 已写成规则（防复发）
- `[LRN-20261006-083]`：**方向/符号结论必须来自"不可被约定翻转"的证据**（结构/几何/物理/解剖）；
  **判据必须配已知正对照**；**更正要"替换"不要"叠加"**；翻转过的结论带"翻转史"
- 已提升到：`MEMORY.md` 坑表第 23 条｜`references/operational-gotchas.md`｜**`SOUL.md`**（方法准则）
  —— **⚠️ 改的是 SOUL.md（我的"灵魂"文件），需告知用户**

### 清理与审计结论
- 7 个文件的横幅已换成"带证据+翻转史"的版本；3 个文件的双横幅已删；
  使用说明 8 处（含 §1「必开」）已改；正文最后一处相反论断已加删除线
- **审计后仍安全**：C2/C3 的应变表是**几何量**（由节点位移算得）⇒ 与约定无关 ✓
- MEMORY 10222/10240 ✅
"""
d = io.open(DAY, encoding="utf-8").read()
if "反思：用户问" in d:
    print("· 日记已有")
else:
    io.open(DAY, "a", encoding="utf-8").write(DADD)
    print(f"日记 +{len(DADD)} 字符")

print("\n" + "=" * 70)
print("最终校验")
print("=" * 70)
import ast
ok = True
CHK = [(os.path.join(W, "MEMORY.md"), 10240), (os.path.join(W, "SOUL.md"), None),
       (os.path.join(W, "references", "operational-gotchas.md"), None),
       (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
       (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
       (os.path.join(R, ".learnings", "ERRORS.md"), None)]
for p, cap in CHK:
    n = len(io.open(p, encoding="utf-8").read().encode("utf-8"))
    f = ""
    if cap:
        f = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    print(f"  {os.path.basename(p):<32} {n:>8} B{f}")
for f in (".learnings/LEARNINGS.md", ".learnings/ERRORS.md"):
    x = io.open(os.path.join(R, f), encoding="utf-8").read()
    ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", x, re.M)
    dup = sorted({i for i in ids if ids.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ids)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup
# 双横幅检查
for rel in ("ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md", "TIE_PORT_TODO.md",
            "ANKLE_FIX_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
            "docs/总体进度.md"):
    t = io.open(os.path.join(R, rel), encoding="utf-8").read()
    n_old = t.count("2026-10-05 重大更正")
    n_new = t.count("翻转史")
    flag = "✓" if (n_old == 0 and n_new == 1) else f"⚠ old={n_old} new={n_new}"
    print(f"  横幅 {rel:<32} {flag}")
    ok &= (n_old == 0 and n_new == 1)
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
