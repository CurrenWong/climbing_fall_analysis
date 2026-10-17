"""记 LRN（旧数字不可跨配置沿用 + 三轮假设全错模式）+ 最终校验。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
GOT = os.path.join(W, "references", "operational-gotchas.md")

t = io.open(LEARN, encoding="utf-8").read()
ids = re.findall(r"\[LRN-(\d{8})-(\d{3})\]", t)
dates = sorted({d for d, _ in ids})
today = "20261006"
nums = sorted(int(n) for d, n in ids if d == today)
print(f"今日 LRN 号段: {nums[:3]}..{nums[-3:] if nums else []}  已用 {len(nums)}")
nxt = (max(nums) + 1) if nums else 1
print(f"⇒ 下一个可用号: LRN-{today}-{nxt:03d}")

ENTRY = f"""

---

## [LRN-{today}-{nxt:03d}] 旧配置的数字不能用在当前配置上 —— 三轮假设全错的根因

**日期**: 2026-10-06
**类别**: correction + methodology
**项目**: climbing_fall_analysis / 踝关节 FE

### 发生了什么
为修"退化软骨单元"，连续提出并实测了 **三个假设，全被自己的数据推翻**：

| 假设 | 实测反驳 |
|---|---|
| 软骨网格有近退化单元 | `tibiotalar_talus__penta6` det(J) **min=1.2245**；全模型 det(J)≤0 的单元 **0 个** |
| 存在超薄单元（t=0.5×gap 的产物） | 厚度 min **0.6828 mm**，<0.1 mm 的 **0 个** |
| 软骨被压穿 ⇒ 单元翻转 | +1.6mm 末态只压 **0.0142 mm**，det(J) 末态 min 0.6116 |

### 根因（两条，都是可复用的）
1. **跨配置沿用旧数字**：我一直在引用「软骨近退化单元 det(J)=2.39e-3」——
   那是**旧软骨网格**（`cart_hex`/早期 `cartilage.py` 路线）的数字，
   而当前是 `cart_patch m=1` 网格。**同一模型的两个配置，网格质量完全不同。**
   ⇒ 规则：**引用任何网格质量数字前，必须说明它属于哪个配置，并在当前配置重测。**
   （同类已发生过：压缩标定曲线在新旧软骨间差 24%，见 §第六轮。）
2. **先动手修、后验证前提**：我连续三轮都在"修"，而没有先问
   **"前提成立吗？网格到底健康吗？"**。第一个诊断（det(J) 扫描）成本仅几秒，
   却能把三轮全部省掉。**⇒ 规则：修之前先测前提；前提为假时，任务是错的，不是方法错。**

### 本轮真正有价值的产物
- **σ_max=104.85 改判为罚因子接触尖峰**（判据与约定无关：E=12 MPa 下 100 MPa
  需约 800% 应变 ⇒ 物理不可能）⇒ **软骨应力口径用 p99**。
- **柔度全在韧带**（逐域 dz 跨度：韧带 0.85~1.65 mm vs 骨 0.01~0.10 mm，关节 0.014 mm）
  ⇒ 状态表 #4「接触刚度 150 N/mm」重新解释为**韧带主导的系统刚度**。

### 与既有规则的关系
这是 `[LRN-20261006-083]`（方向/符号须有不可被约定翻转的证据）的**同族**：
那条管"方向类结论"，这条管"**数值类结论的配置依赖**"。
两条共同的元规则：**结论必须绑定它的证据与配置；换配置就要重测。**
"""
if f"[LRN-{today}-{nxt:03d}]" in t:
    print("· LRN 已有")
else:
    io.open(LEARN, "a", encoding="utf-8").write(ENTRY)
    t2 = io.open(LEARN, encoding="utf-8").read()
    ids2 = re.findall(r"\[LRN-(\d{8})-(\d{3})\]", t2)
    print(f"LEARNINGS +{len(ENTRY)} 字符 / 现在 {len(ids2)} 条 / "
          f"{len(t2.encode('utf-8'))} B")

G = """

## 🔁 「数值类结论的配置依赖」（[LRN-20261006-084]）

**引用任何网格质量 / 标定数字前，必须①说明它属于哪个配置 ②在当前配置重测。**
本项目实测两次跨配置失效：
- 软骨网格质量：旧网格 det(J)=2.39e-3 ⇒ 当前 `cart_patch m=1` 是 **1.22**（健康）
- 压缩标定：新软骨 `CART_NEW+TIE+cart_m1` 与旧格式同位移差 **24%**

**并且：修之前先测前提。** 前提为假时，任务是错的，不是方法错
（曾连提 3 个假设全被推翻，而第一个诊断只花几秒）。
"""
gg = io.open(GOT, encoding="utf-8").read()
if "数值类结论的配置依赖" in gg:
    print("· gotchas 已有")
else:
    io.open(GOT, "a", encoding="utf-8").write(G)
    print(f"operational-gotchas +{len(G)} 字符")

print("\n" + "=" * 74)
print("最终校验")
print("=" * 74)
ok = True
CHK = [(os.path.join(W, "MEMORY.md"), 10240),
       (os.path.join(W, "references", "operational-gotchas.md"), None),
       (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
       (os.path.join(R, "docs", "总体进度.md"), None),
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
    ii = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", x, re.M)
    dup = sorted({i for i in ii if ii.count(i) > 1})
    print(f"  {os.path.basename(f)}: {len(ii)} 条，重号 {dup or '无 ✓'}")
    ok &= not dup
need = [(os.path.join(R, "ANKLE_FIX_TODO.md"),
         ["P2 收口 —— 三个假设全被推翻", "柔度全在韧带", "P2'"]),
        (os.path.join(R, "docs", "总体进度.md"),
         ["骨折判据", "软骨网格健康度"]),
        (os.path.join(W, ".learnings", "LEARNINGS.md"), ["旧配置的数字"])]
for p, ks in need:
    if not os.path.exists(p):
        continue
    tt = io.open(p, encoding="utf-8").read()
    miss = [k for k in ks if k not in tt]
    print(f"  {os.path.basename(p):<32} 缺失 {miss or '无 ✓'}")
    ok &= not miss
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
