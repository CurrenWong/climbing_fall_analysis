"""修复 LRN 真实重号：后追加的三条改号 054/055/056（保住早先的为规范号）。"""
import io
import re

LEARN = r"D:\Project\climbing_fall_analysis\.learnings\LEARNINGS.md"
t = io.open(LEARN, encoding="utf-8").read()
print("改号前各号出现次数：",
      {k: t.count(k) for k in ("LRN-20261004-027", "LRN-20261004-028",
                               "LRN-20261004-029")})
REMAP = {
    "## LRN-20261004-027 — 提速的正解是「把模型做小且做对」，不是「调大迭代上限」":
        "## LRN-20261004-054 — 提速的正解是「把模型做小且做对」，不是「调大迭代上限」",
    "## LRN-20261004-028 — 收敛失败时**必须读迭代日志的 `convergence norms`**":
        "## LRN-20261004-055 — 收敛失败时**必须读迭代日志的 `convergence norms`**",
    "## LRN-20261004-029 — 复现文献数字前，先确认**指标定义**与**运动学方向**":
        "## LRN-20261004-056 — 复现文献数字前，先确认**指标定义**与**运动学方向**",
}
for old, new in REMAP.items():
    c = t.count(old)
    print(f"  {'✓' if c == 1 else '⚠'} ({c}) {old[:56]}")
    if c == 1:
        t = t.replace(old, new)
io.open(LEARN, "w", encoding="utf-8").write(t)

# 复查
ids = re.findall(r"^#+\s*\[?([A-Z]{3}-\d{8}-\d{3})", t, re.M)
dup = sorted({i for i in ids if ids.count(i) > 1})
print(f"\n复查：条目 {len(ids)} 条，真实重号 {dup or '无 ✅'}")
print("改号后出现次数：",
      {k: t.count(k) for k in ("LRN-20261004-027", "LRN-20261004-054",
                               "LRN-20261004-028", "LRN-20261004-055",
                               "LRN-20261004-029", "LRN-20261004-056")})
print(f"文件 {len(t.encode('utf-8'))} B")
