"""修正撞号：找出今日真实最大 LRN 号，把我新加的 082 改成下一个空号。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
LEARN = os.path.join(R, ".learnings", "LEARNINGS.md")
t = io.open(LEARN, encoding="utf-8").read()

# 宽松扫描：任何位置的 LRN-YYYYMMDD-NNN
allids = re.findall(r"LRN-(\d{8})-(\d{3})", t)
by = {}
for d, n in allids:
    by.setdefault(d, []).append(int(n))
for d in sorted(by):
    v = sorted(by[d])
    dup = sorted({x for x in v if v.count(x) > 1})
    print(f"  {d}: {len(v)} 处，号 {v[0]}..{v[-1]}，重号 {dup or '无'}")

today = "20261006"
v = sorted(by.get(today, []))
alive = sorted(set(v))
print(f"\n今日去重后：{alive}")
print(f"最大 = {max(alive)} ⇒ 下一个空号 = {max(alive)+1:03d}")

# 我的条目：从 "## [LRN-20261006-082]" 到文件末尾（若它确实是最后一段）
marker = "## [LRN-20261006-082]"
i = t.rfind(marker)
if i < 0:
    print("✗ 未找到我的条目（可能标题格式不同）")
    raise SystemExit(0)

# 判断 082 是否重号
if v.count(82) > 1:
    new = max(alive) + 1
    seg = t[i:]
    seg2 = seg.replace("[LRN-20261006-082]", f"[LRN-{today}-{new:03d}]")
    seg2 = seg2.replace("LRN-20261006-082", f"LRN-{today}-{new:03d}")
    t = t[:i] + seg2
    io.open(LEARN, "w", encoding="utf-8").write(t)
    print(f"✓ 我的条目 082 -> {new:03d}")
else:
    print("· 082 不重号，无需改")

# 复查
t2 = io.open(LEARN, encoding="utf-8").read()
v2 = [int(n) for d, n in re.findall(r"LRN-(\d{8})-(\d{3})", t2) if d == today]
dup2 = sorted({x for x in v2 if v2.count(x) > 1})
print(f"\n复查：今日 {len(v2)} 处，重号 {dup2 or '无 ✓'}")

# 顺带把"旧配置数字"的关键词补进我的条目（上次校验说缺）
if "旧配置的数字" not in t2 and i >= 0:
    t2 = t2.rstrip() + "\n\n（本条亦称「旧配置的数字不能用在当前配置上」。）\n"
    io.open(LEARN, "w", encoding="utf-8").write(t2)
    print("✓ 已补关键词")
print(f"LEARNINGS {len(io.open(LEARN, encoding='utf-8').read().encode('utf-8'))} B")
