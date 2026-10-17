"""压 MEMORY 最后 205 B：紧凑化第 30 行（386 B 最大行），细节已在 references。"""
import io

MEM = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
m = io.open(MEM, encoding="utf-8").read()
n0 = len(m.encode("utf-8"))
OLD = ('| 21 | **FE 四铁律**（接触/翻转/`.feb` 语法）+ **"某驱动不行先验产物再看报错"** | '
       '`references/fe-modeling-gotchas.md`。① "只靠接触⇒刚体模态"是**不完整诊断**（`[LRN-058]`：'
       '真因常是载荷被丢/初始未咬合）② 接触面须与真实网格面重合且合法 2-流形 ③ tet4 翻转须**奇置换** '
       '④ 优化前先确认它是瓶颈 |')
NEW = ('| 21 | **FE 四铁律** + **"某驱动不行先验产物再看报错"** | '
       '`references/fe-modeling-gotchas.md`（接触面须与真实网格面重合且 2-流形／tet4 翻转须奇置换／'
       '"只靠接触⇒刚体模态"是**不完整诊断**`[LRN-058]`／优化前先确认它是瓶颈）|')
if OLD in m:
    m = m.replace(OLD, NEW)
    print("· 第 30 行已紧凑化")
else:
    print("⚠ 第 30 行未命中")
io.open(MEM, "w", encoding="utf-8").write(m)
n1 = len(m.encode("utf-8"))
print(f"MEMORY.md {n0} -> {n1} B  {'✅ ≤10240' if n1 <= 10240 else f'❌ 超 {n1-10240}'}")
