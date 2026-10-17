"""按行重写 MEMORY 的第 21/23 行（唯一定位），压到 ≤10240。"""
import io

MEM = r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md"
ls = io.open(MEM, encoding="utf-8").read().splitlines()
n0 = sum(len(l.encode("utf-8")) for l in ls) + len(ls)  # 粗略
for i, l in enumerate(ls):
    if l.startswith("| 21 |"):
        ls[i] = ('| 21 | **FE 四铁律** + **"先验产物再看报错"** | `references/fe-modeling-gotchas.md`'
                 '（接触面须=真实网格面且 2-流形／tet4 翻转须奇置换／"只靠接触⇒刚体模态"是'
                 '**不完整诊断**／优化前先确认它是瓶颈）|')
    elif l.startswith("| 23 |"):
        ls[i] = ('| 23 | 🧭 **方向/符号结论须有「不可被约定翻转」的证据**（结构/几何/物理/解剖）'
                 '+ 判据配**已知正对照**；**更正要替换不要叠加** | '
                 '`references/operational-gotchas.md` |')
t = "\n".join(ls) + "\n"
io.open(MEM, "w", encoding="utf-8").write(t)
n = len(t.encode("utf-8"))
print(f"MEMORY.md -> {n} B  {'✅ ≤10240' if n <= 10240 else f'❌ 超 {n-10240}'}")
