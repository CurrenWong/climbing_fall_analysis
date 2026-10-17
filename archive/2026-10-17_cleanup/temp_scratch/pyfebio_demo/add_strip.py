"""加 LIG_MODE=strip：韧带做成**只用真实骨节点**的 quad 条带。

为什么这么做（两条实测教训合起来推出来的）：
  ① truss/BeamDomain（elastic-truss）在真实模型里让 FEBio **读文件就崩**
     （最小实验能过 —— 实验性特性，不值得赌）。
  ② shell 用 THUMS 原始网格能读，但引入了 **137 个「只属韧带」的节点**；
     这些节点的**钻孔转动自由度无约束** ⇒ 刚体模态 ⇒ 首个载荷步就翻 22 个单元。
  ③ ⇒ 把韧带做成**两侧附着节点之间的 quad 条带**：节点全是骨节点（零新增），
     既保留 shell 的行为，又没有自由转动自由度。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① lig_recs 也接受 strip
s = s.replace('if LIG_MODE in ("shell", "truss"):\n    n_lig = int(lg["n_lig"])',
              'if LIG_MODE in ("shell", "truss", "strip"):\n    n_lig = int(lg["n_lig"])')
assert '"strip"' in s

# ② 生成条带
old = """        if LIG_MODE == "shell":
            used.update({int(v) for q in quads for v in q})
        else:
            used.update({int(v) for t in trusses for v in t})"""
assert s.count(old) == 1, f"gen {s.count(old)}"
s = s.replace(old, """        if LIG_MODE == "shell":
            used.update({int(v) for q in quads for v in q})
        else:
            used.update({int(v) for t in trusses for v in t})
        # strip：把 truss 束两两连成 quad（只用已有节点，零新增）
        strips = []
        for k in range(len(trusses) - 1):
            a0, b0 = trusses[k]
            a1, b1 = trusses[k + 1]
            q = [a0, a1, b1, b0]
            if len(set(q)) == 4:
                strips.append(q)
        lig_recs[-1]["strips"] = strips if not lig_recs or True else []""")

# ③ lig_blocks 支持 strip
old2 = """        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["trusses"]]"""
assert s.count(old2) == 1, f"blk {s.count(old2)}"
s = s.replace(old2, """        elif LIG_MODE == "truss":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["trusses"]]
        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["strips"]]""")

# ④ 元素类型
s = s.replace('LIG_ETYPE = "quad4" if LIG_MODE == "shell" else "line2"',
              'LIG_ETYPE = "line2" if LIG_MODE == "truss" else "quad4"')

# ⑤ MeshDomains：strip 用 ShellDomain
old3 = '    if LIG_MODE == "shell":\n        for r in lig_recs:'
assert s.count(old3) == 1
s = s.replace(old3, '    if LIG_MODE in ("shell", "strip"):\n        for r in lig_recs:')

# ⑥ 写出 guard
s = s.replace('if LIG_MODE in ("shell", "truss"):\n    for nm, conns in lig_blocks.items():',
              'if LIG_MODE in ("shell", "truss", "strip"):\n'
              '    for nm, conns in lig_blocks.items():')

# ⑦ 壳退化检查覆盖 strip
s = s.replace('if LIG_MODE == "shell":\n    nbad_shell = 0',
              'if LIG_MODE in ("shell", "strip"):\n    nbad_shell = 0')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
