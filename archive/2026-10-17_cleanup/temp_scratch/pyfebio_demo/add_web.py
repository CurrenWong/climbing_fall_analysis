"""加 LIG_MODE=web：韧带做成**三角网**（只用骨节点）。

每个 A 侧节点连到对侧最近 2 个 B 节点成一个三角；B 侧同理反向。
⇒ 元素数 ≈ 2(|A|+|B|)，连通、无新增节点、无自由转动自由度。
同时修：元素数为 0 的韧带**不写空 <Elements>**（空块会让 FEBio 报
"tag Elements: missing attribute id"）。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① 生成 web
old = '''        # strip：把 truss 束两两连成 quad（只用已有节点，零新增）
        strips = []
        for k in range(len(trusses) - 1):
            a0, b0 = trusses[k]
            a1, b1 = trusses[k + 1]
            q = [a0, a1, b1, b0]
            if len(set(q)) == 4:
                strips.append(q)
        lig_recs[-1]["strips"] = strips if not lig_recs or True else []'''
assert s.count(old) == 1, f"gen {s.count(old)}"
s = s.replace(old, '''        # strip：把 truss 束两两连成 quad（只用已有节点，零新增）
        strips = []
        for k in range(len(trusses) - 1):
            a0, b0 = trusses[k]
            a1, b1 = trusses[k + 1]
            q = [a0, a1, b1, b0]
            if len(set(q)) == 4:
                strips.append(q)
        # web：三角网 —— 每个 A 节点连到对侧最近 2 个 B 节点，反之亦然
        Aset = sorted({t[0] for t in trusses})
        Bset = sorted({t[1] for t in trusses})
        tris = []
        if Aset and Bset:
            PA = xyz_all[[idx_all[a] for a in Aset]]
            PB = xyz_all[[idx_all[b] for b in Bset]]
            D = np.linalg.norm(PA[:, None, :] - PB[None, :, :], axis=2)
            for i, a in enumerate(Aset):
                for j in np.argsort(D[i])[:2]:
                    t3 = [a, Bset[int(j)],
                          Bset[int(np.argsort(D[i])[1])]]
                    if len(set(t3)) == 3 and tuple(sorted(t3)) not in \\
                            {tuple(sorted(x)) for x in tris}:
                        tris.append(t3)
            Dt = D.T
            for j, b in enumerate(Bset):
                for i in np.argsort(Dt[j])[:2]:
                    t3 = [b, Aset[int(i)],
                          Aset[int(np.argsort(Dt[j])[1])]]
                    if len(set(t3)) == 3 and tuple(sorted(t3)) not in \\
                            {tuple(sorted(x)) for x in tris}:
                        tris.append(t3)
        lig_recs[-1]["strips"] = strips
        lig_recs[-1]["tris"] = tris''')

# ② lig_blocks 支持 web
old2 = '''        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["strips"]]'''
assert s.count(old2) == 1, f"blk {s.count(old2)}"
s = s.replace(old2, '''        elif LIG_MODE == "strip":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["strips"]]
        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["tris"]]''')

# ③ 元素类型：web 用 tri3
s = s.replace('LIG_ETYPE = "line2" if LIG_MODE == "truss" else "quad4"',
              'LIG_ETYPE = {"truss": "line2", "web": "tri3"}.get(LIG_MODE, "quad4")')

# ④ 写出时跳过空块
old3 = '''    for nm, conns in lig_blocks.items():
        w(f'    <Elements type="{LIG_ETYPE}" name="lig_{nm}">')'''
assert s.count(old3) == 1
s = s.replace(old3, '''    for nm, conns in lig_blocks.items():
        if not conns:                      # ★ 空 <Elements> 会让 FEBio 报
            continue                        #   "missing attribute id"
        w(f'    <Elements type="{LIG_ETYPE}" name="lig_{nm}">')''')

# ⑤ MeshDomains 跳过空块
old4 = '''    if LIG_MODE in ("shell", "strip"):
        for r in lig_recs:
            w(f'    <ShellDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}">')'''
assert s.count(old4) == 1
s = s.replace(old4, '''    if LIG_MODE in ("shell", "strip", "web"):
        for r in lig_recs:
            if not lig_blocks.get(r["name"]):
                continue
            w(f'    <ShellDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}">')''')

# ⑥ 各处 guard 加 web
s = s.replace('if LIG_MODE in ("shell", "truss", "strip")',
              'if LIG_MODE in ("shell", "truss", "strip", "web")')
s = s.replace('if LIG_MODE in ("shell", "strip")', 'if LIG_MODE in ("shell", "strip", "web")')

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
