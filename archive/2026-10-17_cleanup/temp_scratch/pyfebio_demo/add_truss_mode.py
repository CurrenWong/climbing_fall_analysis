"""给 ankle_lig_feb.py 加 LIG_MODE=truss：按吴恺用杆单元（line2 + BeamDomain
type="elastic-truss"），把每条韧带的截面积**折进材料 E**。

为什么折进 E：FEBio 的 elastic-truss 不暴露截面积参数（实测 <area> 材料标签不认），
设隐式面积为 A_imp，则 k_i = E_i·A_imp/L_i。取 E_i = EA_i · A_cross_i 可得
k_i ∝ EA_i·A_cross_i/L_i —— **相对刚度完全正确，且与 A_imp 无关**。
所以排序结论不受 A_imp 影响；绝对量级再据实说明。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

# ① LIG_MODE 文档 + 韧带材料/元素生成支持 truss
a1 = """lig_recs = []
if LIG_MODE == "shell":"""
assert s.count(a1) == 1
s = s.replace(a1, """lig_recs = []
if LIG_MODE in ("shell", "truss"):""")

a2 = """        lig_recs.append(dict(idx=i, name=nm, pid=pid, quads=quads,
                             mat=mname, thick=th, ea=ea))
        used.update({int(v) for q in quads for v in q})"""
assert s.count(a2) == 1
s = s.replace(a2, """        a_cross = float(lg["area_cross"][i])
        trusses = [[int(x) for x in t] for t in lg["truss"][i]
                   if int(t[0]) > 0 and int(t[1]) > 0]
        # 杆单元：截面积折进 E（E_eff = EA · A_cross），详见文件头注释
        MATERIALS[mname] = ("isotropic elastic",
                            {"density": 1.1e-9,
                             "E": (ea * a_cross) if LIG_MODE == "truss" else (
                                 260.0 if LIG_E_UNIFORM else ea),
                             "v": 0.40})
        lig_recs.append(dict(idx=i, name=nm, pid=pid, quads=quads,
                             mat=mname, thick=th, ea=ea,
                             a_cross=a_cross, trusses=trusses))
        if LIG_MODE == "shell":
            used.update({int(v) for q in quads for v in q})
        else:
            used.update({int(v) for t in trusses for v in t})""")

# ② lig_blocks：truss 用 (a,b) 对
a3 = """lig_blocks = {}
if LIG_MODE == "shell":
    for r in lig_recs:
        lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in q]
                                 for q in r["quads"]]"""
assert s.count(a3) == 1
s = s.replace(a3, """lig_blocks = {}
LIG_ETYPE = "quad4" if LIG_MODE == "shell" else "line2"
if LIG_MODE in ("shell", "truss"):
    for r in lig_recs:
        if LIG_MODE == "shell":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in q]
                                     for q in r["quads"]]
        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["trusses"]]""")

# ③ 韧带壳退化检查：只对 shell 模式
a4 = 'if LIG_MODE == "shell":\n    nbad_shell = 0'
assert s.count(a4) == 1
s = s.replace(a4, 'if LIG_MODE == "shell":\n    nbad_shell = 0')

# ④ 写出：按 LIG_ETYPE
a5 = """    for nm, quads in lig_blocks.items():
        w(f'    <Elements type="quad4" name="lig_{nm}">')
        for c in quads:
            eid += 1
            w(f'      <elem id="{eid}">{",".join(str(x) for x in c)}</elem>')
        w('    </Elements>')"""
assert s.count(a5) == 1
s = s.replace(a5, """    for nm, conns in lig_blocks.items():
        w(f'    <Elements type="{LIG_ETYPE}" name="lig_{nm}">')
        for c in conns:
            eid += 1
            w(f'      <elem id="{eid}">{",".join(str(x) for x in c)}</elem>')
        w('    </Elements>')""")

# ⑤ MeshDomains：truss 用 BeamDomain
a6 = """    for r in lig_recs:
        w(f'    <ShellDomain name="lig_{r["name"]}" '
          f'mat="{r["mat"]}">')
        w(f'      <shell_thickness>{r["thick"]}</shell_thickness>')
        w('    </ShellDomain>')"""
assert s.count(a6) == 1
s = s.replace(a6, """    if LIG_MODE == "shell":
        for r in lig_recs:
            w(f'    <ShellDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}">')
            w(f'      <shell_thickness>{r["thick"]}</shell_thickness>')
            w('    </ShellDomain>')
    else:
        # ★ 实测：BeamDomain type 只认 linear-truss / elastic-truss / linear-beam；
        #   elastic-truss + "isotropic elastic" 可跑通（linear-truss 不收标准材料，
        #   linear-beam 直接崩溃）。截面积折进 E，见文件头。
        for r in lig_recs:
            w(f'    <BeamDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}" type="elastic-truss"/>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
