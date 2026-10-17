"""给 traction 和落位位移加 LOAD_SIGN 开关（默认保持 -1 不改行为）。"""
import io

P = r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py"
s = io.open(P, encoding="utf-8").read()
if "LOAD_SIGN" in s:
    print("已存在 LOAD_SIGN")
else:
    # 1) 环境变量声明（插在 MULTISTEP 声明前）
    anchor = 'MULTISTEP = os.environ.get("MULTISTEP", "0") == "1"'
    decl = ('# ★ LOAD_SIGN：载荷/落位位移相对加载轴 LAX 的符号。\n'
            '#   LAX = 胫骨->距骨（远端）方向；ALIGN_AXIS=1 后 = +z。\n'
            '#   -1 = 沿 -LAX（沿轴负向）；+1 = 沿 +LAX（沿轴向远端）。\n'
            'LOAD_SIGN = float(os.environ.get("LOAD_SIGN", "-1"))\n'
            + anchor)
    assert s.count(anchor) == 1, s.count(anchor)
    s = s.replace(anchor, decl)
    # 2) 单步 traction
    a1 = "w(f'      <traction>{-LAX[0]*tr:.6f},{-LAX[1]*tr:.6f},{-LAX[2]*tr:.6f}'"
    n1 = ('w(f\'      <traction>{LOOKSIGNLAX[0]*tr:.6f},'
          '{LOOKSIGNLAX[1]*tr:.6f},{LOOKSIGNLAX[2]*tr:.6f}\'')
    assert s.count(a1) == 1, ("single", s.count(a1))
    s = s.replace(a1, n1)
    # 3) 多步 step2 traction
    a2 = ("w(f'        <traction>{-LAX[0]*_tr:.6f},{-LAX[1]*_tr:.6f},'\n"
          "              f'{-LAX[2]*_tr:.6f}</traction>')")
    n2 = ("w(f'        <traction>{LOOKSIGNLAX[0]*_tr:.6f},'\n"
          "              f'{LOOKSIGNLAX[1]*_tr:.6f},'\n"
          "              f'{LOOKSIGNLAX[2]*_tr:.6f}</traction>')")
    assert s.count(a2) == 1, ("multi", s.count(a2))
    s = s.replace(a2, n2)
    # 4) 定义 LOOKSIGNLAX：紧跟 axis 计算之后（LAX 定义处之后）
    a3 = "print(f\"加载轴（胫骨->距骨）= {np.round(axis,4)}\")"
    assert s.count(a3) == 1
    s = s.replace(a3, a3 + "\nLOOKSIGNLAX = LOAD_SIGN * LAX"
                  "\nprint(f\"[load] LOAD_SIGN={LOAD_SIGN} ⇒ 载荷方向 = "
                  "{np.round(LOOKSIGNLAX,4)}\")")
    io.open(P, "w", encoding="utf-8").write(s)
    print("已插入 LOAD_SIGN")
import ast
ast.parse(io.open(P, encoding="utf-8").read())
print("AST OK")
for k in ("LOAD_SIGN", "LOOKSIGNLAX"):
    print(f"  {k} 出现 {io.open(P, encoding='utf-8').read().count(k)} 次")
