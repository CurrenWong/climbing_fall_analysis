"""加 DRIVER=dispx：**只**指定顶面轴向位移，横向自由度交给韧带+接触。

为什么这才是对标的正确配置：
  · Wú Kǎi 加的是 600 N **力**，胫骨可自由转动 ⇒ 关节在韧带+接触下自行沉降
  · DRIVER=disp（三自由度全夹持）把胫骨锁死 ⇒ 韧带只有轴向差动
    ⇒ 位移比吴恺小两个数量级、排序失真（实测 ΔL 仅 0.002~0.10 mm）
  · DRIVER=force 又因无横向约束而病态
  ⇒ 折中：轴向指定位移（可控、不病态）+ 横向自由（运动学交给韧带/接触）
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'if DRIVER == "disp":'
assert s.count(old) == 1, s.count(old)
s = s.replace(old, 'if DRIVER in ("disp", "dispx"):')

old2 = """        for _di, _dn in enumerate("xyz"):
            _v = (-DISP_MM if axis[_di] > 0 else DISP_MM) if _di == _kax else 0.0
            w(f'    <bc name="{_nm}_top_{_dn}" type="prescribed displacement" '
              f'node_set="{_ns}">')
            w(f'      <dof>{_dn}</dof><value lc="1">{_v:.6f}</value>'
              '<relative>0</relative>')
            w('    </bc>')"""
assert s.count(old2) == 1
s = s.replace(old2, """        _dofs = "xyz" if DRIVER == "disp" else "xyz"[_kax]
        for _di, _dn in enumerate("xyz"):
            if DRIVER == "dispx" and _di != _kax:
                continue     # 横向自由：交给韧带 + 接触约束
            _v = (-DISP_MM if axis[_di] > 0 else DISP_MM) if _di == _kax else 0.0
            w(f'    <bc name="{_nm}_top_{_dn}" type="prescribed displacement" '
              f'node_set="{_ns}">')
            w(f'      <dof>{_dn}</dof><value lc="1">{_v:.6f}</value>'
              '<relative>0</relative>')
            w('    </bc>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
