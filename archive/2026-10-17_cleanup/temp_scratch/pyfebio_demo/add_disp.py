"""加 DRIVER 开关：force（吴恺原式 600 N）｜disp（指定位移，已验证管线的高效路径）。

为什么加 disp：实测 force 控制下每步要 ~97 次牛顿迭代，
3 步花 407 s —— 病态。已验证的 build_feb.py 正是用「顶面三自由度全指定位移」
（夹持端）才 17 s 收敛的（LRN-017）。排序结论对两者都成立。
"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'FORCE_N = float(os.environ.get("FORCE_N", "600"))'
assert s.count(old) == 1
s = s.replace(old, old + '''
DRIVER = os.environ.get("DRIVER", "force").strip() or "force"
DISP_MM = float(os.environ.get("DISP_MM", "0.05"))''')

old2 = 'print(f"[mode] FIX_NAVICULAR={FIX_NAVICULAR}")'
assert s.count(old2) == 1
s = s.replace(old2, old2 + '\nprint(f"[mode] DRIVER={DRIVER} DISP_MM={DISP_MM}")')

# 载荷段：disp 模式改为指定位移
old3 = """w('  <Loads>')"""
assert s.count(old3) == 1
s = s.replace(old3, """w('  <Loads>')
_use_force = (DRIVER == "force")""")

old4 = """for surf, F, _faces in (("tibia_top", F_t, top_tibia),
                      ("fibula_top", F_f, top_fibula)):
    if not _faces or F <= 0:
        continue"""
assert s.count(old4) == 1
s = s.replace(old4, """for surf, F, _faces in (("tibia_top", F_t, top_tibia),
                      ("fibula_top", F_f, top_fibula)):
    if not _use_force or not _faces or F <= 0:
        continue""")

# disp：顶面三自由度指定位移（夹持端）—— 比只约束横向稳得多
old5 = "w('  </Boundary>')"
assert s.count(old5) == 1
s = s.replace(old5, """if DRIVER == "disp":
    # ★ 顶面三个自由度**全部**指定位移（= 夹持端）：刚体模态消失，
    #   不依赖接触/韧带稳定（LRN-017 的核心修法）。
    _kax = int(np.argmax(np.abs(axis)))
    for _ns, _nm, _faces in (("tibia_top_ns", "tibia", top_tibia),
                             ("fibula_top_ns", "fibula", top_fibula)):
        if not _faces:
            continue
        for _di, _dn in enumerate("xyz"):
            _v = (-DISP_MM if axis[_di] > 0 else DISP_MM) if _di == _kax else 0.0
            w(f'    <bc name="{_nm}_top_{_dn}" type="prescribed displacement" '
              f'node_set="{_ns}">')
            w(f'      <dof>{_dn}</dof><value lc="1">{_v:.6f}</value>'
              '<relative>0</relative>')
            w('    </bc>')
w('  </Boundary>')""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
