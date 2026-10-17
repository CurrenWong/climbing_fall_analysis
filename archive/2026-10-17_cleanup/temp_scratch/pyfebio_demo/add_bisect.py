"""给 ankle_lig_feb.py 加 BONES 环境变量（逗号分隔），用于二分定位负 Jacobian。"""
from pathlib import Path

p = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\ankle_lig_feb.py")
s = p.read_text(encoding="utf-8")

old = 'FORCE_N = float(os.environ.get("FORCE_N", "600"))'
new = ('FORCE_N = float(os.environ.get("FORCE_N", "600"))\n'
       '# 二分定位用：逗号分隔的域键子集，空 = 全部\n'
       'ONLY = [x for x in os.environ.get("BONES", "").split(",") if x]')
assert s.count(old) == 1
s = s.replace(old, new)

old2 = 'for dom, (key, mat) in BONE_PARTS.items():'
new2 = ('for dom, (key, mat) in BONE_PARTS.items():\n'
        '    if ONLY and dom not in ONLY:\n'
        '        continue')
assert s.count(old2) == 1
s = s.replace(old2, new2)

old3 = 'print(f"[mode] MATSET={MATSET} LIG_MODE={LIG_MODE} TIME_STEPS={TIME_STEPS} "\n      f"FIX_TALUS={FIX_TALUS} FORCE={FORCE_N}N")'
assert s.count(old3) == 1
s = s.replace(old3, old3 + '\nprint(f"[mode] BONES={ONLY or \'all\'}")')

# 输出文件名带子集标识，避免互相覆盖
old4 = 'feb = WORK / f"ankle_lig_{MATSET}.feb"'
new4 = ('_tag = MATSET + ("_" + "_".join(ONLY) if ONLY else "")\n'
        'feb = WORK / f"ankle_lig_{_tag}.feb"')
assert s.count(old4) == 1
s = s.replace(old4, new4)
old5 = '(WORK / f"febio_{MATSET}.log")'
new5 = '(WORK / f"febio_{_tag}.log")'
assert s.count(old5) == 1
s = s.replace(old5, new5)

# 子集模式下，若缺 tibia_cort/talus_cort，跳过需要它们的几何计算
old6 = 'tib_pts = XYZ[[i - 1 for i in sorted(dom_nodes["tibia_cort"])]]'
new6 = ('if "tibia_cort" not in dom_nodes or "talus_cort" not in dom_nodes:\n'
        '    print("子集模式：跳过几何/载荷，仅测网格初始化")\n'
        '    tib_pts = XYZ[[i - 1 for i in sorted(dom_nodes[sorted(dom_nodes)[0]])]]\n'
        '    tal_pts = tib_pts\n'
        'else:\n'
        '    tib_pts = XYZ[[i - 1 for i in sorted(dom_nodes["tibia_cort"])]]')
assert s.count(old6) == 1
s = s.replace(old6, new6)
# 修正缩进：tal_pts 原行
old7 = 'tal_pts = XYZ[[i - 1 for i in sorted(dom_nodes["talus_cort"])]]\naxis = tib_pts.mean(0) - tal_pts.mean(0)'
new7 = ('    tal_pts = XYZ[[i - 1 for i in sorted(dom_nodes["talus_cort"])]]\n'
        'axis = tib_pts.mean(0) - tal_pts.mean(0)')
assert s.count(old7) == 1
s = s.replace(old7, new7)

p.write_text(s, encoding="utf-8")
import ast
ast.parse(s)
print("patched + syntax OK")
