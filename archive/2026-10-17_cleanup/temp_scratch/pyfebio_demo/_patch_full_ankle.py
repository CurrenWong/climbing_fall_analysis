"""build_feb.py -> 支持多接触对 + 按 FIX_BONE 固定（行级手术，避开缩进/空白脆弱匹配）。"""
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\scripts\ankle_fe\build_feb.py")
lines = P.read_text(encoding="utf-8").split("\n")


def find(pred, start=0):
    for i in range(start, len(lines)):
        if pred(lines[i]):
            return i
    raise AssertionError("not found")


out = []
i = 0
n_fix = n_sp = n_ct = n_ns = 0
while i < len(lines):
    l = lines[i]

    # ① 固定面：talus_cort -> FIX_BONE
    if l.strip() == 'tal_faces = []':
        j = i
        while 'print(f"距骨外表面' not in lines[j]:
            j += 1
        out += [
            'fix_faces = []',
            'for et in ("hex8", "tet4", "penta6"):',
            '    if (FIX_BONE, et) in elem_blocks:',
            '        fix_faces += surface_faces(FIX_BONE, et)',
            'projT = XYZ[[f[0] - 1 for f, _ in fix_faces]] @ axis',
            'lo = projT.min()',
            'bot_sel = [f for f, p in zip([f for f, _ in fix_faces], projT) '
            'if p < lo + 3.0]',
            'print(f"{FIX_BONE} 外表面 {len(fix_faces):,} 面；远端固定面 '
            '{len(bot_sel)} 面")',
        ]
        n_fix += 1
        i = j + 1
        continue

    # ② Surface 定义里的 talus_bottom 名字
    if l.startswith('all_surfs = {"tibia_top": top_sel, "talus_bottom": bot_sel,'):
        out.append('all_surfs = {"tibia_top": top_sel, '
                   'f"{FIX_BONE}_bottom": bot_sel,')
        n_ns += 1
        i += 1
        continue

    # ③ 固定 BC 的 node_set
    if 'node_set="talus_bottom_ns"' in l:
        out.append(l.replace('talus_bottom_ns", bot_sel',
                             'f"{FIX_BONE}_bottom_ns", bot_sel')
                   .replace('node_set="talus_bottom_ns"',
                            'node_set=f"{FIX_BONE}_bottom_ns"')
                   .replace('fix_talus', 'fix_distal'))
        n_ns += 1
        i += 1
        continue

    # ④ SurfacePair -> 多对
    if l.strip() == 'if CARTS:' and i + 1 < len(lines) and \
            'SurfacePair name="tibiotalar_pair"' in lines[i + 1]:
        j = i
        while "w('    </SurfacePair>')" not in lines[j]:
            j += 1
        out += [
            'for _pn, _pa, _pb in PAIRS:',
            '    w(f\'    <SurfacePair name="{_pn}">\')',
            '    w(f\'      <primary>surf_{_pa}</primary>\')',
            '    w(f\'      <secondary>surf_{_pb}</secondary>\')',
            "    w('    </SurfacePair>')",
        ]
        n_sp += 1
        i = j + 1
        continue

    # ⑤ Contact -> 多对
    if l.strip() == 'if CARTS:' and i + 1 < len(lines) and \
            "<Contact>" in lines[i + 1]:
        j = i
        while "w('    </contact>')" not in lines[j]:
            j += 1
        body = lines[i + 2:j]
        body = [(b[4:] if b.startswith('    ') else b) for b in body]
        out.append('if PAIRS:')
        out.append("    w('  <Contact>')")
        out.append('    for _pn, _pa, _pb in PAIRS:')
        for b in body:
            b = b.replace('_pn[:-5]', '_pn[:-5]')
            b = b.replace('w(\'    <contact name="tibiotalar" '
                          'surface_pair="tibiotalar_pair" \'',
                          'w(f\'    <contact name="{_pn[:-5]}" '
                          'surface_pair="{_pn}" \'')
            out.append('    ' + b)
        out.append("        w('    </contact>')")
        n_ct += 1
        i = j + 1
        continue

    out.append(l)
    i += 1

P.write_text("\n".join(out), encoding="utf-8")
print(f"OK: fix={n_fix} surf={n_ns} pair={n_sp} contact={n_ct}")
