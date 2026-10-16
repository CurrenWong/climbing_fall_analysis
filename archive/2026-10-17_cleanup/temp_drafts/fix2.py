import io

R = r'D:\Project\climbing_fall_analysis\\'
FILES = ['docs\\paper\\论文v2_数学建模版.md', 'docs\\paper\\论文v2_出处索引.md',
         'docs\\methods\\方案更新_v2.md', 'docs\\reports\\最终结论_2026-10-10.md',
         'docs\\reports\\总体进度.md', 'docs\\reports\\线2绝对高度路线_封档建议_2026-10-16.md']

START = '可解的最软点（总刚度'
ENDMK = 'V1_CORRECTION_PER_NODE_K.md`）'
PLAIN = ('可解的最软点（总刚度 875 N/mm = 实测 6.4×）首折 3.05 m；'
         '实测锚点（总 100–175 N/mm）在两网格上均不可解 ⇒ 物理首折 < 3.05 m；'
         '目标 6–9 m ⇒ 差 ≥1.97–2.95×。'
         '⛔ 更正：原「最可辩护配置 ≈4.8–5.1 m／差 1.25–1.87×」作废'
         '（把 per-node k 当物理点；见 `results/opensim_fe/V1_CORRECTION_PER_NODE_K.md`）')

for f in FILES:
    p = R + f
    t = io.open(p, encoding='utf-8').read()
    n = t.count(START)
    fixed = 0
    # 把插入段整段替换为**纯文本版**（去掉嵌套 **，避免与外层粗体冲突）
    while START in t:
        i = t.index(START)
        j = t.index(ENDMK, i) + len(ENDMK)
        t = t[:i] + PLAIN + t[j:]
        fixed += 1
    extra = 0
    # 去掉重复的「实测跖面支撑…不可解；」
    DUP = ('实测跖面支撑 100–175 N/mm 在本模型内**不可解**；'
           '实测跖面支撑 100–175 N/mm 在本模型内**不可解**；')
    if DUP in t:
        t = t.replace(DUP, '实测跖面支撑 100–175 N/mm 在本模型内**不可解**；')
        extra += 1
    # 补上未命中的那处（§5 机制归因段）
    OLD = '；最可辩护配置 $\\approx 4.8$–$5.1$ m vs 目标 6–9 m）；详 §4.6'
    if OLD in t:
        t = t.replace(OLD, '；' + PLAIN + '］；详 §4.6')
        extra += 1
    io.open(p, 'w', encoding='utf-8', newline='\n').write(t)
    print(f'  {f}: 段替换 {fixed} 处，附加修 {extra} 处')

print('\n=== 残留检查 ===')
for f in FILES + ['results\\opensim_fe\\V1_PHYSICAL_ANCHOR_AND_FINAL_VERDICT.md',
                  'results\\opensim_fe\\V1_ABSOLUTE_VALUE_UNDERDETERMINED.md',
                  'results\\opensim_fe\\V1_CORRECTION_PER_NODE_K.md']:
    t = io.open(R + f, encoding='utf-8').read()
    bad = t.count('4.8$–$5.1') + t.count('4.8–5.1')
    ok = '作废' in t or '不可解' in t or '3.05' in t
    print(f'  {"⚠️" if bad else "✅"} {f}: 4.8-5.1 残留 {bad} 处；有更正={ok}')
