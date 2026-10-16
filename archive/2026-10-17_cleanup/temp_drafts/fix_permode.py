import io
import json

R = r'D:\Project\climbing_fall_analysis'

NEW_CLAUSE = ('**可解的最软点（总刚度 **875 N/mm** = 实测 **6.4×**）首折 **3.05 m**；'
              '实测锚点（总 **100–175 N/mm**）在两网格上**均不可解** ⇒ **物理首折 < 3.05 m**；'
              '目标 6–9 m ⇒ **差 ≥ 1.97–2.95×**。**⛔ 更正**：原「最可辩护配置 ≈4.8–5.1 m／'
              '差 1.25–1.87×」**作废**（把 per-node k 当物理点；见 '
              '`results/opensim_fe/V1_CORRECTION_PER_NODE_K.md`）**')

log = []


def do(path, variants, new, tag):
    p = R + '\\' + path
    t = io.open(p, encoding='utf-8').read()
    for v in variants:
        if v in t:
            io.open(p, 'w', encoding='utf-8', newline='\n').write(t.replace(v, new, 1))
            log.append(f'  ✅ {tag}  ← {v[:44]}…')
            return True
    log.append(f'  ❌ {tag}  未命中（试了 {len(variants)} 个变体）')
    return False


# ---- 1) 论文 v2（5 处）----
paper = 'docs/paper/论文v2_数学建模版.md'
do(paper,
   ['最可辩护配置 ≈4.8–5.1 m vs 目标 6–9 m）'],
   '实测跖面支撑 100–175 N/mm 在本模型内**不可解**；' + NEW_CLAUSE + '）', '论文·编辑注')
do(paper,
   ['**最可辩护配置（875 N/mm + $\\delta\\approx 0.1$–0.15）首折 $\\approx 4.8$–$5.1$ m '
    'vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '论文·表5 V1 行')
do(paper,
   ['最可辩护配置 $\\approx 4.8$–$5.1$ m vs 目标 6–9 m；详 §4.6'],
   NEW_CLAUSE + '；详 §4.6', '论文·§5 机制归因段')
do(paper,
   ['**最可辩护配置 $\\approx 4.8$–$5.1$ m**（目标 6–9 m，差 1.25–1.87×）'],
   NEW_CLAUSE, '论文·关键局限(6)')
do(paper,
   ['**最可辩护配置首折 $\\approx 4.8$–$5.1$ m vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '论文·§6 Q3')

# ---- 2) 出处索引 ----
do('docs/paper/论文v2_出处索引.md',
   ['**最可辩护配置 $\\approx 4.8$–$5.1$ m vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '出处索引·V1 行')

# ---- 3) 结论件 C-V1ABS ----
do('docs/reports/最终结论_2026-10-10.md',
   ['**最可辩护配置 ≈4.8–5.1 m vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '结论件·C-V1ABS §二')
do('docs/reports/最终结论_2026-10-10.md',
   ['**最可辩护配置 ≈4.8–5.1 m（目标 6–9 m，差 1.25–1.87×）**'],
   NEW_CLAUSE, '结论件·C-V1ABS（变体）')

# ---- 4) 方案更新_v2 ----
do('docs/methods/方案更新_v2.md',
   ['**最可辩护配置 ≈4.8–5.1 m vs 目标 6–9 m**（差 1.25–1.87×）'],
   NEW_CLAUSE, '方案更新_v2·V1 行')

# ---- 5) 总体进度 §六 第 2 条 ----
do('docs/reports/总体进度.md',
   ['最可辩护 ≈4.8–5.1 m vs 目标 6–9 m（差 1.25–1.87×）'],
   NEW_CLAUSE, '总体进度·§六#2')

# ---- 6) 封档建议件 ----
do('docs/reports/线2绝对高度路线_封档建议_2026-10-16.md',
   ['**最可辩护配置 ≈4.8–5.1 m vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '封档建议件·判决表')
do('docs/reports/线2绝对高度路线_封档建议_2026-10-16.md',
   ['**最可辩护配置 ≈4.8–5.1 m vs 目标 6–9 m（差 1.25–1.87×）**'],
   NEW_CLAUSE, '封档建议件·判决表(2)')

# ---- 7) V1_PHYSICAL_ANCHOR_AND_FINAL_VERDICT.md（顶部加 banner + §2）----
p = R + r'\results\opensim_fe\V1_PHYSICAL_ANCHOR_AND_FINAL_VERDICT.md'
t = io.open(p, encoding='utf-8').read()
BAN = ('\n> ⛔ **本件 §2 的「最物理可辩护配置 ≈4.8–5.1 m」已作废（2026-10-16 自查更正）**：\n'
       '> 该数在 `v1_phys_corner.py` 里把 cl=3.0 的 **per-node `k=5`** 当作物理点，'
       '而 cl=3.0 有 **479** 个跖面节点 ⇒ 总刚度 **2395 N/mm = 实测的 17.4×**，**不是物理点**。\n'
       '> **现行认值**：实测锚点（总 100–175 N/mm）**在两网格上均不可解**（负 Jacobian 随软度\n'
       '> 由 5 → 3646 → 6349，机械性欠约束）；**可解最软点（粗网格，总 875 N/mm = 实测 6.4×）\n'
       '> ⇒ 首折 3.05 m**，且首折随支撑刚度**单调递增** ⇒ **物理首折 < 3.05 m** ⇒\n'
       '> **目标 6–9 m 差 ≥1.97–2.95×**。详见 `V1_CORRECTION_PER_NODE_K.md`。\n'
       '> （§1 的粗网格扫描、§3 的两旋钮耦合、§4.3 的封档结论**均不受影响**。）\n')
t = t.replace('---\n\n## §1', BAN + '\n---\n\n## §1', 1)
io.open(p, 'w', encoding='utf-8', newline='\n').write(t)
log.append('  ✅ PHYSICAL_ANCHOR banner')

# ---- 8) PHYSICAL_CORNER.json 标作废 ----
p = R + r'\results\opensim_fe\V1_PHYSICAL_CORNER.json'
d = json.loads(io.open(p, encoding='utf-8').read())
d = {"_SUPERSEDED": ("⛔ 作废（2026-10-16 自查）：本件把 cl=3.0 的 **per-node k=5** 当作「物理点」，"
                     "但 cl=3.0 有 479 个跖面节点 ⇒ 总刚度 2395 N/mm = 实测的 17.4×。"
                     "现行认值见 V1_CORRECTION_PER_NODE_K.md：实测锚点在两网格均不可解；"
                     "可解最软点（粗网格总 875 N/mm）首折 3.05 m ⇒ 物理首折 < 3.05 m。"),
     "_see": "results/opensim_fe/V1_CORRECTION_PER_NODE_K.md", **d}
io.open(p, 'w', encoding='utf-8', newline='\n').write(
    json.dumps(d, ensure_ascii=False, indent=2))
log.append('  ✅ PHYSICAL_CORNER.json 标作废')

print('\n'.join(log))
