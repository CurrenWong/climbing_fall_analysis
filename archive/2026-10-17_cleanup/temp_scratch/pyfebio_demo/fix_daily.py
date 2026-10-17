"""修 daily 索引：纠正 §十四 的过时结论 + 补 §十五。"""
from pathlib import Path

D = Path(r"C:\Users\curre\.copaw\workspaces\default\memory\2026-10-04.md")
t = D.read_text(encoding="utf-8")

OLD14 = ("- `2026-10-04/ankle-fe-solve-and-handover.md` → §十四 tie 发散结案（误诊）："
         "真根因=**关节两软骨面不共形 ⇒ 关节未闭合 ⇒ F 全为罚因子伪影**"
         "（PENALTY 1e-4→1 ⇒ 0.83→419.56 N，AUGLAG 也不收敛）⇒ `[LRN-054]`，`436.8N` 失效")
NEW14 = ("- `2026-10-04/ankle-fe-solve-and-handover.md` → §十四 tie 发散结案（误诊）："
         "tie 无问题；`|u|max` 10× = 400mm 胫骨悬臂弯曲（刚体拟合 0.17° 占 79%）"
         "⇒ `[LRN-054]` ⚠️ **该条的后两条结论已被 §十五/`[LRN-056]` 撤回**")

NEW15 = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十五 **自我更正**（`[LRN-20261005-056]`）："
         "撤回「关节不共形」（节点稀疏漏检；稠密面积加权 t=0 接触 31mm²、加载后 **344mm²**≈生理 300mm²；"
         "点在-tet 互穿仅 6 节点）与「F 是罚因子伪影」（`k=F/δ` 随罚趋近物理值 **49/151/451/880** N/mm "
         "vs Anderson 估计 **1059**）；**真堵点=`PENALTY≥2` 全负 Jacobian ⇒ F 无法验收敛，真值可能高 1.5~2×**"
         "⇒ `436.8N` 偏低")

if OLD14 in t:
    t = t.replace(OLD14, NEW14, 1)
    print("✓ §十四 索引已纠正")
else:
    print("⚠ 未匹配 §十四 索引")
if "§十五 **自我更正**" not in t:
    # 插在 §十六 之前
    k = t.find("→ §十六 踝子线下一步")
    if k > 0:
        j = t.rfind("\n", 0, k)
        t = t[:j] + NEW15 + t[j:]
        print("✓ 已补 §十五 索引")
    else:
        t = t.rstrip() + NEW15 + "\n"
        print("✓ 已补 §十五 索引（尾部）")
D.write_text(t, encoding="utf-8")
print(f"daily note -> {D.stat().st_size:,} B")
tt = D.read_text(encoding="utf-8")
for k in ("§十四 tie 发散结案", "§十五 **自我更正**", "§十六 踝子线下一步"):
    print(f"  [{'OK' if k in tt else 'XX'}] {k}")
