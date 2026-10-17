"""固化 tied-elastic 这一轮的结论：写 LRN + 更新 ANKLE_FIX_TODO / CONV_TODO / MEMORY。"""
import io
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
LRN = ROOT / ".learnings/LEARNINGS.md"

RECS = [
    ("LRN-20261004-046", "knowledge_gap",
     "FEBio `<Surface>` 必须列出**单元真实 facet**，否则静默丢弃",
     "给 tied-elastic 写 primary/secondary 面时，我按 npz 里的四边形照抄，"
     "结果 FEBio 报 `tiesurf_* has 61 invalid facets`（93 面里 61 个无效）并"
     "**静默丢弃**，tie 形同虚设。",
     "93 面里 61 个无效 = 正好是 hex8 来源的 quad 数。软骨在 build_feb 里已全部"
     "拆成 **tet4**，它的 facet 只能是**三角形**；四边形根本不是网格的面。",
     "面必须从**已装配的 elem_blocks** 反查真实 facet（`cart_outer()` 用的是同一招），"
     "不能信任上游文件的四边形。修正后 invalid 从 61 → 1。"
     "**这与 `cart_outer` 当年的坑同源，说明它是个反复出现的一类错误，不是一次性事故。**"),

    ("LRN-20261004-047", "best_practice",
     "mortar/tie 的 secondary 面**不完整**会伪装成「模型变硬 2.2 倍」",
     "我用「与软骨底面坐标重合」挑 secondary 面，tibia 匹配 63/93、talus 166/200，"
     "缺口让部分 primary 节点无处投影，被 `search_radius` 拉到**远处** facet 上锁死，"
     "支反力从 436.8 N 虚假升到 **709 → 826 → 969 N**（正好 ~2.2×）。",
     "面缺口 ⇒ 局部软骨底面被绑到远处骨面 ⇒ 虚假刚度。",
     "改用「骨外表面中靠近软骨底面的**完整局部片**」（117/118 面，饱和后不随 TOL 变），"
     "支反力立刻回到 **404.7 N**（−7.4%），且罚参数单调收敛。"
     "**排查判据：支反力随 penalty 单调上升背离真值 = 面有洞；"
     "从下方收敛 = mortar 正常。**"),

    ("LRN-20261004-048", "best_practice",
     "验证 mortar 必须做「同加载状态」对比，且必须看 `time` 不看步数",
     "FEBio 的 time_stepper 会自适应放大步长：`TIME_STEPS=20` 最终只报"
     "「14 步完成」，但 xplt 里 `time=1.0` —— **确实是满载荷**。"
     "只看「14 步 vs 13 步」会以为两边状态不同而误判不可比。",
     "步数 ≠ 加载比例。",
     "对比一律读 xplt 的 `state attrs['time']`，取 **t=1.0** 那一帧；"
     "罚参数扫描时也只看 t=1.0，避免拿中间态互比。"),

    ("LRN-20261004-049", "best_practice",
     "参数扫描必须把「建网格」和「跑求解」解耦",
     "每轮重建 .feb（含 Jacobian 自检 26,788 单元 + 写 1.9 MB）约 15-20 s，"
     "而扫参数根本不动网格 ⇒ 迭代极慢（单轮 20-350 s）。",
     "把求解器参数写进 .feb 后再跑，等于每轮白付一次建网格成本。",
     "`temp/pyfebio_demo/fasttie.py`：只对已写好的 .feb 做**文本替换**再调 FEBio，"
     "并把 `max_retries=0` 让它一失败就退。单轮 2-40 s（快 5-10×）。"
     "**并被这一招直接救回：** 早前 12 例扫描 0/12 通过，快扫后立刻暴露"
     "「penalty 1/2/4 都能过」—— 说明失败是 secondary 面造成的，不是罚参数。"),

    ("LRN-20261004-050", "negative_result",
     "软骨面内红色加密被母单元扭曲挡住（30% 子单元翻负）",
     "想通过「软骨 1→4 面内加密 + tied-elastic 绑回粗骨」提升接触分辨率。"
     "加密本身正确（wedge6 加密 0 翻负、L=0 体积偏差 0.0000%），但 hex8 加密恒有"
     "约 30% 子单元翻负（69/244、97/300、296/976），且 4 个子位置分布均匀"
     "（20/19/16/14）⇒ **不是我的绕向 bug**。",
     "母单元本身严重扭曲：同一面上 offset 方向从「主要 −x」摆到「主要 −z」"
     "（60-90° 变化），因为 `CART_OFFSET=toward`（指向对侧骨）在关节边缘处方向不良定义。"
     "母单元 |V|min 只有 0.0127 / 0.00066。",
     "**先修几何再加密。** 已实测两条路："
     "① `CART_OFFSET=normal` 单独用 → |V|min 0.0063（更差）；"
     "② `CART_OFFSET=normal + ART_ALIGN=0.5` → |V|min **0.429（34×）**，"
     "但关节面节点 73→33（砍掉 55% 的边缘脊线节点）。"
     "③ 下一步：在**曲面参数域**内做重划分/光顺（而不是在 3D 偏移空间加密），"
     "或对边缘节点单独用「切平面法向」而不是顶点法向。"),
]

s = LRN.read_text(encoding="utf-8") if LRN.exists() else "# LEARNINGS\n"
add = ["\n\n---\n"]
for rid, cat, title, what, why, how in RECS:
    if rid in s:
        add.append(f"（{rid} 已存在，跳过）\n")
        continue
    add.append(f"\n## {rid} [{cat}] {title}\n\n"
               f"- **现象**：{what}\n- **根因**：{why}\n"
               f"- **修复/下一步**：{how}\n- **日期**：2026-10-04\n")
LRN.write_text(s + "".join(add), encoding="utf-8")
print(f"LEARNINGS.md -> {LRN.stat().st_size:,} B")
for rid, *_ in RECS:
    print(f"  [{'OK' if rid in LRN.read_text(encoding='utf-8') else 'XX'}] {rid}")
