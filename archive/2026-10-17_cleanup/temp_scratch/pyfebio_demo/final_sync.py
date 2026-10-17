"""追加当日会话记忆 §十一（tied-elastic），并做最终只读校验。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
MEMF = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"

SEC = """

---

## 十一、tied-elastic 通路：机制打通并验证（2026-10-04 下午）

### 起因
P1d2 的死结是「加密骨 ⇒ 骨 hex8 负 Jacobian」。手册 §3.14.5/§3.14.6 指出
`tied` / `tied-elastic` 专门用来**连接两个非共形网格**，于是设想去掉骨的加密、
只加密软骨，再用 tie 把细软骨底面绑回粗骨。

### 做了什么
- 新增 `scripts/ankle_fe/cart_refine.py`：把软骨底面节点**复制**成软骨专属节点
  （`CART_DUP=1`），再对底面/外面同步做**面内红色加密**（`CART_L=N`，支持
  hex8→4 / tri→3）。输出 `cartilage_mesh_t{N}d.npz`，含 `_tie_faces` / `_bone_faces`。
- `scripts/ankle_fe/build_feb.py` 加 `TIE=1` 分支（**默认关**，不影响既有路径）。
- 快扫工具 `temp/pyfebio_demo/fasttie.py`（只改已写好的 .feb 再跑 FEBio）
  与 `steps_time.py`（读每个 state 的 time + 轴向支反力）。

### 结果：罚参数单调从下方收敛到共形解

| penalty | R @ t=1.0 | 偏差 |
|---|---|---|
| 1 | 404.672 | −7.4% |
| 2 | 420.930 | −3.6% |
| **4** | **432.497** | **−0.99%** |
| 8+ | 首步即崩 | — |

**L=0 时几何与共形模型完全相同**，所以 436.817 N 就是真值 ⇒ 收敛到 −0.99%
即证明 mortar 正确。甜点 `penalty=4`。

### 五个坑（LRN-046~050）

1. **046 tie 面必须是单元真实 facet** —— 软骨已整体拆成 tet4，它的 facet 只能是
   **三角形**；照抄 npz 的四边形 ⇒ FEBio 报 `61 invalid facets` 并**静默丢弃**。
2. **047 secondary 面有洞 ⇒ 伪装成「模型硬 2.2×」** —— 只挑「与底面坐标重合」的
   骨面时 tibia 只匹配 63/93、talus 166/200；缺口让 primary 节点被 `search_radius`
   拉到远处 facet 上锁死，支反力虚假升到 709 → 826 → 969 N。
   **判据：随 penalty 单调上升背离真值 = 有洞；从下方收敛 = 正常。**
   改用「骨外面中靠近软骨底面的完整局部片」（117/118 面，饱和后不随 TOL 变）。
3. **048 必须同 `t=1.0` 比较** —— 步长自适应，`TIME_STEPS=20` 只报「14 步完成」，
   但 xplt 里 `time=1.0`。只看步数会误判两个模型不可比。
4. **049 参数扫描必须与建网格解耦** —— 重建 .feb 一轮 15-20 s，扫参数却不动网格。
   快扫后单轮 2-40 s，并**直接救回结论**：早前 12 例 0/12 通过，快扫后立刻暴露
   「penalty 1/2/4 都能过」⇒ 失败源于面，不是罚参数。
5. **050 面内红加密本身仍被母单元扭曲挡住**（负面结果）—— 加密实现正确
   （wedge6 3 个子单元 0 翻负、L=0 体积偏差 **+0.0000%**、wedge 加密后 |V|min
   改善 4×），但 **hex8 恒有 ~30% 子单元翻负**（69/244、97/300、296/976），
   且 4 个子位置分布均匀（20/19/16/14）⇒ **不是绕向 bug**，是母单元严重扭曲：
   同一面上 offset 方向从「主要 −x」摆到「主要 −z」（60–90°），因为
   `CART_OFFSET=toward`（指向对侧骨）在关节边缘处方向不良定义；母单元
   |V|min 只有 0.0127 / 0.00066。

### 下一步（P1d2b）
**先修几何再加密。** 已实测：
- `CART_OFFSET=normal` 单独用 → |V|min 0.0063（**更差**）
- `normal + ART_ALIGN=0.5` → |V|min **0.429（34×）**，但关节面节点 73→33（砍 55%）

建议方向：在**曲面参数域**内做重划分/光顺（而不是在 3D 偏移空间加密），
或对边缘节点单独用**切平面法向**替代顶点法向。

### 口径纠正
MEMORY.md 是 CRLF，`len(text.encode())` 比 `stat().st_size` **每行少 1 字节**
（146 行 ⇒ 少 146）。**10 KB 上限必须按磁盘口径判**，否则会误以为达标
（本轮就是这样先误判了 10,189 达标、实际 10,334 超限）。顺带删掉一行重复的
「详见」指针。
"""

if MEMF.exists():
    t = MEMF.read_text(encoding="utf-8")
    if "十一、tied-elastic" not in t:
        MEMF.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"memory note -> {MEMF.stat().st_size:,} B")
else:
    print("MISS", MEMF)

# ---- 最终只读校验 ----
print("\n" + "=" * 70)
M = WS / "MEMORY.md"
sz = M.stat().st_size
print(f"MEMORY.md disk = {sz:,} / 10240  余量 {10240 - sz:+,}  "
      f"{'OK' if sz <= 10240 else '超限'}")
s = M.read_text(encoding="utf-8")
for k, lbl in [("tied-elastic", "tied-elastic"), ("432.4", "432.4"),
               ("cart_refine", "cart_refine"), ("LRN-046", "LRN-046"),
               ("CONV_TODO", "CONV 指针"), ("ANKLE_FIX_TODO", "修正计划")]:
    print(f"  [{'OK' if k in s else 'XX'}] {lbl}")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", s)))
dead = [p for p in refs
        if not ((WS / p).exists() or (ROOT / p).exists())
        and not p.startswith(("memory/", "references/", "digest/", "docs/"))
        and not list(ROOT.glob(f"**/{Path(p).name}"))]
print(f"  引用 {len(refs)} / 死链 {len(dead)} {dead}")

print("\n--- 本轮文档 ---")
for k, p in [("FE_PIPELINE_TODO.md", ROOT / "FE_PIPELINE_TODO.md"),
             ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
             ("CONV_TODO.md", ROOT / "CONV_TODO.md"),
             ("LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md"),
             ("memory note", MEMF)]:
    t = p.read_text(encoding="utf-8") if p.exists() else ""
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>8,}  "
          f"{k:24s} tied-elastic={'Y' if 'tied-elastic' in t else 'n'}")
