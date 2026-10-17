"""最终校验 + 更新 ANKLE_LOCAL_TODO 为"更正后的路线3"。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

(ROOT / "ANKLE_LOCAL_TODO.md").write_text("""# ANKLE_LOCAL — 路线3：关节局部 + 骨刚性（对标 Anderson）

## ⚠️ 优先级已下调（2026-10-05，见 `[LRN-20261005-056]`）

**发现：关节几何没问题、`F` 也不是伪影。** 所以路线3**不是救火**，它现在的主要收益是：
① 去掉 400mm 长骨悬臂弯曲（让指定位移真正压缩关节）② 对齐 Anderson 设定便于对标。

**当前第一优先级改成了**：**提高接触求解器的罚因子上限**（`PENALTY ≥ 2` 全负 Jacobian ⇒
`F` 无法验收敛，真值可能高 1.5~2×）。详见 `[LRN-20261005-056]`。

## 现状（实测）
| 项 | 值 |
|---|---|
| 关节刚度 `k = F/δ` | 49 / 151 / 451 / **880** N/mm（罚 1e-3→1） |
| Anderson 独立估计 `E·A/t` | **1059 N/mm**（p=1 差 −17%） |
| 加载后接触面积 | **344 mm²** ≈ Anderson 生理 **300 mm²** ✓ |
| 互穿节点 | **6**（几何正常） |
| `PENALTY` 上限 | **1**（2/4/8/16 全负 Jacobian） |

## Anderson 配方（据 `docs/Anderson双论文_关键信息.md`）
| 项 | 值 |
|---|---|
| 软骨 E / ν | **12 MPa / 0.42** |
| 摩擦 | Coulomb **μ = 0.01** |
| 软骨厚度 | 1.5 mm（2006）→ **1.7 mm**（2007） |
| 网格 | **全六面体**，收敛需 **~50,000 单元** |
| 主从比 | 胫骨=master，单元尺寸比 **1:3** |
| 距骨约束 | 弱弹簧 **500 N/mm**（2006）/ A-P 弹簧 **50 N/mm**（2007） |
| 预载流程 | 轻接触 → **200 N 落位**（松开内翻/内旋约束）→ 正式加载 |
| 载荷 | 2006: 10→2800 N + 屈伸；2007: **600 N** |
| 骨 | **刚性**（reference nodes） |
| 韧带 | 不建 |

## 计划（更正后）
- [x] 读 Anderson 建模/边界/材料节 ✓
- [x] 诊断：关节不闭合？→ **关：几何正常**（互穿 6 节点、接触 344mm²≈生理 300mm²）
- [x] 诊断：F 是伪影？→ **不是**（k 随罚趋近 1059）
- [ ] **【新 P0】提高接触罚因子上限**（`symmetric_stiffness`/`node_reloc`/步长/`max_refs`），
      或把 AUGLAG 增广调到真收敛 ⇒ 目标 `k → 1059`
- [ ] 实现 `TRUNC_MM`（只取关节局部）
- [ ] 实现 `RIGID_BONES`（骨刚性，reference nodes）
- [ ] 距骨弱弹簧 500 N/mm + 解除内翻/内旋（替代现在的 talus 底面全固定）
- [ ] 软骨 E=12/ν=0.42 ✓（已有）+ 摩擦 μ=0.01
- [ ] 200 N 落位 → 600 N 正式加载
- [ ] 验证 `k` 是否落到 1059 附近
""", encoding="utf-8")
print(f"ANKLE_LOCAL_TODO.md -> {(ROOT/'ANKLE_LOCAL_TODO.md').stat().st_size:,} B")

print("\n" + "=" * 76)
print("最终校验")
print("=" * 76)
lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
ids = re.findall(r"\[LRN-(\d{8}-\d{3})\]", lt)
ANK = ["TIE_DIVERGE_TODO.md", "CART_HEX_TODO.md", "ANKLE_FIX_TODO.md",
       "ANKLE_LOCAL_TODO.md", "docs/总体进度.md"]
CHK = [
    ("LRN-20261005-056 已入库", "LRN-20261005-056" in lt),
    ("LRN-056 含撤回①", "节点稀疏 ⇒ 接触斑落在节点之间" in lt),
    ("LRN-056 含撤回②", "没看增量倍率" in lt),
    ("LRN-056 含 k 表 880/1059", "880.0" in lt and "1059 N/mm" in lt),
    ("LRN-056 含罚上限表", "| rc | **0** | 1 | 1 | 1 | 1 |" in lt),
    ("LRN-056 含撞号教训", "早已被攀岩跌落线占用" in lt),
    ("LRN-056 记中间错误", "法向定向不一致造成的伪影" in lt),
    ("LRN-054 仍在档", "LRN-20261005-054" in lt),
    ("MEMORY: 非伪影", "F 未收敛但非伪影" in ms),
    ("MEMORY: 引用 056", "LRN-054/056" in ms),
    ("TIE_DIVERGE 已修", "LRN-20261005-056" in
     (ROOT / "TIE_DIVERGE_TODO.md").read_text(encoding="utf-8")),
    ("ANKLE_FIX 引用 056", "LRN-20261005-056" in
     (ROOT / "ANKLE_FIX_TODO.md").read_text(encoding="utf-8")),
    ("CART_HEX 引用 056", "LRN-20261005-056" in
     (ROOT / "CART_HEX_TODO.md").read_text(encoding="utf-8")),
    ("总体进度 引用 056", "LRN-056" in
     (ROOT / "docs/总体进度.md").read_text(encoding="utf-8")),
    ("ANKLE_LOCAL 已更正", "LRN-20261005-056" in
     (ROOT / "ANKLE_LOCAL_TODO.md").read_text(encoding="utf-8")),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")
dup = sorted({x for x in ids if ids.count(x) > 1})
print(f"\n  LEARNINGS: {LT_SIZE if False else (ROOT/'.learnings/LEARNINGS.md').stat().st_size:,} B，"
      f"条目 {len(ids)}，重复 ID {len(dup)} 个（历史遗留）：{dup}")
sz = M.stat().st_size
print(f"  MEMORY.md {sz:,}/10240 余量 {10240-sz:+,} {'✓ OK' if sz<=10240 else '✗ 超'}")
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
print(f"  session note {SESS.stat().st_size:,} B  "
      f"[§十五 {'Y' if '十五、自我更正' in SESS.read_text(encoding='utf-8') else 'n'}]")
