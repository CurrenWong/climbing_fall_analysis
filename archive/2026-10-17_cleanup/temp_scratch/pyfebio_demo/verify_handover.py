"""压缩上下文前的交接完整性验证。"""
import re
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")

entry = (P / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
head = entry[entry.index("# 🧭 接手入口"):entry.index("## C. 解出 THUMS")]

print("=== 入口节引用的每个路径 ===")
paths = sorted(set(re.findall(r"`([\w\-./\\\\]+\.(?:py|md|json|npz|zip|k|fe?b))`", head)))
bad = []
for r in paths:
    cands = [P / r, WS / r, WS / "references" / r, WS / "skills" / r]
    hit = next((c for c in cands if c.exists()), None)
    print(f"  {'OK  ' if hit else 'MISS'} {r}")
    if hit is None:
        bad.append(r)
print(f"  -> {len(paths)} 个，{len(bad)} 个缺失\n")

print("=== 状态一致性（不应再有'进行中/阻塞'的过期说法）===")
stale = [("FE_PIPELINE_TODO.md", "A ⚠️ 进行中"),
         ("FE_PIPELINE_TODO.md", "卡在 Jacobian"),
         ("ANKLE_FE_TODO.md", "阻塞")]
for f, s in stale:
    t = (P / f).read_text(encoding="utf-8")
    hit = s in t
    # ANKLE_FE 的"阻塞"允许出现在 superseded 提示之外吗 -> 全部列出
    print(f"  {'❌ 仍存在' if hit else '✅ 已清除'}  {f}: {s!r}")
d = (P / "docs/踝关节FE建模交接.md").read_text(encoding="utf-8")
print(f"  {'❌ 仍存在' if '13.4 A 的剩余阻塞' in d else '✅ 已清除'}  "
      f"交接文档: '13.4 A 的剩余阻塞'")

print("\n=== 三份文档的当前状态声明 ===")
for f in ("FE_PIPELINE_TODO.md", "ANKLE_FE_TODO.md", "docs/踝关节FE建模交接.md"):
    t = (P / f).read_text(encoding="utf-8")
    st = [l.strip() for l in t.split("\n")[:20]
          if "状态" in l or "已收敛" in l or "已停止" in l]
    print(f"  {f}:")
    for s in st[:2]:
        print(f"      {s[:110]}")

print("\n=== 关键产物存在性 ===")
for r in ("temp/thums/thums_lowerlimb.npz", "temp/thums/bone_materials.json",
          "temp/thums/cartilage.json", "temp/thums/cartilage_mesh.npz",
          "temp/thums/febio/ankle_tibiotalar.feb",
          "temp/thums/febio/ankle_tibiotalar.log",
          "temp/thums/febio/ankle_tibiotalar.h5",
          "temp/pyfebio_demo/a6_reaction.py",
          "temp/pyfebio_demo/diag_facets.py",
          ".learnings/LEARNINGS.md",
          str(WS / "references/fe-modeling-gotchas.md")):
    p = P / r
    if not p.exists():
        p = WS / r
    print(f"  {'OK ' if p.exists() else 'MISS'} {p.stat().st_size if p.exists() else 0:>10,}  {r}")

