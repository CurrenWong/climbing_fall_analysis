"""最终**只读**校验（不再改文件）。"""
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
L = 10240
s = M.read_text(encoding="utf-8")
sz = len(s.encode())
print(f"MEMORY.md = {sz:,} / {L}   余量 {L - sz:+,}   {'OK' if sz <= L else '超限'}")

print("\n--- 踝关节段 ---")
L2 = s.split("\n")
for i, ln in enumerate(L2):
    if "tied-elastic 验证通过" in ln:
        print("\n".join("  " + x for x in L2[i:i + 5]))
        break

print("\n--- 死链 ---")
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", s)))
dead = [p for p in refs
        if not ((WS / p).exists() or (ROOT / p).exists())
        and not p.startswith(("memory/", "references/", "digest/", "docs/"))
        and not list(ROOT.glob(f"**/{Path(p).name}"))]
print(f"  {len(refs)} 引用 / {len(dead)} 死链 {dead}")

print("\n--- 本轮 artifact ---")
for k, p in [("cart_refine.py", ROOT / "scripts/ankle_fe/cart_refine.py"),
             ("build_feb.py(TIE)", ROOT / "scripts/ankle_fe/build_feb.py"),
             ("cartilage_mesh_t0d.npz", ROOT / "temp/thums/cartilage_mesh_t0d.npz"),
             ("cartilage_mesh_t1d.npz", ROOT / "temp/thums/cartilage_mesh_t1d.npz"),
             ("fasttie.py", ROOT / "temp/pyfebio_demo/fasttie.py"),
             ("steps_time.py", ROOT / "temp/pyfebio_demo/steps_time.py"),
             ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
             ("LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md")]:
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {k}")

lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
print("\n--- LRN 046~050 ---")
for n in range(46, 51):
    k = f"LRN-20261004-0{n}"
    print(f"  [{'OK' if k in lrn else 'XX'}] {k}")
