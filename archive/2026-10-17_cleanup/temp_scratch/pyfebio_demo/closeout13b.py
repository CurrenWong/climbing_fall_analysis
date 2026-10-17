"""压缩前收口 ②：更新 docs/总体进度.md 的子线段 + 全量校验。"""
import re
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"

# ── ① 总览：子线表格补两行 ──
P = ROOT / "docs/总体进度.md"
s = P.read_text(encoding="utf-8")
OLD = ("| **P1d2b** 软骨面内红加密 | 🚧 | **hex8 侧恒 28% 子单元翻负**（母单元扭曲）；几何需先修。`[LRN-046~051]` |")
NEW = ("| **P1d2b** 软骨面内红加密 | ⚪ | **已改用重分网路线取代**（细分/平滑/改方向三条捷径实测排除）`[LRN-052]` |\n"
       "| **新** 关节底面重分网 `cart_patch.py` | 🟡 | **可建可解**：m=1..4 全 `rc=0`，底面格 1119→**5424**(m=2)→10848(m=3)；"
       "但 **tie 界面成新瓶颈**（`|u|max` 超指定 10×、F 不收敛）`[LRN-053]` |")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    P.write_text(s, encoding="utf-8")
print(f"docs/总体进度.md -> {P.stat().st_size:,} B  "
      f"[{'OK' if 'cart_patch.py' in P.read_text(encoding='utf-8') else 'XX'}]")

# ── ② 全量校验 ──
print("\n" + "=" * 70)
print("产物清单")
print("=" * 70)
ART = [
    ("docs/总体进度.md", ROOT / "docs/总体进度.md"),
    ("scripts/ankle_fe/cart_patch.py", ROOT / "scripts/ankle_fe/cart_patch.py"),
    ("scripts/ankle_fe/cart_hex.py", ROOT / "scripts/ankle_fe/cart_hex.py"),
    ("scripts/ankle_fe/build_feb.py", ROOT / "scripts/ankle_fe/build_feb.py"),
    ("CART_HEX_TODO.md", ROOT / "CART_HEX_TODO.md"),
    ("FE_PIPELINE_TODO.md", ROOT / "FE_PIPELINE_TODO.md"),
    ("ANKLE_FIX_TODO.md", ROOT / "ANKLE_FIX_TODO.md"),
    ("CONV_TODO.md", ROOT / "CONV_TODO.md"),
    ("LEARNINGS.md", ROOT / ".learnings/LEARNINGS.md"),
    ("m1.npz", ROOT / "temp/thums/cartilage_mesh_patch_m1.npz"),
    ("m2.npz", ROOT / "temp/thums/cartilage_mesh_patch_m2.npz"),
    ("m3.npz", ROOT / "temp/thums/cartilage_mesh_patch_m3.npz"),
    ("MEMORY.md", M),
]
for k, p in ART:
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  {k}")

print("\n" + "=" * 70)
print("关键结论是否在档")
print("=" * 70)
lt = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
ms = M.read_text(encoding="utf-8")
ov = (ROOT / "docs/总体进度.md").read_text(encoding="utf-8")
td = (ROOT / "CART_HEX_TODO.md").read_text(encoding="utf-8")
CHK = [
    ("LRN-046~053 全在档", all(f"LRN-20261004-0{n}" in lt for n in range(46, 54))),
    ("LRN-052 三条捷径证伪", "H3（被证伪）" in lt or "被证伪" in lt),
    ("LRN-053 对角线铁律", "同源同对角" in lt),
    ("LRN-053 AUGLAG 判据", "1/(11·TIME_STEPS)" in lt),
    ("LRN-053 |u|max 10×", "0.1472" in lt),
    ("MEMORY: cart_patch", "cart_patch.py" in ms),
    ("MEMORY: 对角铁律", "0-2" in ms),
    ("MEMORY: AUGLAG", "AUGLAG" in ms),
    ("总览: cart_patch 行", "cart_patch.py" in ov),
    ("TODO: 开关表", "CP_DIAG" in td and "AUGLAG" in td),
    ("TODO: 验收判据", "|u|max" in td),
]
for k, v in CHK:
    print(f"  [{'OK' if v else 'XX'}] {k}")

print("\n" + "=" * 70)
print("MEMORY.md 尺寸（磁盘口径）+ 死链")
print("=" * 70)
sz = M.stat().st_size
refs = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb))`", ms)))
dead = [p for p in refs
        if not ((WS / p).exists() or (ROOT / p).exists())
        and not p.startswith(("memory/", "references/", "digest/", "docs/"))
        and not list(ROOT.glob(f"**/{Path(p).name}"))]
print(f"  {sz:,} / 10240  余量 {10240 - sz:+,}  {'OK' if sz <= 10240 else '超限'}")
print(f"  引用 {len(refs)} / 死链 {len(dead)} {dead}")

print("\n" + "=" * 70)
print("git 可追踪性")
print("=" * 70)
for f in ("docs/总体进度.md", "CART_HEX_TODO.md", "scripts/ankle_fe/cart_patch.py"):
    r = subprocess.run(["git", "check-ignore", str(ROOT / f)], cwd=ROOT,
                       capture_output=True, text=True)
    print(f"  {f:38s} {'不进 git' if r.returncode == 0 else '进 git'}")

print("\n" + "=" * 70)
print("记忆文件")
print("=" * 70)
for k, p, tok in (("session note", WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md", "十三"),
                  ("daily note", WS / "memory/2026-10-04.md", "十三 软骨重分网")):
    t = p.read_text(encoding="utf-8") if p.exists() else ""
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>9,}  "
          f"{k:14s} 含§十三={'Y' if tok in t else 'n'}")
