"""压缩前收口：MEMORY.md 踝关节段重写 + 新增 gotcha 22，并做全量校验。"""
import json
import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240

s = M.read_text(encoding="utf-8")

# ── ① 踝关节段重写（行 100–111）──
old_sec = """- **接手入口**：`FE_PIPELINE_TODO.md`「🧭 接手入口」节；**刚度问题**→`CONV_TODO.md`
- ⚠️ **刚度有 2 倍量级不确定度（未归因）**：同几何/材料下 tet4 **259.3 N** vs
  penta6 **116.9 N**。五假设已排除（格式/锁死/几何/剔除/接触面积），
  最尖锐=**接触面定向**。⇒ 600 N 档**应力/面积双靶可用**，
  **载荷-位移不可信**
- **管线**：C/B/A 全 ✅（`rc=0`/16.6 s，436.8 N，软骨 p99 2.58/4.27
  → 外推 600 N ≈3.4 MPa 对吴恺 3.74）
- **韧带轮**（→`ANKLE_LIG_TODO.md`）：4 条已定位（胫跟/跟腓/下胫腓/胫舟）；
  EA 510/510/**75.4**、厚 2.0；`TRUNC_MM=60` 截断胫腓骨
  （**160 倍**，5.4 s）
- ❌ **排序=负结果**：ΔL 0.001~0.08 vs 吴恺 1.92~5.24 mm；胫跟颠倒
- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十"""

new_sec = """- **入口三件**：接手=`FE_PIPELINE_TODO.md`「🧭 接手入口」｜
  修正计划=`ANKLE_FIX_TODO.md`｜**真值**=`docs/Anderson双论文_关键信息.md`
- ⚠️ **刚度 2 倍不确定度**：tet4 **259.3** vs penta6 **116.9** N；
  **材料已排除**（4 组参数比值恒 2.22~2.35）⇒ 剩离散/接触
- ⭐ **根因取向 = 接触区欠分辨**：关节面仅 **73 节点**（~2.5mm 间距）；
  加密到 692 节点后软骨体积涨 **2.3×** ⇒ 粗网格低估面积
- **Anderson 真值**：all-hex｜距骨**弱弹簧 500~50 N/mm + 解除内外翻内旋**｜
  骨视为刚体｜**不建韧带**｜软骨 E=12 ν=0.42、厚 **1.5~1.7mm**｜精度 **10~15%**
- **韧带线已降级**：Tochigi 实测站姿下韧带本就松弛（CFL 0.2~2.1%）⇒
  我的 ΔL≈0 **与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **与实验不符**
- **详见**：`memory/2026-10-03/ankle-fe-modeling.md` §六–§十"""

assert old_sec in s, "踝关节段未找到"
s = s.replace(old_sec, new_sec, 1)

# ── ② 新增 gotcha 22（插在 gotcha 21 之前）──
anchor21 = "| 21 | **FE 接触建模四铁律"
row22 = ("| 22 | **`build_feb.py` ≠ `ankle_lig_feb.py`**（前者定死、**无开关**）；"
         "**cmd `set` 污染整条命令行** ⇒ 控制组必须无条件 `env.pop` | "
         "`ANKLE_FIX_TODO.md` + `[LRN-20261004-042/045]` |\n")
if "| 22 |" not in s:
    s = s.replace(anchor21, row22 + anchor21, 1)

# ── ③ 落盘 + 裁剪到达标 ──
M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"改后 disk = {sz:,} / {LIMIT}  (余量 {LIMIT - sz:+,})")
for a, b in (("  **材料已排除**（4 组参数比值恒 2.22~2.35）⇒ 剩离散/接触",
              "  **材料已排除**（比值恒 2.22~2.35）⇒ 剩离散/接触"),
             ("- ⭐ **根因取向 = 接触区欠分辨**：关节面仅 **73 节点**（~2.5mm 间距）；\n"
              "  加密到 692 节点后软骨体积涨 **2.3×** ⇒ 粗网格低估面积",
              "- ⭐ **根因取向=接触区欠分辨**：仅 **73 节点**(~2.5mm)；加密到 692 后\n"
              "  软骨体积涨 **2.3×** ⇒ 粗网格低估面积"),
             ("- **Anderson 真值**：all-hex｜距骨**弱弹簧 500~50 N/mm + 解除内外翻内旋**｜\n"
              "  骨视为刚体｜**不建韧带**｜软骨 E=12 ν=0.42、厚 **1.5~1.7mm**｜精度 **10~15%**",
              "- **Anderson 真值**：all-hex｜距骨**弱弹簧 500~50 N/mm + 解除内外翻内旋**｜\n"
              "  骨刚性｜**不建韧带**｜软骨 E=12 ν=0.42、厚 **1.5~1.7mm**｜精度 **10~15%**")):
    if sz > LIMIT and a in s:
        s = s.replace(a, b, 1)
        M.write_text(s, encoding="utf-8")
        sz = M.stat().st_size
        print(f"  裁剪 -> {sz:,} (余量 {LIMIT - sz:+,})")
print(f"最终 disk = {sz:,} / {LIMIT}  {'OK' if sz <= LIMIT else '超限'}")

# ── ④ 全量校验 ──
print("\n" + "=" * 74)
docs = {
    "docs/Anderson双论文_关键信息.md": ROOT / "docs/Anderson双论文_关键信息.md",
    "ANKLE_FIX_TODO.md": ROOT / "ANKLE_FIX_TODO.md",
    "CONV_TODO.md": ROOT / "CONV_TODO.md",
    "ANKLE_600N_TODO.md": ROOT / "ANKLE_600N_TODO.md",
    "ANKLE_LIG_TODO.md": ROOT / "ANKLE_LIG_TODO.md",
    "FE_PIPELINE_TODO.md": ROOT / "FE_PIPELINE_TODO.md",
    "docs/吴恺论文_重解析关键信息.md": ROOT / "docs/吴恺论文_重解析关键信息.md",
    ".learnings/LEARNINGS.md": ROOT / ".learnings/LEARNINGS.md",
    "refine3.py": ROOT / "temp/pyfebio_demo/refine3.py",
    "cmp_ref.py": ROOT / "temp/pyfebio_demo/cmp_ref.py",
    "norm_probe.py": ROOT / "temp/pyfebio_demo/norm_probe.py",
    "thums_lowerlimb_r1.npz": ROOT / "temp/thums/thums_lowerlimb_r1.npz",
    "r1/cartilage.json": ROOT / "temp/thums/r1/cartilage.json",
}
miss = [k for k, p in docs.items() if not p.exists()]
for k, p in docs.items():
    print(f"  [{'OK ' if p.exists() else 'MISS'}] {p.stat().st_size if p.exists() else 0:>7,}  {k}")

lrn = (ROOT / ".learnings/LEARNINGS.md").read_text(encoding="utf-8")
print("\n本轮 LRN：")
for n in (39, 40, 41, 42, 43, 44, 45):
    k = f"LRN-20261004-0{n}"
    print(f"  [{'OK' if k in lrn else 'XX'}] {k}")

print("\nMEMORY.md 死链：")
paths = sorted(set(re.findall(r"`([\w\-./]+\.(?:md|py|npz|json|feb|k))`", s)))
dead = [p for p in paths
        if not ((WS / p).exists() or (ROOT / p).exists()
                or p.startswith(("memory/", "references/", "digest/", "docs/")))]
for p in dead:
    print(f"  [DEAD] {p}")
print(f"  {len(paths)} 引用 / {len(dead)} 死链")
print(f"\nALL = {'PASS' if not miss and not dead else 'CHECK'}")
