"""MEMORY.md 裁到 ≤10240；并把脚本引用写成可解析路径。"""
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
M = WS / "MEMORY.md"
LIMIT = 10240

s = M.read_text(encoding="utf-8")

# ① 脚本引用改成可解析的相对路径（同时略微缩字节）
s = s.replace("**`build_feb.py` ≠ `ankle_lig_feb.py`**（前者定死、**无开关**）；",
              "**`scripts/ankle_fe/build_feb.py` ≠ `.../ankle_lig_feb.py`**（前者无开关）；")

# ② 逐项裁剪到达标
CUTS = [
    ("| 22 | **`scripts/ankle_fe/build_feb.py` ≠ `.../ankle_lig_feb.py`**（前者无开关）；"
     "**cmd `set` 污染整条命令行** ⇒ 控制组必须无条件 `env.pop` | "
     "`ANKLE_FIX_TODO.md` + `[LRN-20261004-042/045]` |",
     "| 22 | **`build_feb.py` ≠ `ankle_lig_feb.py`**（前者无开关）；**cmd `set` 污染整条"
     "命令行** ⇒ 控制组必须无条件 `env.pop` | `[LRN-041~045]` |"),
    ("- **入口三件**：接手=`FE_PIPELINE_TODO.md`「🧭 接手入口」｜\n"
     "  修正计划=`ANKLE_FIX_TODO.md`｜**真值**=`docs/Anderson双论文_关键信息.md`",
     "- **入口**：接手=`FE_PIPELINE_TODO.md`「🧭 接手入口」｜计划=`ANKLE_FIX_TODO.md`｜"
     "真值=`docs/Anderson双论文_关键信息.md`"),
    ("- **韧带线已降级**：Tochigi 实测站姿下韧带本就松弛（CFL 0.2~2.1%）⇒\n"
     "  我的 ΔL≈0 **与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **与实验不符**",
     "- **韧带线已降级**：Tochigi 实测站姿下韧带本就松弛（CFL 0.2~2.1%）⇒ 我的 ΔL≈0 "
     "**与实验一致**；吴恺 1.92~5.24mm 隐含 8~21% **与实验不符**"),
    ("- **Anderson 真值**：all-hex｜距骨**弱弹簧 500~50 N/mm + 解除内外翻内旋**｜\n"
     "  骨刚性｜**不建韧带**｜软骨 E=12 ν=0.42、厚 **1.5~1.7mm**｜精度 **10~15%**",
     "- **Anderson 真值**：all-hex｜距骨**弱弹簧 +解除内外翻内旋**｜骨刚性｜**不建韧带**｜"
     "软骨 E=12 ν=0.42、厚 **1.5~1.7mm**｜精度 **10~15%**"),
]
for a, b in CUTS:
    if a in s:
        s = s.replace(a, b, 1)
    M.write_text(s, encoding="utf-8")
    sz = M.stat().st_size
    print(f"  -> {sz:,} (余量 {LIMIT - sz:+,})")

sz = M.stat().st_size
if sz > LIMIT:      # 还超则继续砍 gotcha 22 的尾部说明
    s = s.replace("**cmd `set` 污染整条命令行** ⇒ 控制组必须无条件 `env.pop` | `[LRN-041~045]` |",
                  "**cmd `set` 污染整条命令行** | `[LRN-042/045]` |")
    M.write_text(s, encoding="utf-8")
    sz = M.stat().st_size
    print(f"  二次裁剪 -> {sz:,} (余量 {LIMIT - sz:+,})")
print(f"\n最终 disk = {sz:,} / {LIMIT}  {'OK' if sz <= LIMIT else '超限'}")

# ③ 复核两个脚本真实位置（消"死链"误报）
for nm in ("build_feb.py", "ankle_lig_feb.py"):
    for cand in (ROOT / "scripts/ankle_fe" / nm, ROOT / nm):
        if cand.exists():
            print(f"  {nm} -> {cand.relative_to(ROOT)}  ({cand.stat().st_size:,} B)")
            break
