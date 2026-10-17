import re
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
t = (ROOT / "FE_PIPELINE_TODO.md").read_text(encoding="utf-8")
head = t[t.index("# 🧭 接手入口"):t.index("## C. 解出 THUMS")]

print("=== 裸名出现在哪 ===")
for name in ("bone_materials.json", "thums_lowerlimb.npz"):
    for i, l in enumerate(head.split("\n"), 1):
        for h in re.findall(r"`([^`]+)`", l):
            if h == name:
                print(f"  行{i}: {l.strip()[:120]}")

print("\n=== 三个收尾检查（逐个单跑，避免输出被吞）===")
for label, ok in (
    ("入口标注已实跑验证", "本节命令已实跑验证" in t),
    ("过时状态已清除",
     "A ⚠️ 进行中" not in t and "卡在 Jacobian" not in t),
    ("交接文档互指入口",
     "先读 `FE_PIPELINE_TODO.md` 的「🧭 接手入口」节" in
     (ROOT / "docs/踝关节FE建模交接.md").read_text(encoding="utf-8")),
):
    print(f"  {label}: {ok}")
