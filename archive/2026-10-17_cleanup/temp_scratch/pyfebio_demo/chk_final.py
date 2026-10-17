"""最终内容抽查（跨文档关键论点是否在位）。"""
import io

W = r"C:\Users\curre\.copaw\workspaces\default"
R = r"D:\Project\climbing_fall_analysis"

CHECKS = [
    (R + r"\ANKLE_ROADMAP_TODO.md",
     ["「519 N 上限」的机制", "距骨低端刚体转动", "615.4 N @ +2.5 mm",
      "A4 首轮结果", "cart_m3", "SUBTALAR_LAUGON=PENALTY", "A1+A2 已打通"]),
    (W + r"\MEMORY.md",
     ["519 N」是收敛代价", "615.4N@+2.5mm", "SUBTALAR_LAUGON=PENALTY", "+14.2%",
      "cart_m1_fix`/`cart_m3` 已否证"]),
    (W + r"\memory\2026-10-06.md",
     ["「519 N 上限」机制查清", "A4 首轮", "距骨转动"]),
]
ok = True
for f, keys in CHECKS:
    t = io.open(f, encoding="utf-8").read()
    miss = [k for k in keys if k not in t]
    nm = f.replace("\\", "/").split("/")[-1]
    print(f"{nm:<26} 缺失 {miss or '无 ✓'}")
    ok &= not miss

# 旧错句是否残留
print("\n--- 旧/错结论残留检查 ---")
STALE = [
    (W + r"\MEMORY.md", "已否证（2026-10-06 实测：接入后 F **+0.0000%**"),
    (R + r"\ANKLE_ROADMAP_TODO.md", "接入距下关节的真实效应 = +0.0000%"),
    (W + r"\MEMORY.md", "单接关节 **零效应**"),
]
for f, s in STALE:
    t = io.open(f, encoding="utf-8").read()
    nm = f.replace("\\", "/").split("/")[-1]
    print(f"{nm:<26} {'⚠ 仍在：' + s[:40] if s in t else '已清 ✓'}")
    ok &= s not in t
print("\n结论：", "PASS ✅" if ok else "有缺 ⚠")
