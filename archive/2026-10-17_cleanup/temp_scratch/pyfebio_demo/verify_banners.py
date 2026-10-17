"""精确复查横幅状态（用旧横幅的特征句，而不是日期串）。"""
import io
import os

R = r"D:\Project\climbing_fall_analysis"
OLD_MARK = "加载方向符号一直是反的"        # 旧横幅特征句
NEW_MARK = "本注记写于 2026-10-06"          # 新横幅特征句
FILES = ["ANKLE_FIX_TODO.md", "ANKLE_GOAL_TODO.md", "ANKLE_600N_TODO.md",
         "TIE_PORT_TODO.md", "ANKLE_LIG_TODO.md", "FE_PIPELINE_TODO.md",
         "docs/总体进度.md"]
ok = True
for rel in FILES:
    p = os.path.join(R, rel)
    t = io.open(p, encoding="utf-8").read()
    a, b = t.count(OLD_MARK), t.count(NEW_MARK)
    good = (a == 0 and b == 1)
    ok &= good
    print(f"  {rel:<32} 旧错句 {a}  新横幅 {b}  {'✓' if good else '⚠'}")
    if a:
        for i, l in enumerate(t.splitlines()):
            if OLD_MARK in l:
                print(f"       行 {i+1}: {l.strip()[:120]}")

# 使用说明里是否还有"必开"或错误方向
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
d = io.open(D, encoding="utf-8").read()
print(f"\n使用说明：含 '必开' = {'是 ⚠' if '必开' in d else '否 ✓'}；"
      f"含 'DISP_MM>0 才是压实' = {'是 ✓' if 'DISP_MM > 0` 才是压实' in d else '否 ⚠'}")
ok &= ("必开" not in d)
print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
