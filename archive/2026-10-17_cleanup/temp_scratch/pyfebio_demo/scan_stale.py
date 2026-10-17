"""扫描所有文档里**已作废的拉伸侧结论/数字**，列出需要更正的位置。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
PAT = {
    "592.6": r"592\.6",
    "209.55": r"209\.55",
    "187.77": r"187\.77",
    "187.90": r"187\.90",
    "202.29": r"202\.29",
    "235N上限": r"235\s*N|≤\s*235|上限 ~235",
    "力控制病态": r"力控制病态|力控制不可用|力控制.{0,6}失败",
    "关节初始张开": r"初始张开|中位间隙 0\.97|0\.973",
    "到不了600": r"600\s*N\s*到不了|到不了",
    "拉伸": r"拉伸",
    "LOAD_SIGN": r"LOAD_SIGN",
}
SKIP = (".venv", "temp", ".git", "node_modules", "__pycache__")
hits = {}
for root, dirs, files in os.walk(R):
    dirs[:] = [d for d in dirs if d not in SKIP]
    for f in files:
        if not f.endswith((".md", ".txt")):
            continue
        p = os.path.join(root, f)
        try:
            t = io.open(p, encoding="utf-8", errors="replace").read()
        except OSError:
            continue
        for k, rx in PAT.items():
            n = len(re.findall(rx, t))
            if n:
                rel = os.path.relpath(p, R)
                hits.setdefault(rel, {})[k] = n
for rel in sorted(hits):
    print(f"\n{rel}")
    for k, v in hits[rel].items():
        print(f"   {k:<14} x{v}")
