from pathlib import Path

t = Path(r"temp/pyfebio_demo/_cr1fix.txt").read_text(encoding="utf-8", errors="replace")
KEEP = ("[", "单元", "负体积", "体积", "底面节点", "写出", "软骨节点")
for ln in t.split("\n"):
    if any(k in ln for k in KEEP):
        print(ln)
