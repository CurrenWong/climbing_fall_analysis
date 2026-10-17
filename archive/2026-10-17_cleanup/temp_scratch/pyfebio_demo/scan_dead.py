"""死代码扫描：同一 block 内 return/continue/break/raise 之后仍有语句。"""
import ast
import io
import sys
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
targets = sys.argv[1:] or ["scripts/ankle_fe/ankle_lig_feb.py"]
KINDS = (ast.FunctionDef, ast.AsyncFunctionDef, ast.For, ast.While,
         ast.If, ast.With, ast.AsyncWith, ast.Try, ast.Module)
TERM = (ast.Return, ast.Continue, ast.Break, ast.Raise)

total = 0
for t in targets:
    p = ROOT / t
    if not p.exists():
        print(f"  (skip, 不存在) {t}")
        continue
    src = io.open(p, encoding="utf-8").read()
    lines = src.split("\n")
    tree = ast.parse(src)
    hits = []
    for node in ast.walk(tree):
        if not isinstance(node, KINDS):
            continue
        body = getattr(node, "body", None)
        if not body:
            continue
        for i, st in enumerate(body[:-1]):
            if isinstance(st, TERM):
                nxt = body[i + 1]
                hits.append((st.lineno, type(st).__name__,
                             nxt.lineno, lines[nxt.lineno - 1].strip()[:72]))
    if hits:
        print(f"\n=== {t} ===")
        for a, k, b, txt in hits:
            print(f"  !! line {a:>4} {k:8} -> line {b:>4}: {txt}")
    total += len(hits)
print(f"\n合计 {total} 处不可达代码")
