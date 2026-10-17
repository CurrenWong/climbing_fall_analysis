"""内省 pyfebio 的类名/标签 —— pyfebio 由 FEBio schema 生成，是最权威的本地依据。"""
import inspect
import pyfebio

print("pyfebio 顶层:", [n for n in dir(pyfebio) if not n.startswith("_")])

for mod in ("loads", "loaddata", "boundary", "mesh", "meshdomains"):
    m = getattr(pyfebio, mod, None)
    if m is None:
        print(f"\n[{mod}] 不存在")
        continue
    names = [n for n in dir(m) if not n.startswith("_")]
    print(f"\n[{mod}] {names}")

# 找带 xml_tag / tag 的类
print("\n=== loads 模块内的类与标签 ===")
try:
    import pyfebio.loads as L
    for n in dir(L):
        if n.startswith("_"):
            continue
        o = getattr(L, n)
        if inspect.isclass(o):
            tag = getattr(o, "xml_tag", None) or getattr(o, "__xml_tag__", None)
            fields = getattr(o, "model_fields", None) or getattr(o, "__fields__", None)
            print(f"  {n:34s} tag={tag}  fields={list(fields)[:8] if fields else None}")
except Exception as exc:
    print("  导入失败:", type(exc).__name__, exc)
