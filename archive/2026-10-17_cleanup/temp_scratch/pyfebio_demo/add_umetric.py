"""给读数加第二种口径：韧带节点的**绝对位移** |u|（吴恺"位移（拉伸长度）"有歧义）。"""
from pathlib import Path
import ast

p = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\read_lig.py")
s = p.read_text(encoding="utf-8")

# 在每条韧带的统计里加 |u|
old = """    ds = np.array(ds)
    res.append((float(ds.mean()), nm, len(ds), float(ds.max()), float(eas[i])))
    print(f"{nm:34s} {len(ds):>4} {ds.mean():>9.4f} {ds.max():>9.4f} "
          f"{eas[i]:>8.1f}")"""
assert s.count(old) == 1
s = s.replace(old, """    ds = np.array(ds)
    # 口径二：该韧带**全部节点的位移模**（吴恺"位移"可能指这个）
    nd = sorted({int(x) - 1 for pr in pairs[i][:npair[i]] for x in pr})
    umag = float(np.linalg.norm(D[nd], axis=1).mean()) if nd else 0.0
    umax = float(np.linalg.norm(D[nd], axis=1).max()) if nd else 0.0
    res.append((float(ds.mean()), nm, len(ds), float(ds.max()),
                float(eas[i]), umag, umax))
    print(f"{nm:34s} {len(ds):>4} {ds.mean():>9.4f} {ds.max():>9.4f} "
          f"{eas[i]:>8.1f} {umag:>9.4f} {umax:>9.4f}")""")

old2 = '''print(f"{'韧带':34s} {'纤维':>4} {'平均ΔL':>9} {'最大ΔL':>9} {'EA(MPa)':>8}")
print("-" * 78)'''
assert s.count(old2) == 1
s = s.replace(old2, '''print(f"{'韧带':34s} {'纤维':>4} {'平均ΔL':>9} {'最大ΔL':>9} "
      f"{'EA':>7} {'|u|均':>9} {'|u|max':>9}")
print("-" * 96)''')

old3 = """res.sort(reverse=True)
for k, (m, nm, n, mx, ea) in enumerate(res, 1):
    print(f"  {k:2d}. {nm:34s} {m:8.4f} mm   (EA={ea:.1f})")"""
assert s.count(old3) == 1
s = s.replace(old3, """res.sort(reverse=True)
print("\\n【口径一 ΔL = 拉伸长度变化】")
for k, (m, nm, n, mx, ea, ug, ux) in enumerate(res, 1):
    print(f"  {k:2d}. {nm:34s} {m:8.4f} mm   (EA={ea:.1f})")
print("\\n【口径二 |u| = 韧带节点绝对位移均值】")
for k, (_, nm, _, _, ea, ug, ux) in enumerate(
        sorted(res, key=lambda t: -t[5]), 1):
    print(f"  {k:2d}. {nm:34s} {ug:8.4f} mm   (EA={ea:.1f})")""")

old4 = """sub = [(m, WU[nm]) for m, nm, *_ in res if nm in WU]"""
assert s.count(old4) == 1
s = s.replace(old4, """sub = [(m, WU[nm]) for m, nm, *_ in res if nm in WU]
usub = [(ug, WU[nm]) for _, nm, _, _, _, ug, _ in sorted(
    res, key=lambda t: -t[5]) if nm in WU]""")

old5 = """if sub:
    print("\\n吴恺要的四条（本模型 vs 吴恺 5.24/5.04/2.05/1.92）：")
    print("  吴恺:   胫跟 > 胫舟 > 跟腓 > 下胫腓前")
    print("  本模型: " + " > ".join(nm for _, nm in sub))
    for m, nm in sub:
        print(f"    {nm:8s} {m:8.4f} mm")"""
assert s.count(old5) == 1
s = s.replace(old5, """if sub:
    print("\\n吴恺要的四条（吴恺: 胫跟 5.24 > 胫舟 5.04 > 跟腓 2.05 > 下胫腓前 1.92）")
    print("  口径一 ΔL    : " + " > ".join(nm for _, nm in sub))
    print("  口径二 |u|   : " + " > ".join(nm for _, nm in usub))
    for m, nm in sub:
        print(f"    {nm:8s} ΔL={m:8.4f} mm")""")

old6 = """if react is not None:
    R = np.atleast_2d(react)
    a = r["axis"]
    kax = int(np.argmax(np.abs(a)))
    print(f"\\n支反力沿加载轴分量之和 = {R[-1][:, kax].sum():.3f}")"""
assert s.count(old6) == 1
s = s.replace(old6, """if react is not None:
    R = np.asarray(react, dtype=float).reshape(-1, 3)
    kax = int(np.argmax(np.abs(r["axis"])))
    print(f"\\n支反力沿加载轴分量之和 = {R[:, kax].sum():.3f}")""")

p.write_text(s, encoding="utf-8")
ast.parse(s)
print("patched + syntax OK")
