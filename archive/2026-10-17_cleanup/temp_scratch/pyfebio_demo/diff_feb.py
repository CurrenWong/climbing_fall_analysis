"""结构化对比：原可跑模型 vs 新脚本的 Mesh/MeshDomains 段。"""
import re
from collections import Counter
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio")
L = W.parent / "febio_lig"

orig = (W / "ankle_tibiotalar.feb").read_text(encoding="latin-1")
newf = sorted(L.glob("ankle_lig_WK_tibia*.feb"))
print("新文件:", [p.name for p in newf])
if not newf:
    raise SystemExit("无新文件")
new = newf[-1].read_text(encoding="latin-1")


def stat(t, tag):
    print(f"\n=== {tag} ===")
    print(f"  总行数 {len(t.splitlines()):,}")
    # 元素块
    for m in re.finditer(r'<Elements type="(\w+)" name="([^"]+)">', t):
        pass
    blocks = re.findall(r'<Elements type="(\w+)" name="([^"]+)">', t)
    print(f"  Elements 块 {len(blocks)}: {Counter(blocks).most_common(6)}")
    etypes = Counter(b[0] for b in blocks)
    print(f"  类型分布: {dict(etypes)}")
    # 每块元素数
    big = sorted(((len(re.findall(r"<elem ", seg)), nm)
                  for nm, seg in
                  [(b[1], t.split(f'name="{b[1]}">')[1].split("</Elements>")[0])
                   for b in blocks]), reverse=True)
    print(f"  各块元素数 top6: {big[:6]}")
    # 计数
    print(f"  elem 总数 {len(re.findall(r'<elem ', t)):,}")
    print(f"  node 总数 {len(re.findall(r'<node ', t)):,}")
    print(f"  SolidDomain {len(re.findall(r'<SolidDomain', t))}")
    print(f"  Surface {len(re.findall(r'<Surface ', t))}")
    print(f"  NodeSet {len(re.findall(r'<NodeSet ', t))}")
    print(f"  Contact {len(re.findall(r'<contact ', t))}")
    print(f"  bc 数 {len(re.findall(r'<bc ', t))}")
    print(f"  surface_load {len(re.findall(r'<surface_load', t))}")
    doms = re.findall(r'<SolidDomain name="([^"]+)" mat="([^"]+)"', t)
    print(f"  域: {doms}")


stat(orig, "原模型 ankle_tibiotalar.feb")
stat(new, f"新脚本 {newf[-1].name}")
