"""P2c（重写）：+1.6mm 末态下软骨单元还剩多少 det(J)。"""
import re
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"

with h5py.File(H5, "r") as f:
    nd = f["meshes/0/nodes"][:]
    nid = np.asarray(nd["id"], int)
    Xref = np.column_stack([nd["x"], nd["y"], nd["z"]]).astype(float)
    print(f"参考坐标 {Xref.shape}  节点 id {nid.min()}..{nid.max()}")

    sids = sorted(f["states"].keys(), key=int)
    last = f["states"][sids[-1]]
    g = last["node_data"]["displacement"]
    raw = list(g.keys())
    num = sorted([k for k in raw if k.isdigit()], key=int)
    print(f"位移 group：{len(raw)} 项；数字键 {len(num)}；非数字 {[k for k in raw if not k.isdigit()]}")
    if num:
        D = np.array([np.asarray(g[k][:], float).ravel()[:3] for k in num])
    else:
        D = np.asarray(g["all"][:], float)
    print(f"位移 {D.shape}")

Xdef = Xref + D
print(f"位移范数: max {np.linalg.norm(D,axis=1).max():.4f} mm  "
      f"均值 {np.linalg.norm(D,axis=1).mean():.4f} mm")

# 序号映射：feb node id -> 行号
pos = {int(v): i for i, v in enumerate(nid)}

txt = FEB.read_text(encoding="latin-1", errors="replace")
DOMS = {}
for m in re.finditer(r'<Elements\s+type="(penta6|hex8)"[^>]*?name="([^"]*)"'
                     r'[^>]*>(.*?)</Elements>', txt, re.S):
    nm, kind, body = m.group(2), m.group(1), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 4])
print(f"域: {list(DOMS.keys())}")


def detp(P, kind="penta6"):
    if kind == "penta6":
        r_, s_, t_ = 1/3, 1/3, 1/2
        dN = np.array([[-(1-t_), -(1-t_), -(1-r_-s_)], [(1-t_), 0, -r_],
                       [0, (1-t_), -s_], [-t_, -t_, (1-r_-s_)], [t_, 0, r_],
                       [0, t_, s_]])
    else:
        raise NotImplementedError
    return np.linalg.det(dN.T @ P)


def thick(P, kind="penta6"):
    bot = P[[0, 1, 2]].mean(0) if kind == "penta6" else P[[0, 3, 4, 7]].mean(0)
    top = P[[3, 4, 5]].mean(0) if kind == "penta6" else P[[1, 2, 6, 5]].mean(0)
    return float(np.linalg.norm(top - bot))


print("\n" + "=" * 96)
print("软骨域：参考态 vs 末态(+1.6mm) 的 det(J) 与厚度")
print("=" * 96)
print(f"{'域':<30}{'detJ参考min':>13}{'detJ末态min':>13}{'末态p1':>12}"
      f"{'<0.5':>7}{'<0.2':>7}")
for nm, (kind, conn) in DOMS.items():
    if kind != "penta6" or "tibiotalar" not in nm:
        continue
    d0, d1, t0, t1 = [], [], [], []
    for c in conn:
        P0 = Xref[[pos[n] for n in c]]
        P1 = Xdef[[pos[n] for n in c]]
        d0.append(detp(P0))
        d1.append(detp(P1))
        t0.append(thick(P0))
        t1.append(thick(P1))
    d0, d1, t0, t1 = map(np.array, (d0, d1, t0, t1))
    print(f"{nm:<30}{d0.min():>13.4e}{d1.min():>13.4e}"
          f"{np.percentile(d1, 1):>12.4e}{(d1<0.5).sum():>7}{(d1<0.2).sum():>7}")
    print(f"    厚度 参考 min {t0.min():.4f} / 中位 {np.median(t0):.4f} mm")
    print(f"    厚度 末态 min {t1.min():.4f} / 中位 {np.median(t1):.4f} mm")
    comp = t0 - t1
    print(f"    压缩量 min {comp.min():.4f} / 中位 {np.median(comp):.4f} / "
          f"max {comp.max():.4f} mm")
    print(f"    ⇒ 最薄压剩 {t1.min()/t0.min()*100:.1f}%；"
          f"压缩量/厚度 最大的单元 = {(comp/np.maximum(t0,1e-9)).max()*100:.1f}%")

print("\n" + "=" * 96)
print("👉 判定")
print("=" * 96)
print("  若末态 det(J) 接近 0 或为负 ⇒ 「软骨被压穿」机制成立 ⇒ 根因=模型柔度太大（问题#4）")
print("  若末态 det(J) 仍健康 ⇒ 翻转发生在别处（骨/韧带）")
