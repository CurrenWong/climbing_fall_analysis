"""P2c 决定性测量：+1.6mm 末态下软骨单元还剩多少 det(J)（用末态位移算）。"""
import re
import sys
from pathlib import Path

import numpy as np
import h5py

R = Path(r"D:\Project\climbing_fall_analysis")
H5 = R / "temp/thums/febio_lig/ankle_lig_WK_stress.h5"
FEB = R / "temp/thums/febio_lig/ankle_lig_WK.feb"

# ── h5 结构探查 ──
with h5py.File(H5, "r") as f:
    print("顶层:", list(f.keys()))
    st = f["states"]
    sids = sorted(st.keys(), key=int)
    last = st[sids[-1]]
    print(f"末 state={sids[-1]} 子项: {list(last.keys())}")
    for k in last:
        try:
            print(f"   {k}: {list(last[k].keys())}")
        except AttributeError:
            pass
    if "meshes" in f:
        print("meshes/0:", list(f["meshes/0"].keys()))
        nd = f["meshes/0/nodes"]
        print(f"  nodes: Dataset {getattr(nd, 'shape', '?')} {getattr(nd,'dtype','')}")

    # 找位移
    disp = None
    for cand in ("node_data", "nodes"):
        if cand in last and "displacement" in last[cand]:
            disp = np.asarray(last[cand]["displacement"][:])
            print(f"\n✓ 位移来自 states/<id>/{cand}/displacement {disp.shape}")
            break
    if disp is None:
        for k in last:
            if hasattr(last[k], "keys") and "displacement" in last[k]:
                disp = np.asarray(last[k]["displacement"][:])
                print(f"\n✓ 位移来自 states/<id>/{k}/displacement {disp.shape}")
                break

    # 参考坐标
    XYZ = {}
    for m in re.finditer(r'<node id="(\d+)">\s*([-\d.eE+]+)\s*,\s*([-\d.eE+]+)\s*,'
                         r'\s*([-\d.eE+]+)', FEB.read_text(encoding="latin-1")):
        XYZ[int(m.group(1))] = np.array([float(m.group(2)), float(m.group(3)),
                                         float(m.group(4))])
    ids = sorted(XYZ)
    X0 = np.array([XYZ[i] for i in ids])
    print(f"✓ 参考坐标 {X0.shape}（节点 {min(ids)}..{max(ids)}）")

    if disp is not None and len(disp) == len(X0):
        Xd = X0 + disp
    else:
        print(f"⚠ 位移 {None if disp is None else disp.shape} 与节点 {X0.shape} "
              f"对不上，改用 h5 网格坐标")
        Xd = None

# ── 逐单元 det(J)：参考 vs 末态 ──
txt = FEB.read_text(encoding="latin-1", errors="replace")
DOMS = {}
for m in re.finditer(r'<Elements\s+type="(penta6|hex8)"[^>]*?name="([^"]*)"'
                     r'[^>]*>(.*?)</Elements>', txt, re.S):
    kind, nm, body = m.group(2), m.group(1), m.group(3)
    conn = [[int(x) for x in c.split(",")] for c in
            re.findall(r'<elem[^>]*>\s*([\d,\s]+?)\s*</elem>', body)]
    DOMS[nm] = (kind, [c for c in conn if len(c) >= 4])


def detp(P):
    r_, s_, t_ = 1/3, 1/3, 1/2
    dN = np.array([[-(1-t_), -(1-t_), -(1-r_-s_)], [(1-t_), 0, -r_],
                   [0, (1-t_), -s_], [-t_, -t_, (1-r_-s_)], [t_, 0, r_],
                   [0, t_, s_]])
    return np.linalg.det(dN.T @ P)


print("\n" + "=" * 92)
print("软骨域 det(J)：参考态 vs 末态（+1.6mm）")
print("=" * 92)
print(f"{'域':<30}{'参考 min':>12}{'末态 min':>12}{'末态 p1':>12}"
      f"{'<0.3':>7}{'<0.1':>7}{'<0.01':>7}")
for nm, (kind, conn) in DOMS.items():
    if kind != "penta6" or "tibiotalar" not in nm:
        continue
    d0, d1 = [], []
    for c in conn:
        P0 = np.array([XYZ[n] for n in c], float)
        d0.append(detp(P0))
        if Xd is not None:
            P1 = np.array([Xd[n - min(ids)] for n in c], float)
            d1.append(detp(P1))
    d0 = np.array(d0)
    d1 = np.array(d1)
    print(f"{nm:<30}{d0.min():>12.4e}{d1.min():>12.4e}"
          f"{np.percentile(d1,1):>12.4e}{(d1<0.3).sum():>7}{(d1<0.1).sum():>7}"
          f"{(d1<0.01).sum():>7}")

print("\n" + "=" * 92)
print("关节闭合量（末态 vs 参考）")
print("=" * 92)
for nm, (kind, conn) in DOMS.items():
    if kind != "penta6" or "tibiotalar" not in nm:
        continue
    t0, t1 = [], []
    for c in conn:
        P0 = np.array([XYZ[n] for n in c], float)
        bot0, top0 = P0[[0, 1, 2]].mean(0), P0[[3, 4, 5]].mean(0)
        t0.append(np.linalg.norm(top0 - bot0))
        if Xd is not None:
            P1 = np.array([Xd[n - min(ids)] for n in c], float)
            bot1, top1 = P1[[0, 1, 2]].mean(0), P1[[3, 4, 5]].mean(0)
            t1.append(np.linalg.norm(top1 - bot1))
    t0, t1 = np.array(t0), np.array(t1)
    print(f"  {nm}")
    print(f"    厚度 参考: min {t0.min():.4f} 中位 {np.median(t0):.4f} mm")
    print(f"    厚度 末态: min {t1.min():.4f} 中位 {np.median(t1):.4f} mm")
    print(f"    ⇒ 压缩量: min {(t0-t1).min():.4f} 中位 {np.median(t0-t1):.4f} mm")
    print(f"    ⇒ 最薄的压剩 {t1.min()/t0.min()*100:.1f}%")
