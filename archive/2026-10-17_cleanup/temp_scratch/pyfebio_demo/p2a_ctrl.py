"""正对照：合成正规楔形 penta6 ⇒ det(J) 应处处为正且恒定。"""
import numpy as np

dN = lambda r, s, t: np.array([  # noqa: E731
    [-(1-t), -(1-t), -(1-r-s)],
    [(1-t), 0, -r],
    [0, (1-t), -s],
    [-t, -t, (1-r-s)],
    [t, 0, r],
    [0, t, s],
])

PTS = [(1/6, 1/6, 1/2), (2/3, 1/6, 1/2), (1/6, 2/3, 1/2),
       (1/3, 1/3, 0.0), (1/3, 1/3, 1.0), (1/3, 1/3, 1/2),
       (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
       (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0),
       (0.5, 0.0, 0.0), (0.0, 0.5, 0.0), (0.5, 0.0, 1.0), (0.0, 0.5, 1.0),
       (1/3, 1/3, 0.25), (1/3, 1/3, 0.75)]

print("=" * 92)
print("正对照 1：正规楔形（底三角 (0,0,0),(1,0,0),(0,1,0)，高 1）")
print("=" * 92)
P = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
              [0., 0., 1.], [1., 0., 1.], [0., 1., 1.]])
ds = [(p, np.linalg.det(dN(*p).T @ P)) for p in PTS]
for p, d in ds:
    print(f"  r,s,t={tuple(round(x,3) for x in p)!s:<22} det(J) = {d:>10.5f}")
print(f"  ⇒ 全部为正？ {all(d > 0 for _, d in ds)}   范围 "
      f"[{min(d for _,d in ds):.5f}, {max(d for _,d in ds):.5f}]")

print("\n" + "=" * 92)
print("正对照 2：把同一楔形的顶面整体平移 0.5（仍是合法楔形）")
print("=" * 92)
P2 = P + np.array([0.5, 0.0, 0.0])
ds2 = [(p, np.linalg.det(dN(*p).T @ P2)) for p in PTS]
print(f"  ⇒ 全部为正？ {all(d > 0 for _, d in ds2)}   范围 "
      f"[{min(d for _,d in ds2):.5f}, {max(d for _,d in ds2):.5f}]")

print("\n" + "=" * 92)
print("正对照 3：故意造一个「几乎压扁」的楔形（顶面压到 z=0.01）")
print("=" * 92)
P3 = P.copy()
P3[[3, 4, 5], 2] = 0.01
ds3 = [(p, np.linalg.det(dN(*p).T @ P3)) for p in PTS]
print(f"  ⇒ 全部为正？ {all(d > 0 for _, d in ds3)}   范围 "
      f"[{min(d for _,d in ds3):.6f}, {max(d for _,d in ds3):.6f}]")
print("  ⇒ 若合成压扁楔形也报负值 ⇒ 我的度量点超出有效域，需修")
print("  ⇒ 若压扁楔形全正、而实测楔形有负 ⇒ 实测单元**真的翻转**")

print("\n" + "=" * 92)
print("结论用途")
print("=" * 92)
print("  正对照成立 ⇒ 我的 Gauss 度量可信 ⇒ 软骨里确有近奇异/翻转采样点")
print("  ⇒ 上轮「网格健康、P2 删除」的结论**必须撤回**；P2 原命题成立。")
