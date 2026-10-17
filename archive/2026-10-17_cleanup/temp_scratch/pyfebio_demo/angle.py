import numpy as np

a = np.array([0.6766, 0.1363, 0.7237])
a = a / np.linalg.norm(a)
for g, nm in ((-np.array([0, 0, 1.0]), "global -Z (重力向)"),
              (np.array([0, 0, 1.0]), "global +Z"),
              (-np.array([0, 1.0, 0]), "global -Y"),
              (np.array([1.0, 0, 0]), "global +X"),
              (-np.array([1.0, 0, 0]), "global -X"),
              (np.array([0, 1.0, 0]), "global +Y")):
    g = g / np.linalg.norm(g)
    ang = np.degrees(np.arccos(np.clip(float(a @ g), -1, 1)))
    print(f"  胫骨轴 与 {nm:18s} 夹角 = {ang:6.2f} deg")
print(f"\n  胫骨轴分量: {np.round(a, 4)}")
