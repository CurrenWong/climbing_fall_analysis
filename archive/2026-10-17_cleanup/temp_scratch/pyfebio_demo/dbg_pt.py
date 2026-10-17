import sys
import numpy as np
sys.path.insert(0, r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo")
from pt_tri import dist_point_to_tris

T = np.random.rand(5, 3, 3)
p = np.random.rand(3)
print("T", T.shape, "p", p.shape)
print("dist", dist_point_to_tris(p, T))
