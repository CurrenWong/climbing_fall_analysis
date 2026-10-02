"""只跑 visualize_b2.py 里的部分图（分批执行，避免一次跑太久）。"""
import sys, pathlib, time
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import visualize_b2 as V

which = sys.argv[1] if len(sys.argv) > 1 else "all"

if which in ("pad_matrix", "all"):
    print("--- 峰值压力热力图：软垫 ---", flush=True)
    V.fig_pressure_matrix(on_pad=True)
if which in ("hard_matrix", "all"):
    print("--- 峰值压力热力图：硬地面 ---", flush=True)
    V.fig_pressure_matrix(on_pad=False)
if which in ("rope", "all"):
    print("--- 绳索图 ---", flush=True)
    V.fig_rope_constitutive()
    V.fig_rope_devices()
print("ALL DONE", flush=True)
