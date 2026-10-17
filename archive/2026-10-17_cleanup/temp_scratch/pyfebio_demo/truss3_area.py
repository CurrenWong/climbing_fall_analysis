"""测出 elastic-truss 的截面积约定：给定位移 -> 读支反力 -> 反推 A。"""
import subprocess
from pathlib import Path
import numpy as np

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

# E=260, L=10 mm, 指定位移 1.0 mm -> 理论 F = E*A*dl/L = 26*A
E, L, DL = 260.0, 10.0, 1.0

src = (WORK / "et_iso.feb").read_text(encoding="latin-1")
for tag, extra in (("baseline", ""),
                   ("mdata_area5", """
  <MeshData>
    <ElementData name="area" elem_set="tr">
      <elem id="1">5.0</elem>
    </ElementData>
  </MeshData>""")):
    t = src.replace("  <LoadData>", extra + "\n  <LoadData>") if extra else src
    f = WORK / f"area_{tag}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    ok = "T E R M I N A T I O N" in o and "NORMAL" in o.replace(" ", "")
    print(f"\n=== {tag} rc={r.returncode} "
          f"{'OK' if 'NORMAL' in o.replace(' ', '') else 'FAIL'} ===")
    if "unrecognized" in o or "ERROR" in o:
        for l in o.split("\n"):
            if "ERROR" in l or "unrecognized" in l or "invalid" in l:
                print("   ", l.strip()[:140])
    # 读支反力
    h5 = f.with_suffix(".xplt")
    if not h5.exists():
        h5 = WORK / f"{f.stem}.xplt"
    try:
        import pyfebio
        d = pyfebio.xplt.to_hdf5(str(h5))
        import h5py
        with h5py.File(d, "r") as h:
            def walk(g, pre=""):
                for k in g:
                    if isinstance(g[k], h5py.Group):
                        walk(g[k], pre + "/" + k)
                    else:
                        yield pre + "/" + k, g[k]
            keys = [k for k, _ in walk(h)]
            rk = [k for k in keys if "reaction" in k.lower()]
            for k in rk[:2]:
                a = h[k][()]
                print(f"    {k}: shape={getattr(a,'shape',None)}")
                if a.size:
                    ax = np.atleast_2d(a)
                    print(f"      末行 = {ax[-1]}")
                    fz = float(np.atleast_1d(ax[-1])[-1])
                    print(f"      |F|≈{abs(fz):.4f} N  => 隐含 A = {abs(fz)/(E*DL/L):.4f} mm²")
    except Exception as e:
        print("    读 xplt 失败:", type(e).__name__, e)
