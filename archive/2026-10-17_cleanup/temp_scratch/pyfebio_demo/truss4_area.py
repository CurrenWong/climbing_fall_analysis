"""修正版：测 elastic-truss 截面积，并试 MeshData 的写法。"""
import subprocess
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, L, DL = 260.0, 10.0, 1.0

src = (WORK / "et_iso.feb").read_text(encoding="latin-1")

CASES = {
    "baseline": src,
    "elemdata": src.replace("  <LoadData>", """  <MeshData>
    <ElementData name="area" elem_set="tr">
      <elem id="1">5.0</elem>
    </ElementData>
  </MeshData>
  <LoadData>"""),
    "bd_area": src.replace('<BeamDomain name="tr" mat="lig" type="elastic-truss"/>',
                           '<BeamDomain name="tr" mat="lig" type="elastic-truss">'
                           '<area>5.0</area></BeamDomain>'),
}

for tag, t in CASES.items():
    f = WORK / f"a_{tag}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    normal = "NORMAL" in o.replace(" ", "")
    print(f"\n=== {tag} rc={r.returncode} {'OK' if normal else 'FAIL'} ===")
    if not normal:
        for l in o.split("\n"):
            if "ERROR" in l.upper() or "unrecognized" in l or "invalid" in l:
                print("   ", l.strip()[:140])
        continue
    xp = WORK / f"a_{tag}.xplt"
    if not xp.exists():
        print("    无 xplt")
        continue
    dst = WORK / f"a_{tag}.h5"
    X.to_hdf5(str(xp), str(dst))
    with h5py.File(dst, "r") as h:
        found = []

        def walk(g, pre=""):
            for k in g:
                if isinstance(g[k], h5py.Group):
                    walk(g[k], pre + "/" + k)
                else:
                    found.append(pre + "/" + k)
        walk(h)
        for k in found:
            if "reaction" in k.lower():
                a = h[k][()]
                ax = np.atleast_2d(a)
                z = float(np.atleast_1d(ax[-1])[-1])
                print(f"    {k} 末行={ax[-1]}")
                print(f"    |F_z|≈{abs(z):.4f} N  =>  隐含面积 A = "
                      f"{abs(z)/(E*DL/L):.4f} mm²")
                break
