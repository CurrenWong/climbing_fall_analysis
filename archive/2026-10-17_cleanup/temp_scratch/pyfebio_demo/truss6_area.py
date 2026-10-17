"""打补丁读 truss 的支反力 -> 实测面积约定；并测 MeshData 面积映射。"""
import subprocess
import textwrap
import inspect
from pathlib import Path
import numpy as np
import h5py
import pyfebio.xplt as X

# ── 容错补丁：domain / surface / nodeset 缺 name 时用 id 兜底 ──
src = textwrap.dedent(inspect.getsource(X.parse_mesh))
src = src.replace('domain["name"]', 'domain.get("name", str(domain["id"][0]))')
src = src.replace('surface["name"]', 'surface.get("name", str(surface["id"][0]))')
ns = X.__dict__
exec(compile(src, "<patched_parse_mesh>", "exec"), ns)
X.parse_mesh = ns["parse_mesh"]
print("已打补丁 parse_mesh")

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, L, DL = 260.0, 10.0, 1.0

src_feb = (WORK / "et_iso.feb").read_text(encoding="latin-1")

CASES = {
    "baseline": src_feb,
    "elemdata": src_feb.replace("  <LoadData>", """  <MeshData>
    <ElementData name="area" elem_set="tr">
      <elem id="1">5.0</elem>
    </ElementData>
  </MeshData>
  <LoadData>"""),
}

for tag, t in CASES.items():
    f = WORK / f"b_{tag}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    normal = "NORMAL" in o.replace(" ", "")
    print(f"\n=== {tag} rc={r.returncode} {'OK' if normal else 'FAIL'} ===")
    if not normal:
        for l in o.split("\n"):
            if "ERROR" in l.upper() or "unrecognized" in l.lower() \
               or "invalid" in l.lower() or "must" in l.lower():
                print("   ", l.strip()[:150])
        continue
    xp = WORK / f"b_{tag}.xplt"
    dst = WORK / f"b_{tag}.h5"
    try:
        X.to_hdf5(str(xp), str(dst))
    except Exception as e:
        print("   to_hdf5 仍失败:", type(e).__name__, e)
        continue
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
                A_imp = abs(z) / (E * DL / L)
                print(f"    {k} shape={ax.shape}")
                print(f"    末行={ax[-1]}")
                print(f"    |F_z| = {abs(z):.4f} N  =>  隐含截面积 A = {A_imp:.4f} mm²")
                break
