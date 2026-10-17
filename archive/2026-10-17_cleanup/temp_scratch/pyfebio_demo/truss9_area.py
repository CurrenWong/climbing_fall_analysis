"""补全 pyfebio.parse_mesh 的 name 兜底，读 truss 支反力 -> 实测截面积。"""
import subprocess, textwrap, inspect
from pathlib import Path
import numpy as np, h5py
import pyfebio.xplt as X

src = textwrap.dedent(inspect.getsource(X.parse_mesh))
for old, new in (
    ("domain.get(\"name\", str(domain[\"id\"][0]))",
     "domain.get(\"name\", str(domain[\"id\"][0]))"),
    ("{domain['name']}", "{domain.get('name', str(domain['id'][0]))}"),
    ("{surface['name']}", "{surface.get('name', str(surface['id'][0]))}"),
    ('if nodeset["name"] == "":', 'if nodeset.get("name", "") == "":'),
    ("{nodeset['name']}", "{nodeset.get('name', str(nodeset['id'][0]))}"),
    ("{elementset['name']}", "{elementset.get('name', str(elementset['id'][0]))}"),
):
    src = src.replace(old, new)
assert "domain['name']" not in src and "surface['name']" not in src
ns = X.__dict__
exec(compile(src, "<patched>", "exec"), ns)
X.parse_mesh = ns["parse_mesh"]
print("patch applied")

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, L, DL = 260.0, 10.0, 1.0
src_feb = (WORK / "et_iso.feb").read_text(encoding="latin-1")

CASES = {"baseline": src_feb}
# 换个长度，验证 F ∝ 1/L（同时排除"面积恰好=1"的巧合）
CASES["L20"] = src_feb.replace('<node id="2">0,0,-10</node>',
                               '<node id="2">0,0,-20</node>')
CASES["E520"] = src_feb.replace("<E>260</E>", "<E>520</E>")

for tag, t in CASES.items():
    f = WORK / f"c_{tag}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    if "NORMAL" not in o.replace(" ", ""):
        print(f"{tag}: FAILED")
        continue
    xp, dst = WORK / f"c_{tag}.xplt", WORK / f"c_{tag}.h5"
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
                ax = np.atleast_2d(h[k][()])
                z = float(np.atleast_1d(ax[-1])[-1])
                Lc = 20.0 if tag == "L20" else 10.0
                Ec = 520.0 if tag == "E520" else 260.0
                A = abs(z) * Lc / (Ec * DL)
                print(f"  {tag:9s} |F|={abs(z):9.4f} N  ->  A = {A:.4f} mm²")
                break
