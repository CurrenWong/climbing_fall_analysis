import textwrap, inspect, traceback
from pathlib import Path
import pyfebio.xplt as X

src = textwrap.dedent(inspect.getsource(X.parse_mesh))
src = src.replace('domain["name"]', 'domain.get("name", str(domain["id"][0]))')
src = src.replace('surface["name"]', 'surface.get("name", str(surface["id"][0]))')
ns = X.__dict__
exec(compile(src, "<p>", "exec"), ns)
X.parse_mesh = ns["parse_mesh"]

W = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
try:
    X.to_hdf5(str(W / "b_baseline.xplt"), str(W / "b_baseline.h5"))
    print("OK")
except Exception:
    tb = traceback.format_exc()
    print(tb)
    print("\n=== parse_state 里的 ['name'] 用法 ===")
    for name in ("parse_state", "_parse_domain_section", "parse_mesh"):
        fn = getattr(X, name, None)
        if fn is None:
            continue
        s = textwrap.dedent(inspect.getsource(fn))
        for i, l in enumerate(s.split("\n"), 1):
            if '"name"' in l or "[" + "'" + "name" + "'" + "]" in l:
                print(f"  {name}:{i}: {l.strip()[:120]}")
