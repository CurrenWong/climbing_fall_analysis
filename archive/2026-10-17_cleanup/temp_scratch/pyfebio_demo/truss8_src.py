import textwrap, inspect
import pyfebio.xplt as X

src = textwrap.dedent(inspect.getsource(X.parse_mesh))
src = src.replace('domain["name"]', 'domain.get("name", str(domain["id"][0]))')
src = src.replace('surface["name"]', 'surface.get("name", str(surface["id"][0]))')
for i, l in enumerate(src.split("\n"), 1):
    mark = "  <<<<" if '"name"' in l or "'name'" in l else ""
    print(f"{i:3d}| {l}{mark}")
