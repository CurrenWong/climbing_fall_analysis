"""用 FEBio 的 <logfile><node_data> 写 ASCII，实测 truss 的截面积约定。"""
import subprocess
from pathlib import Path

WORK = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
E, L, DL = 260.0, 10.0, 1.0

base = (WORK / "et_iso.feb").read_text(encoding="latin-1")
# 把 <Output> 换成带 logfile 的版本
i0 = base.index("  <Output>")
i1 = base.index("</Output>") + len("</Output>")
OUT = """  <Output>
    <plotfile type="febio"><var type="displacement"/><var type="reaction forces"/></plotfile>
    <logfile>
      <node_data data="reaction forces" file="react.txt" node_set="pull"/>
      <node_data data="displacement" file="disp.txt" node_set="pull"/>
    </logfile>
  </Output>"""

CASES = {
    "L10_E260": (10.0, 260.0),
    "L20_E260": (20.0, 260.0),
    "L10_E520": (10.0, 520.0),
}
for tag, (Lc, Ec) in CASES.items():
    t = base[:i0] + OUT + base[i1:]
    t = t.replace('<node id="2">0,0,-10</node>', f'<node id="2">0,0,{-Lc}</node>')
    t = t.replace("<E>260</E>", f"<E>{Ec}</E>")
    f = WORK / f"d_{tag}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(WORK),
                       capture_output=True, text=True, timeout=300)
    o = (r.stdout or "") + (r.stderr or "")
    print(f"\n=== {tag} rc={r.returncode} "
          f"{'OK' if 'NORMAL' in o.replace(' ','') else 'FAIL'} ===")
    for fn in ("react.txt", "disp.txt"):
        p = WORK / fn
        if p.exists():
            txt = p.read_text(errors="ignore").strip()
            last = txt.split("\n")[-1] if txt else "(空)"
            print(f"   {fn}: {last[:120]}")
        else:
            print(f"   {fn}: 不存在")
    if "NORMAL" not in o.replace(" ", ""):
        for l in o.split("\n"):
            if "ERROR" in l.upper() or "unrecognized" in l.lower() \
               or "invalid" in l.lower():
                print("   ", l.strip()[:140])
