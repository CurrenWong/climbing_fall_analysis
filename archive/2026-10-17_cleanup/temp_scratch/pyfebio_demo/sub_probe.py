"""拿崩溃的 truss .feb 做减法，逐项定位崩溃触发点。"""
import re
import subprocess
from pathlib import Path

W = Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_lig")
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")
PROBE = W.parent.parent / "pyfebio_demo" / "sub_probe"
PROBE.mkdir(parents=True, exist_ok=True)

src = (W / "ankle_lig_WK.feb").read_text(encoding="latin-1")
print(f"源 .feb {len(src):,} 字符")


def cut(t, tag):
    """删掉顶层 <tag>...</tag> 段。"""
    i0 = t.find(f"  <{tag}")
    if i0 < 0:
        return t
    i1 = t.find(f"</{tag}>", i0)
    if i1 < 0:
        # 自闭合或无闭合：删到下一个顶层标签
        i1 = t.find("\n  <", i0 + 3)
    else:
        i1 += len(f"</{tag}>") + 1
    return t[:i0] + t[i1:]


VARIANTS = {
    "00_asis": src,
    "01_no_loads": cut(src, "Loads"),
    "02_no_contact": cut(src, "Contact"),
    "03_no_loads_contact": cut(cut(src, "Loads"), "Contact"),
    "04_no_truss": re.sub(r'\s*<Elements type="line2".*?</Elements>', "",
                          src, flags=re.S),
    "05_no_cart_surf": re.sub(r'\s*<Surface name="surf_.*?</Surface>', "",
                              src, flags=re.S),
}

for name, t in VARIANTS.items():
    f = PROBE / f"{name}.feb"
    f.write_text(t, encoding="latin-1")
    r = subprocess.run([str(FEBIO), "-i", str(f)], cwd=str(PROBE),
                       capture_output=True, text=True, timeout=600)
    o = (r.stdout or "") + (r.stderr or "")
    crashed = r.returncode == 3221225477
    read_ok = "SUCCESS" in o
    print(f"\n=== {name}: rc={r.returncode} "
          f"{'READ-OK' if read_ok else 'read-FAIL'} "
          f"{'**CRASH**' if crashed else ''} ===")
    if not read_ok and not crashed:
        for l in [x.strip() for x in o.split("\n") if x.strip()]:
            if "ERROR" in l.upper() or "invalid" in l.lower() \
               or "unrecognized" in l.lower() or "syntax" in l.lower():
                print("   ", l[:130])
    elif read_ok:
        neg = [l.strip() for l in o.split("\n") if "negative" in l.lower()]
        print("   ", neg[0][:90] if neg else "有解算输出")
