"""直接看 _parse_domain_section 对 beam 域返回了什么。"""
import inspect, textwrap, traceback
from pathlib import Path
import pyfebio.xplt as X

W = Path(r"D:\Project\climbing_fall_analysis\temp\pyfebio_demo\truss_probe")
raw = (W / "c_baseline.xplt").read_bytes()
print("xplt bytes:", len(raw))

# 手工遍历 block，找到 domain section 并打印其内容
buf = raw[4:]
i = 0
while i < len(buf) - 8:
    try:
        tag, offset = X.parse_prefix(buf[i:i + 8])
    except Exception as e:
        print("prefix 解析停:", e)
        break
    child = buf[i + 8:i + 8 + offset]
    nm = X.TAG_LUT[tag].name
    print(f"\n@ {i}: {nm} (offset {offset}, child {len(child)}B)")
    if nm == "PLT_DOMAIN_SECTION":
        doms = X._parse_domain_section(child)
        print(f"  域数 {len(doms)}")
        for d in doms:
            print("   keys:", list(d.keys()))
            print("   ", {k: (v if not hasattr(v, 'shape') else f'arr{v.shape}')
                          for k, v in d.items()})
    if nm == "PLT_ELEMENTSET_SECTION":
        try:
            es = X.parse_elementset_section(child)
            for e in es:
                print("  elementset keys:", list(e.keys()), e.get("name"))
        except Exception as ex:
            print("  ", ex)
    i += 8 + offset
