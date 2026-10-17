import re, numpy as np
FEB = r"temp/thums/febio_lig/ankle_lig_WK.feb"
txt = open(FEB, errors="replace").read()
for nm in ["tibia_top_ns", "fibula_top_ns", "talus_low_ns", "calcaneus_bottom_ns", "navicular_all_ns"]:
    m = re.search(r'<NodeSet name="%s">(.*?)</NodeSet>' % re.escape(nm), txt, re.S)
    body = m.group(1)
    toks = [x for x in re.split(r"[,\s]+", body.strip()) if x]
    vals = [int(x) for x in toks]
    big = [v for v in vals if v > 4029]
    print("%-22s len(body)=%6d  ntok=%4d  nuniq=%4d  n>4029=%4d  start=%s  end=%s" % (
        nm, len(body), len(toks), len(set(vals)), len(big), toks[:3], toks[-3:]))
    if big:
        print("      sample >4029:", big[:12])
print()
# what does the raw text look like around tibia_top_ns
i = txt.find('<NodeSet name="tibia_top_ns">')
print("--- raw 400 chars ---")
print(repr(txt[i:i + 400]))
print("--- around its close ---")
j = txt.find("</NodeSet>", i)
print(repr(txt[j - 200:j + 200]))
