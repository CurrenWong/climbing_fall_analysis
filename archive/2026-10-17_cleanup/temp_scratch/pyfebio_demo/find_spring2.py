"""在 User Manual 里搜 'spring' 并打印上下文（找弹簧边界条件语法）。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
PAT = re.compile(r"spring", re.I)

found = []
with pdfplumber.open(PDF) as pdf:
    print(f"[PDF {len(pdf.pages)} 页]")
    for i, pg in enumerate(pdf.pages):
        t = pg.extract_text() or ""
        if PAT.search(t):
            # 只留带 "bc" / "boundary" / "discrete" 的页（更可能是语法页）
            low = t.lower()
            tag = []
            if 'type="spring"' in t:
                tag.append('TYPE=SPRING')
            if "spring_stiffness" in t:
                tag.append("STIFFNESS")
            if "<bc" in low and "spring" in low:
                tag.append("BC+SPRING")
            if "discrete" in low:
                tag.append("DISCRETE")
            found.append((i + 1, ",".join(tag) or "-", t))

print(f"含 spring 的页 {len(found)}: {[(p, t) for p, t, _ in found]}")
for p, tag, t in found:
    if "TYPE=SPRING" in tag or "STIFFNESS" in tag or ("BC+SPRING" in tag):
        print(f"\n{'='*74}\n### 第 {p} 页  [{tag}]\n{'='*74}")
        # 只打印含 spring 的段落附近
        for m in PAT.finditer(t):
            a, b = max(0, m.start() - 700), min(len(t), m.end() + 700)
            print("..." + t[a:b].replace("\n", " | ") + "...")
            print("-" * 60)
        break
