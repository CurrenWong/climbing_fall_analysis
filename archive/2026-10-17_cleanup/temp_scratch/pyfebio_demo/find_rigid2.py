"""查：刚体材料声明 + MeshDomains 里刚体怎么写 + 完整最小例子。"""
import re

import pdfplumber

PDF = r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf"
with pdfplumber.open(PDF) as pdf:
    pages = [(i + 1, pg.extract_text() or "") for i, pg in enumerate(pdf.pages)]

print("=== 含 'rigid body' 且含 'Material' 的页 ===")
for p, t in pages:
    if "rigid body" in t.lower() and ("material" in t.lower()):
        print(f"  第 {p} 页")

print("\n=== MeshDomains 段里 rigid 的写法（搜 <rigid_domain / rigid_domain / domain=）===")
for k in ("rigid_domain", "<SolidDomain", "domain_type", "MeshDomains", "rigid>",
          "type=\"rigid\""):
    hit = [p for p, t in pages if k in t]
    print(f"  {k!r:16} -> 页 {hit[:14]}")

print("\n=== 第 73 页（材料定义附近）===")
for p, t in pages:
    if p == 73:
        print(t[:3000])
