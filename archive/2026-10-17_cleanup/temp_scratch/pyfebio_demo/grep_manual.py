"""从本地 FEBio User Manual PDF 提取与 Surface / surface_load 相关的说明。"""
import re
from pathlib import Path

PDF = Path(r"D:\Program\FEBioStudio\doc\FEBio_User_Manual.pdf")
print(f"存在: {PDF.exists()}  {PDF.stat().st_size/1024/1024:.1f} MB")

text = None
for mod in ("pypdf", "PyPDF2", "pdfplumber", "fitz"):
    try:
        if mod == "fitz":
            import fitz
            d = fitz.open(str(PDF))
            text = "\n".join(p.get_text() for p in d)
        elif mod == "pdfplumber":
            import pdfplumber
            with pdfplumber.open(str(PDF)) as pdf:
                text = "\n".join((p.extract_text() or "") for p in pdf.pages)
        else:
            m = __import__(mod)
            R = m.PdfReader(str(PDF))
            text = "\n".join((p.extract_text() or "") for p in R.pages)
        print(f"用 {mod} 提取成功，{len(text):,} 字符")
        break
    except Exception as exc:
        print(f"  {mod} 失败: {type(exc).__name__}")

if text is None:
    raise SystemExit("没有可用的 PDF 库")

Path(r"D:\Project\climbing_fall_analysis\temp\thums\febio_manual.txt").write_text(
    text, encoding="utf-8")

for kw in ("<Surface", "surface_load", "<Elements", "SurfacePair",
           "surface name", "<NodeSet"):
    print(f"\n===== 含 {kw!r} 的上下文 =====")
    n = 0
    for m in re.finditer(re.escape(kw), text):
        s = max(0, m.start() - 260)
        snippet = text[s:m.start() + 320].replace("\n", " | ")
        print(f"  ...{snippet}...")
        n += 1
        if n >= 3:
            break
    if n == 0:
        print("  （未找到）")
