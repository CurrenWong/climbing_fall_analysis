"""从 Anderson 2006（Intra-articular...）md 里抽 §2.2 原文语境。"""
import re
from pathlib import Path

D = Path(r"D:\Project\climbing_fall_analysis\paper")
F = D / "MinerU_markdown_Intra-articular_Contact_Stress_Distributions_at_th_2106695331008294912.md"
t = F.read_text(encoding="utf-8", errors="replace")
# 去掉 markdown 里的图片/公式噪声，保留文字
t = re.sub(r"\n{3,}", "\n\n", t)

for key in ("weakly restrained", "500N/mm", "linear spring", "rigid",
            "inversion-eversion", "seating", "pre-loading", "internal-external"):
    for m in re.finditer(re.escape(key), t, re.I):
        a, b = max(0, m.start() - 1100), min(len(t), m.end() + 1100)
        print(f"\n{'#'*76}\n### 关键词 {key!r}\n{'#'*76}")
        print(t[a:b].replace("\n", " "))
