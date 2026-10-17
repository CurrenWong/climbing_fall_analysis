import io, json
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
t = io.open(P / "docs/踝关节FE建模交接.md", encoding="utf-8").read()
need = ["## 0.", "## 2.", "### 2.3", "## 3.", "## 4.", "## 5.", "## 6.", "## 7.",
        "## 8.", "## 9.", "## 10.", "## 11.", "## 12.",
        "### 12.5", "### 12.6", "12.5 ✅", "12.6 决策"]
miss = [n for n in need if n not in t]
print("交接文档 缺失:", miss or "无 ✅", f"({len(t.encode('utf-8'))} bytes)")

lrn = io.open(P / ".learnings/LEARNINGS.md", encoding="utf-8").read()
print("LEARNINGS 含 LRN-014:", "LRN-20261003-014" in lrn, f"({len(lrn.encode('utf-8'))} bytes)")

for f in ("temp/thums/thums_lowerlimb.npz", "temp/thums/summary.json",
          "temp/thums/verify.json", "scripts/ankle_fe/thums_extract.py",
          "scripts/ankle_fe/thums_verify.py", "THUMS_EXTRACT_TODO.md"):
    p = P / f
    print(f"  {'✅' if p.exists() else '❌'} {p.stat().st_size if p.exists() else 0:>10,}  {f}")

s = json.loads((P / "temp/thums/summary.json").read_text(encoding="utf-8"))
print(f"\nsummary: 总节点 {s['n_nodes_total']:,}, 保留 {s['n_nodes_kept']:,}, "
      f"部件 {len(s['parts'])}, 解析失败 {s['parse_fails']}")
