"""收官 2：压 climbing_fall_analysis 段。"""
from pathlib import Path

P = Path(r"C:\Users\curre\.copaw\workspaces\default\MEMORY.md")
t = P.read_text(encoding="utf-8")
before = len(t.encode("utf-8"))
done = []


def sub(old, new, tag):
    global t
    assert t.count(old) == 1, f"[{tag}] 命中 {t.count(old)}"
    t = t.replace(old, new)
    done.append(tag)


sub("- **权威文档**：`阶段总结.md`（v1 结论 + 经验 §3.1–3.14）·\n"
    "  `docs/方案更新_v2.md`（结构性诊断 + P0–P5 路线）· `.learnings/`",
    "- **权威文档**：`阶段总结.md`（v1 结论 + 经验）· `docs/方案更新_v2.md`（P0–P5 路线）· `.learnings/`",
    "docs")

sub("- **核心诊断**：模型**无转动自由度**，结构上算不出**踝扭伤**（占 26–28%，\n"
    "  两个独立人群一致）——「跑得对」和「答得上」是两件事",
    "- **核心诊断**：模型**无转动自由度**，结构上算不出**踝扭伤**（占 26–28%）"
    "——「跑得对」≠「答得上」", "diag")

sub("- **P1 的硬结果**：真实抱石场景下**一例骨折都报不出来**（1000 次 MC，最高利用率 0.725）\n"
    "  ⇒ 损伤部位分布是**空的**。覆盖度 100% → 无旋转 32.6% → 结论可信 13.8% → 能产出的损伤类型 19.3%",
    "- **P1 硬结果**：真实抱石场景**一例骨折都报不出来**（1000 次 MC，最高利用率 0.725）"
    "⇒ 损伤部位分布是**空的**。覆盖度 100% → 无旋转 32.6% → 结论可信 13.8% → 可产出损伤类型 19.3%",
    "p1")

sub("- **可复现结果**（`DRIVER=disp DISP_MM=0.1472`）：`rc=0`/`NORMAL TERMINATION`/16.6 s，\n"
    "  支反力 **436.8 N**，软骨 von Mises p99 **2.58 / 4.27 MPa** → 外推 600 N 胫侧 ≈3.4 MPa，\n"
    "  **与吴恺的 3.74 MPa 吻合**",
    "- **可复现结果**（`DRIVER=disp DISP_MM=0.1472`）：`rc=0`/`NORMAL TERMINATION`/16.6 s，"
    "支反力 **436.8 N**，软骨 p99 **2.58 / 4.27 MPa** → 外推 600 N 胫侧 ≈3.4 MPa，"
    "**与吴恺 3.74 MPa 吻合**", "res")

P.write_text(t, encoding="utf-8")
after = len(t.encode("utf-8"))
print(f"已应用 {len(done)} 处: {done}")
print(f"MEMORY.md {before:,} -> {after:,} bytes "
      f"({'OK, 余量 ' + str(10240 - after) if after <= 10240 else 'FAIL, 差 ' + str(after-10240)})")
