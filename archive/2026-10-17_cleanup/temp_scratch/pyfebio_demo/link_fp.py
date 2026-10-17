"""找 FE_PIPELINE_TODO.md 里的「接手入口」并插入论文重解析指针。"""
from pathlib import Path

FP = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
f = FP.read_text(encoding="utf-8")
L = f.split("\n")

hits = [i for i, l in enumerate(L) if "接手入口" in l or "\U0001f9ed" in l]
print("命中行：", [i + 1 for i in hits])
for i in hits:
    print(f"  {i+1:4d}| {L[i][:120]}")

if hits:
    k = hits[0]
    INS = ("> 🆕 **2026-10-04 论文重解析（MinerU）**：新增可算靶子 **接触面积 335.5 mm²**"
           "（我的模型 250.9~456.9 mm² 覆盖它，已对上）；§3 讨论给出完整机理"
           "**外踝长且靠后 ⇒ 踝轴斜 ⇒ 负重时小腿外翻 + 内旋 0°~0.5° ⇒ "
           "外侧(跟腓)<内侧、胫跟>胫舟**；**韧带位移判定不可复现**（差 30~40 倍）；"
           "**新阻塞：模型到不了 600 N**（薄软骨单元挤反，det(J) min 0.0105）"
           "→ `docs/吴恺论文_重解析关键信息.md`")
    L.insert(k, INS)
    L.insert(k + 1, "")
    FP.write_text("\n".join(L), encoding="utf-8")
    print(f"\n已插入 -> {FP.stat().st_size:,} bytes")
else:
    print("未找到，改为插在「## 一句话」节后")
    j = next(i for i, l in enumerate(L) if l.startswith("## 三件必须知道的事"))
    L.insert(j, "")
    L.insert(j, "> 🆕 **2026-10-04 论文重解析**：新靶子接触面积 **335.5 mm²**（已对上）、"
                "§3 给出**外翻+内旋**机理、韧带位移**判定不可复现**、"
                "新阻塞**到不了 600 N** → `docs/吴恺论文_重解析关键信息.md`")
    FP.write_text("\n".join(L), encoding="utf-8")
    print(f"已插入 -> {FP.stat().st_size:,} bytes")
