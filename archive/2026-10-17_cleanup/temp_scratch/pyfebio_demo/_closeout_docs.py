"""压缩上下文前的文档收口：标注实跑验证 + 给旧文档加 superseded 指路 + 交接文档互指。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis")
done = []


def patch(path, subs):
    t = path.read_text(encoding="utf-8")
    b = len(t.encode("utf-8"))
    for old, new, tag in subs:
        assert t.count(old) == 1, f"[{path.name}][{tag}] 命中 {t.count(old)}"
        t = t.replace(old, new)
        done.append(tag)
    path.write_text(t, encoding="utf-8")
    print(f"  {path.name} {b:,} -> {len(t.encode('utf-8')):,} bytes")


# ① FE_PIPELINE_TODO：标注入口命令已实跑验证
patch(P / "FE_PIPELINE_TODO.md", [(
    "> ⚠️ cmd 的 `set X=Y &&` **行尾不能带空格**（否则 40013，见 MEMORY gotcha #5）。",
    "> ✅ **本节命令已实跑验证（2026-10-04）**：4 步全部 rc=0，FEBio `rc=0`/17.3 s，\n"
    "> 读回支反力 **436.817 N**、软骨 p99 **2.5799 MPa** —— 与本文件各处记录**逐字一致**。\n"
    "> ⚠️ cmd 的 `set X=Y &&` **行尾不能带空格**（否则 40013，见 MEMORY gotcha #5）。",
    "verified")])

# ② ANKLE_FE_TODO：标为已取代，指路到 THUMS 线
patch(P / "ANKLE_FE_TODO.md", [(
    "# ANKLE_FE — 按吴恺2012 论文建踝关节有限元模型\n",
    "# ANKLE_FE — 按吴恺2012 论文建踝关节有限元模型\n\n"
    "> ## ⚠️ 本文件是**历史文档，IITD 路线已停止**（2026-10-03）\n"
    "> **接手请看 `FE_PIPELINE_TODO.md` 的「🧭 接手入口」节** —— 几何源已换成\n"
    "> **THUMS V7.1**（`paper/AM50_V71_Occupant.zip` 内 219 MB `.k`），\n"
    "> **FEBio 求解已收敛并通过吴恺量级验证**（436.8 N / 软骨 p99 2.58–4.27 MPa）。\n"
    "> 本文件记录的是被取代的 IITD 路线（IITD 软骨 STL 是缺陷几何、\n"
    "> **无独立腓骨**、其 OSIM 在原生 OpenSim 4.6 **加载失败**），仅作对照留档。\n"
    "> 完整历史见 `docs/踝关节FE建模交接.md`。\n",
    "superseded")])

# ③ 交接文档 §13 顶部指回入口节
patch(P / "docs/踝关节FE建模交接.md", [(
    "## 13. 材料 / 软骨 / FEBio 管线（2026-10-03，见 `FE_PIPELINE_TODO.md`）",
    "## 13. 材料 / 软骨 / FEBio 管线（2026-10-03，见 `FE_PIPELINE_TODO.md`）\n\n"
    "> **压缩上下文后接手 → 先读 `FE_PIPELINE_TODO.md` 的「🧭 接手入口」节**\n"
    "> （状态 / 一键复现 / 下一步排序 / 负面结果清单 / 环境开关 / 别再踩的坑）。\n"
    "> 本节 §13.1–13.7 是该入口的展开细节。",
    "crossref")])

print(f"\n共 {len(done)} 处: {done}")
