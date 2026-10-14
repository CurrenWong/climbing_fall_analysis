# -*- coding: utf-8 -*-
"""J · 冻结 BC2 数据保全（合规形式）：定位冻结期代码提交 + 生成 sha256 清单。

背景：`temp/opensim_fe/fracture_matrix_bc2/`（51.6 MB / 83 文件）是 4 个沿用部位
（calcaneus_r / R_HIPBONE / L3 / C5）的**唯一数据副本**，且**已实测不可复现**（差 6–22%）。
但 `.gitignore:72` 明确把 `*.feb *.xplt *.hdf5 *.npz` 与 `temp/` 归为
「regenerable simulation artefacts (safe to delete)」，仓库**对外发布**只收
code + docs + 可再生 results ⇒ **不得入库**（会违反发布政策）。

⇒ 本件做的是把「静默丢失」变成「可检测丢失 + 可恢复」：
  1. 定位冻结期代码提交（可 `git checkout` 回退 ⇒ 恢复可再生性）
  2. 生成逐文件 sha256 清单（任何副本都能验证逐字节一致）
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "opensim_fe"
SRC = ROOT / "temp" / "opensim_fe" / "fracture_matrix_bc2"

# 与 A2「沿用冻结的 4 部位」判定直接相关的模块
MODULES = [
    "scripts/opensim_fe/bone_feb.py",
    "scripts/opensim_fe/fracture_matrix.py",
    "scripts/opensim_fe/rerun_matrix_bc2.py",
    "scripts/opensim_fe/paper_alignment.py",
    "scripts/opensim_fe/s4_alignment_bc2.py",
]


def git(*a: str) -> str:
    return subprocess.run(["git", *a], cwd=ROOT, capture_output=True,
                          text=True, encoding="utf-8", errors="replace").stdout.strip()


def sha256(p: pathlib.Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def main() -> int:
    # ---- 1) 冻结期代码提交 ---------------------------------------------------
    hist: dict[str, list[dict]] = {}
    for m in MODULES:
        rows = []
        for ln in git("log", "--format=%H|%ad|%s", "--date=short", "--", m).splitlines():
            h, d, s = ln.split("|", 2)
            rows.append({"commit": h[:12], "date": d, "subject": s})
        hist[m] = rows[:6]
    # 冻结登记日 = 2026-10-09（`FROZEN_BC2` 记忆）；取该日及之前的最后一次改动
    FREEZE = "2026-10-09"
    at_freeze, after = {}, []
    for m, rows in hist.items():
        pick = None
        for r in rows:
            if r["date"] <= FREEZE:
                pick = r
                break
        at_freeze[m] = pick or (rows[-1] if rows else None)
        if pick:
            for r in rows:
                if r["date"] > FREEZE:
                    after.append({**r, "module": m, "date_vs_freeze": r["date"]})

    # ---- 2) 逐文件 sha256 清单 ------------------------------------------------
    files = sorted(p for p in SRC.rglob("*") if p.is_file())
    man = []
    for p in files:
        st = p.stat()
        man.append({
            "path": str(p.relative_to(ROOT)).replace("\\", "/"),
            "size_b": st.st_size,
            "mtime": dt.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
            "sha256": sha256(p),
        })
    total = sum(m["size_b"] for m in man)

    doc = {
        "what": "冻结 BC2 逐骨 FE（FEBio 输入 + 线性自检 + xplt/hdf5 结果 + 日志）",
        "why_frozen_not_in_git": (
            "`.gitignore:72` 明确把 *.feb/*.xplt/*.hdf5/*.npz 与 temp/ 归为"
            "「Regenerable simulation artefacts (safe to delete — reproduced by scripts)」，"
            "仓库对外发布只收 code + docs + 可再生 results ⇒ 入库会违反该政策。"
            "但本目录**已实测不可复现**（同条件重跑差 6–22%）。"
        ),
        "why_it_matters": (
            "A2 有 4 个部位（calcaneus_r / R_HIPBONE / L3 / C5）的判定**直接沿用本目录的 h5**；"
            "本目录是它们的唯一副本。"
        ),
        "src_dir_local_only": str(SRC.relative_to(ROOT)).replace("\\", "/"),
        "total_files": len(man),
        "total_size_b": total,
        "total_size_mb": round(total / 1e6, 2),
        "freeze_registered_date": FREEZE,
        "code_at_freeze": at_freeze,
        "code_changed_after_freeze": after,
        "files": man,
    }
    outj = RES / "FROZEN_BC2_MANIFEST.json"
    outj.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")

    # ---- 3) 人读版 ------------------------------------------------------------
    L = [
        "# 冻结 BC2 数据保全 —— 清单与恢复路径",
        "",
        f"> **日期**：2026-10-13 ｜ **机器可读**：`FROZEN_BC2_MANIFEST.json`",
        f"> **本地路径（不入库）**：`{doc['src_dir_local_only']}/`",
        f"> **体量**：{doc['total_files']} 文件 / **{doc['total_size_mb']} MB**",
        "",
        "---",
        "",
        "## 1. 为什么它不在 git 里（**政策 vs 实测的冲突**）",
        "",
        f"> `{doc['src_dir_local_only']}/` 是 4 个沿用部位判定的唯一数据副本。",
        "> 但 `.gitignore:72` 明确把 `*.feb`/`*.xplt`/`*.hdf5`/`*.npz` 与 `temp/` 归为",
        "> 「**Regenerable simulation artefacts (safe to delete — reproduced by scripts)**」，",
        "> 理由是仓库**对外发布**、只收 code + docs + **可再生** results。",
        "",
        "> **⇒ 不可入库（会违反该发布政策）；但「可再生」这个前提已被实测证伪**（同条件重跑差 **6–22%**）。",
        "> **⇒ 本件的作用：把「静默丢失」变成「可检测丢失 + 可恢复」。**",
        "",
        "## 2. 为什么会不可复现（**代码漂移**）",
        "",
        f"冻结登记日 **{FREEZE}**；此后以下模块被改动 ⇒ 用当前代码重跑得到的是**另一套 FE**：",
        "",
        "| 模块 | 冻结期最后改动 |",
        "|---|---|",
    ]
    for m, r in at_freeze.items():
        L.append(f"| `{m}` | {r['commit']} ({r['date']}) |" if r else f"| `{m}` | — |")
    if after:
        L += ["", "**冻结后被改动的提交**（恢复时须回退这些）：", ""]
        L += [f"- `{a['module']}` `{a['commit']}` ({a['date']}) — {a['subject']}"
              for a in after]
    L += [
        "",
        "## 3. 恢复路径（**任一条即可**）",
        "",
        "1. **回到冻结期代码**：上表提交哈希 ⇒ `git checkout <commit> -- <modules>` ⇒ 重跑",
        "2. **逐字节校验副本**：本清单含**逐文件 sha256** ⇒ 任何外部副本都能验证一致",
        f"3. **降级使用**：若已丢失 ⇒ A2 那 4 个部位**退回「不可复现数据」**状态，须登记",
        "",
        "## 4. 建议（**待拍板**）",
        "",
        f"- 把 `{doc['src_dir_local_only']}/` 另存到**仓库外**（或私有仓库），本清单入 git 作校验基准",
        "- 恢复后，**重跑并与本清单逐文件比对**，差异本身是有价值的证据（可定位是哪一段代码漂移造成）",
        "",
        "## 5. 溯源",
        "",
        "| 项 | 来源 |",
        "|---|---|",
        "| 不可复现实测 | `FE_BC_PROBE.md`（`V0_control` 未复现冻结值，差 6–22%） |",
        "| 代码漂移核查 | `git log` on 冻结期相关模块 |",
        "| 「沿用 4 骨」的判据 | `A2_RECOMPUTE_SCOPE.md` |",
        "| A2 新计数 | `A2_RECOMPUTE_RESULT.md` |",
        "",
    ]
    outm = RES / "FROZEN_BC2_MANIFEST.md"
    outm.write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"wrote {outj.name} ({outj.stat().st_size} B)")
    print(f"wrote {outm.name} ({outm.stat().st_size} B)")
    print(f"\n{len(man)} files, {total/1e6:.2f} MB hashed")
    print("\ncode at freeze:")
    for m, r in at_freeze.items():
        print(f"  {m}: {r['commit'] if r else '-'} ({r['date'] if r else '-'})")
    print(f"\ncommits AFTER freeze: {len(after)}")
    for a in after:
        print(f"  {a['date']} {a['commit']} {a['module']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())