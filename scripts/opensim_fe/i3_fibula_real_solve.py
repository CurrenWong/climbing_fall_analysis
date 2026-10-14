# -*- coding: utf-8 -*-
"""I″ · `share = 0.2241` 的**真解** FE（`long95` 支撑带）—— I′ 的定案步骤。

**不修改 G9 的预注册常量**：`g9_fibula_fe_support.py::SHARES` 原样保留在文件里，
本驱动在运行时**只追加**一档 `(0.100, 0.2241)`，并**重跑 0.10** 作为
「同条件可复现性」检查（I′ 的正对照已用它差 0.6%）。

**为什么只放 0.10 与 0.2241 两档**：G9 的线性自检是相对最低档取 `ratio_median`，
两档即直接给出 **`0.2241 / 0.10 = 2.241×`** 这一关键量 —— 检验「σ(h) 在已验证的 2× 之外是否仍线性」。

用法：`.venv\\Scripts\\python.exe scripts\\opensim_fe\\i3_fibula_real_solve.py`
"""
from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import g9_fibula_fe_support as g9  # noqa: E402

# ★ 纯增量：保留 G9 原两档 + 追加待测档（不是替换）
BASE = [("share_0.100", 0.10)]


def key_of(s: float) -> str:
    return f"share_{s:.4f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shares", default="0.2241",
                    help="逗号分隔的待测 share（均追加到 0.10 基准之后）")
    ap.add_argument("--bands", default="long95")
    ap.add_argument("--suffix", default="_share02241")
    a = ap.parse_args()

    shares = [float(x) for x in a.shares.split(",") if x.strip()]
    g9.SHARES = BASE + [(key_of(s), s) for s in shares]

    print("原 G9.SHARES =", [("share_0.050", 0.05), ("share_0.100", 0.10)])
    print("本件 SHARES  =", g9.SHARES, "（0.10 保留作基准/复现检查）")
    print("BAND_LONG95  =", g9.BAND_LONG95)
    print(f"LIN_TOL = {g9.LIN_TOL}；FIBULA_WINDOW = {g9.FIBULA_WINDOW}\n")
    rc = g9.main([
        "--bands", a.bands,
        "--primary-band", a.bands.split(",")[0].strip(),
        "--out-suffix", a.suffix,
    ])
    print(f"\ng9.main rc={rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())