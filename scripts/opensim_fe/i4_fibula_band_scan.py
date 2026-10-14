# -*- coding: utf-8 -*-
"""I4 · 把 `share = 0.2241`（M1 定案值）推入可解区 —— **只扫支撑带**。

**为什么只扫支撑带**：③「加密 load step」已被实测证伪 ——
`fracture_matrix.py:936`：「`time_steps=1`（模块默认）——实测步数 >1 会让
FEBio 对本 NodalForce 设置误加载而发散，故不采用」。
本 BC 用端部 NodalForce 加载 ⇒ 步数不可调。

**带的设计依据**：F 路线实测「侧向覆盖 ≥80% 自由长 ⇒ 夹持 DOF 无关」。
`long95=[(0.31,0.69)]` 只约束 38% 全长 ⇒ 覆盖自由长 (62%) 的 **61%**，未达 80%。
⇒ 追加 `wide50=[(0.25,0.75)]`（约束 50% 全长 ⇒ 覆盖自由长 **100%**）与
`wide70=[(0.27,0.73)]`（约束 46% ⇒ 覆盖 **85%**）。

**不修改 G9 的预注册常量**：`BAND_ALL` / `SHARES` 均在运行时注入，
产物走 `--out-suffix`，既有预注册产物不动。
"""
from __future__ import annotations

import argparse
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "opensim_fe"))
sys.path.insert(0, str(ROOT / "src"))

import g9_fibula_fe_support as g9  # noqa: E402

# 带宽 → 约束全长比例 / 覆盖自由长比例（用于报告）
BANDS = [
    ("long95", [(0.31, 0.69)], "G9 口径（基线，0.18 通 / 0.20 不通）"),
    ("long85", [(0.33, 0.67)], "G9 候选，略加约束"),
    ("wide70", [(0.27, 0.73)], "覆盖自由长 85%（≥F 路线 80% 门槛）"),
    ("wide50", [(0.25, 0.75)], "覆盖自由长 100%（对齐 G9 主设定 three 的外包络）"),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--share", default="0.2241")
    ap.add_argument("--suffix", default="_bandscan")
    a = ap.parse_args()
    s = float(a.share)

    g9.SHARES = [("share_0.100", 0.10), (f"share_{s:.4f}", s)]
    g9.BAND_ALL = dict(g9.BAND_ALL)
    for name, band, _ in BANDS:
        g9.BAND_ALL[name] = band

    print("SHARES =", g9.SHARES)
    for name, band, why in BANDS:
        c = sum(hi - lo for lo, hi in band)
        free = 1.0 - c
        print(f"  {name:<8} {band}  约束 {c:.0%} 全长  ⇒ 覆盖自由长 "
              f"{c/free:.0%}   {why}")

    rc = g9.main(["--bands", ",".join(n for n, _, _ in BANDS),
                  "--primary-band", BANDS[-1][0],
                  "--out-suffix", a.suffix])
    print(f"\ng9.main rc={rc}")
    return rc


if __name__ == "__main__":
    raise SystemExit(main())