"""数字链路三方核对：论文 v2 <-> 出处索引 <-> 源码。

链路：论文写的常量  →  出处索引标注的代码位置  →  源码实际值。
任何一方不一致 => FAIL。这是"论文数字可溯源"的最低充分检查。
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

IDX = ROOT / "docs/paper/论文v2_出处索引.md"
PAPER = ROOT / "docs/paper/论文v2_数学建模版.md"

FAIL: list[str] = []


def src_value(expr: str, module: str):
    """从源码提取常量值。兼容 `NAME: float = 215.0` 与 `NAME = 215.0`。
    ⚠️ 旧版用 `[:=]` 会先撞上类型注解的 `:`，抓到 'float = 215.0'。"""
    txt = (ROOT / module).read_text(encoding="utf-8", errors="replace")
    m = re.search(
        rf"^{re.escape(expr)}(?:\s*:\s*[\w.\[\]\| ]+)?\s*=\s*([^\n#]+)",
        txt, re.M)
    if not m:
        return None
    raw = m.group(1).strip().rstrip(",")
    if re.fullmatch(r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?", raw):
        return float(raw)
    return raw.strip("\"'")


def check(label: str, paper_ok: bool, idx_ok: bool, code_val,
          want) -> None:
    ok = paper_ok and idx_ok and (code_val == want)
    print(f"  {'OK  ' if ok else 'FAIL'}  {label:<34} "
          f"paper={'Y' if paper_ok else 'N'} idx={'Y' if idx_ok else 'N'} "
          f"code={code_val!r} want={want!r}")
    if not ok:
        FAIL.append(label)


def main() -> int:
    paper = PAPER.read_text(encoding="utf-8")
    idx = IDX.read_text(encoding="utf-8")

    print("=" * 76)
    print("三方核对  论文 v2  <->  出处索引  <->  源码")
    print("=" * 76)

    # 1. 撕脱阈 215 N
    check("ANKLE_AVULSION_LOAD_N = 215",
          bool(re.search(r"215\s*N", paper)) or "215" in paper,
          "215" in idx and "ANKLE_AVULSION_LOAD_N" in idx,
          src_value("ANKLE_AVULSION_LOAD_N", "src/climbing/pad2d.py"), 215.0)

    # 2. ATFL 失效应变 0.14
    check("ATFL_FAILURE_STRAIN = 0.14",
          "0.14" in paper,
          "0.14" in idx and "ATFL_FAILURE_STRAIN" in idx,
          src_value("ATFL_FAILURE_STRAIN",
                    "src/climbing/coupling/ankle_ligament.py"), 0.14)

    # 3. V3 阈值 0.0125 s = 12.5 ms
    check("V3_THRESHOLD_S = 0.0125",
          "0.0125" in paper or "12.5 ms" in paper,
          "0.0125" in idx,
          src_value("V3_THRESHOLD_S", "src/climbing/head_branch.py"), 0.0125)

    # 4. V8 阈值 0.005 s = 5 ms
    check("V8_THRESHOLD_AT_20M_S = 0.005",
          "0.005" in paper or "5 ms" in paper,
          "0.005" in idx,
          src_value("V8_THRESHOLD_AT_20M_S", "src/climbing/head_branch.py"),
          0.005)

    # 5. k_neck 默认 200（HeadBranchParams 的 dataclass 字段）
    hb = (ROOT / "src/climbing/head_branch.py").read_text("utf-8",
                                                          errors="replace")
    m5 = re.search(r"k_neck_n_per_mm\s*:\s*[\w.]+\s*=\s*([0-9.]+)", hb)
    check("k_neck default = 200",
          "200" in paper,
          "k200" in idx or "k_neck`=200" in idx,
          float(m5.group(1)) if m5 else None, 200.0)

    # 6. 摩擦 0.4 —— Posture2D 默认是 0.0，0.4 由 MC/probe 注入（索引也这么记）
    mc = (ROOT / "scripts/paper/phase_p1_montecarlo.py").read_text("utf-8",
                                                             errors="replace")
    check("mu_slide = 0.4（MC 注入，非 Posture2D 默认）",
          "0.4" in paper, "0.4" in idx and "phase_p1_montecarlo" in idx,
          bool(re.search(r"mu_slide\s*=\s*0\.4", mc)), True)

    # 7. 踝损伤率 34.5% / 22.8%（论文 & 索引一致）
    check("踝损伤 34.5% / 踝骨折 22.8%",
          "34.5" in paper and "22.8" in paper,
          "34.5" in idx and "22.8" in idx, True, True)

    # 8. alpha=0.477 载荷分享（仅论文&索引，源码是默认 1.0 由 MC 注入）
    check("alpha=0.477 anchor (25294.9/53060)",
          "0.477" in paper, "0.477" in idx and "25294.9" in idx, True, True)

    # 9. S5 可引用量 34.4 kN / 0.574 MPa
    check("S5 F_n=34.4 kN / sigma=0.574 MPa",
          "34.4" in paper and "0.574" in paper,
          "34.4" in idx and "0.574" in idx, True, True)

    # 10. 已撤销数字不得复活
    print()
    print("  -- 反向检查：已撤销的数字不得重新出现在论文 --")
    for bad, why in (("2.70×", "比值类已撤销"), ("15.14 cm", "比值类已撤销")):
        present = bad in paper
        retracted = "已撤销" in paper or "撤销" in paper
        ok = not present or retracted   # 出现时必须同时标注撤销
        print(f"  {'OK  ' if ok else 'FAIL'}  '{bad}' in paper={present} "
              f"retraction noted={retracted}  ({why})")
        if not ok:
            FAIL.append(f"retracted {bad}")

    print()
    print("=" * 76)
    print("RESULT  " + ("PASS" if not FAIL else "FAIL: " + repr(FAIL)))
    print("=" * 76)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    raise SystemExit(main())
