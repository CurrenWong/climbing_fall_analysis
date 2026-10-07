"""B1 力控找根 —— 外层二分，把「位移驱动」的踝 FE 变成「给定力求解」。

背景（`ANKLE_ROADMAP_TODO.md` 阶段 B）：
    B1 二分法找根 —— 压缩曲线单调 0→∞ ⇒ 3~4 次命中 | 能按给定力（200 N / 500 N）求解
    B2 双步序：step1 落位/软组织预载 → step2 加载到目标力

为什么是**外层二分**而不是 `ankle_lig_feb.py` 的原生 `DRIVER=force`：
    原生力控在 t=0 因关节初始张开（中位间隙 ~1 mm、接触仅 ~31 mm²）而病态
    （`ankle_lig_feb.py:1408-1411`），需要 `MULTISTEP=1` 落位协议才稳。
    位移驱动 + 外层二分**不改变求解器**，只在包线内换点 ⇒ 更稳、且能量化
    「哪一段位移对应哪一段力」。

硬约束（`docs/踝子线_scope决策_2026-10-07.md` §三）：
    - 目标力必须落在 **≤519 N** 包线内（A1+A2 量程上界）；超出直接拒绝，不硬跑。
    - **不做** route C（TALUS_K × DISP 扫描）/ route B（换软骨公式）。

用法：
    python scripts/ankle_fe/force_root.py F_TARGET=400
    python scripts/ankle_fe/force_root.py F_TARGET=500 TOL=0.02 MAXITER=6 \\
        FIX_TALUS=transverse CART_NEW=1 TIE=1 CART_DIR=temp/thums/cart_m1 POSE_IE=0
    # 其余 KEY=VAL 原样透传给 acceptance.py（= ankle_lig_feb.py 的 env 接口）
"""
from __future__ import annotations

import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
PY = ROOT / ".venv/Scripts/python.exe"
ACC = ROOT / "scripts/ankle_fe/acceptance.py"

# ── 已知单调压缩曲线（`ANKLE_FIX_TODO.md` 真·压缩曲线，基线配置）──────
#   DISP_MM -> 两顶面合力 (N)。用作**初始括号**，不是最终答案。
KNOWN = ((0.25, 124.0), (0.50, 230.7), (1.00, 418.1), (1.20, 518.88))

#: A1+A2 量程上界（`ANKLE_ROADMAP_TODO.md` §2：+1.2 mm ⇒ 518.88 N）。
#: 超出此值 = 无解（负 Jacobian），不是"错误的大数值"。
F_ENVELOPE_N = 519.0
D_ENVELOPE_MM = 1.25          # 略放宽给二分留余量；再往上必崩

F_RE = re.compile(r"投影到加载轴 = (-?[\d.eE+]+) N")
FMAG_RE = re.compile(r"\|F\|=([\d.eE+]+)")
NORMAL = "NORMAL TERMINATION"


def _interp_disp(f_target: float) -> float:
    """按已知曲线线性反插一个 DISP_MM 初值。"""
    for (d0, f0), (d1, f1) in zip(KNOWN, KNOWN[1:]):
        if f0 <= f_target <= f1:
            return d0 + (d1 - d0) * (f_target - f0) / (f1 - f0)
    return KNOWN[-1][0] if f_target > KNOWN[-1][1] else KNOWN[0][0]


RC_RE = re.compile(r"febio=(-?\d+)\s+febio用时=([\d.]+)s")


def run_point(disp_mm: float, passthrough: dict) -> tuple[float, bool, str]:
    """跑一次 acceptance.py（位移驱动），返回 (|F_ax|, ok, tail)。

    ok = NORMAL 且 FEBio **真的跑过**（rc=0 且用时 ≥ 1 s）。
    ★ 防陈旧读数：acceptance.py 失败时会**保留上一次**的 log/reac 文件，
      直接解析会拿到旧值（实测：配置被吞时每次迭代都返回同一个 401.67 N）。
    """
    args = [str(PY), "-u", str(ACC), f"DISP_MM={disp_mm:.4f}", "DRIVER=dispx"]
    args += [f"{k}={v}" for k, v in passthrough.items()]
    r = subprocess.run(args, cwd=str(ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=3600)
    out = (r.stdout or "") + (r.stderr or "")
    m = F_RE.search(out)
    f_ax = abs(float(m.group(1))) if m else float("nan")
    normal = NORMAL in out
    rc_m = RC_RE.search(out)
    febio_rc = int(rc_m.group(1)) if rc_m else None
    febio_t = float(rc_m.group(2)) if rc_m else 0.0
    fresh = (febio_rc == 0 and febio_t >= 1.0)
    ok = normal and fresh
    tail = " | ".join(out.strip().splitlines()[-3:])[:200]
    if not fresh:
        tail = (f"febio_rc={febio_rc} febio用时={febio_t}s ⇒ 判为未真跑/失败 || " + tail)
    return f_ax, ok, tail


def main() -> int:
    raw = {a.split("=", 1)[0]: a.split("=", 1)[1]
           for a in sys.argv[1:] if "=" in a and not a.startswith("--")}
    if "F_TARGET" not in raw:
        print(__doc__)
        return 2
    f_target = float(raw.pop("F_TARGET"))
    tol = float(raw.pop("TOL", "0.02"))
    max_iter = int(raw.pop("MAXITER", "6"))

    print("=" * 74)
    print(f"B1 力控找根   F_TARGET={f_target:.2f} N   tol={tol:.1%}   max_iter={max_iter}")
    print("=" * 74)
    if f_target > F_ENVELOPE_N:
        print(f"  ✗ 拒绝：目标 {f_target:.1f} N 超出包线 ≤{F_ENVELOPE_N:.0f} N")
        print("    （A1+A2 配置 +1.2 mm ⇒ 518.88 N；再往上是负 Jacobian 无解）")
        print("    → 见 docs/踝子线_scope决策_2026-10-07.md §三")
        return 1

    d_guess = _interp_disp(f_target)
    lo, hi = max(0.02, d_guess * 0.60), min(D_ENVELOPE_MM, d_guess * 1.45)
    print(f"  初值（已知曲线反插）= {d_guess:.4f} mm   括号 [{lo:.4f}, {hi:.4f}]")

    history: list[tuple[float, float, bool]] = []
    for it in range(1, max_iter + 1):
        t0 = time.time()
        f_ax, ok, tail = run_point(d_guess, raw)
        el = time.time() - t0
        rel = abs(f_ax - f_target) / f_target if f_ax == f_ax else float("nan")
        history.append((d_guess, f_ax, ok))
        flag = "✓" if ok else "✗ 未真跑/非 NORMAL"
        print(f"  [{it}] DISP_MM={d_guess:.4f}  F={f_ax:.2f} N  rel={rel:.2%}  "
              f"{flag}  ({el:.1f}s)")
        if not ok:
            print(f"       tail: {tail}")
            print("       ✗ 该点未产生可信新解 ⇒ 中止（禁止拿陈旧读数凑结果）")
            return 1
        if f_ax != f_ax:                    # 解析失败
            print("       ✗ 无法从 acceptance.py 输出解析『投影到加载轴』")
            return 1
        if rel <= tol:
            print(f"\n  ✅ 收敛：DISP_MM={d_guess:.4f} mm ⇒ F={f_ax:.2f} N"
                  f"（目标 {f_target:.2f}，偏差 {rel:.2%} ≤ {tol:.1%}）")
            print(f"     迭代 {it} 次；曲线点："
                  + "  ".join(f"{d:.3f}→{f:.1f}" for d, f, _ in history))
            return 0
        # 单调递增 ⇒ 二分收紧括号
        if f_ax < f_target:
            lo = d_guess
        else:
            hi = d_guess
        d_guess = 0.5 * (lo + hi)
        if hi - lo < 1e-3:
            print(f"\n  ⚠ 括号已收至 {hi-lo:.2e} mm 仍未达 tol；"
                  f"末次 DISP_MM={d_guess:.4f} ⇒ F={f_ax:.2f} N")
            return 1

    print(f"\n  ✗ {max_iter} 次内未达 tol={tol:.1%}；末次 F={history[-1][1]:.2f} N")
    print(f"     曲线点：" + "  ".join(f"{d:.3f}→{f:.1f}" for d, f, _ in history))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
