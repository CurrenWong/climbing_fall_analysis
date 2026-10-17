"""A1 零回归检查（只生成 .feb，不跑 FEBio）：
1) 用 **git HEAD 原版** ankle_lig_feb.py 生成基线 .feb
2) 用 **改后版** + SUBTALAR=0 生成对照 .feb
3) 逐字节 diff ⇒ 证明「关掉 subtalar 时输出与改前完全一致」
用法: python temp/pyfebio_demo/a1_gen_compare.py
"""
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
WORK = ROOT / "temp/thums/febio_lig"
TMP = ROOT / "temp/pyfebio_demo"
PY = ROOT / ".venv/Scripts/python.exe"
NEW = ROOT / "scripts/ankle_fe/ankle_lig_feb.py"
ORIG = TMP / "_orig_ankle_lig_feb.py"
FEB = WORK / "ankle_lig_WK.feb"

ENV = dict(
    DRIVER="dispx", DISP_MM="+1.6", DET_MIN="0.0", TRUNC_MM="60",
    AXIS_SRC="full", ALIGN_AXIS="1", FIX_TALUS="transverse",
    CART_NEW="1", TIE="1", CART_DIR="temp/thums/cart_m1", POSE_IE="0",
    LOG_DATA="1", NO_RUN="1",
)


def gen(script: Path, extra: dict, out_copy: Path) -> Path:
    env = os.environ.copy()
    env.update(ENV)
    env.update(extra)
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run([str(PY), "-u", str(script)], cwd=str(ROOT), env=env,
                       capture_output=True, text=True, encoding="utf-8",
                       errors="replace", timeout=900)
    print(f"--- {script.name} {extra} rc={r.returncode}")
    for l in (r.stdout or "").splitlines():
        if l.startswith(("[subtalar]", "[cart]", "  接触面", "[tie]", "翻转修正")):
            print("   ", l[:150])
    if r.returncode != 0:
        print((r.stdout or "")[-2000:])
        print((r.stderr or "")[-2000:])
        raise SystemExit(f"生成失败 rc={r.returncode}")
    if not FEB.exists():
        raise SystemExit(f"没有生成 {FEB}")
    shutil.copy(FEB, out_copy)
    return out_copy


def main():
    a = gen(ORIG, {}, TMP / "_cmp_orig.feb")
    b = gen(NEW, {"SUBTALAR": "0"}, TMP / "_cmp_off.feb")
    ta = a.read_text(encoding="latin-1").splitlines()
    tb = b.read_text(encoding="latin-1").splitlines()
    if ta == tb:
        print(f"\n✅ 零回归：SUBTALAR=0 时输出与原版**逐字节一致**（{len(ta)} 行）")
    else:
        n = sum(1 for x, y in zip(ta, tb) if x != y) + abs(len(ta) - len(tb))
        print(f"\n🚨 有差异：{len(ta)} vs {len(tb)} 行，不一致行 ≈ {n}")
        shown = 0
        for i, (x, y) in enumerate(zip(ta, tb), 1):
            if x != y and shown < 20:
                print(f"  L{i}\n   - {x[:120]}\n   + {y[:120]}")
                shown += 1
        for i in range(min(len(ta), len(tb)), max(len(ta), len(tb))):
            print(f"  L{i+1} 仅一侧有: {(ta + tb)[i if i < len(ta) else len(ta) + i - len(tb)][:120]}")
            shown += 1
            if shown > 30:
                break


if __name__ == "__main__":
    main()
