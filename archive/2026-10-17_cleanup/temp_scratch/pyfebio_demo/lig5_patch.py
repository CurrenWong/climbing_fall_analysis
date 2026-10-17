"""L2: 扩充 thums_extract.py 的 TARGETS —— 加舟/骰/楔骨 + 足部韧带。"""
import shutil
from pathlib import Path

ROOT = Path(r"D:\Project\climbing_fall_analysis")
T = ROOT / "scripts/ankle_fe/thums_extract.py"
NPZ = ROOT / "temp/thums/thums_lowerlimb.npz"

# 备份（改提取器会重生成 npz，现有已验证管线依赖它）
bak = NPZ.with_suffix(".npz.bak-preLIG")
if not bak.exists():
    shutil.copy2(NPZ, bak)
    print(f"backup -> {bak.name} ({bak.stat().st_size:,} bytes)")
else:
    print(f"backup exists: {bak.name}")

t = T.read_text(encoding="utf-8")
anchor = """    81100100: ("R_LIG_COLLATERAL_FIBULAR_TISSUES", "solid"),
    81100200: ("R_LIG_COLLATERAL_TIBIAL_TISSUES", "solid"),
}"""
assert t.count(anchor) == 1, f"anchor hit {t.count(anchor)}"

add = """    81100100: ("R_LIG_COLLATERAL_FIBULAR_TISSUES", "solid"),
    81100200: ("R_LIG_COLLATERAL_TIBIAL_TISSUES", "solid"),
    # ── 远端跗骨（★ 2026-10-04 补：胫舟韧带需要舟骨锚定）────────────
    #    实测：跗骨 SPON=solid(hex/tet) / CORT=shell（与距跟骨的 CORT=solid 不同）
    81001400: ("R_NAVICULAR_SPON", "solid"),
    81001401: ("R_NAVICULAR_CORT", "shell"),
    81001500: ("R_CUBOID_SPON", "solid"),
    81001501: ("R_CUBOID_CORT", "shell"),
    81001600: ("R_MEDIAL_CUNEIFORM_SPON", "solid"),
    81001601: ("R_MEDIAL_CUNEIFORM_CORT", "shell"),
    81001700: ("R_INTERMED_CUNEIFORM_SPON", "solid"),
    81001701: ("R_INTERMED_CUNEIFORM_CORT", "shell"),
    81001800: ("R_LATERAL_CUNEIFORM_SPON", "solid"),
    81001801: ("R_LATERAL_CUNEIFORM_CORT", "shell"),
    # ── 足部韧带（★ 补：胫舟 / 弹簧韧带可能在其中）──────────────────
    81101301: ("L_LIG_FOOT_1", "shell"),
    81101302: ("L_LIG_FOOT_2", "shell"),
    81101303: ("L_LIG_FOOT_3", "shell"),
    81101304: ("L_LIG_FOOT_4", "shell"),
    81101305: ("L_LIG_FOOT_5", "shell"),
    81101306: ("L_LIG_FOOT_6", "shell"),
    81101307: ("L_LIG_FOOT_7", "shell"),
    81101308: ("L_LIG_FOOT_8", "shell"),
    81101408: ("L_LIG_FOOT_10", "shell"),
    81101409: ("L_LIG_FOOT_11", "shell"),
    81101410: ("L_LIG_CAPITIS_FIBULAE", "shell"),
}"""
T.write_text(t.replace(anchor, add), encoding="utf-8")

import ast
ast.parse(T.read_text(encoding="utf-8"))
n = sum(1 for l in T.read_text(encoding="utf-8").split("\n")
        if l.strip().startswith("8110") or l.strip().startswith("8100"))
print(f"thums_extract.py TARGETS 现在约 {n} 条；语法 OK")
