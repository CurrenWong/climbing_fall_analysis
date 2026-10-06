"""踝关节 + 韧带 FE 模型（对标吴恺 2012 的验证条件）。

吴恺设置的原文（paper/_text/吴恺_足踝FE.txt）：
  · 材料：骨 7300 MPa/0.3，软骨 10 MPa/0.4，韧带 260 MPa/0.4
  · 「韧带使用只承受拉应力单轴梁单元来模拟」
  · 边界：胫、腓骨下端**上截面 600 N 垂直压缩**；**固定跟骨**；**约束距骨**
  · 结果：韧带位移 胫跟 5.24 > 胫舟 5.04 > 跟腓 2.05 > 下胫腓前 1.92 mm

本脚本的两套材料（MATSET 环境变量）：
  WK   = 吴恺原文材料（默认）—— 用于直接对标他的位移数字
  THUMS= THUMS 材料（骨 18000 等）—— 与已验证的 436.8 N 线一致

模型构成：tibia(+spon+marrow) / fibula(+spon) / talus(+spon)
          / calcaneus(+spon) / navicular(+cort+spon)
        + 胫距软骨(已生成) + 距下软骨(A1 接入) + 10 组韧带壳

★ 距下关节接触（A1，2026-10-06 接入）：npz 里 `subtalar_talus`/`subtalar_calcaneus`
  软骨早已生成却从未接进 `.feb` ⇒ 距骨下方只有韧带剪切（伪影 44.8%）。
  现按 `SUBTALAR` 环境变量接入（默认 auto = npz 有新格式键就接）。
  几何预检（`temp/pyfebio_demo/a1_subtalar_gap.py`，只测不改）：
    subtalar 外表面间隙中位 **0.52/0.55 mm**、法向对面（dot 中位 −0.96）
    ⇒ 比已在跑的 tibiotalar 对（1.00/1.16 mm、−0.91）**更贴合**。
  ⚠️ 仍存在少量初始穿透（>0.1 mm 的面片 17/584 与 17/688，最大 −0.56/−0.79 mm）
    —— 与正对照 tibiotalar 同量级（11/174、最大 −1.23 mm）⇒ 一并如实报告。
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

OUT = Path(r"D:\Project\climbing_fall_analysis\temp\thums")
WORK = OUT / "febio_lig"
WORK.mkdir(parents=True, exist_ok=True)
FEBIO = Path(r"D:\Program\FEBioStudio\bin\febio4.exe")

MATSET = os.environ.get("MATSET", "WK").strip() or "WK"
TIME_STEPS = int(os.environ.get("TIME_STEPS", "20"))
FIX_TALUS = os.environ.get("FIX_TALUS", "transverse").strip()
FIB_SHARE = float(os.environ.get("FIB_SHARE", "-1"))   # <0 = 按面积分配
LIG_MODE = os.environ.get("LIG_MODE", "shell").strip()  # shell | none
FORCE_N = float(os.environ.get("FORCE_N", "600"))
DRIVER = os.environ.get("DRIVER", "force").strip() or "force"
DISP_MM = float(os.environ.get("DISP_MM", "0.05"))
# 二分定位用：逗号分隔的域键子集，空 = 全部
ONLY = [x for x in os.environ.get("BONES", "").split(",") if x]
NO_CONTACT = os.environ.get("NO_CONTACT", "0") == "1"
# 实测：模型在 '不写这个 BC' 的情况下就收敛（26 s）——原假设「韧带专属节点
# 转动自由度导致刚体模态」是**错的**，真凶是力控制（见 DRIVER）。默认关。
NO_ROTLOCK = os.environ.get("NO_ROTLOCK", "1") == "1"
# ★ 舟骨只靠 胫舟(→胫骨) + 弹簧(→跟骨) 两条软壳约束 ⇒ 近奇异（LRN-016/017 铁律）。
#   真实足部它被大量未建模结构固定，这里固定它是合理的建模选择。
FIX_NAVICULAR = os.environ.get("FIX_NAVICULAR", "1") == "1"
TRUNC_MM = float(os.environ.get("TRUNC_MM", "0"))
# LOADDIR=axis（默认，沿胫骨长轴）｜z（沿全局 Z）—— 验证"垂直"的定义
LOADDIR = os.environ.get("LOADDIR", "axis").strip() or "axis"
# ★ LOAD_SIGN：载荷/落位位移相对加载轴 LAX 的符号。
#   LAX = 胫骨->距骨（远端）方向；ALIGN_AXIS=1 后 = +z。
#   -1 = 沿 -LAX；+1 = 沿 +LAX（沿轴向远端）。默认 -1 保持既有行为。
LOAD_SIGN = float(os.environ.get("LOAD_SIGN", "-1"))
# force 模式下顶面横向是否约束：none=不加（吴恺配置，允许外翻/内旋）
TOPFIX = os.environ.get("TOPFIX", "lateral").strip() or "lateral"

# ★ 文本读数出口（2026-10-05）：把节点位移/反力写到独立 txt，供 acceptance 判据使用。
#   手册 §3.19.1 + 附录 E.2.1：node_data 的 data 支持 ux;uy;uz / Rx;Ry;Rz。
#   默认关（0）⇒ 不改变既有行为。
LOG_DATA = os.environ.get("LOG_DATA", "0") == "1"

# ★ A1 载荷路径证据（2026-10-06）：logfile 的 `face_data` 读关节接触面的
#   `contact gap` / `contact pressure`（手册 §E.2.2）。**默认 0** ⇒ 不动既有
#   验收；要拿「距下接触力」读数时显式 LOG_CONTACT=1（属附加诊断跑）。
LOG_CONTACT = os.environ.get("LOG_CONTACT", "0") == "1"

# ★ 加载轴来源：trunc（默认，截断后质心差 = 既有行为）| full（截断前全网格质心差）
#   动机：TRUNC_MM 截断会移动质心 ⇒ 加载轴被"顺带"转掉 ~12°（实测）。见 [LRN-20261005-060]
AXIS_SRC = os.environ.get("AXIS_SRC", "trunc").strip() or "trunc"

# ★ ALIGN_AXIS=1（2026-10-05）：把整个模型旋转到「加载轴 = +z」，再让 dispx 直接
#   prescribe z。**根治 `dispx` 只按全局单轴（argmax）近似轴向的语义缺陷**。
#   实测该缺陷会让同一 0.2 mm 的位移给出 87 N vs 159 N（1.8 倍）的力 ⇒ 必须消除。
#   注意：开启后**所有输出位移/坐标都在旋转后的坐标系**里（投影轴=+z，便于比较）。
ALIGN_AXIS = os.environ.get("ALIGN_AXIS", "0") == "1"

# ★ 姿势输入（2026-10-05，为「判不同姿势」而加）：把**足部组**（距骨/跟骨/舟骨
#   + 距侧软骨）整体绕踝关节中心**刚性旋转**。这是关节角的正确运动学表达：
#   骨的相对朝向改变，跨关节的韧带随之被拉伸 ⇒ 姿势成为输入。
#   POSE_DF：背屈(+) / 跖屈(−)，绕"内外侧轴"
#   POSE_IE：内翻(+) / 外翻(−)，绕"前后轴"
#   轴与中心都会**打印出来**（近似必须有名字 —— [LRN-20261005-060] 的规则）。
POSE_DF = float(os.environ.get("POSE_DF", "0") or 0)
POSE_IE = float(os.environ.get("POSE_IE", "0") or 0)

# ★ 姿势的**加载方向参数化**（2026-10-05，**推荐用法**）：不改几何，只把加载/驱动方向
#   相对关节倾斜。理由：不同姿势对关节的差别，本质是**压缩+剪切/力矩配置不同**
#   （T22 的"倾斜地面"正是此思路）。相比"刚性旋转足部组"，它
#     ① 不破坏"相切"的关节面（旋转必然在局部造干涉 ⇒ t=0 就翻单元）
#     ② 全程仍是**位移控制** ⇒ 落在已验证的稳定区
#   POSE_TILT    ：绕**内外侧轴**倾斜（矢状面）⇒ 压缩 + 前后剪切（背屈/跖屈类）
#   POSE_TILT_ML ：绕**前后轴**倾斜（冠状面）⇒ 压缩 + 内外剪切（内外翻类）
POSE_TILT = float(os.environ.get("POSE_TILT", "0") or 0)
POSE_TILT_ML = float(os.environ.get("POSE_TILT_ML", "0") or 0)

# ★ 接触参数可测（2026-10-05）：实测模型只有 ~[0.1,0.3]mm 轴向的窄稳定窗，
#   根因=关节 t=0 未咬合（中位间隙 0.97mm）⇒ 接触病态。`node_reloc` 正是 FEBio 为
#   「初始张开的接触」提供的**落位**机制（把从节点投影到主面）⇒ 值得实测。
#   默认值与既有行为逐字一致。
NODE_RELOC = os.environ.get("NODE_RELOC", "0").strip() or "0"
LAUGON = os.environ.get("LAUGON", "PENALTY").strip() or "PENALTY"
PENALTY = os.environ.get("PENALTY", "1").strip() or "1"
# 软骨单元类型：tet4（旧，拆楔形）/ prism（penta6+hex8，抗压缩好得多）
CART_ELEM = os.environ.get("CART_ELEM", "tet4").strip() or "tet4"
# 软骨材料（锁死判定：线性 tet4 的体积锁死随 ν→0.5 放大）
CART_E = float(os.environ.get("CART_E", "10.0"))
CART_NU = float(os.environ.get("CART_NU", "0.40"))

# ★ TIE（2026-10-05）：骨 ↔ 软骨底面用 `tied-elastic` 连接(`[LRN-051]` 手册 §3.14.6)。
#   仅在使用"新格式软骨"(cart_patch 输出，底面脱离骨节点)时必须开。
#   默认 0 ⇒ 旧格式(底面与骨共享节点)行为逐字不变。
TIE = os.environ.get("TIE", "0") == "1"
TIE_PRIM: dict = {}       # j -> [[全局节点号]...] 软骨底面三角（tie primary）
TIE_PEN = float(os.environ.get("TIE_PENALTY", "1.0"))
# auto_penalty=1 会让 <penalty> 失效（FEBio 自动定刚度）⇒ 罚因子扫描时必须关掉
TIE_AUTO = os.environ.get("TIE_AUTO_PENALTY", "1")
TIE_TOL = float(os.environ.get("TIE_TOL", "1.0"))
TIE_STOL = float(os.environ.get("TIE_SEARCH_TOL", "0.01"))
TIE_SRAD = float(os.environ.get("TIE_SEARCH_RADIUS", "1.0"))
TIE_SEC_TOL = float(os.environ.get("TIE_SEC_TOL", "0.05"))
# ── 路线 A 环境变量（`[LRN-20261005-067/068]`）—— 早于 Material/MeshDomains 段使用 ──
TALUS_RIGID = os.environ.get("TALUS_RIGID", "0") == "1"
TALUS_SPRING_K = float(os.environ.get("TALUS_SPRING_K", "500.0"))   # 接地平动弹簧
TAL_ROTLOCK = os.environ.get("TAL_ROTLOCK", "1") == "1"   # 1=固定转轴(落位) 0=解除(duty cycle)
RIGID_TALUS_MAT = "RIGIDtalus"
RIGID_GROUND_MAT = "RIGIDground"
# ★ 刚体材料的**ID**（追加在 MATERIALS 之后；rigid_connector 的 body_a/body_b 要 ID）
RIGID_TALUS_ID = None       # 在 Material 段写出后回填
RIGID_GROUND_ID = None
# 收敛容差/迭代上限（env 可调，便于调参）
ETOL = os.environ.get("ETOL", "0.001").strip() or "0.001"
DTOL = os.environ.get("DTOL", "0.01").strip() or "0.01"
MAX_REFS = os.environ.get("MAX_REFS", "50").strip() or "50"
MAX_UPS = os.environ.get("MAX_UPS", "10").strip() or "10"
LIG_SUB = [int(x) for x in os.environ.get("LIG_SUB", "").split(",") if x]

print(f"[mode] MATSET={MATSET} LIG_MODE={LIG_MODE} TIME_STEPS={TIME_STEPS} "
      f"FIX_TALUS={FIX_TALUS} FORCE={FORCE_N}N")
print(f"[mode] BONES={ONLY or 'all'}")
print(f"[mode] NO_CONTACT={NO_CONTACT} LIG_SUB={LIG_SUB or 'all'}")
print(f"[mode] FIX_NAVICULAR={FIX_NAVICULAR}")
print(f"[mode] DRIVER={DRIVER} DISP_MM={DISP_MM}")
print(f'[mode] TRUNC_MM={TRUNC_MM}')

# ── 骨部件：域键 -> (npz键, 材料名) ──
BONE_PARTS = {
    "tibia_cort": ("81000700_solid", "mat_tibia_cort"),
    "tibia_spon": ("81000600_solid", "mat_long_spon"),
    "tibia_marrow": ("81000601_solid", "mat_marrow"),
    "fibula_cort": ("81000900_solid", "mat_fibula_cort"),
    "fibula_spon": ("81000800_solid", "mat_long_spon"),
    "talus_cort": ("81001100_solid", "mat_tarsal_cort"),
    "talus_spon": ("81001000_solid", "mat_tarsal_spon"),
    "calcaneus_cort": ("81001300_solid", "mat_tarsal_cort"),
    "calcaneus_spon": ("81001200_solid", "mat_tarsal_spon"),
    # ★ navicular_cort(81001401) 是**壳**(quad4)，不能当 solids 写进 SolidDomain
    #   —— 当作 tet4 会立刻报 Negative jacobian。用 navicular 的 solid 部分即可。
    "navicular_spon": ("81001400_solid", "mat_tarsal_spon"),
}

if MATSET == "WK":
    MATERIALS = {
        "mat_tibia_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_fibula_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_tarsal_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_tarsal_cort_shell": ("isotropic elastic", {"density": 2.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_long_spon": ("isotropic elastic", {"density": 1.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_tarsal_spon": ("isotropic elastic", {"density": 1.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_marrow": ("isotropic elastic", {"density": 1.0e-9, "E": 7300.0, "v": 0.30}),
        "mat_cartilage": ("isotropic elastic", {"density": 1.1e-9, "E": CART_E, "v": CART_NU}),
        "mat_ligament": ("isotropic elastic", {"density": 1.1e-9, "E": 260.0, "v": 0.40}),
    }
    LIG_E_UNIFORM = True
else:
    MATERIALS = {
        "mat_tibia_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 18000.0, "v": 0.30}),
        "mat_fibula_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 18500.0, "v": 0.30}),
        "mat_tarsal_cort": ("isotropic elastic", {"density": 2.0e-9, "E": 15000.0, "v": 0.30}),
        "mat_tarsal_cort_shell": ("isotropic elastic", {"density": 2.0e-9, "E": 15000.0, "v": 0.30}),
        "mat_long_spon": ("isotropic elastic", {"density": 8.6e-10, "E": 160.0, "v": 0.45}),
        "mat_tarsal_spon": ("isotropic elastic", {"density": 1.0e-9, "E": 73.4, "v": 0.45}),
        "mat_marrow": ("isotropic elastic", {"density": 1.0e-9, "E": 12.0, "v": 0.499}),
        "mat_cartilage": ("isotropic elastic", {"density": 1.1e-9, "E": CART_E, "v": CART_NU}),
    }
    LIG_E_UNIFORM = False

# ── 读数据 ──（★ THUMS_NPZ / CART_DIR 可指向加密后的网格；默认与原行为一致）
_NPZ = os.environ.get("THUMS_NPZ", "").strip()
_CART_DIR = os.environ.get("CART_DIR", "").strip()
_NPZ_PATH = Path(_NPZ) if _NPZ else (OUT / "thums_lowerlimb.npz")
_CART_PATH = Path(_CART_DIR) if _CART_DIR else OUT
print(f"[data] npz   = {_NPZ_PATH}")
print(f"[data] carts = {_CART_PATH}")
db = np.load(_NPZ_PATH, allow_pickle=True)
nid_all, xyz_all = db["node_ids"], db["node_xyz"]
idx_all = {int(v): i for i, v in enumerate(nid_all)}
cg = np.load(_CART_PATH / os.environ.get("CART_MESH", "cartilage_mesh.npz"))
lg = np.load(OUT / "ligaments.npz", allow_pickle=True)

WEDGE_T = [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
HEX_T = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
         (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]

mesh_elems = []   # (dom, etype, mat, conn)
used = set()
for dom, (key, mat) in BONE_PARTS.items():
    if ONLY and dom not in ONLY:
        continue
    for row in db[key]:
        conn = []
        for x in row:
            xi = int(x)
            if xi not in (0, -1) and xi not in conn:
                conn.append(xi)
        u = len(conn)
        if u == 8:
            mesh_elems.append((dom, "hex8", mat, conn)); used.update(conn)
        elif u == 4:
            mesh_elems.append((dom, "tet4", mat, conn)); used.update(conn)
        elif u == 6:
            for t in WEDGE_T:
                tet = [conn[z] for z in t]
                mesh_elems.append((dom, "tet4", mat, tet)); used.update(tet)

# 软骨（胫距，已生成）—— 两种格式（2026-10-05）
#   旧 = `{j}_elems` 行内是「底面节点 + 负哨兵」，外表面由 `{j}_offsets` 提供；
#        **底面与骨共享节点** ⇒ 不需要 TIE。实测含 291 个翻转单元 / 22 个被剔除（留洞）。
#   新 = `cart_patch.py` 输出：`{j}_nodes`(全部节点，含底面) + `{j}_elems`(6 节点直接连接)
#        + `{j}_is_offset` + `{j}_tie_faces`。**底面脱离骨节点 ⇒ 必须 TIE=1**。
#        其价值：底面翘曲恒 0、无翻转单元（治「~235 N 上限」）。
CARTS = ["tibiotalar_tibia", "tibiotalar_talus"]
# ★ A1（2026-10-06）：距下关节接入开关。默认 auto = npz 含**新格式**键
#   （`subtalar_*_nodes`，cart_patch 输出）时才接 —— 旧格式 `_offsets` 没在
#   本管线里验证过，不接（打提示）。SUBTALAR=0 强制关（做基线/消融对照）。
SUBTALAR = (os.environ.get("SUBTALAR", "auto").strip() or "auto").lower()
_HAS_SUB = ("subtalar_talus_nodes" in cg.files
            and "subtalar_calcaneus_nodes" in cg.files)
if SUBTALAR == "0":
    pass
elif SUBTALAR in ("1", "auto", "true", "on") and _HAS_SUB:
    CARTS += ["subtalar_talus", "subtalar_calcaneus"]
elif SUBTALAR in ("1", "true", "on"):
    raise SystemExit("[subtalar] SUBTALAR=1 但 npz 缺 subtalar_*_nodes（新格式）键")
print(f"[subtalar] 开关={SUBTALAR} 新格式键={_HAS_SUB} "
      f"⇒ {'接入' if 'subtalar_talus' in CARTS else '不接'}"
      f"（CARTS={CARTS}）")
CART_NEW = (os.environ.get("CART_NEW", "").strip() == "1"
            or "tibiotalar_tibia_nodes" in cg.files)
CN_ISOFF: dict = {}


def _cnidx(x):
    """`{j}_elems` / `{j}_tie_faces` 的索引约定（实测）：
    **负索引** -k ⇒ `{j}_nodes` 的第 k 行（0-based = k-1）。
    例：`-1` ⇒ 第 0 行；`-5` ⇒ 第 4 行。值域 [-N, -1]，无 0、无正数。
    """
    x = int(x)
    return (-x - 1) if x < 0 else (x - 1)


# ★ penta6 雅可比多点采样（实测判据，与 FEBio 的 "N negative jacobians" 逐字吻合）
_PENTA_PTS = [(1.0, 0.0, -1.0), (0.0, 1.0, -1.0), (0.0, 0.0, -1.0),
              (1.0, 0.0, 1.0), (0.0, 1.0, 1.0), (0.0, 0.0, 1.0),
              (1 / 3, 1 / 3, 0.0), (0.25, 0.25, -0.5), (0.25, 0.25, 0.5),
              (0.5, 0.25, -0.5), (0.25, 0.5, -0.5), (0.5, 0.25, 0.5),
              (0.25, 0.5, 0.5), (0.4, 0.4, 0.0)]


def _penta_dN(r, s, t):
    """∂N/∂(r,s,t)，形状 (6,3)。节点 0..2 底面、3..5 顶面。"""
    L = (1 - r - s, r, s)
    dL = ((-1.0, -1.0), (1.0, 0.0), (0.0, 1.0))
    G = np.zeros((6, 3))
    for i in range(3):
        G[i] = ((1 - t) / 2 * dL[i][0], (1 - t) / 2 * dL[i][1], -L[i] / 2)
        G[i + 3] = ((1 + t) / 2 * dL[i][0], (1 + t) / 2 * dL[i][1], L[i] / 2)
    return G


def _min_detJ_penta6(N, idx):
    """该楔形在 14 个采样点上的最小 det(J)。<=0 ⇒ FEBio 会判负雅可比。"""
    X = N[list(idx)]
    worst = np.inf
    for (r, s, t) in _PENTA_PTS:
        d = float(np.linalg.det(X.T @ _penta_dN(r, s, t)))
        if d < worst:
            worst = d
    return worst


# ★ 18 种保持楔形拓扑的排列 = 底面 3 循环 × 顶面 3 循环 × (底面/顶面是否整体对调)
_PENTA_PERM = []
for _bt in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
    for _tp in ((3, 4, 5), (4, 5, 3), (5, 3, 4)):
        _PENTA_PERM.append(tuple(_bt) + tuple(_tp))
        _PENTA_PERM.append(tuple(_tp) + tuple(_bt))

if CART_NEW:
    print("[cart] ★ 新格式（{j}_nodes + {j}_elems + is_offset）"
          f" ⇒ TIE={'1' if TIE else '0'}")
    if not TIE:
        print("  ⚠ 新格式软骨底面**已脱离骨节点** ⇒ 不开 TIE=1 界面会脱开（结果不可信）")
    for j in CARTS:
        CN_ISOFF[j] = np.asarray(cg[f"{j}_is_offset"], bool)
        if j == CARTS[0]:
            print("  [cart] 索引约定：**负索引** -k ⇒ {j}_nodes 的第 k 行（0-based k-1）")
        _ne = 0
        # ★ 单元类型：默认 **penta6**（直接写楔形，保住网格本身的单元形状）。
        #   曾经默认拆成 3 个 tet4 —— 那会在薄层里造出细长 sliver，自检当场
        #   剔除 219/3222 个（留洞）⇒ FEBio t=0 报负 Jacobian。实测。
        _new_et = os.environ.get("CART_NEW_ELEM", "penta6").strip() or "penta6"
        _nd = np.asarray(cg[f"{j}_nodes"], float)
        _nfix = 0
        for row in cg[f"{j}_elems"]:
            _raw = [int(x) for x in row if int(x) != 0]
            conn = [f"cn:{j}:{_cnidx(x)}" for x in _raw]
            if len(conn) == 6:
                # ★ 修"负雅可比"楔形（实测判据与 FEBio 逐字吻合）：
                #   FEBio 报的 "N negative jacobians" == 用**多点采样**(6 顶点 + 8 内部点)
                #   算 penta6 det(J) 为负的单元数（实测 20+42=62，减已修 16 = **46**，
                #   与 FEBio 输出逐字一致）。只在**质心**采样会因正负贡献抵消而漏判。
                #   修法：在 18 种保持楔形拓扑的排列里搜一个让所有采样点 det(J)>0 的
                #   （纯重新编号：几何不变、网格仍共形、不留洞）。
                _idx = [_cnidx(x) for x in _raw]
                if _min_detJ_penta6(_nd, _idx) <= 0.0:
                    for _pm in _PENTA_PERM:
                        _cand = [_idx[z] for z in _pm]
                        if _min_detJ_penta6(_nd, _cand) > 0.0:
                            _idx = _cand
                            _nfix += 1
                            break
                if _new_et == "penta6":
                    mesh_elems.append((j, "penta6", "mat_cartilage",
                                       [f"cn:{j}:{z}" for z in _idx]))
                    _ne += 1
                else:
                    _cc = [f"cn:{j}:{z}" for z in _idx]
                    for t in WEDGE_T:
                        mesh_elems.append((j, "tet4", "mat_cartilage",
                                           [_cc[z] for z in t]))
                    _ne += 3
            elif len(conn) == 4:
                for t in HEX_T:
                    mesh_elems.append((j, "tet4", "mat_cartilage",
                                       [conn[z] for z in t]))
                _ne += 6
            elif len(conn) == 8:
                mesh_elems.append((j, "hex8", "mat_cartilage", conn))
                _ne += 1
        print(f"[cart] {j}: {len(cg[f'{j}_elems'])} 输入 -> {_ne} 单元 / "
              f"{len(cg[f'{j}_nodes'])} 节点（外表面 {int(CN_ISOFF[j].sum())}）"
              f"｜扭转修正 {_nfix}")
else:
    for j in CARTS:
        for row in cg[f"{j}_elems"]:
            vals = [int(x) for x in row if int(x) != 0]
            if not vals or vals[0] <= 0:
                continue
            m = len(vals) // 2
            base = vals[:m]
            offs_ = [f"off:{j}:{-x - 1}" for x in vals[m:]]
            full = base + offs_
            if CART_ELEM == "hex8" and m == 4:
                # ★ 四边形底 → **单个 hex8**（Anderson 2006 的配方：沿法向挤出的砖层）
                #   前提：配合 CART_OFFSET=normal（顶点法向 + 反向翻转），
                #   否则偏移方向不垂直于面 ⇒ hex 歪斜、被闸门剔除。
                mesh_elems.append((j, "hex8", "mat_cartilage", full))
            elif CART_ELEM in ("hex8", "prism", "prism15") and m in (3, 4):
                # ★ 三角形底 → penta6（6 节点楔形）：薄片拆 tet 必为细长 sliver
                #   ⇒ 一压就翻负。四边形底也走这里（hex8 模式除外）。
                tris = ((0, 1, 2),) if m == 3 else ((0, 1, 2), (0, 2, 3))
                for tri in tris:
                    conn = ([full[z] for z in tri]
                            + [full[z + m] for z in tri])
                    mesh_elems.append((j, "penta6", "mat_cartilage", conn))
            else:
                for t in (HEX_T if m == 4 else WEDGE_T):
                    mesh_elems.append((j, "tet4", "mat_cartilage",
                                       [full[z] for z in t]))
            used.update(base)

AXIS_FULL = None      # ★ 截断前（全网格）长轴 —— 见下方 AXIS_SRC
if TRUNC_MM > 0:
    _tal_pts = []
    for _c in mesh_elems:
        if _c[0] == "talus_cort":
            _tal_pts += [xyz_all[idx_all[x]] for x in _c[3]
                         if not isinstance(x, str)]
    _ref = np.mean(_tal_pts, axis=0) if _tal_pts else np.zeros(3)
    _tib_pts = []
    for _c in mesh_elems:
        if _c[0] == "tibia_cort":
            _tib_pts += [xyz_all[idx_all[x]] for x in _c[3]
                         if not isinstance(x, str)]
    _axis = np.mean(_tib_pts, axis=0) - _ref if _tib_pts else np.array(
        [0.0, 0.0, 1.0])
    _axis = _axis / np.linalg.norm(_axis)
    AXIS_FULL = _axis.copy()
    _new = []
    _drop = 0
    for _c in mesh_elems:
        _dom = _c[0]
        if _dom.startswith(("tibia_", "fibula_")):
            _pts = np.array([xyz_all[idx_all[x]] for x in _c[3]
                             if not isinstance(x, str)])
            if float((_pts.mean(0) - _ref) @ _axis) > TRUNC_MM:
                _drop += 1
                continue
        _new.append(_c)
    mesh_elems = _new
    used = {x for _c in mesh_elems for x in _c[3]
            if not isinstance(x, str)}
    print(f"[trunc] 截断到 {TRUNC_MM} mm 内：剔除 {_drop:,} 单元，"
          f"剩 {len(mesh_elems):,} 单元 / {len(used):,} 节点")

print(f"骨+软骨：单元 {len(mesh_elems):,}，节点 {len(used):,}")

# ── 韧带壳 ──
lig_recs = []
if LIG_MODE in ("shell", "truss", "strip", "web"):
    n_lig = int(lg["n_lig"])
    for i in range(n_lig):
        if LIG_SUB and (i + 1) not in LIG_SUB:
            continue
        nm = str(lg["names"][i])
        pid = int(lg["pids"][i])
        ea = float(lg["ea"][i])
        th = float(lg["thick"][i])
        key = f"{pid:08d}_shell"
        quads = [[int(v) for v in q] for q in db[key]]
        mname = f"mat_lig_{nm}"
        if LIG_E_UNIFORM:
            MATERIALS[mname] = ("isotropic elastic",
                                {"density": 1.1e-9, "E": 260.0, "v": 0.40})
        else:
            MATERIALS[mname] = ("isotropic elastic",
                                {"density": 1.1e-9, "E": ea, "v": 0.40})
        a_cross = float(lg["area_cross"][i])
        trusses = [[int(x) for x in t] for t in lg["truss"][i]
                   if int(t[0]) > 0 and int(t[1]) > 0]
        # 杆单元：截面积折进 E（E_eff = EA · A_cross），详见文件头注释
        MATERIALS[mname] = ("isotropic elastic",
                            {"density": 1.1e-9,
                             "E": (ea * a_cross) if LIG_MODE == "truss" else (
                                 260.0 if LIG_E_UNIFORM else ea),
                             "v": 0.40})
        lig_recs.append(dict(idx=i, name=nm, pid=pid, quads=quads,
                             mat=mname, thick=th, ea=ea,
                             a_cross=a_cross, trusses=trusses))
        if LIG_MODE == "shell":
            used.update({int(v) for q in quads for v in q})
        else:
            used.update({int(v) for t in trusses for v in t})
        # strip：把 truss 束两两连成 quad（只用已有节点，零新增）
        strips = []
        for k in range(len(trusses) - 1):
            a0, b0 = trusses[k]
            a1, b1 = trusses[k + 1]
            q = [a0, a1, b1, b0]
            if len(set(q)) == 4:
                strips.append(q)
        # web：三角网 —— 每个 A 节点连到对侧最近 2 个 B 节点，反之亦然
        Aset = sorted({t[0] for t in trusses})
        Bset = sorted({t[1] for t in trusses})
        tris = []
        if Aset and len(Bset) >= 2:
            PA = xyz_all[[idx_all[a] for a in Aset]]
            PB = xyz_all[[idx_all[b] for b in Bset]]
            D = np.linalg.norm(PA[:, None, :] - PB[None, :, :], axis=2)
            for i, a in enumerate(Aset):
                for j in np.argsort(D[i])[:2]:
                    t3 = [a, Bset[int(j)],
                          Bset[int(np.argsort(D[i])[1])]]
                    if len(set(t3)) == 3 and tuple(sorted(t3)) not in \
                            {tuple(sorted(x)) for x in tris}:
                        tris.append(t3)
        if Bset and len(Aset) >= 2:
            if not (Aset and len(Bset) >= 2):
                PB = xyz_all[[idx_all[b] for b in Bset]]
                PA = xyz_all[[idx_all[a] for a in Aset]]
                D = np.linalg.norm(PA[:, None, :] - PB[None, :, :], axis=2)
            Dt = D.T
            for j, b in enumerate(Bset):
                for i in np.argsort(Dt[j])[:2]:
                    t3 = [b, Aset[int(i)],
                          Aset[int(np.argsort(Dt[j])[1])]]
                    if len(set(t3)) == 3 and tuple(sorted(t3)) not in \
                            {tuple(sorted(x)) for x in tris}:
                        tris.append(t3)
        lig_recs[-1]["strips"] = strips
        lig_recs[-1]["tris"] = tris
    print(f"韧带：{len(lig_recs)} 条，"
          f"单元 {sum(len(r['quads']) for r in lig_recs)}，"
          f"并入节点后总计 {len(used):,}")

# ── 紧凑编号 ──
bone_ids = sorted(used)
node_id_of = {v: i + 1 for i, v in enumerate(bone_ids)}
coords = [tuple(xyz_all[idx_all[v]]) for v in bone_ids]
for j in CARTS:
    if CART_NEW:
        # 新格式：底面+外表面**全部**是新节点（底面已脱离骨节点 ⇒ 靠 TIE 连接）
        for k, p in enumerate(cg[f"{j}_nodes"]):
            node_id_of[f"cn:{j}:{k}"] = len(coords) + 1
            coords.append(tuple(p))
        for row in cg[f"{j}_tie_faces"]:
            TIE_PRIM.setdefault(j, []).append(
                [node_id_of[f"cn:{j}:{_cnidx(x)}"] for x in row if int(x) != 0])
    else:
        for k, p in enumerate(cg[f"{j}_offsets"]):
            node_id_of[f"off:{j}:{k}"] = len(coords) + 1
            coords.append(tuple(p))
XYZ = np.asarray(coords, float)
# 软骨外表面（offset）节点的 FE 编号集合 —— 用于识别接触面
OFFSET_IDS = set()
for j in CARTS:
    if CART_NEW:
        OFFSET_IDS |= {node_id_of[f"cn:{j}:{k}"]
                       for k in range(len(cg[f"{j}_nodes"])) if CN_ISOFF[j][k]}
    else:
        OFFSET_IDS |= {node_id_of[f"off:{j}:{k}"]
                       for k in range(len(cg[f"{j}_offsets"]))}
print(f"其中 offset（软骨外表面）节点 {len(OFFSET_IDS):,}")
print(f"合并后节点 {len(coords):,}")

# ── 软骨 p 升级：penta6 -> penta15（9 条边中点，**全局共享**）──
if CART_ELEM == "prism15":
    _id2xyz = {i + 1: np.asarray(coords[i], float)
               for i in range(len(coords))}

    def _xyz(key):
        return _id2xyz[node_id_of[key]]

    # VTK 二次楔形边序：6..14 = mid(0,1) mid(1,2) mid(2,0) mid(3,4) mid(4,5)
    #                             mid(5,3) mid(0,3) mid(1,4) mid(2,5)
    P15_EDGES = [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5),
                 (5, 3), (0, 3), (1, 4), (2, 5)]
    _new, _nmid = [], 0
    for _dom, _et, _mat, _conn in mesh_elems:
        if _et == "penta6" and _mat == "mat_cartilage":
            _mids = []
            for _a, _b in P15_EDGES:
                _ka, _kb = _conn[_a], _conn[_b]
                _key = "mid:" + "|".join(sorted([str(_ka), str(_kb)]))
                if _key not in node_id_of:
                    node_id_of[_key] = len(coords) + 1
                    coords.append(tuple((_xyz(_ka) + _xyz(_kb)) / 2.0))
                    _id2xyz[node_id_of[_key]] = np.asarray(coords[-1], float)
                    _nmid += 1
                _mids.append(_key)
            _new.append((_dom, "penta15", _mat, list(_conn) + _mids))
        else:
            _new.append((_dom, _et, _mat, _conn))
    mesh_elems = _new
    XYZ = np.asarray(coords, float)
    print(f"[prism15] 新增边中点 {_nmid:,} 个；软骨改 penta15；"
          f"节点总数 {len(coords):,}")

# 域 -> 节点集合（用 FE 编号）
dom_nodes = defaultdict(set)
for dom, et, mat, conn in mesh_elems:
    for c in conn:
        dom_nodes[dom].add(node_id_of[c])

# ── 几何：长轴 / 顶面 / 底面 ──
if "tibia_cort" not in dom_nodes or "talus_cort" not in dom_nodes:
    print("子集模式：跳过几何/载荷，仅测网格初始化")
    tib_pts = XYZ[[i - 1 for i in sorted(dom_nodes[sorted(dom_nodes)[0]])]]
    tal_pts = tib_pts
else:
    tib_pts = XYZ[[i - 1 for i in sorted(dom_nodes["tibia_cort"])]]
    tal_pts = XYZ[[i - 1 for i in sorted(dom_nodes["talus_cort"])]]
axis = tib_pts.mean(0) - tal_pts.mean(0)
axis = axis / np.linalg.norm(axis)
# ★ AXIS_SRC（2026-10-05 新增，默认 trunc = 保持既有行为）
#   背景：软网格截断 TRUNC_MM 会改变**截断后**质心 ⇒ 加载轴随之偏转（实测 TRUNC_MM=60
#   时与全网格长轴差 ~12°）。截断本身用的是全网格轴（正确），所以此处给出选择与对照打印。
if AXIS_SRC == "full" and AXIS_FULL is not None:
    _ang = np.degrees(np.arccos(np.clip(float(abs(axis @ AXIS_FULL)), -1, 1)))
    print(f"[axis] 截断后轴={np.round(axis,4)} 全网格轴={np.round(AXIS_FULL,4)} "
          f"夹角={_ang:.2f}deg ⇒ 采用 **全网格轴**（AXIS_SRC=full）")
    axis = AXIS_FULL.copy()
elif AXIS_FULL is not None:
    _ang = np.degrees(np.arccos(np.clip(float(abs(axis @ AXIS_FULL)), -1, 1)))
    print(f"[axis] 截断后轴={np.round(axis,4)} 全网格轴={np.round(AXIS_FULL,4)} "
          f"夹角={_ang:.2f}deg ⇒ 采用截断后轴（默认；置 AXIS_SRC=full 可用全网格轴）")
print(f"加载轴（胫骨->距骨）= {np.round(axis,4)}")
_kax_show = int(np.argmax(np.abs(axis)))
print(f"[axis] argmax(|axis|) = {_kax_show} ⇒ DRIVER=dispx 将 prescribe 'xyz'[{_kax_show}]"
      f" = {'xyz'[_kax_show]}（力控制则沿该轴加载）")

# ── ALIGN_AXIS：把模型旋转到「加载轴 = +z」 ──
#   这样 dispx 的“按全局单轴”就**恰好等于真轴向**，消除 argmax 近似带来的 1.8 倍歧义。
if ALIGN_AXIS:
    _a = axis / np.linalg.norm(axis)
    _z = np.array([0.0, 0.0, 1.0])
    _v = np.cross(_a, _z)
    _c = float(_a @ _z)
    if np.linalg.norm(_v) < 1e-12:
        _R = np.eye(3) if _c > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        _vx = np.array([[0.0, -_v[2], _v[1]],
                        [_v[2], 0.0, -_v[0]],
                        [-_v[1], _v[0], 0.0]])
        _R = np.eye(3) + _vx + _vx @ _vx * ((1.0 - _c) / float(_v @ _v))
    for _i in range(len(coords)):
        coords[_i] = tuple((_R @ np.asarray(coords[_i], float)).tolist())
    XYZ = np.asarray(coords, float)
    print(f"[align] ALIGN_AXIS=1：模型已旋转，加载轴 -> +z（旋转角 "
          f"{np.degrees(np.arccos(np.clip(_c, -1, 1))):.2f}deg）；"
          f"输出坐标/位移均在旋转后坐标系")
    axis = _z.copy()
    print(f"加载轴（胫骨->距骨）= {np.round(axis,4)}  [已对齐 +z]")

# ── 姿势输入：把足部组绕踝关节中心刚性旋转 ──
#   组 = 距骨 / 跟骨 / 舟骨（皮质+松质） + 距侧软骨 tibiotalar_talus
#   （胫侧软骨、胫腓骨、韧带近端留在原处 ⇒ 跨关节韧带被拉伸 = 生理的）
#   轴：AP ≈ 跟骨质心 → 舟骨质心（足的长轴）；ML = AP × 加载轴（与 AP 正交的"内外侧轴"）
if POSE_DF or POSE_IE or POSE_TILT or POSE_TILT_ML:
    _FOOT = [d for d in dom_nodes if d.startswith(("talus_", "calcaneus_", "navicular_"))] \
        + [c for c in CARTS if c.endswith("_talus")]
    # ★ A1 连带修正：距下关节的**两侧**软骨都长在足部组上（距骨 + 跟骨）⇒ 姿势
    #   旋转必须带上 `subtalar_calcaneus`，否则只转距下距骨侧会把这对接触面撕开
    #   （tibiotalar 只有距骨侧在足部组 = 正确，因为胫侧留在原处）。
    if "subtalar_calcaneus" in CARTS:
        _FOOT.append("subtalar_calcaneus")
    _foot_ids = sorted({i for d in _FOOT for i in dom_nodes[d]})
    _cart_pt = np.array([XYZ[i - 1] for i in sorted(dom_nodes.get(
        "tibiotalar_talus", _foot_ids))])
    _ctr = _cart_pt.mean(0) if len(_cart_pt) else XYZ[[i - 1 for i in _foot_ids]].mean(0)
    _heel = XYZ[[i - 1 for i in sorted(dom_nodes.get("calcaneus_cort", _foot_ids))]].mean(0)
    _mid = XYZ[[i - 1 for i in sorted(dom_nodes.get("navicular_spon", _foot_ids))]].mean(0)
    _ap = _mid - _heel
    _ap = _ap / max(np.linalg.norm(_ap), 1e-12)
    # ★ 内外侧轴 = AP × 加载轴（⊥ 矢状面 ⇒ 屈伸绕它）。**不要**再做第二次 cross：
    #   那会把轴转进 AP–加载轴平面（实测得到与加载轴夹角仅 ~16° 的错误轴）。
    _ml = np.cross(_ap, axis)
    _ml = _ml / max(np.linalg.norm(_ml), 1e-12)
    print(f"[pose] 自检 |ML·加载轴| = {abs(float(_ml @ axis)):.4f}（应≈0）"
          f"   AP 与加载轴夹角 = {np.degrees(np.arccos(np.clip(abs(float(_ap @ axis)), -1, 1))):.1f}deg")

    def _rot(u, deg):
        u = np.asarray(u, float)
        u = u / max(np.linalg.norm(u), 1e-12)
        th = np.radians(deg)
        K = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
        return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    _R = _rot(_ml, POSE_DF) @ _rot(_ap, POSE_IE)
    print(f"[pose] 足部组 {len(_foot_ids)} 节点 / {len(_FOOT)} 域")
    print(f"[pose] 旋转中心 = {np.round(_ctr,3)} mm")
    print(f"[pose] 内外侧轴 ML = {np.round(_ml,4)}   前后轴 AP = {np.round(_ap,4)}")
    print(f"[pose] POSE_DF(背屈+/-跖屈) = {POSE_DF} deg   POSE_IE(内翻+/-外翻) = {POSE_IE} deg")

    # ── 关节面间距：用软骨**外表面**（off: 节点）的最小距离做度量 ──
    _off = {}

    def _outer(j):
        # ★ `[ERR-20261005-010]`：原实现只认**旧格式** `{j}_offsets`，
        #   而 `CART_NEW`（cart_patch）用 `{j}_nodes` + `{j}_is_offset`（键名 `cn:{j}:{k}`）
        #   ⇒ 只要开 POSE_TILT 就会 `KeyError: '..._offsets is not a file in the archive'`
        #   ⇒ **姿势功能自新软骨格式起一直是坏的**（姿势矩阵曾因此读到旧 .feb 而"无差异"）。
        if j not in _off:
            if CART_NEW:
                _off[j] = [node_id_of[f"cn:{j}:{k}"]
                           for k in range(len(cg[f"{j}_nodes"])) if CN_ISOFF[j][k]]
            else:
                _off[j] = [i for k in range(len(cg[f"{j}_offsets"]))
                           for i in [node_id_of.get(f"off:{j}:{k}")]
                           if i and (i - 1) < len(coords)]
        return _off[j]

    def _gap_at(XYZ_):
        """给定坐标数组时的关节面最小间距。"""
        A = XYZ_[[i - 1 for i in _outer("tibiotalar_tibia")]]
        B = XYZ_[[i - 1 for i in _outer("tibiotalar_talus")]]
        if not len(A) or not len(B):
            return float("nan")
        return float(np.sqrt(((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)).min())

    def _gap():
        return _gap_at(XYZ)

    _d0 = _gap()
    for _i in _foot_ids:
        coords[_i - 1] = tuple((_ctr + _R @ (np.asarray(coords[_i - 1], float) - _ctr)).tolist())
    XYZ = np.asarray(coords, float)
    _d1 = _gap()
    # ★★ 重新落位：沿加载轴平移足部组，把关节面**尽量贴近**（Anderson 的"minimal
    #   contact"）。`[ERR-20261005-011]` 三轮试错记录：
    #     ① `_shift = _d1 - _d0` 沿 +axis 平移 —— 符号反了（越推越开 0.3625→0.4416）
    #     ② 割线法 —— 不收敛（残差 +0.33），因为"最近点对"随平移**换点**、gap(s) 非光滑
    #     ③ 括入+二分 —— 更差（0.4664），因为 **gap(s) 是 V 形不是单调**，二分找目标值无解
    #   ⇒ 正解：对 s 做**一维极小化**（粗扫 + 细化）。旋转后两软骨面已不再相切，
    #     单纯平移**不一定**能回到转前的 0.0211 ⇒ 目标是"可达的最小间距"而非某个指定值。
    _base = {_i: np.asarray(coords[_i - 1], float).copy() for _i in _foot_ids}

    def _xyz_at(s):
        X = np.asarray(coords, float).copy()
        for _i in _foot_ids:
            X[_i - 1] = _base[_i] + s * axis
        return X

    _shift_tot = 0.0
    if _d0 == _d0 and _d1 == _d1:
        # ★ 平移量**必须有界**：实测无界极小化会给出 +2.04 mm 的荒谬平移
        #   （关节面被滑到无关位形、韧带被预拉伸）⇒ 限制 |s| ≤ 0.6 mm，
        #   与 Anderson「小幅调整到 minimal contact」一致。`[ERR-20261005-011]`
        _LIM = float(os.environ.get("POSE_SEAT_LIM", "0.6"))
        _best_s, _best_g = 0.0, _d1
        _s = -_LIM
        while _s <= _LIM + 1e-9:                 # 粗扫 0.01 mm
            _g = _gap_at(_xyz_at(_s))
            if _g == _g and _g < _best_g:
                _best_s, _best_g = _s, _g
            _s += 0.01
        _s = _best_s - 0.01                      # 细化 0.001 mm
        while _s <= _best_s + 0.0101:
            _g = _gap_at(_xyz_at(_s))
            if _g == _g and _g < _best_g:
                _best_s, _best_g = _s, _g
            _s += 0.001
        _shift_tot = _best_s
    for _i in _foot_ids:
        coords[_i - 1] = tuple((_base[_i] + _shift_tot * axis).tolist())
    XYZ = np.asarray(coords, float)
    _d2 = _gap()
    print(f"[pose] 关节面最小间距：转前 {_d0:.4f} → 转后 {_d1:.4f} → "
          f"**落位后 {_d2:.4f} mm**（平移 {_shift_tot:+.4f} mm 沿加载轴，"
          f"一维极小化）")
    if _d2 == _d2 and _d2 > 0.15:
        print(f"[pose] ⚠ 旋转后可达的最小间距仍 {_d2:.4f} mm ⇒ 驱动位移须 > 该值"
              f"（否则关节在加载前是张开的，会出负 Jacobian）")
    print(f"[pose] 已应用（输出坐标/位移在旋转后的姿势下）")
LAX = axis.copy()
if LOADDIR == "z":
    _z = np.array([0.0, 0.0, 1.0])
    LAX = _z if float(axis @ _z) > 0 else -_z
    print(f"[loaddir] 沿全局 Z 加载，与胫骨轴夹角 "
          f"{np.degrees(np.arccos(abs(float(axis @ _z)))):.1f} deg")
# ★ 载荷/落位方向 = LOAD_SIGN × LAX（LAX 在此处才最终确定，故定义放这里）
LOOKSIGNLAX = LOAD_SIGN * LAX
_SIGN_TXT = "压缩(推向关节)" if LOAD_SIGN > 0 else "★拉伸(把关节拉开)★"
print(f"[load] LOAD_SIGN={LOAD_SIGN:+.0f} ⇒ 载荷方向 = {np.round(LOOKSIGNLAX,4)}"
      f"（LAX={np.round(LOOKSIGNLAX/LOAD_SIGN if LOAD_SIGN else 1,4)}）")
print(f"[load] 该方向在关节处表现为：{_SIGN_TXT}")
if LOAD_SIGN < 0:
    # `[ERR-20261005-009]`：`-LAX` 把关节**拉开**（拉伸）⇒ 能被韧带撑住的量级只有 ~10-200 N，
    # 历史上所有"力控制崩在 ~12 N"的怪现象都源于此（含被当成成功的 592.61 N）。
    print("[load] ⚠️ 警告：当前是**拉伸**方向 ⇒ 反力由韧带承担、量级虚低。"
          "要压实关节请设 LOAD_SIGN=+1（压缩侧 −1.0 mm ⇒ −608.82 N，全程 NORMAL）。")


def face_of(dom, etype, key):
    """取某域的极值面片（key=max 取最高，min 取最低）。"""
    faces = []
    for dom2, et, mat, conn in mesh_elems:
        if dom2 != dom or et != etype:
            continue
        for f in ([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)] if et == "tet4"
                  else [(0, 1, 2, 3)]):
            faces.append([conn[z] for z in f])
    ids = [node_id_of[c] for f in faces for c in f]
    pr = XYZ[[i - 1 for i in ids]] @ axis
    lim = pr.max() if key == "max" else pr.min()
    span = pr.max() - pr.min()
    sel = [f for f in faces
           if abs(np.mean(XYZ[[node_id_of[c] - 1 for c in f]] @ axis) - lim)
           < 0.02 * span]
    return sel


FACES_TET = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
FACES_HEX = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]


TET_FACES = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]
FACES_TET = TET_FACES
MIN_FACET_AREA = 1e-8


def facet_ok(f):
    if len(set(f)) != len(f):
        return False
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return np.linalg.norm(n) / 2.0 >= MIN_FACET_AREA


def surface_faces(dom, etype):
    """该域的**真外表面**：朝外定向后，只保留出现 1 次的面。

    ★ 这两点都不能省：
      ① 既遍历单元的**所有**面再按「出现次数==1」筛 —— 只看极值会把
         **内部面**当成载荷面（实测胫骨 2415 个六面体却选出 1551 个"顶面"）。
      ② 朝向必须朝外（法向背离单元质心），否则 contact/load 的正负号会反。
    """
    cnt, info = {}, {}
    for ids in blocks[(dom, etype)]:
        P = XYZ[[i - 1 for i in ids]]
        cen = P.mean(0)
        for f in (FACES_TET if etype == "tet4" else FACES_HEX):
            fv = [ids[z] for z in f]
            pts = XYZ[[i - 1 for i in fv]]
            n = np.zeros(3)
            for i in range(len(pts)):
                n += np.cross(pts[i], pts[(i + 1) % len(pts)])
            if np.dot(n, pts.mean(0) - cen) < 0:
                n = -n
                fv = fv[::-1]
            key = tuple(sorted(fv))
            if not facet_ok(fv):
                continue
            info[key] = tuple(fv)
            cnt[key] = cnt.get(key, 0) + 1
    return [info[k] for k, c in cnt.items() if c == 1]


def _face_normal(f):
    P = XYZ[[i - 1 for i in f]]
    n = np.zeros(3)
    for i in range(len(P)):
        n += np.cross(P[i], P[(i + 1) % len(P)])
    return n / 2.0


def extreme_faces(dom, which, band=25.0, cos_min=0.6):
    """该域**朝向加载轴**的端面（顶 or 底）。

    ★ 不能只用「沿轴投影带」挑：加载轴是斜的时候（THUMS 里胫骨轴 ≈
      [0.68,0.14,0.72]），水平截面的投影会散开好几 mm ⇒ 3mm 带只切到平台一角
      （实测 56.8 mm²，真值 ~1000 mm²）。**按面法向与轴对齐**才是几何上正确的判据。
    """
    fs = []
    for et in ("hex8", "tet4"):
        if (dom, et) in blocks:
            fs += surface_faces(dom, et)
    if not fs:
        return []
    sgn = 1.0 if which == "max" else -1.0
    aligned = [f for f in fs if float(_face_normal(f) @ axis) * sgn > cos_min]
    if not aligned:
        return []
    proj = XYZ[[f[0] - 1 for f in aligned]] @ axis
    lim = proj.max() if which == "max" else proj.min()
    if which == "max":
        return [f for f, p in zip(aligned, proj) if p > lim - band]
    return [f for f, p in zip(aligned, proj) if p < lim + band]




def area_proj(faces):
    """面的**投影面积**（沿加载轴 `axis`）。

    ⚠️ 2026-10-06 修 `[ERR-20261006-005]`：原 quad4 分支用
    `0.5*|(P1-P0)×(P3-P2)·axis|`（"对边之差"的叉积）。对**平面**四边形它碰巧
    等于法向，但骨的六面体外表面是**翘曲四边形**（非平面），该向量根本不是法向
    ⇒ 投影面积被严重低估。实测（胫骨顶 114 面，`p2a_repro.py`）：

        原式  47.75 mm²   ← 模型一直用的值
        正解 555.52 mm²   ← 投影面积（拆 2 三角）
        真实 2148.31 mm²  ← 未投影

    ⇒ 差 11.6×（投影）/ 45×（真实）。影响：`build_feb.py` 力控的
    `traction = FORCE_N / top_area` 除数错 11.6× ⇒ 力控实验曾整体加载过量。
    **位移控制不受影响**（BC 是 prescribe 位移，与面积无关）。

    正解：quad4 拆成 2 个三角形，各自取 `0.5*|(b-a)×(c-a)·axis|`。
    """
    a = 0.0
    for n in faces:
        P = XYZ[[c - 1 for c in n]]
        if len(P) == 3:
            a += 0.5 * abs(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), axis))
        else:
            a += 0.5 * abs(np.dot(np.cross(P[1] - P[0], P[2] - P[0]), axis))
            a += 0.5 * abs(np.dot(np.cross(P[2] - P[0], P[3] - P[0]), axis))
    return a



# ── 组装 .feb ──
lines = []
w = lines.append
w('<?xml version="1.0" encoding="ISO-8859-1"?>')
w('<febio_spec version="4.0">')
w('  <Module type="solid"/>')
w('  <Globals><Constants><T>0</T><R>0</R><Fc>0</Fc></Constants></Globals>')
w('  <Material>')
for i, (name, (mt, pr)) in enumerate(MATERIALS.items(), 1):
    w(f'    <material id="{i}" name="{name}" type="{mt}">')
    for k, v in pr.items():
        w(f'      <{k}>{v}</{k}>')
    w('    </material>')
if TALUS_RIGID:
    # ★ 刚体材料：E/ν 对结果无影响，但 contact 的 auto-penalty 会用到 ⇒ 给合理值
    _nmat = len(MATERIALS)
    _nmat += 1
    RIGID_TALUS_ID = _nmat
    w(f'    <material id="{_nmat}" name="{RIGID_TALUS_MAT}" type="rigid body">')
    w('      <density>1e-9</density>')
    w('      <E>12000</E><v>0.42</v>')
    w('    </material>')
    _nmat += 1
    RIGID_GROUND_ID = _nmat
    w(f'    <material id="{_nmat}" name="{RIGID_GROUND_MAT}" type="rigid body">')
    w('      <density>1e-9</density>')
    w('      <E>12000</E><v>0.42</v>')
    w('    </material>')
w('  </Material>')

# 元素块：按 (dom, etype) 分组
blocks = defaultdict(list)
for dom, et, mat, conn in mesh_elems:
    blocks[(dom, et)].append([node_id_of[c] for c in conn])
# 韧带块
lig_blocks = {}
LIG_ETYPE = {"truss": "line2", "web": "tri3"}.get(LIG_MODE, "quad4")
if LIG_MODE in ("shell", "truss", "strip", "web"):
    for r in lig_recs:
        if LIG_MODE == "shell":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in q]
                                     for q in r["quads"]]
        elif LIG_MODE == "truss":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["trusses"]]
        elif LIG_MODE == "strip":
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["strips"]]
        else:
            lig_blocks[r["name"]] = [[node_id_of[int(v)] for v in t]
                                     for t in r["tris"]]

# 韧带专属节点（仅被壳用、不被任何实体用）—— shell 模式下它们的**转动
# 自由度无约束** ⇒ 刚体模态（实测：首个载荷步 22 个单元翻负）。
LIG_ONLY_NS = []
if LIG_MODE == "shell":
    _solid_nodes = set()
    for _k, _v in blocks.items():
        for _c in _v:
            _solid_nodes.update(_c)
    _lig_nodes = {x for _v in lig_blocks.values() for _c in _v for x in _c}
    LIG_ONLY_NS = sorted(_lig_nodes - _solid_nodes)
    print(f"韧带专属节点（转动自由度原本自由）: {len(LIG_ONLY_NS)}")
    if NO_ROTLOCK:
        LIG_ONLY_NS = []


# ── ★ Jacobian 自检 + 朝向修正 + 剔除 ────────────────────────────────
# 这一步在 build_feb.py 里就有，第一版新脚本**漏搬**了 ⇒ FEBio 直接报
# "Negative jacobian detected during mesh initialization"。
# 教训（LRN-20261003-018）：**复用已验证管线的结构时，不能只搬"看得见"的部分**
# —— 这类"清理/修正"环节最容易被漏掉，而它恰恰是能不能跑起来的关键。
GP2 = 0.577350269189626


def detJ_hex8(P):
    """P:(8,3) 在 8 个 Gauss 点算 det(J)，返回最小值。"""
    xi = np.array([[-GP2, -GP2, -GP2], [GP2, -GP2, -GP2], [GP2, GP2, -GP2],
                   [-GP2, GP2, -GP2], [-GP2, -GP2, GP2], [GP2, -GP2, GP2],
                   [GP2, GP2, GP2], [-GP2, GP2, GP2]])
    ref = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
    worst = 1e30
    for sg in xi:
        dN = np.zeros((8, 3))
        for i in range(8):
            a, b, c = ref[i]
            dN[i, 0] = a * (1 + b * sg[1]) * (1 + c * sg[2]) / 8
            dN[i, 1] = b * (1 + a * sg[0]) * (1 + c * sg[2]) / 8
            dN[i, 2] = c * (1 + a * sg[0]) * (1 + b * sg[1]) / 8
        worst = min(worst, float(np.linalg.det(P.T @ dN)))
    return worst


def detJ_tet4(P):
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0],
                                         P[3] - P[0]])))


def detJ_penta6(P):
    """楔形 Jacobian 最小值：在 (r,s,t) 采样网格上取 min det(J)。

    形状函数 N0=(1-r-s)(1-t), N1=r(1-t), N2=s(1-t),
             N3=(1-r-s)t,     N4=r t,     N5=s t
    J = Σ_i P_i ⊗ [dN/dr, dN/ds, dN/dt]

    ★ 只在底面形心算会漏掉：三角**内部**或沿 t 的退化会让 FEBio 在积分点
      报 "Negative jacobian during mesh initialization"，而单点检查算出来是正的。

    ⚠️ 2026-10-06 修 `[LRN-20261006-085]`：原 `RS` 只到 **2/3**
    ⇒ 三角的两个顶点 `(r,s)=(1,0)` 与 `(0,1)` **从未被采样**
    （只有 `(0,0)` 被采到）⇒ 漏掉「顶点附近局部翻转」的楔形单元。
    实测（`temp/pyfebio_demo/p2_locate.py`，软骨 tibiotalar 两域）：

        RS 到 2/3（原）：min = **2.3914e-03**（看起来"只是小"，静默通过）
        RS 到 1.0（现）：min = **−7.8900e+00**（11 个单元 ≤0、12 个 ≤1e-2）

    这些单元全在关节接触面/tie 面上，正是 σ_max=104.85 尖峰的来源之一。
    **判据 = 采样点 + 阈值 + 方向，三者缺一不可。**

    ⚠️ 2026-10-06 第二修：把完整采样**降级为诊断用**（`detJ_penta6_full`），
    本函数**恢复原 RS**（只到 2/3）。原因：用完整采样驱动"剔除"会多剔 11 个
    **接触面上的**单元 ⇒ **FEBio rc=1（在关节面打洞）**，实测。
    ⇒ **检测与处置必须分离**：完整采样负责**报警**，剔除仍走原判据，
      真正的修复走 `cart_patch.py` 的局部重铺（见 ANKLE_FIX_TODO §第十一轮）。
    """
    RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3)
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    worst = None
    for r in RS:
        for sq in RS:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in TS:
                d = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                J = P.T @ d                      # (3,3)
                det = float(np.linalg.det(J))
                if worst is None or det < worst:
                    worst = det
    return worst


def sv_penta6(P):
    """楔形**带符号体积**（依赖节点顺序，用于朝向判定/修正）。

    J = [P1-P0, P2-P0, dt]，dt = (P3+P4+P5-P0-P1-P2)/3
    ★ 交换上下三角（3 次对换 = 奇置换）会正确变号。
      采样版 detJ_penta6 是按形状函数求和、**不能**可靠地靠调换节点来变号，
      所以朝向判定必须另用本函数 —— 这正是此前负朝向楔形「修不动、只能剔」
      的原因（剔除率 tet4 2% vs penta6 14%）。
    """
    dt = (P[3] + P[4] + P[5] - P[0] - P[1] - P[2]) / 3.0
    return float(np.linalg.det(np.stack([P[1] - P[0], P[2] - P[0], dt])))


def detJ_penta15(P):
    """penta15 的角点落在 penta6 位置上 ⇒ 用其 6 角点做检查（代理）。"""
    return detJ_penta6(P[:6])


def detJ_penta6_full(P):
    """楔形 det(J) 的**完整采样**版 —— **只用于诊断报警，绝不驱动剔除**
    （`[LRN-20261006-085]`）。

    与原版唯一区别：`RS` 含 **1.0** ⇒ 覆盖三角的两个"远端"顶点
    `(r,s)=(1,0)` / `(0,1)`。原采样只到 2/3，漏掉它们。

    实测（软骨 tibiotalar 两域，`temp/pyfebio_demo/p2_locate.py`）：
        原版 min = **2.3914e-03**（"只是小"，静默通过）
        完整版 min = **−7.8900e+00**（11 个 ≤0、12 个 ≤1e-2）
    这些单元**全在关节接触面/tie 面上** ⇒ 只能局部重铺，**不能剔除**
    （实测剔除 42 个后 FEBio **rc=1**：在关节面打洞）。
    """
    RS = (0.0, 1 / 6, 1 / 3, 1 / 2, 2 / 3, 1.0)
    TS = (0.0, 0.25, 0.5, 0.75, 1.0)
    worst = None
    for r in RS:
        for sq in RS:
            if r + sq > 1.0 + 1e-9:
                continue
            for t in TS:
                d = np.array([
                    [-(1 - t), -(1 - t), -(1 - r - sq)],
                    [1 - t, 0.0, -r],
                    [0.0, 1 - t, -sq],
                    [-t, -t, (1 - r - sq)],
                    [t, 0.0, r],
                    [0.0, t, sq],
                ])
                det = float(np.linalg.det(P.T @ d))
                if worst is None or det < worst:
                    worst = det
    return worst


# ── ★ 完整采样诊断（报警用，不驱动剔除）──
if os.environ.get("DET_REPORT", "1") == "1":
    _nneg = _nsm = 0
    _wmin = None
    for _key in blocks:
        _dom, _et = _key
        if _et not in ("penta6", "penta15"):
            continue
        for _ids in blocks[_key]:
            _w = detJ_penta6_full(XYZ[[i - 1 for i in _ids]])
            _wmin = _w if _wmin is None else min(_wmin, _w)
            if _w <= 0:
                _nneg += 1
            elif _w <= 1e-2:
                _nsm += 1
    print(f"  [detJ自检·完整采样] worst-det(J) min = {_wmin:.4e}；"
          f"≤0 的单元 {_nneg} 个、0<x≤1e-2 的 {_nsm} 个")
    if _nneg or _nsm:
        print(f"  [detJ自检·完整采样] 🚨 有 {_nneg + _nsm} 个楔形在顶点附近翻转"
              f" ⇒ 局部应力不可信（σ_max 尖峰来源）；"
              f"**不可剔除**（在关节接触面上）⇒ 走 cart_patch.py 局部重铺")


FN = {"hex8": detJ_hex8, "tet4": detJ_tet4, "penta6": detJ_penta6, "penta15": detJ_penta15}
DET_MIN = float(os.environ.get("DET_MIN", "0.5"))
# ★ 不能只 >0：近退化薄片会最先翻负。但**闸门太严会把承力单元也挖掉**
#   （实测 DET_MIN=0.5 时软骨楔形被剔 14+29 个，接触面 150/198→138/171，
#     同位移下支反力从 259 N 掉到 117 N）⇒ 可按需放宽。
#   FEBio 自身只要求积分点 det>0，故放宽到 0.02~0.1 通常仍能初始化。
bad_total = fixed_total = 0
print("\n=== Jacobian 自检 + 朝向修正（det(J) 必须 > 0）===")
for key in list(blocks):
    dom, et = key
    fn = FN[et]
    keep, nbad, dmins = [], 0, []
    for ids in blocks[key]:
        P = XYZ[[i - 1 for i in ids]]

        def _ok(PP):
            """朝向 + 质量双判据。"""
            if et in ("penta6", "penta15"):
                if sv_penta6(PP[:6]) <= 0:
                    return None
            return FN[et](PP)

        d = _ok(P)
        if d is not None and d > DET_MIN:
            keep.append(ids); dmins.append(d); continue
        # ★ 朝向修正必须用**奇置换**：
        #   tet4 = 交换两个节点（1 次对换）
        #   penta6 = 交换上下三角（3 次对换）
        #   ✗ 整体倒序是偶置换，行列式符号不变
        if et == "tet4":
            alt = [ids[1], ids[0]] + ids[2:]
        elif et == "penta6":
            alt = ids[3:] + ids[:3]
        elif et == "penta15":
            alt = None      # 15 节点的奇置换会把中点挪到别的棱（几何变了）
        else:
            alt = ids[len(ids) // 2:] + ids[:len(ids) // 2]
        d2 = _ok(XYZ[[i - 1 for i in alt]]) if alt is not None else None
        if d2 is not None and d2 > DET_MIN:
            keep.append(alt); dmins.append(d2); fixed_total += 1
        else:
            nbad += 1
    blocks[key] = keep
    bad_total += nbad
    mn = min(dmins) if dmins else float("nan")
    print(f"  {'OK ' if not nbad else '!! '} {dom:22s} {et:6s} "
          f"保留 {len(keep):>6,} 剔除 {nbad:>4,}  det(J) min={mn:.4e}")
print(f"  翻转修正 {fixed_total} 个；仍无法修正而剔除 {bad_total} 个")

# 韧带壳：quad4 的退化检查（绕向是任意的，只查退化，不翻转）
if LIG_MODE in ("shell", "strip", "web"):
    nbad_shell = 0
    for nm in list(lig_blocks):
        keep = []
        for ids in lig_blocks[nm]:
            P = XYZ[[i - 1 for i in ids]]
            if len(P) == 3:
                d = np.cross(P[1] - P[0], P[2] - P[0])
            else:
                d = np.cross(P[2] - P[0], P[3] - P[1])
            if np.linalg.norm(d) > 1e-6:
                keep.append(ids)
            else:
                nbad_shell += 1
        lig_blocks[nm] = keep
    print(f"  韧带壳退化剔除 {nbad_shell} 个")

# ── 软骨接触面（★ 必须取 tet 的**真实外表面**，不能用 npz 原始四边形）──
# 这两条都是实测换来的：
#   ① 接触面与网格实际面不重合 ⇒ 接触**完全不检测**，但 rc=0、无警告、只是应力全 0
#   ② 反向：面必须朝外定向（法向背离对顶点）
CART_FACES = {
    # penta15 用角点面即可：中点落在棱上，几何几乎一致（接触只需 tri3）
    "penta15": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "tet4": TET_FACES,
    "penta6": [(0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)],
    "hex8": FACES_HEX,
}


def cart_outer(j):
    """软骨**外表面**（全部由 offset 节点构成的面），朝外定向。

    ★ 支持 tet4 / penta6 / hex8 —— 换单元类型后必须同步换面表，
      否则接触面与网格真实面不重合 ⇒ 接触完全不检测（rc=0 无警告、应力全 0）。
    """
    out, seen = [], set()
    for et in ("tet4", "penta6", "penta15", "hex8"):
        for ids in blocks.get((j, et), []):
            for f in CART_FACES[et]:
                fv = [ids[z] for z in f]
                if not all(x in OFFSET_IDS for x in fv):
                    continue
                key = tuple(sorted(fv))
                if key in seen or not facet_ok(fv):
                    continue
                seen.add(key)
                rest = [ids[z] for z in range(len(ids)) if z not in f]
                cen = XYZ[[i - 1 for i in rest]].mean(0)
                pts = XYZ[[i - 1 for i in fv]]
                n = np.zeros(3)
                for i in range(len(pts)):
                    n += np.cross(pts[i], pts[(i + 1) % len(pts)])
                if np.dot(n, pts.mean(0) - cen) < 0:
                    fv = fv[::-1]
                out.append(fv)
    return out


contact_surfs = {j: cart_outer(j) for j in CARTS}

# ── TIE：骨 ↔ 软骨底面（`tied-elastic`，手册 §3.14.6，连接非共形网格）──
#   primary   = 软骨**底面**三角（直接取新格式的 `{j}_tie_faces`，不经反查）
#   secondary = 该关节「A 骨」外表面中**靠近软骨底面**的局部片（KD-tree within TOL）
#   ★ 用整块骨外表面太贵且易误抓远处 facet；只取"与底面重合的 facet"又会缺 ⇒
#     折中：质心距任一底面节点 < TIE_SEC_TOL 的那些面（两端要各取一次）
TIE_SEC: dict = {}
if TIE:
    from scipy.spatial import cKDTree as _KD
    _BONE_OF = {"tibiotalar_tibia": "tibia_cort",
                "tibiotalar_talus": "talus_cort",
                "subtalar_talus": "talus_cort",
                "subtalar_calcaneus": "calcaneus_cort"}
    for j in CARTS:
        # ★ primary 必须来自**存活的**软骨单元（反查真实 facet），不能照抄 npz 的
        #   `{j}_tie_faces`：那些面可能属于被 Jacobian 闸门**剔除**的单元 ⇒
        #   其余软骨在该处与骨脱开 ⇒ 自由飞 ⇒ t=0 报 11 个负 Jacobian（实测）。
        #   build_feb.py 同样做反查，原因一致。
        _bot = {node_id_of[f"cn:{j}:{k}"]
                for k in range(len(cg[f"{j}_nodes"])) if not CN_ISOFF[j][k]}
        prim, _seen = [], set()
        # ★ 底面 facet 表：tet4 的 4 个三角 / penta6 的 5 个面（底面 = (0,1,2)）。
        #   不能只扫 tet4 —— 默认单元已改 penta6，只扫 tet4 会**静默取到空集**，
        #   tie 被跳过、软骨脱开，而日志只是轻描淡写一句"跳过"（实测踩过）。
        for _et, _flist in (("tet4", TET_FACES), ("penta6", CART_FACES["penta6"])):
            for ids in blocks.get((j, _et), []):
                for f in _flist:
                    fv = [ids[z] for z in f]
                    if not all(x in _bot for x in fv):
                        continue
                    key = tuple(sorted(fv))
                    if key in _seen or not facet_ok(fv):
                        continue
                    _seen.add(key)
                    _rest = [ids[z] for z in range(len(ids)) if z not in f]
                    _cen = XYZ[[i - 1 for i in _rest]].mean(0)
                    _pts = XYZ[[i - 1 for i in fv]]
                    _n = np.zeros(3)
                    for i in range(len(_pts)):
                        _n += np.cross(_pts[i], _pts[(i + 1) % len(_pts)])
                    if np.dot(_n, _pts.mean(0) - _cen) < 0:
                        fv = fv[::-1]
                    prim.append(fv)
        if not prim:
            print(f"  [tie] {j}: 存活单元里找不到底面三角，跳过")
            continue
        _pids = sorted({x for f in prim for x in f})
        _tree = _KD(XYZ[[i - 1 for i in _pids]])
        A_bone = _BONE_OF.get(j)
        sec, _seen = [], set()
        for et in ("hex8", "tet4", "penta6", "penta15"):
            if (A_bone, et) not in blocks:
                continue
            for f in surface_faces(A_bone, et):
                key = tuple(sorted(f))
                if key in _seen:
                    continue
                d, _ = _tree.query(XYZ[[i - 1 for i in f]])
                if float(d.min()) <= TIE_SEC_TOL:
                    _seen.add(key)
                    sec.append(f)
        if not prim or not sec:
            print(f"  [tie] {j}: primary {len(prim)} / secondary {len(sec)}，跳过")
            continue
        TIE_SEC[j] = sec
        print(f"  [tie] {j}: primary(软骨底面) {len(prim)} 面 / "
              f"secondary(骨关节面, 坐标匹配≤{TIE_SEC_TOL}) {len(sec)} 面")

for j, v in contact_surfs.items():
    print(f"  接触面 {j}: {len(v)} 面")

# ── 面的选取（必须在 blocks 建好之后）──
top_tibia = extreme_faces("tibia_cort", "max")
top_fibula = extreme_faces("fibula_cort", "max")
bot_calc = extreme_faces("calcaneus_cort", "min")
print(f"外表面顶/底：胫骨顶 {len(top_tibia)} / 腓骨顶 {len(top_fibula)} / "
      f"跟骨底 {len(bot_calc)}")

A_tib, A_fib = area_proj(top_tibia), area_proj(top_fibula)


def _area_real(faces):
    """未投影的真实面积（quad 拆 2 三角）—— 只用于自检。"""
    a = 0.0
    for n in faces:
        P = XYZ[[c - 1 for c in n]]
        if len(P) == 3:
            a += 0.5 * np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
        else:
            a += 0.5 * np.linalg.norm(np.cross(P[1]-P[0], P[2]-P[0]))
            a += 0.5 * np.linalg.norm(np.cross(P[2]-P[0], P[3]-P[0]))
    return a


# ── ★ 面积自检（`[ERR-20261006-005]`）：投影面积必须与真实面积相容 ──
#   投影 ≤ 真实（投影是真实面元的轴向分量）；若投影/真实 过小，说明面的
#   法向算错或选到了侧向面 —— 这曾让 top_area 低估 11.6×，力控整体过量加载。
_R_tib = _area_real(top_tibia)
_R_fib = _area_real(top_fibula)
for _nm, _ap, _ar in (("胫骨", A_tib, _R_tib), ("腓骨", A_fib, _R_fib)):
    if _ar > 0:
        _r = _ap / _ar
        _flag = "✅" if 0.20 <= _r <= 1.05 else "🚨"
        print(f"  [area自检] {_nm}顶：投影 {_ap:.1f} / 真实 {_ar:.1f} mm² "
              f"=> 比 {_r:.3f} {_flag}")
        if not (0.20 <= _r <= 1.05):
            print(f"  [area自检] 🚨 投影/真实 = {_r:.3f} 越界（应 0.2~1.0）"
                  f" ⇒ 面积不可信，禁止用于力控 traction！")
share = A_fib / (A_tib + A_fib) if FIB_SHARE < 0 else FIB_SHARE
print(f"顶面投影面积：胫骨 {A_tib:.1f} mm²，腓骨 {A_fib:.1f} mm² "
      f"=> 腓骨分担 {share*100:.1f}% 载荷")
tal_faces = extreme_faces("talus_cort", "min")

# ── ★ 距骨「弱弹簧」支撑（Anderson 2006 §2.2 的思路；`[LRN-20261005-065]`）────────
#  现状(`FIX_TALUS=transverse`)：距骨底面把**横向分量硬置零** ⇒ 两端横向硬约束
#  ⇒ 胫距两轴被强制平行 ⇒ **内外翻/内旋在运动学上做不到**
#  （注意：实体节点只有 3 个平动自由度，**没有"旋转自由度"可锁**；
#   卡住转动的是"相对转动所需的横向位移被禁"这件事本身）。
#  `TALUS_SUPPORT=spring`：把横向硬约束换成**有限刚度弹簧**
#  = 接地哑节点 + 2 节点线性弹簧 ⇒ **既消刚体模态奇异、又保留姿势自由度**
#    ⇒ 姿势可以由力/力矩驱动而**产生**（而不是像 `POSE_TILT` 那样只能当输入）。
#
#  语法已**实测**（User Manual §3.6.7 `DiscreteSet` + §3.16.1.1 `linear spring`）：
#    <DiscreteSet name="s"><delem>i,g</delem></DiscreteSet>        （Mesh 段内）
#    <Discrete><discrete_material id="1" type="linear spring"><E>k</E></...>
#              <discrete dmat="1" discrete_set="s"/></Discrete>     （顶层）
#  最小算例校验：Δk=99、δ=0.1 ⇒ ΔF 实测 **104.70 N** = 理论 105 N（21 步求和）✓
#  另两个实测语法点：`<NodeSet name="x">1,2,3</NodeSet>`（**逗号列表**，不是 <node> 子元素）；
#  `<zero displacement>` 用 `<x_dof>1</x_dof>`，`<prescribed displacement>` 用
#  `<dof>x</dof><value lc="1">…</value>`（`lc` 指向 load curve ⇒ 需要 LoadData 段）。
TALUS_SUPPORT = os.environ.get("TALUS_SUPPORT", "fix").strip() or "fix"
TALUS_K = float(os.environ.get("TALUS_K", "500.0"))          # 每个横向方向的合力刚度
TALUS_SPRING_L = float(os.environ.get("TALUS_SPRING_L", "5.0"))   # 哑节点臂长（定向用）
_tal_dofs = ("xyz" if FIX_TALUS == "all" else
             "".join(d for d, v in zip("xyz", np.abs(axis))
                     if v < np.abs(axis).max() - 1e-6)) or "xz"
TAL_GROUND: dict = {}          # d -> (ground node ids, [(talus_node, ground_node), ...])
TAL_FIXD: dict = {}            # d -> 是否真的建了弹簧
if TALUS_SUPPORT == "spring":
    _tal_ids = sorted({c for f in tal_faces for c in f}) if tal_faces else []
    if not _tal_ids:
        print("⚠ [talus-spring] 距骨底面节点为空 ⇒ 退化为 fix")
    else:
        _knode = TALUS_K / max(len(_tal_ids), 1)
        for _d in _tal_dofs:
            _e = np.zeros(3)
            _e["xyz".index(_d)] = 1.0
            _gids, _dl = [], []
            for _ni in _tal_ids:
                _p = np.asarray(coords[_ni - 1], float)
                coords.append(tuple((_p + TALUS_SPRING_L * _e).tolist()))
                _g = len(coords)
                _gids.append(_g)
                _dl.append((_ni, _g))
            TAL_GROUND[_d] = (_gids, _dl)
            TAL_FIXD[_d] = True
        print(f"[talus-spring] {len(_tal_ids)} 个距骨节点 × 方向 {_tal_dofs}"
              f" ⇒ {sum(len(v[1]) for v in TAL_GROUND.values())} 根弹簧；"
              f"单根 k={_knode:.4f} N/mm（合计 {TALUS_K} N/mm/方向）；"
              f"臂长 {TALUS_SPRING_L} mm ⇒ 哑节点 {len(coords) - len(_tal_ids)*len(_tal_dofs) + 1}"
              f"..{len(coords)}")
        print("[talus-spring] ⚠ 横向硬约束已**移除**（改为弹簧）⇒ 内外翻/内旋自由度已释放")

# ★★ 多步落位（Anderson §2.2 的落位协议；`[LRN-20261005-070]`）
#   为什么需要：`DRIVER=force` 失败的根因是**关节初始张开**（t=0 中位间隙 ~1 mm、
#   接触仅 ~31 mm²）⇒ 加载时胫骨与距骨只靠弹簧 + 韧带相连 ⇒ 近似机构 ⇒ 不收敛。
#   Anderson 的解法：step1 把胫骨带到**最小接触** → step2 才加载。
#   FEBio 规则（手册 §3.17 / §6.1）：**主干 `<Boundary>` 贯穿所有 step，step 内的只在该 step 生效**；
#   **多步时 `<Control>` 不再放文件顶部**，改放文件底部 `<Step>` 里。
# ★ LOAD_SIGN 已在文件顶部声明（须早于 LAX 使用处）
MULTISTEP = os.environ.get("MULTISTEP", "0") == "1"
SEAT_MM = float(os.environ.get("SEAT_MM", "0.50"))          # step1 落位位移（轴向 mm）
TIME_STEPS_SEAT = int(os.environ.get("TIME_STEPS_SEAT", "10"))
# ★ SEAT_HOLD=1：step2 里**继续保持**落位位移（不撤约束）
#   诊断用：区分"撤约束后弹性回弹导致关节张开"与"载荷本身问题"。
SEAT_HOLD = os.environ.get("SEAT_HOLD", "0") == "1"

# ── ★★ 路线 A：距骨做成 FEBio **刚体**（Anderson §2.2 的忠实复现；`[LRN-20261005-067/068]`）──
#  为什么必须刚体：Anderson 的"固定/解除**绕内外翻轴与内旋轴**的转动"，
#  在纯可变形实体网格里**无法表达**（实体节点只有 3 个平动 DOF，见 `[LRN-20261005-065]`）；
#  他是靠 `assumed rigid ... surfaces` + `reference nodes`（刚体 6 DOF，含 3 转动）做到的。
#  实测通过的原语（最小算例已验证）：
#    <material type="rigid body">            / <SolidDomain type="rigid-solid">
#    <bc type="rigid" node_set="…"><rb>名字</rb></bc>
#    <Rigid><rigid_bc type="rigid_fixed"> Rx_dof..Rw_dof  （★把 R*_dof 写 0 = 解除该转轴）
#           <rigid_connector type="rigid spring"> body_a/body_b/insertion_a/insertion_b/k
#           <rigid_load type="rigid_force"> <dof>Rz</dof><value lc="1">…</value>
#  坑：rigid_force 用 <dof>+<value>（不是 <force>）；刚体"接地"必须造哑刚体；
#      <plotfile> 没有 rigid displacement/rotation 变量 ⇒ 读节点位移。
#  FEBio 4 行为：刚体 DOF **初始全自由** ⇒ 必须自己加约束。
TAL_RIGID_NS: list = []
GROUND_IDS: list = []
if TALUS_RIGID:
    _tn = set()
    for _d in ("talus_cort", "talus_spon", "talus_marrow"):
        _tn |= set(dom_nodes.get(_d, ()))
    TAL_RIGID_NS = sorted(_tn)
    # 地面：2 个哑节点（挂在全约束刚体上），位置取距骨几何中心上方/下方远处
    if TAL_RIGID_NS:
        _c = np.asarray([coords[i - 1] for i in TAL_RIGID_NS], float).mean(0)
    else:
        _c = np.zeros(3)
    for _k in (1.0, -1.0):
        coords.append(tuple((_c + _k * 50.0 * axis).tolist()))
        GROUND_IDS.append(len(coords))
    print(f"[talus-rigid] ★ 路线 A：距骨 {len(TAL_RIGID_NS)} 个节点 → 刚体 `{RIGID_TALUS_MAT}`；"
          f"地面哑节点 {GROUND_IDS}；接地弹簧 k={TALUS_SPRING_K} N/mm")
    if TAL_GROUND:
        print("[talus-rigid] ⚠ 与 TALUS_SUPPORT=spring 同时开：刚体模式下节点级弹簧无意义，"
              "已自动忽略 TALUS_SUPPORT")
        TAL_GROUND = {}
        TAL_FIXD = {}


w('  <Mesh>')
w('    <Nodes name="all">')
for i, p in enumerate(coords, 1):
    w(f'      <node id="{i}">{p[0]:.6f},{p[1]:.6f},{p[2]:.6f}</node>')
w('    </Nodes>')
eid = 0
for (dom, et), v in blocks.items():
    sel = [c for c in v if False] or v
    w(f'    <Elements type="{et}" name="{dom}__{et}">')
    for c in v:
        eid += 1
        w(f'      <elem id="{eid}">{",".join(str(x) for x in c)}</elem>')
    w('    </Elements>')
LIG_EID0 = eid
if LIG_MODE in ("shell", "truss", "strip", "web"):
    for nm, conns in lig_blocks.items():
        if not conns:                      # ★ 空 <Elements> 会让 FEBio 报
            continue                        #   "missing attribute id"
        w(f'    <Elements type="{LIG_ETYPE}" name="lig_{nm}">')
        for c in conns:
            eid += 1
            w(f'      <elem id="{eid}">{",".join(str(x) for x in c)}</elem>')
        w('    </Elements>')
LIG_EID1 = eid


def ns_of(faces):
    return sorted({c for f in faces for c in f})


if top_tibia:
    w(f'    <NodeSet name="tibia_top_ns">'
  f'{",".join(str(x) for x in ns_of(top_tibia))}</NodeSet>')
if top_fibula:
    w(f'    <NodeSet name="fibula_top_ns">'
  f'{",".join(str(x) for x in ns_of(top_fibula))}</NodeSet>')
if bot_calc:
    w(f'    <NodeSet name="calcaneus_bottom_ns">'
  f'{",".join(str(x) for x in ns_of(bot_calc))}</NodeSet>')
# 距骨：取其沿轴最低的一批节点做横向约束

w(f'    <NodeSet name="talus_low_ns">'
  f'{",".join(str(x) for x in ns_of(tal_faces))}</NodeSet>')

if LIG_ONLY_NS:
    w(f'    <NodeSet name="lig_only_ns">'
      f'{",".join(str(x) for x in LIG_ONLY_NS)}</NodeSet>')
if FIX_NAVICULAR and "navicular_spon" in dom_nodes:
    _nav = sorted(dom_nodes["navicular_spon"])
    w(f'    <NodeSet name="navicular_all_ns">'
      f'{",".join(str(x) for x in _nav)}</NodeSet>')
if top_tibia:
    w('    <Surface name="tibia_top">')
    k4 = [f for f in top_tibia if len(f) == 4]
    k3 = [f for f in top_tibia if len(f) == 3]
    for n, f in enumerate(k4, 1):
        w(f'      <quad4 id="{n}">{",".join(str(c) for c in f)}</quad4>')
    for n, f in enumerate(k3, 1):
        w(f'      <tri3 id="{n}">{",".join(str(c) for c in f)}</tri3>')
    w('    </Surface>')
if top_fibula:
    w('    <Surface name="fibula_top">')
    k4 = [f for f in top_fibula if len(f) == 4]
    k3 = [f for f in top_fibula if len(f) == 3]
    for n, f in enumerate(k4, 1):
        w(f'      <quad4 id="{n}">{",".join(str(c) for c in f)}</quad4>')
    for n, f in enumerate(k3, 1):
        w(f'      <tri3 id="{n}">{",".join(str(c) for c in f)}</tri3>')
    w('    </Surface>')
for j, sel in contact_surfs.items():
    w(f'    <Surface name="surf_{j}">')
    for n, f in enumerate(sel, 1):
        w(f'      <tri3 id="{n}">{",".join(str(x) for x in f)}</tri3>')
    w('    </Surface>')
# ★ A1：关节接触对统一定义（与 build_feb.py 的 PAIRS 同序：(pair, primary, secondary)）
JOINT_PAIRS = [("tibiotalar_pair", "tibiotalar_tibia", "tibiotalar_talus")]
if "subtalar_talus" in CARTS and "subtalar_calcaneus" in CARTS:
    JOINT_PAIRS.append(("subtalar_pair", "subtalar_talus", "subtalar_calcaneus"))
for _pn, _pa, _pb in JOINT_PAIRS:
    w(f'    <SurfacePair name="{_pn}">')
    w(f'      <primary>surf_{_pa}</primary>')
    w(f'      <secondary>surf_{_pb}</secondary>')
    w('    </SurfacePair>')
for j, sec in TIE_SEC.items():
    w(f'    <Surface name="tiesurf_{j}">')
    for n, f in enumerate(TIE_PRIM[j], 1):
        w(f'      <tri3 id="{n}">{",".join(str(x) for x in f)}</tri3>')
    w('    </Surface>')
    w(f'    <Surface name="bonesurf_{j}">')
    for n, f in enumerate(sec, 1):
        _tg = {3: "tri3", 4: "quad4"}.get(len(f), "tri3")
        w(f'      <{_tg} id="{n}">{",".join(str(x) for x in f)}</{_tg}>')
    w('    </Surface>')
    w(f'    <SurfacePair name="tiepair_{j}">')
    w(f'      <primary>tiesurf_{j}</primary>')
    w(f'      <secondary>bonesurf_{j}</secondary>')
    w('    </SurfacePair>')
if TAL_GROUND:
    _allg = sorted({_g for _gd, _dl in TAL_GROUND.values() for _g in _gd})
    w(f'    <NodeSet name="tal_ground_ns">'
      f'{",".join(str(g) for g in _allg)}</NodeSet>')
    for _d, (_gids, _dl) in TAL_GROUND.items():
        w(f'    <DiscreteSet name="tal_spring_{_d}">')
        for _a, _b in _dl:
            w(f'      <delem>{_a},{_b}</delem>')
        w('    </DiscreteSet>')
if TALUS_RIGID:
    w(f'    <NodeSet name="talus_rigid_ns">'
      f'{",".join(str(v) for v in TAL_RIGID_NS)}</NodeSet>')
    w(f'    <NodeSet name="ground_rigid_ns">'
      f'{",".join(str(v) for v in GROUND_IDS)}</NodeSet>')
w('  </Mesh>')
# ★ <Discrete> 必须写在 <MeshDomains> **之后**（实测：放在 </Mesh> 与 <MeshDomains>
#   之间时 FEBio 报 `tag "discrete" : invalid value for attribute`）。
#   小节顺序依 User Manual §3.1 的推荐序。

w('  <MeshDomains>')
_dommat = {d: m for d, (k, m) in BONE_PARTS.items()}
_dommat.update({j: "mat_cartilage" for j in CARTS})
_TAL_DOMS = {"talus_cort", "talus_spon", "talus_marrow"}
for dom, et in blocks:
    if TALUS_RIGID and dom in _TAL_DOMS:
        # ★ Anderson 的 "assumed rigid ... surfaces" ⇒ 距骨域为刚性
        w(f'    <SolidDomain name="{dom}__{et}" type="rigid-solid" '
          f'mat="{RIGID_TALUS_MAT}"/>')
    else:
        w(f'    <SolidDomain name="{dom}__{et}" mat="{_dommat[dom]}"/>')
if LIG_MODE == "shell":
    if LIG_MODE in ("shell", "strip", "web"):
        for r in lig_recs:
            if not lig_blocks.get(r["name"]):
                continue
            w(f'    <ShellDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}">')
            w(f'      <shell_thickness>{r["thick"]}</shell_thickness>')
            w('    </ShellDomain>')
    else:
        # ★ 实测：BeamDomain type 只认 linear-truss / elastic-truss / linear-beam；
        #   elastic-truss + "isotropic elastic" 可跑通（linear-truss 不收标准材料，
        #   linear-beam 直接崩溃）。截面积折进 E，见文件头。
        for r in lig_recs:
            w(f'    <BeamDomain name="lig_{r["name"]}" '
              f'mat="{r["mat"]}" type="elastic-truss"/>')
w('  </MeshDomains>')
# ★★ 路线 A 的 <Rigid> 段移到 </Boundary> **之后**（实测：放在 MeshDomains 之后会在
#    FEBio 初始化时报 `Nonlinear constraint 1 (<unnamed>) failed to initialize`；
#    最小算例跑通的顺序是 …</Boundary> → <Rigid> → <Output>）。
_RIGID_LINES: list = []
if TALUS_RIGID:
    _R = _RIGID_LINES.append
    _ax = int(np.argmax(np.abs(axis)))          # 轴向（ALIGN_AXIS 后应为 2=z）
    # (a) 地面：全约束
    _R('  <Rigid>')
    _R('    <rigid_bc type="rigid_fixed">')
    _R(f'      <rb>{RIGID_GROUND_MAT}</rb>')
    _R('      <Rx_dof>1</Rx_dof><Ry_dof>1</Ry_dof><Rz_dof>1</Rz_dof>')
    _R('      <Ru_dof>1</Ru_dof><Rv_dof>1</Rv_dof><Rw_dof>1</Rw_dof>')
    _R('    </rigid_bc>')
    # (b) 距骨：沿轴(垂直于地面)平动约束；横向平动由弹簧提供 ⇒ 不硬约束
    #     ★ 转动 Ru/Rv/Rw 由 TAL_ROTLOCK 控制：1=固定（落位期），0=解除（duty cycle）
    _R('    <rigid_bc type="rigid_fixed">')
    _R(f'      <rb>{RIGID_TALUS_MAT}</rb>')
    for _di, _dn in enumerate("xyz"):
        if _di == _ax:
            _R(f'      <R{_dn}_dof>1</R{_dn}_dof>')     # 垂直地面位移 = 约束
    _R(f'      <Ru_dof>{1 if TAL_ROTLOCK else 0}</Ru_dof>')
    _R(f'      <Rv_dof>{1 if TAL_ROTLOCK else 0}</Rv_dof>')
    _R(f'      <Rw_dof>{1 if TAL_ROTLOCK else 0}</Rw_dof>')
    _R('    </rigid_bc>')
    # (c) 接地弱弹簧（translational）：挂在距骨几何中心 ↔ 地面点
    #     ★★ body_a/body_b 是**材料 ID**（不是名字！）。真模型里刚体材料追加在最后，
    #     不能像最小算例那样写死 1/2 —— 写错会引用到可变形材料 ⇒ FEBio 报
    #     `Nonlinear constraint 1 (<unnamed>) failed to initialize`（实测踩过）。
    _tc = np.asarray([coords[i - 1] for i in TAL_RIGID_NS], float).mean(0)
    _gc = np.asarray([coords[i - 1] for i in GROUND_IDS[0:1]], float)[0]
    _R('    <rigid_connector type="rigid spring">')
    _R(f'      <body_a>{RIGID_TALUS_ID}</body_a><body_b>{RIGID_GROUND_ID}</body_b>')
    _R(f'      <insertion_a>{_tc[0]:.4f},{_tc[1]:.4f},{_tc[2]:.4f}</insertion_a>')
    _R(f'      <insertion_b>{_gc[0]:.4f},{_gc[1]:.4f},{_gc[2]:.4f}</insertion_b>')
    _R(f'      <k>{TALUS_SPRING_K}</k>')
    _R('    </rigid_connector>')
    _R('  </Rigid>')
    print(f"[talus-rigid] 接地弹簧：距骨中心 {np.round(_tc,2)} ↔ 地面 {np.round(_gc,2)}，"
          f"k={TALUS_SPRING_K}；转动锁定={TAL_ROTLOCK}；"
          f"body_a=mat#{RIGID_TALUS_ID} body_b=mat#{RIGID_GROUND_ID}")
# ★ 距骨弱弹簧的 <Discrete> 段（必须在 MeshDomains 之后；语法见上方 TAL_GROUND 注释）
if TAL_GROUND:
    w('  <Discrete>')
    # 每个方向一个材料，刚度 = TALUS_K / 该方向节点数（⇒ 方向合力 = TALUS_K）
    for _n, (_d, (_gids, _dl)) in enumerate(TAL_GROUND.items(), 1):
        w(f'    <discrete_material id="{_n}" type="linear spring">')
        w(f'      <E>{TALUS_K / max(len(_dl), 1):.6f}</E>')
        w('    </discrete_material>')
    for _n, _d in enumerate(TAL_GROUND, 1):
        w(f'    <discrete dmat="{_n}" discrete_set="tal_spring_{_d}"/>')
    w('  </Discrete>')

w('  <LoadData>')
w('    <load_controller id="1" type="loadcurve">')
w('      <interpolate>LINEAR</interpolate><extend>CONSTANT</extend>')
w('      <points><pt>0.0,0.0</pt><pt>1.0,1.0</pt></points>')
w('    </load_controller>')
if MULTISTEP:
    # ★★ 第二条曲线：**step1(t∈[0,1]) 恒为 0，step2(t∈[1,2]) 爬升到 1**。
    #   没有它的话 <traction> 无 lc ⇒ FEBio 把**满载荷立刻加上**
    #   ⇒ 实测 step2 第一个增量就报 48 个负 Jacobian（t=1.00091、载荷≈满值）。
    w('    <load_controller id="2" type="loadcurve">')
    w('      <interpolate>LINEAR</interpolate><extend>CONSTANT</extend>')
    w('      <points><pt>0.0,0.0</pt><pt>1.0,0.0</pt><pt>2.0,1.0</pt></points>')
    w('    </load_controller>')
w('  </LoadData>')
w('  <Loads>')
_use_force = (DRIVER == "force") and not MULTISTEP   # 多步时载荷移到 step2
F_t = FORCE_N * (1 - share)
F_f = FORCE_N * share
for surf, F, _faces in (("tibia_top", F_t, top_tibia),
                      ("fibula_top", F_f, top_fibula)):
    if not _use_force or not _faces or F <= 0:
        continue
    # ★ BUGFIX(2026-10-05, [ERR-20261005-007])：原先以下 7 行缩进在 if 体内，
    #   紧跟 continue ⇒ 全为死代码 ⇒ DRIVER=force 时 <Loads> 恒为空 ⇒
    #   模型无载荷。此前"力控制病态(刚体模态)"的诊断为误诊。已 dedent 修正。
    # traction = 力 / 面积，方向沿 -axis（向下压）
    a = area_proj(top_tibia if surf == "tibia_top" else top_fibula)
    tr = F / a if a > 0 else 0.0
    w(f'    <surface_load type="traction" surface="{surf}">')
    # ★ 单步路径同样要显式斜坡（否则满载荷瞬加；曲线 1 = (0,0)->(1,1)）
    w('      <scale lc="1">1.0</scale>')
    w(f'      <traction>{LOOKSIGNLAX[0]*tr:.6f},{LOOKSIGNLAX[1]*tr:.6f},{LOOKSIGNLAX[2]*tr:.6f}'
      f'</traction>')
    w('    </surface_load>')
w('  </Loads>')

w('  <Boundary>')
if bot_calc:
    w('    <bc name="fix_calcaneus" type="zero displacement" '
      'node_set="calcaneus_bottom_ns">')
    w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
    w('    </bc>')
# 距骨约束：transverse = 只约束横向；all = 全约束
# ★ 若开了弱弹簧(TAL_GROUND)，则**不再**对距骨底面做横向硬约束 ——
#   改由接地弹簧提供有限刚度（消奇异 + 保留内外翻/内旋自由度）。
if TAL_GROUND:
    _allg = sorted({_g for _gd, _dl in TAL_GROUND.values() for _g in _gd})
    print(f"[talus-spring] 距骨横向硬约束**已跳过**；改为固定 {len(_allg)} 个哑节点")
    w('    <bc name="fix_ground" type="zero displacement" '
      'node_set="tal_ground_ns">')
    w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
    w('    </bc>')
tal_dofs = ("xyz" if FIX_TALUS == "all" else
            "".join(d for d, v in zip("xyz", np.abs(axis))
                    if v < np.abs(axis).max() - 1e-6))
if not tal_dofs:
    tal_dofs = "xz"
if not TAL_GROUND:
    if TALUS_RIGID:
        # ★ 刚体模式下不能再用节点级 BC（这些节点已无自由度）⇒ 用 rigid BC 挂到刚体
        w(f'    <bc name="rigid_talus" type="rigid" node_set="talus_rigid_ns">')
        w(f'      <rb>{RIGID_TALUS_MAT}</rb>')
        w('    </bc>')
        w(f'    <bc name="rigid_ground" type="rigid" node_set="ground_rigid_ns">')
        w(f'      <rb>{RIGID_GROUND_MAT}</rb>')
        w('    </bc>')
    else:
        w(f'    <bc name="fix_talus" type="zero displacement" node_set="talus_low_ns">')
        for d in tal_dofs:
            w(f'      <{d}_dof>1</{d}_dof>')
        w('    </bc>')
# 胫骨/腓骨顶面：约束横向，轴向自由（载荷驱动）
kax = int(np.argmax(np.abs(axis)))
for ns, nm, _faces in (("tibia_top_ns", "tibia", top_tibia),
                       ("fibula_top_ns", "fibula", top_fibula)):
    if not _faces:      # 空 NodeSet 不能写 BC（FEBio: invalid value for node_set）
        continue
    if DRIVER == "force" and TOPFIX == "none":
        continue        # ★ 顶面完全不约束：靠接触 + 韧带定运动学（允许外翻/内旋）
    w(f'    <bc name="{nm}_lateral" type="zero displacement" node_set="{ns}">')
    for di, dn in enumerate("xyz"):
        if di != kax:
            w(f'      <{dn}_dof>1</{dn}_dof>')
    w('    </bc>')
if LIG_ONLY_NS:
    # ★ FEBio 壳的附加自由度叫 back-displacement，dof 名是 sx/sy/sz，
    #   不是 <x_rotation> / <u_rotation>（那两个标签 FEBio 不认）。
    w('    <bc name="lig_back_lock" type="zero shell displacement" '
      f'node_set="lig_only_ns">')
    for _d in ("sx", "sy", "sz"):
        w(f'      <dof>{_d}</dof><value lc="1">0</value><relative>0</relative>')
    w('    </bc>')
    w('    <bc name="fix_navicular" type="zero displacement" '
      'node_set="navicular_all_ns">')
    w('      <x_dof>1</x_dof><y_dof>1</y_dof><z_dof>1</z_dof>')
    w('    </bc>')
if DRIVER in ("disp", "dispx") and not MULTISTEP:   # 多步时顶面位移移到 step1
    # ★ 顶面三个自由度**全部**指定位移（= 夹持端）：刚体模态消失，
    #   不依赖接触/韧带稳定（LRN-017 的核心修法）。
    _kax = int(np.argmax(np.abs(LAX)))
    # ★ 姿势=加载方向倾斜（见文件头 POSE_TILT 说明）：把驱动方向相对加载轴倾斜，
    #   等效于给关节施加"压缩 + 剪切/力矩"，全程仍是位移控制（落在稳定区）。
    if POSE_TILT or POSE_TILT_ML:
        _vec = -DISP_MM * (LAX - np.tan(np.radians(POSE_TILT)) * _ml
                           - np.tan(np.radians(POSE_TILT_ML)) * _ap)
        print(f"[pose] 加载方向倾斜：POSE_TILT={POSE_TILT}deg(内外侧) "
              f"POSE_TILT_ML={POSE_TILT_ML}deg(前后) ⇒ 驱动位移 = {np.round(_vec,5)} mm "
              f"(轴向分量 {_vec[_kax]:.5f})")
        for _ns, _nm, _faces in (("tibia_top_ns", "tibia", top_tibia),
                                 ("fibula_top_ns", "fibula", top_fibula)):
            if not _faces:
                continue
            for _di, _dn in enumerate("xyz"):
                w(f'    <bc name="{_nm}_top_{_dn}" type="prescribed displacement" '
                  f'node_set="{_ns}">')
                w(f'      <dof>{_dn}</dof><value lc="1">{_vec[_di]:.6f}</value>'
                  '<relative>0</relative>')
                w('    </bc>')
    for _ns, _nm, _faces in (("tibia_top_ns", "tibia", top_tibia),
                             ("fibula_top_ns", "fibula", top_fibula)):
        if not _faces:
            continue
        if POSE_TILT or POSE_TILT_ML:
            continue      # 倾斜模式已在上方写完全部 3 个 dof（避免重复 BC）
        _dofs = "xyz" if DRIVER == "disp" else "xyz"[_kax]
        for _di, _dn in enumerate("xyz"):
            if DRIVER == "dispx" and _di != _kax:
                continue     # 横向自由：交给韧带 + 接触约束
            _v = (-DISP_MM if LAX[_di] > 0 else DISP_MM) if _di == _kax else 0.0
            w(f'    <bc name="{_nm}_top_{_dn}" type="prescribed displacement" '
              f'node_set="{_ns}">')
            w(f'      <dof>{_dn}</dof><value lc="1">{_v:.6f}</value>'
              '<relative>0</relative>')
            w('    </bc>')
w('  </Boundary>')
# ★★ 路线 A：<Rigid> 段在这里刷出（**必须在 </Boundary> 之后**，与已跑通的最小算例同序）
lines.extend(_RIGID_LINES)

if (not NO_CONTACT) or TIE_SEC:
  w('  <Contact>')
  # ★ 顺序与 build_feb.py 对齐：**tie 先写，关节接触后写**（实测差异点之一）
  for j in TIE_SEC:
    w(f'    <contact name="tie_{j}" surface_pair="tiepair_{j}" '
      'type="tied-elastic">')
    w(f'      <penalty>{TIE_PEN}</penalty>')
    w(f'      <auto_penalty>{TIE_AUTO}</auto_penalty>')
    w('      <two_pass>0</two_pass>')
    # ★ 必须写 <laugon>PENALTY</laugon>（与 build_feb.py 一致）：漏掉时
    #   实测 t=0 即报负 Jacobian，且与位移大小无关（0.001mm 也失败）
    w('      <laugon>PENALTY</laugon>')
    w(f'      <tolerance>{TIE_TOL}</tolerance>')
    w('      <symmetric_stiffness>0</symmetric_stiffness>')
    w(f'      <search_tol>{TIE_STOL}</search_tol>')
    w(f'      <search_radius>{TIE_SRAD}</search_radius>')
    w('    </contact>')
  if not NO_CONTACT:
    # ★ A1：关节接触对统一写出。
    #   · tibiotalar 的参数与旧实现**逐字一致**（LAUGON 原样透传、auto_penalty
    #     同样的精确比较）⇒ 基线 .feb 不变、零回归。
    #   · subtalar 默认 AUGLAG（细网格接触 PENALTY 会失败的项目铁律），
    #     `SUBTALAR_LAUGON` 可覆盖。
    for _pn, _pa, _pb in JOINT_PAIRS:
        _nm = _pn[:-5]
        if _nm == "subtalar":
            _lg = (os.environ.get("SUBTALAR_LAUGON", "AUGLAG").strip()
                   or "AUGLAG").upper()
        else:
            _lg = LAUGON
        w(f'    <contact name="{_nm}" surface_pair="{_pn}" '
          'type="sliding-elastic">')
        w(f'      <laugon>{_lg}</laugon>')
        w(f'      <penalty>{PENALTY}</penalty>')
        w(f'      <auto_penalty>{0 if _lg == "AUGLAG" else 1}</auto_penalty>')
        w('      <two_pass>0</two_pass>')
        w(f'      <node_reloc>{NODE_RELOC}</node_reloc>')
        w('      <symmetric_stiffness>0</symmetric_stiffness>')
        w('      <tolerance>0.02</tolerance>')
        w('    </contact>')
  w('  </Contact>')
_TAG_EARLY = MATSET + ("_" + "_".join(ONLY) if ONLY else "")
w('  <Output>')
w('    <plotfile type="febio"><var type="displacement"/>'
  '<var type="stress"/><var type="reaction forces"/></plotfile>')
if LOG_DATA:
    # 手册 §3.19.1：file= 写到独立文本；delim= 定分隔符；item 列表省略=全部节点
    w(f'    <logfile file="log_{_TAG_EARLY}.txt">')
    w(f'      <node_data data="ux;uy;uz" name="disp" '
      f'file="disp_{_TAG_EARLY}.txt" delim=","/>')
    for _ns, _nm in (("tibia_top_ns", "tibia"), ("fibula_top_ns", "fibula")):
        _faces = top_tibia if _nm == "tibia" else top_fibula
        if not _faces:
            continue
        _ids = ",".join(str(x) for x in ns_of(_faces))
        w(f'      <node_data data="Rx;Ry;Rz" name="reac_{_nm}" '
          f'file="reac_{_nm}_{_TAG_EARLY}.txt" delim=",">{_ids}</node_data>')
    w(f'      <node_data data="Rx;Ry;Rz" name="reac_talus" '
      f'file="reac_talus_{_TAG_EARLY}.txt" delim=",">'
      f'{",".join(str(x) for x in ns_of(tal_faces))}</node_data>')
    if LOG_CONTACT and not NO_CONTACT:
        # 逐面片接触间隙/压力 ⇒ × 面片面积 可积分出接触力（A1 判据③）
        _csurfs = [x for x in os.environ.get("LOG_CONTACT_SURF", "").split(",")
                   if x] or [f"surf_{_pb}" for _pn, _pa, _pb in JOINT_PAIRS]
        for _sn in _csurfs:
            w(f'      <face_data data="contact gap;contact pressure" '
              f'surface="{_sn}" file="cface_{_sn}_{_TAG_EARLY}.txt" delim=","/>')
    w('    </logfile>')
w('  </Output>')
if not MULTISTEP:
    w('  <Control>')
    w('    <analysis>STATIC</analysis><time_steps>{}</time_steps>'
      f'<step_size>{1.0/TIME_STEPS:.6f}</step_size>'.format(TIME_STEPS))
    w('    <solver type="solid">')
    # ★★ `symmetric_stiffness=0`（**必须是 <solver> 的子标签**，放 Control 直下会
    #   `tag "symmetric_stiffness" : unrecognized tag`）—— 刚体连接件的硬要求，
    #   FEBio 自己在日志里警告："Rigid connectors require non-symmetric stiffness matrix."
    #   用默认对称化矩阵 ⇒ 弹簧约束力算错 ⇒ 实测 step2 立刻 47 个负 Jacobian。
    #   （顺带也是 sliding-elastic 接触推荐的设置。）
    w('      <symmetric_stiffness>0</symmetric_stiffness>')
    w('      <linear_solver type="pardiso"/>')
    w(f'      <max_refs>{MAX_REFS}</max_refs>')
    w('      <diverge_reform>1</diverge_reform>')
    w('      <reform_each_time_step>1</reform_each_time_step>')
    w(f'      <dtol>{DTOL}</dtol>')
    w(f'      <etol>{ETOL}</etol>')
    w('      <rtol>0</rtol>')
    w('      <lstol>0.9</lstol>')
    w('      <min_residual>1e-20</min_residual>')
    w('      <rhoi>-2</rhoi>')
    w(f'      <qn_method type="BFGS"><max_ups>{MAX_UPS}</max_ups></qn_method>')
    w('    </solver>')
    w('    <time_stepper>')
    w('      <dtmin>0.0001</dtmin><dtmax>0.05</dtmax>')
    w('      <max_retries>10</max_retries><opt_iter>10</opt_iter>')
    w('    </time_stepper>')
    w('  </Control>')
# ★★★ 多步落位：`<Step>` 放在文件底部（手册 §6.1：「Control 不再在顶部，Steps 加到文件底部」）
#   step1 = 位移落位（把关节从 ~1 mm 初始间隙压到"最小接触"，对齐 Anderson）
#   step2 = 力控制加载（step1 的位移 BC 只在该 step 生效 ⇒ 此步轴向自由，力接管）
_RIGID_STEP1 = [ln for ln in _RIGID_LINES]      # 落位期沿用（转轴锁定按 TAL_ROTLOCK）
if MULTISTEP:
    def _ctl(ts, extra=""):
        w('      <Control>')
        w(f'        <analysis>STATIC</analysis><time_steps>{ts}</time_steps>'
          f'<step_size>{1.0/ts:.6f}</step_size>')
        w('        <solver type="solid">')
        w('          <symmetric_stiffness>0</symmetric_stiffness>')
        w('          <linear_solver type="pardiso"/>')
        w(f'          <max_refs>{MAX_REFS}</max_refs>')
        w('          <diverge_reform>1</diverge_reform>')
        w('          <reform_each_time_step>1</reform_each_time_step>')
        w(f'          <dtol>{DTOL}</dtol>')
        w(f'          <etol>{ETOL}</etol>')
        w('          <rtol>0</rtol><lstol>0.9</lstol>')
        w('          <min_residual>1e-20</min_residual><rhoi>-2</rhoi>')
        w(f'          <qn_method type="BFGS"><max_ups>{MAX_UPS}</max_ups></qn_method>')
        w('        </solver>')
        w('        <time_stepper>')
        w('          <dtmin>0.0001</dtmin><dtmax>0.05</dtmax>')
        w('          <max_retries>10</max_retries><opt_iter>10</opt_iter>')
        w('        </time_stepper>')
        w('      </Control>')

    def _topdisp():
        """step1 的落位位移：**只 prescribe 轴向那一个 dof**（压向距骨）。
        ★ 不能把横向也写 0：全局 <Boundary> 已有横向 zero-displacement，
          同一节点同一 dof 两条 BC 重复施加 ⇒ 实测 step1"收敛"却在 step2 起点
          留下 42 个翻负单元（`dispx` 只写轴向，所以干净）。"""
        for _nm, _ns, _faces in (("tibia_top_ns", "tibia", top_tibia),
                                 ("fibula_top_ns", "fibula", top_fibula)):
            if not _faces:
                continue
            # ★★ 落位方向必须 = **压缩方向**。
            #   实测（`[LRN-20261005-072]`）：`dispx` 下 DISP_MM=-1.6 ⇒ 约束位移 +1.6(+z)
            #   ⇒ NORMAL 且 Rz=-1017 N（压缩）；而 DISP_MM=+1.6 ⇒ 约束位移 -1.6
            #   ⇒ 不收敛且 Rz=+555 N（**拉伸**，把关节拉开）。
            #   ⇒ 压缩方向 = +LAX（当 LAX[kax] > 0）。原实现用 -SEAT_MM ⇒ **落位在拉关节**。
            _dn = "xyz"[kax]
            _v = (SEAT_MM if LAX[kax] > 0 else -SEAT_MM)
            w(f'      <bc name="seat_{_nm}" '
              f'type="prescribed displacement" node_set="{_nm}">')
            w(f'        <dof>{_dn}</dof><value lc="1">{_v:.6f}</value>'
              '<relative>0</relative>')
            w('      </bc>')

    def _step_loads():
        """step2 的力控制：沿加载轴施加 FORCE_N（按顶面投影面积分到胫/腓骨）"""
        _Ft, _Ff = FORCE_N * (1 - share), FORCE_N * share
        for _surf, _F, _faces in (("tibia_top", _Ft, top_tibia),
                                  ("fibula_top", _Ff, top_fibula)):
            if not _faces or _F <= 0:
                continue
            _a = area_proj(top_tibia if _surf == "tibia_top" else top_fibula)
            _tr = _F / _a if _a > 0 else 0.0
            w(f'      <surface_load type="traction" surface="{_surf}">')
            # ★★★ 载荷曲线必须挂在 <scale> 上，**不能挂在 <traction> 上**
            #   （手册 §3.13.2.2 原文："{scale lc=1}1.0{/scale}" + "{traction}0,0,1{/traction}"，
            #    "An optional load curve can be defined for the **scale** element"）。
            #   挂错位置会被**静默忽略** ⇒ 满载荷在第一增量瞬加 ⇒ 42 个负 Jacobian。
            #   `[ERR-20261005-008]`
            w('        <scale lc="2">1.0</scale>')
            w(f'        <traction>{LOOKSIGNLAX[0]*_tr:.6f},'
              f'{LOOKSIGNLAX[1]*_tr:.6f},'
              f'{LOOKSIGNLAX[2]*_tr:.6f}</traction>')
            w('      </surface_load>')

    print(f"[multistep] step1 位移落位 SEAT_MM={SEAT_MM}（{TIME_STEPS_SEAT} 步）"
          f" → step2 力控制 FORCE_N={FORCE_N} N（{TIME_STEPS} 步）")
    w('  <Step>')
    w('    <step id="1">')
    _ctl(TIME_STEPS_SEAT)
    w('      <Boundary>')
    _topdisp()
    w('      </Boundary>')
    w('    </step>')
    w('    <step id="2">')
    _ctl(TIME_STEPS)
    if SEAT_HOLD:
        w('      <Boundary>')
        _topdisp()
        w('      </Boundary>')
    w('      <Loads>')
    _step_loads()
    w('      </Loads>')
    w('    </step>')
    w('  </Step>')
w('</febio_spec>')

_tag = MATSET + ("_" + "_".join(ONLY) if ONLY else "")
# ── 导出读数所需数据（供 read_lig.py 用）──
_dump = {}
_dump["n_node"] = np.array(len(coords))
_fe = np.zeros(len(coords) + 1, dtype=np.int64)
for _k, _v in node_id_of.items():
    if isinstance(_k, int):
        _fe[_v] = _k
_dump["fe_to_thums"] = _fe
_dump["axis"] = axis
_names, _eas, _pairs, _npair = [], [], [], []
for r in lig_recs:
    pairs = r.get("trusses") or []
    fp = [[node_id_of[int(a)], node_id_of[int(b)]] for a, b in pairs]
    _names.append(r["name"])
    _eas.append(r["ea"])
    _pairs.append(fp)
    _npair.append(len(fp))
_mx = max(_npair) if _npair else 0
_TR = -np.ones((len(_names), max(_mx, 1), 2), dtype=np.int64)
for i, fp in enumerate(_pairs):
    for j, pr in enumerate(fp):
        _TR[i, j] = pr
_dump["lig_names"] = np.array(_names)
_dump["lig_ea"] = np.array(_eas)
_dump["lig_pairs"] = _TR
_dump["lig_npair"] = np.array(_npair)
np.savez(WORK / f"lig_read_{_tag}.npz", **_dump)
print(f"读数数据 -> {WORK / f'lig_read_{_tag}.npz'}")

feb = WORK / f"ankle_lig_{_tag}.feb"
feb.write_text("\n".join(lines), encoding="latin-1")
print(f"\n写出 {feb}（{feb.stat().st_size/1024:.0f} KB）"
      f"  单元：骨+软骨 {LIG_EID0}，韧带 {LIG_EID1-LIG_EID0}，共 {LIG_EID1}")

# ★ 只建不解（默认 0 = 照旧）。
#   用途：几何/姿势检查（`check_c2.py`）只需参考构型，不需要求解 ⇒ 省掉 FEBio 时间。
#   注意：`NO_RUN=1` 时不写 FEBio 日志，`acceptance.py` 会读到上一次的旧日志 ⇒
#   不要在验收流程里用它。
if os.environ.get("NO_RUN", "0") == "1":
    print("[no-run] NO_RUN=1 ⇒ 只写出 .feb，不调用 FEBio")
    sys.exit(0)

if not FEBIO.exists():
    print(f"❌ 找不到 FEBio: {FEBIO}")
    sys.exit(1)
t0 = time.time()
r = subprocess.run([str(FEBIO), "-i", str(feb)], cwd=str(WORK),
                   capture_output=True, text=True, timeout=3600)
print(f"\nFEBio rc={r.returncode}  用时 {time.time()-t0:.1f}s")
out = (r.stdout or "") + (r.stderr or "")
(WORK / f"febio_{_tag}.log").write_text(out, encoding="utf-8")
for l in out.strip().split("\n")[-45:]:
    print("  " + l)
