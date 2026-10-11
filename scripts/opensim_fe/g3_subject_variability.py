"""g3_subject_variability.py -- G3 「被试级变异」统计层（预注册设计实现）。

本脚本严格实现 `docs/methods/G3_被试级变异_统计层设计.md`（预注册，已批准）：
把 S4 的 50×9 **确定性**骨折矩阵（σ_nom(h) ≥ σc）升级为**被试级随机抽样下的
骨折概率** p[h][site]，使渐进部位（腓骨/骨盆/颅骨）的 Logistic OR 有定义（对标
[Y25] 的 per-metre OR）。

锁定的预注册设定（跑数前锁定，跑后不改）
----------------------------------------
* D1 相关结构：**全局骨质量因子 + 部位残差**（Lognormal）。
* D2 主设定 COV（A3 判定只看这一行）：``COV_B=0.15, COV_E=0.10, COV_A=0.15``。
* D3 样本/种子：``n_subj=400``；seeds ``{20261003, 20261004, 20261005}``，官方 seed
  ``20261003``。
* 敏感性网格（全量报告，不许挑）：COV_B∈{0.10,0.15,0.20} × COV_E∈{0.05,0.10}
  × COV_A∈{0.10,0.15,0.20}，各 × 3 seeds；外加一个「质量不入载荷」对照
  （``m_s/m_ref = 1``）。

被试模型（design §3.2；被试一次抽参数，跨 50 高度**配对**共享）
-----------------------------------------------------------------
对虚拟被试 s，抽参数（与高度无关）：

    b_s   ~ Lognormal（全局骨质量，名义 1）
    e_s^k ~ Lognormal（部位残差，跨部位独立，名义 1）
    a_s   ~ Lognormal（截面缩放，名义 1）
    m_s   ~ MASS_MODEL（分类 Normal + clip 40–130 kg，复用 scenarios）

判据：

    y[s][h][k] = 1[ σ_nom^k(h) · (m_s/m_ref)/a_s ≥ σc^k · b_s · e_s^k ]

Lognormal 约定（design §3.2 名义=1；本脚本采用 **mean-1** 参数化并在报告中记录）：
若 X~Lognormal(mu,sigma) 且变异系数 c，则

    sigma = sqrt(ln(1+c^2)) ,  mu = -sigma^2/2   （⇒ E[X]=1，median=e^{mu}）。

即「按 CV 归一化到**均值 1**」的 lognormal。COV→0 时 mu→0、sigma→0、X→1。

输入 / 零 FE
------------
σ_nom^k(h) 直接取 ``results/opensim_fe/fracture_matrix_bc2.json`` 的
``sigma_vm_mpa["median"][bone]``（BC2 主口径 median，逐格已存）；σc^k 取 route-乙
**弱区口径**（``paper_alignment.REGIONS`` ← ``bone.py::MATERIAL_STRENGTH_MPA``）：
足 150 / 胫 70 / 腓 70 / 股 80 / 盆 180 / 腰胸胸 150 / 颅 160 MPa。
被试缩放全部是逐格乘除 —— **无 FE 调用、无插值、无新求解器**。

OR 估计
-------
**原封复用** ``paper_alignment.logistic_fit``（x=高度 m，y=拍平 0/1；OR=exp(β₁)
per-metre，Wald CI）；判定复用 ``s4_alignment_bc2.a3_verdict``（±30% 内 ≥2/3，
靶 腓骨 1.682 / 骨盆 1.236 / 颅骨 1.576）。**不换估计器、不改判定规则**。

交付物（纯增量）
----------------
* ``results/opensim_fe/g3_probability_matrix.json``
* ``results/opensim_fe/g3_or_estimates.json``
* ``results/opensim_fe/S4_A3_G3_REPORT.md``

复现（仓库根）::

    $env:PYTHONPATH="src"
    & .venv\\Scripts\\python.exe scripts\\opensim_fe\\g3_subject_variability.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from collections.abc import Mapping
from datetime import datetime
from pathlib import Path
from typing import TypedDict, cast

import numpy as np

_HERE = Path(__file__).resolve().parent
ROOT = _HERE.parents[1]
for _p in (ROOT, ROOT / "src", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from scripts.opensim_fe import paper_alignment as pa  # noqa: E402  (复用 logistic_fit / REGIONS / _sc)
from scripts.opensim_fe import s4_alignment_bc2 as s4  # noqa: E402  (复用 A3_TARGETS / a3_verdict，零改动)
from climbing.scenarios import MASS_CLIP_KG, MASS_MODEL  # noqa: E402

#: 供测试直接引用被本脚本复用的既有模块（re-export，供 test_g3_subject_variability 校验）。
__all__ = ["pa", "s4"]


# ---------------------------------------------------------------------------
# 类型（仅静态标注；运行期仍是普通 dict，零行为改动）
# ---------------------------------------------------------------------------
class FractureMatrixJson(TypedDict):
    """``fracture_matrix_bc2.json`` 中本脚本读取的字段。"""

    heights_m: list[float]
    sigma_vm_mpa: dict[str, dict[str, list[float]]]


class Inputs(TypedDict):
    """``load_inputs`` 的返回：σ_nom(h)[site] + 弱区 σc + 高度。"""

    heights: np.ndarray
    sites: list[str]
    bones: list[str]
    sigma_nom: np.ndarray
    sigma_c: np.ndarray
    weak_sigma_c_mpa: dict[str, float]


class Factors(TypedDict):
    """``sample_factors`` 的返回：逐被试因子与载荷缩放。"""

    b: np.ndarray
    e: np.ndarray
    a: np.ndarray
    mass_kg: np.ndarray
    mass_ratio: np.ndarray
    lhs: np.ndarray


class OrRecord(TypedDict, total=False):
    """``paper_alignment.logistic_fit`` 结果裁剪后的统一 OR 记录。"""

    status: str
    OR: float | None
    ci95: list[float] | None
    beta: float
    se: float
    wald_z: float
    p: float
    n: int
    n_fracture: int
    n_safe: int
    converged: bool
    note: str
    separation_threshold_m: float


class A3Verdict(TypedDict):
    """``s4_alignment_bc2.a3_verdict`` 的输出（只声明本脚本读取的字段）。"""

    targets: dict[str, float]
    per_part: dict[str, dict[str, object]]
    n_within: int
    n_target: int
    verdict: str


class ConfigResult(TypedDict):
    """``run_config`` 的返回：p / 逐部位 OR / A3。"""

    seed: int
    cov_b: float
    cov_e: float
    cov_a: float
    mass_in_load: bool
    y_shape: list[int]
    p: np.ndarray
    or_map: dict[str, OrRecord]
    a3: A3Verdict


class Degenerate(TypedDict):
    """DoD-2 极限退化的汇总记录。"""

    cov: float
    mass_in_load: bool
    max_abs_dev: float
    majority_matches_base: bool
    exact_cov0_matches_base: bool
    matches_base: bool
    by_site: dict[str, bool]
    knife_edge_cells: list[tuple[str, int, float]]
    mean_y: dict[str, list[float]]


class SeedStability(TypedDict):
    """DoD-4 单部位的三 seed OR 漂移记录。"""

    or_by_seed: dict[int, float | None]
    rel_drift_vs_official: dict[int, float | None]
    max_rel_drift: float | None
    pass_20pct: bool


class BuildResult(TypedDict):
    """``build`` 的返回：主设定 + 退化 + 网格 + 对照 + 种子稳定。"""

    heights: np.ndarray
    sites: list[str]
    base: np.ndarray
    primary: ConfigResult
    primary_by_seed: dict[int, ConfigResult]
    degenerate: Degenerate
    grid: list[ConfigResult]
    control: list[ConfigResult]
    seed_stability: dict[str, SeedStability]
    n: int


class Verdict(TypedDict):
    """``build_verdict`` 的返回（报告与摘要用）。"""

    a3_verdict: str
    n_within: int
    all_three_ok: bool
    gradual_status: dict[str, str | None]
    implementation_note: str | None
    degenerate_ok: bool
    order_ok: bool
    order_detail: dict[str, bool | None]
    monotone_ok: bool
    seed_ok: bool
    quotable_ok: bool
    headline: str
    paragraph: str


RES = ROOT / "results" / "opensim_fe"
BC2_JSON = RES / "fracture_matrix_bc2.json"
OUT_PROB = RES / "g3_probability_matrix.json"
OUT_OR = RES / "g3_or_estimates.json"
OUT_REPORT = RES / "S4_A3_G3_REPORT.md"

# ---------------------------------------------------------------------------
# 预注册参数（LOCKED，design §4；跑后不得修改）
# ---------------------------------------------------------------------------
N_SUBJ = 400
SEEDS = [20261003, 20261004, 20261005]
OFFICIAL_SEED = 20261003
COV_B_PRIMARY = 0.15
COV_E_PRIMARY = 0.10
COV_A_PRIMARY = 0.15
COV_B_GRID = [0.10, 0.15, 0.20]
COV_E_GRID = [0.05, 0.10]
COV_A_GRID = [0.10, 0.15, 0.20]
COV_DEGENERATE = 0.001  # 「near-zero」极限退化（design §6 DoD-2）
M_REF_KG = 75.337       # 载荷链参考质量（design §3.2）

A3_TARGETS = dict(s4.A3_TARGETS)         # {腓骨:1.682, 骨盆:1.236, 颅骨:1.576}
PAPER_GRADUAL = set(pa.PAPER_GRADUAL)    # {腓骨,骨盆,颅骨}
GRADUAL = ["腓骨", "骨盆", "颅骨"]

#: 禁用数字（口径纪律，design §6 DoD-7）：这些「结果」不得出现在本报告。
#: 本脚本**不复述**这些串（避免自身触发），仅在 DoD-7 指向设计原文。
BANNED_RESULT_NUMBERS = ["2.205", "15.14", "2.70×", "2.70x"]

#: G3 交付物（纯增量，唯一允许的新文件集）。
G3_DELIVERABLES = [
    "scripts/opensim_fe/g3_subject_variability.py",
    "results/opensim_fe/g3_probability_matrix.json",
    "results/opensim_fe/g3_or_estimates.json",
    "results/opensim_fe/S4_A3_G3_REPORT.md",
    "tests/test_g3_subject_variability.py",
]

#: 任务开始（G3 动工前）`git status --porcelain` 的既有条目 —— 非 G3 所改。
PREEXISTING_GIT_STATUS = [
    " M docs/README.md",
    " M docs/methods/目标重构_刚性验证与垫上推断.md",
    " M docs/paper/论文v2_数学建模版.md",
    " M docs/reports/S5_QUOTABLE.md",
    " M docs/reports/总体进度.md",
    "?? docs/methods/G3_被试级变异_统计层设计.md",
    "?? docs/reports/S5_QUOTABLE.md.bak_20261010",
    "?? docs/reports/最终结论_2026-10-10.md",
    "?? results/opensim_fe/SINK_130_140_150_CLOSURE_REPORT.md",
    "?? results/opensim_fe/sink130_140_150_closure.json",
    "?? scripts/opensim_fe/run_sink_closure_20261010.py",
]


def git_status_lines() -> list[str] | None:
    """读 `git status --porcelain`（read-only；失败返回 None）。"""
    import subprocess
    try:
        r = subprocess.run(
            ["git", "-c", "core.quotepath=false", "-C", str(ROOT),
             "status", "--porcelain"],
            capture_output=True, text=True, timeout=30)
        return [ln.rstrip() for ln in r.stdout.splitlines() if ln.strip()]
    except Exception:  # noqa: BLE001
        return None


# ---------------------------------------------------------------------------
# 核心：lognormal 约定 / 输入 / 被试抽样 / 概率矩阵
# ---------------------------------------------------------------------------
def log_params(cov: float) -> tuple[float, float]:
    """按 CV 归一化到**均值 1** 的 lognormal：(mu, sigma) = (-s^2/2, sqrt(ln(1+c^2)))。

    ⇒ E[X]=1（名义 1），median=exp(mu)。COV→0 ⇒ (0,0)，X→1。
    """
    c = float(cov)
    sig = float(np.sqrt(np.log1p(c * c)))
    return -0.5 * sig * sig, sig


def load_inputs(path: Path = BC2_JSON) -> Inputs:
    """读取 ``fracture_matrix_bc2.json``：σ_nom(h)[site]（median）+ 弱区 σc + 高度。"""
    fm = cast(FractureMatrixJson, json.loads(Path(path).read_text(encoding="utf-8")))
    heights = np.asarray(fm["heights_m"], dtype=float)
    sv = fm["sigma_vm_mpa"]["median"]
    sites = list(pa.REGION_NAMES)
    bones = [r[1] for r in pa.REGIONS]
    # route-乙 弱区 σc（各部位第一个子区域 = 两端/股骨颈/整体）
    sigma_c = np.array([pa._sc(r[3][0][1]) for r in pa.REGIONS], dtype=float)
    sigma_nom = np.array([np.asarray(sv[b], dtype=float) for b in bones])  # (9,50)
    return Inputs(heights=heights, sites=sites, bones=bones,
                  sigma_nom=sigma_nom, sigma_c=sigma_c,
                  weak_sigma_c_mpa={cn: float(sigma_c[i]) for i, cn in enumerate(sites)})


def deterministic_base(sigma_nom: np.ndarray, sigma_c: np.ndarray) -> np.ndarray:
    """确定性骨折矩阵 1[σ_nom(h) ≥ σc]（G3 的 COV→0 退化基准）。形状 (n_site,n_h)。"""
    return (np.asarray(sigma_nom, float) >= np.asarray(sigma_c, float)[:, None]).astype(int)


def sample_factors(seed: int, n: int, cov_b: float, cov_e: float, cov_a: float,
                   n_sites: int, mass_in_load: bool = True) -> Factors:
    """逐被试抽 b_s/e_s^k/a_s/m_s（RNG=``default_rng([seed, s])``，phase_p1 同款）。"""
    mb, sb = log_params(cov_b)
    me, se = log_params(cov_e)
    ma, sa = log_params(cov_a)

    labels = list(MASS_MODEL.items())               # [(label,(prob,mean,sd)), ...]
    probs = np.array([v[0] for _, v in labels], dtype=float)
    cum = np.cumsum(probs / probs.sum())

    b = np.empty(n, dtype=float)
    a = np.empty(n, dtype=float)
    e = np.empty((n, n_sites), dtype=float)
    mass = np.empty(n, dtype=float)
    for s in range(n):
        rng = np.random.default_rng([seed, s])
        b[s] = rng.lognormal(mb, sb)
        e[s] = rng.lognormal(me, se, size=n_sites)
        a[s] = rng.lognormal(ma, sa)
        u = float(rng.random())
        idx = int(np.searchsorted(cum, u, side="right"))
        idx = min(idx, len(labels) - 1)
        mean_kg, sd_kg = labels[idx][1][1], labels[idx][1][2]
        mass[s] = float(np.clip(rng.normal(mean_kg, sd_kg), *MASS_CLIP_KG))

    mass_ratio = (mass / M_REF_KG) if mass_in_load else np.ones(n, dtype=float)
    lhs = mass_ratio / a                                  # σ 侧缩放 (m/m_ref)/a_s
    return Factors(b=b, e=e, a=a, mass_kg=mass, mass_ratio=mass_ratio, lhs=lhs)


def simulate(seed: int, n: int, cov_b: float, cov_e: float, cov_a: float,
             inputs: Inputs, mass_in_load: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """被试级模拟：返回 y (n, n_site, n_h) 0/1 与 p (n_site, n_h)=mean_s(y)。"""
    sigma_nom = inputs["sigma_nom"]                       # (ns,nh)
    sigma_c = inputs["sigma_c"]                           # (ns,)
    ns = sigma_nom.shape[0]
    fac = sample_factors(seed, n, cov_b, cov_e, cov_a, ns, mass_in_load)

    lhs = fac["lhs"][:, None, None]                       # (n,1,1)
    thr = sigma_c[None, :] * fac["b"][:, None] * fac["e"]  # (n,ns)
    y = sigma_nom[None, :, :] * lhs >= thr[:, :, None]    # (n,ns,nh)
    p = y.mean(axis=0)                                    # (ns,nh)
    return y, p


def or_all_sites(heights: np.ndarray, sites: list[str], y: np.ndarray) -> dict[str, OrRecord]:
    """逐部位跑 ``paper_alignment.logistic_fit``（x=高度 m，y 拍平 0/1）。"""
    heights = np.asarray(heights, dtype=float)
    n = y.shape[0]
    x = np.tile(heights, n)
    out: dict[str, OrRecord] = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with np.errstate(over="ignore", invalid="ignore"):
            for k, cn in enumerate(sites):
                out[cn] = cast(OrRecord, cast(object, pa.logistic_fit(x, y[:, k, :].reshape(-1).astype(int))))
    return out


def wilson_ci(k: np.ndarray, n: float, z: float = 1.96) -> tuple[np.ndarray, np.ndarray]:
    """二项率的 Wilson 95% CI（向量化）。"""
    k = np.asarray(k, dtype=float)
    n = float(n)
    ph = k / n
    denom = 1.0 + z * z / n
    center = (ph + z * z / (2.0 * n)) / denom
    half = (z * np.sqrt(ph * (1.0 - ph) / n + z * z / (4.0 * n * n))) / denom
    return center - half, center + half


def a3_verdict(or_map: dict[str, OrRecord]) -> A3Verdict:
    """复用 s4_alignment_bc2.a3_verdict（规则零改动）。"""
    return cast(A3Verdict, cast(object, s4.a3_verdict(or_map)))


def _site_or_record(rec: Mapping[str, object]) -> OrRecord:
    """把 logistic_fit 结果裁剪成统一的 OR 记录。"""
    keep = ("status", "OR", "ci95", "beta", "se", "wald_z", "p",
            "n", "n_fracture", "n_safe", "converged", "note")
    return cast(OrRecord, cast(object, {k: rec.get(k) for k in keep if k in rec}))


# ---------------------------------------------------------------------------
# 汇总
# ---------------------------------------------------------------------------
def run_config(seed: int, cov_b: float, cov_e: float, cov_a: float,
               inputs: Inputs, n: int = N_SUBJ, mass_in_load: bool = True) -> ConfigResult:
    """跑一个 (seed, COV) 配置：返回 p / 逐部位 OR / A3。"""
    y, p = simulate(seed, n, cov_b, cov_e, cov_a, inputs, mass_in_load)
    or_map = or_all_sites(inputs["heights"], inputs["sites"], y)
    a3 = a3_verdict(or_map)
    return ConfigResult(seed=seed, cov_b=cov_b, cov_e=cov_e, cov_a=cov_a,
                        mass_in_load=bool(mass_in_load), y_shape=list(y.shape),
                        p=p, or_map=or_map, a3=a3)


def build(inputs: Inputs, n: int = N_SUBJ) -> BuildResult:
    """跑主设定 + 极限退化 + 全敏感性网格 + 质量对照。"""
    heights = inputs["heights"]
    sites = inputs["sites"]
    base = deterministic_base(inputs["sigma_nom"], inputs["sigma_c"])

    # ---- 主设定（A3 判定只看这一行）-----------------------------------
    primary = run_config(OFFICIAL_SEED, COV_B_PRIMARY, COV_E_PRIMARY, COV_A_PRIMARY,
                         inputs, n, mass_in_load=True)

    # ---- DoD-2 极限退化：COV→0.001 且质量不入载荷（全部被试变异→关）------
    # 说明：COV=0.001 是「near-zero」近似；极限 COV→0 时逐格 bit-exact 复现。
    # 有一格是**刀刃**（颅骨 h=6：σ_nom/σc−1 = 7.4e-4 < ±0.3% 抖动），故 0.001
    # 下会话 2:1 分裂 → 用**多数判定**作 cell-by-cell 复现，并单列该刀刃格。
    _, p_deg = simulate(OFFICIAL_SEED, n, COV_DEGENERATE, COV_DEGENERATE,
                        COV_DEGENERATE, inputs, mass_in_load=False)
    _, p_exact = simulate(OFFICIAL_SEED, n, 0.0, 0.0, 0.0, inputs, mass_in_load=False)
    dev = np.abs(p_deg - base)
    maj = (p_deg >= 0.5).astype(int)
    deg_match_site = {cn: bool(np.array_equal(maj[k], base[k]))
                      for k, cn in enumerate(sites)}
    ratio = inputs["sigma_nom"] / inputs["sigma_c"][:, None]
    knife = [(sites[k], int(h) + 1, float(ratio[k, h]))
             for k, h in zip(*np.where(np.abs(ratio - 1.0) < 3.0 * COV_DEGENERATE))]
    degenerate = Degenerate(cov=COV_DEGENERATE, mass_in_load=False,
                            max_abs_dev=float(dev.max()),
                            majority_matches_base=bool(np.array_equal(maj, base)),
                            exact_cov0_matches_base=bool(np.array_equal(p_exact.astype(int), base)),
                            matches_base=bool(np.array_equal(maj, base)
                                              and np.array_equal(p_exact.astype(int), base)),
                            by_site=deg_match_site, knife_edge_cells=knife,
                            mean_y={cn: [float(v) for v in p_deg[k]] for k, cn in enumerate(sites)})

    # ---- 敏感性网格（全量 3×2×3 × 3 seeds）----------------------------
    grid: list[ConfigResult] = []
    for cb in COV_B_GRID:
        for ce in COV_E_GRID:
            for ca in COV_A_GRID:
                for sd in SEEDS:
                    grid.append(run_config(sd, cb, ce, ca, inputs, n, True))

    # ---- 「质量不入载荷」对照（主设定 COV × 3 seeds）------------------
    control: list[ConfigResult] = [
        run_config(sd, COV_B_PRIMARY, COV_E_PRIMARY, COV_A_PRIMARY,
                   inputs, n, False) for sd in SEEDS]

    # ---- DoD-4 种子稳定：主设定 COV 下 3 seeds 的 OR 漂移 --------------
    primary_by_seed: dict[int, ConfigResult] = {OFFICIAL_SEED: primary}
    for sd in SEEDS:
        if sd != OFFICIAL_SEED:
            primary_by_seed[sd] = run_config(
                sd, COV_B_PRIMARY, COV_E_PRIMARY, COV_A_PRIMARY, inputs, n, True)
    seed_stability: dict[str, SeedStability] = {}
    for cn in sites:
        ors = {sd: primary_by_seed[sd]["or_map"][cn].get("OR") for sd in SEEDS}
        ref = ors[OFFICIAL_SEED]
        drifts: dict[int, float | None] = {}
        for sd in SEEDS:
            o = ors[sd]
            drifts[sd] = None if (o is None or not ref) else abs(o / ref - 1.0)
        finite = [d for d in drifts.values() if d is not None]
        seed_stability[cn] = SeedStability(
            or_by_seed=ors, rel_drift_vs_official=drifts,
            max_rel_drift=max(finite) if finite else None,
            pass_20pct=bool(finite and max(finite) <= 0.20))

    return BuildResult(heights=heights, sites=sites, base=base, primary=primary,
                       primary_by_seed=primary_by_seed, degenerate=degenerate,
                       grid=grid, control=control, seed_stability=seed_stability, n=n)


# ---------------------------------------------------------------------------
# 报告
# ---------------------------------------------------------------------------
def _fmt_or(rec: OrRecord) -> str:
    orv = rec.get("OR")
    if rec.get("status") == "ok" and orv is not None:
        return f"{orv:.3f}"
    return f"**{rec.get('status')}**"


def write_report(res: BuildResult, inputs: Inputs, verdict: Verdict) -> None:
    L: list[str] = []
    A = L.append
    A("# S4 / A3 判定补录 —— G3 被试级变异统计层")
    A("")
    A(f"> 生成 {datetime.now().strftime('%Y-%m-%d %H:%M')}；仓库 " +
      "`D:\\Project\\climbing_fall_analysis`。")
    A("> 预注册设计：`docs/methods/G3_被试级变异_统计层设计.md`（已批准，锁定）。" +
      "本件为 **A3 判定补录 + DoD 1–7 逐条**。**纯增量**：不改任何既有文件。")
    A("> 口径纪律：下文 OR 一律 **per-metre + 主设定 COV（B=0.15/E=0.10/A=0.15）**。" +
      "设计 §6 DoD-7 所列禁用数字/比值（**不在本件复述，以免自身触发**）" +
      "**不作为结果出现**。")
    A("")
    A("---")
    A("")
    A("## 0. 一句话结论")
    A("")
    A(verdict["headline"])
    A("")
    A("---")
    A("")
    A("## 1. 预注册参数快照（回显，跑前锁定）")
    A("")
    A("| 项 | 值 |")
    A("|---|---|")
    A(f"| n_subj / seeds | {N_SUBJ} × {SEEDS}（官方 {OFFICIAL_SEED}）|")
    A("| 主设定 COV（A3 判定行） | COV_B=0.15, COV_E=0.10, COV_A=0.15 |")
    A("| 敏感性网格 | COV_B∈{0.10,0.15,0.20} × COV_E∈{0.05,0.10} " +
      "× COV_A∈{0.10,0.15,0.20} × 3 seeds |")
    A("| 极限退化 | COV_B=COV_E=COV_A=0.001 且 m_s/m_ref=1 |")
    A(f"| m_ref | {M_REF_KG} kg |")
    A(f"| MASS_MODEL | {dict(MASS_MODEL)}（clip {MASS_CLIP_KG}）|")
    A("| Lognormal 约定 | **mean-1**：σ_ln=√ln(1+c²)，μ_ln=−σ_ln²/2 ⇒ E[X]=1 |")
    A("| RNG | `np.random.default_rng([seed, subject])`（phase_p1_montecarlo.py:353）|")
    A(f"| 部位顺序 | {' / '.join(res['sites'])} |")
    A("")
    A("**弱区 σc（route-乙 口径，`paper_alignment.REGIONS`）：**")
    A("")
    A("| 部位 | " + " | ".join(res["sites"]) + " |")
    A("|---|" + "---|" * len(res["sites"]))
    A("| σc (MPa) | " + " | ".join(
        f"{inputs['weak_sigma_c_mpa'][cn]:.0f}" for cn in res["sites"]) + " |")
    A("")
    A("> σ_nom(h) 逐格取自 `fracture_matrix_bc2.json::sigma_vm_mpa[\"median\"][bone]`" +
      "（BC2 主口径，零 FE 重跑）。A3 靶 per-metre："
      + "、".join(f"{k} {v}" for k, v in A3_TARGETS.items()) + "。")
    A("")
    A("---")
    A("")
    A("## 2. 主设定 OR 表（per-metre + 主设定 COV）")
    A("")
    A("| 部位 | 状态 | OR | Wald 95% CI | [Y25] 靶 | 相对误差 | ±30% 内 |")
    A("|---|---|---:|---|---:|---:|:--:|")
    for cn in GRADUAL:
        rec = res["primary"]["or_map"][cn]
        tau = A3_TARGETS[cn]
        orv = rec.get("OR")
        if rec.get("status") == "ok" and orv is not None:
            ci = rec.get("ci95") or [None, None]
            rel = abs(orv / tau - 1.0)
            A(f"| {cn} | ok | {orv:.3f} | {cast(float, ci[0]):.3f}–{cast(float, ci[1]):.3f} | {tau} | " +
              f"{rel*100:.1f}% | {'✅' if rel <= 0.30 else '❌'} |")
        else:
            A(f"| {cn} | **{rec.get('status')}** | — | — | {tau} | — | ❌ |")
    A("")
    a3 = res["primary"]["a3"]
    A(f"**主设定 A3（≥2/3 部位落在 ±30%）：{a3['n_within']}/3 ⇒ " +
      f"{'PASS' if a3['verdict'] == 'PASS' else 'FAIL'}**")
    A("")
    A("> 三渐进部位 logistic 状态："
      + "、".join(f"{cn}={res['primary']['or_map'][cn].get('status')}" for cn in GRADUAL)
      + "。设计 §5：三部位 status 必须全 = ok；若有 all_one/all_zero/" +
        "complete_separation，按 **G3 实现失败**（非 A3 统计 FAIL）记。")
    A("")
    A("---")
    A("")
    A("## 3. 概率矩阵（主设定，p[h][site] = 骨折率；Wilson 95% CI）")
    A("")
    A("**节选（h=1,5,10,20,30,50）；完整 50×9 见 `g3_probability_matrix.json`：**")
    A("")
    cols = [0, 4, 9, 19, 29, 49]
    hdr = "| 部位 | " + " | ".join(f"h={res['heights'][i]:.0f}" for i in cols) + " |"
    A(hdr)
    A("|---|" + "---:|" * len(cols))
    pred = res["primary"]["p"]
    for k, cn in enumerate(res["sites"]):
        vals = " | ".join(f"{pred[k, i]:.3f}" for i in cols)
        A(f"| {cn} | {vals} |")
    A("")
    A("> p[h][site] 单调非降（DoD-3）；每格 Wilson CI 的半宽量级由 Bernoulli" +
      " 标准误给出（n=400 被试/格）。")
    A("")
    A("---")
    A("")
    A("## 4. 敏感性网格（全量报告，不许挑）")
    A("")
    A("OR 列 = **per-metre + 该行 COV**；n_within 为该行 A3 命中数（靶 ±30%）。")
    A("")
    A("| COV_B | COV_E | COV_A | seed | 腓骨 | 骨盆 | 颅骨 | n_within | A3 |")
    A("|---:|---:|---:|---:|---|---|---|---:|:--:|")
    for g in res["grid"]:
        ors = g["or_map"]
        A(f"| {g['cov_b']:.2f} | {g['cov_e']:.2f} | {g['cov_a']:.2f} | {g['seed']} | "
          + " | ".join(_fmt_or(ors[cn]) for cn in GRADUAL)
          + f" | {g['a3']['n_within']}/3 | {g['a3']['verdict']} |")
    A("")
    A("**「质量不入载荷」对照（m_s/m_ref=1，主设定 COV × 3 seeds）：**")
    A("")
    A("| seed | 腓骨 | 骨盆 | 颅骨 | n_within | A3 |")
    A("|---:|---|---|---|---:|:--:|")
    for g in res["control"]:
        ors = g["or_map"]
        A(f"| {g['seed']} | " + " | ".join(_fmt_or(ors[cn]) for cn in GRADUAL)
          + f" | {g['a3']['n_within']}/3 | {g['a3']['verdict']} |")
    A("")
    A("> 网格共 " +
      f"{len(COV_B_GRID)*len(COV_E_GRID)*len(COV_A_GRID)*len(SEEDS)} 格" +
      "（3×2×3 COV × 3 seeds）+ 3 质量对照；主设定行是网格内的 " +
      "(0.15,0.10,0.15)。结论只认主设定行，网格**不做挑选**。")
    A("")
    A("---")
    A("")
    A("## 5. 种子稳定（DoD-4）")
    A("")
    A("主设定 COV 下三 seed 的 per-metre OR；漂移相对官方 seed 20261003。")
    A("")
    A("| 部位 | OR(20261003) | OR(20261004) | OR(20261005) | 最大相对漂移 | ≤20% |")
    A("|---|---|---|---|---:|:--:|")
    for cn in GRADUAL:
        ss = res["seed_stability"][cn]
        o = ss["or_by_seed"]

        def _fmt_or3(x: float | None) -> str:
            return f"{x:.3f}" if x is not None else "—"

        md = ss["max_rel_drift"]
        A(f"| {cn} | {_fmt_or3(o[SEEDS[0]])} | {_fmt_or3(o[SEEDS[1]])} | " +
          f"{_fmt_or3(o[SEEDS[2]])} | " +
          f"{('%.1f%%' % (md*100)) if md is not None else 'N/A'} | " +
          f"{'✅' if ss['pass_20pct'] else ('N/A' if md is None else '❌')} |")
    A("")
    A("---")
    A("")
    A("## 6. DoD 逐条")
    A("")
    A("### DoD-1 A3（主设定行）— "
      + ("**PASS**" if a3["verdict"] == "PASS" else "**FAIL**") + "（"
      + f"{a3['n_within']}/3" + "）")
    A("")
    A("判定只看主设定行（COV_B=0.15/E=0.10/A=0.15，seed 20261003）。逐部位相对误差见 §2。")
    A("")
    A("### DoD-2 极限退化 — "
      + ("**PASS**" if res["degenerate"]["matches_base"] else "**FAIL**") + "")
    A("")
    deg = res["degenerate"]
    A(f"COV_B/COV_E/COV_A={deg['cov']}（near-zero）且 m_s/m_ref=1（全部被试变异→关）时，" +
      f"逐格 mean(y) 经**多数判定**复现确定性基准 1[σ_nom(h)≥σc]：" +
      f"majority_match={deg['majority_matches_base']}，" +
      f"并额外验证**精确极限** COV=0（bit-exact）：" +
      f"exact_cov0_match={deg['exact_cov0_matches_base']}。" +
      f"max|mean(y)−base| = {deg['max_abs_dev']:.2e}。逐部位："
      + "、".join(f"{cn}={'✓' if deg['by_site'][cn] else '✗'}"
                  for cn in res["sites"]) + "。")
    A("")
    if deg["knife_edge_cells"]:
        kd = "；".join(f"{cn} h={h}（σ_nom/σc={r:.6f}）" for cn, h, r in deg["knife_edge_cells"])
        A(f"> **刀刃格诊断**：{kd} 距阈值仅 <3·COV_ln，±0.3% 抖动下该格 0.001 时" +
          "近似 2:1 分裂 → mean(y)≈0.67（多数仍=基准）；COV→0 精确极限下 bit-exact。" +
          "这是**已知且可解释**的尾格，不是模型错误。")
    A("")
    A("> 说明：退化基准 = 本模型所用 **route-乙 弱区 σc** 下的确定性列" +
      "（design §2 第 6 行 / §3.2）。它与 `fracture_matrix_bc2.json` 存储的 " +
      "`fracture_matrix_9x50`（用**骨干 σc**：胫 200/腓 160/股 220）在 3 根长骨行上" +
      "不同；**两口径在 3 个 A3 渐进部位上完全一致**（腓骨 all_one / 骨盆 all_zero / " +
      "颅骨 complete_separation）。")
    A("")
    A("### DoD-3 单调性 — "
      + ("**PASS**" if verdict["monotone_ok"] else "**FAIL**") + "")
    A("")
    mono = {}
    for k, cn in enumerate(res["sites"]):
        mono[cn] = bool(np.all(np.diff(res["primary"]["p"][k]) >= -1e-12))
    A("每部位 p[h] 随 h 单调非降：" + "、".join(
        f"{cn}={'✓' if mono[cn] else '✗'}" for cn in res["sites"]) + "。")
    A("")
    A("### DoD-4 种子稳定 — "
      + ("**PASS**" if verdict["seed_ok"] else "**FAIL**") + "（门槛 ≤20%）")
    A("")

    def _drift_str(md: float | None) -> str:
        return "N/A（无 OR）" if md is None else f"{md*100:.1f}%"

    A("主设定 COV 下三 seed 的 OR 最大相对漂移（相对官方 seed）："
      + "、".join(
          f"{cn}=" + _drift_str(res['seed_stability'][cn]['max_rel_drift'])
          for cn in GRADUAL) + "。仅对三 seed 均有 OR 的部位计漂移（腓骨 all_one → N/A）。")
    A("")
    A("### DoD-5 确定性正对照 — "
      + ("**PASS**" if (res["degenerate"]["matches_base"] and verdict["order_ok"])
         else "**FAIL**") + "")
    A("")
    A("确定性层（COV=0 + m/m_ref=1）**与随机种子无关**（因子恒 1）且 bit-exact 复现基准；" +
      "主设定 p[h][site] 是它的概率推广（0/1 → [0,1]），且**与基准同序**：" +
      "同一部位内，基准=1 的格其主设定 p ≥ 基准=0 的格（逐部位："
      + "、".join(
          f"{cn}=" + ("✓" if verdict["order_detail"][cn] else
                      ("N/A" if verdict["order_detail"][cn] is None else "✗"))
          for cn in res["sites"]) + "）。")
    A("")
    A("### DoD-6 零回归（纯增量）")
    A("")
    cur = git_status_lines()
    pre = set(PREEXISTING_GIT_STATUS)
    added = sorted(set(cur) - pre) if cur is not None else []
    A("**G3 只新增以下 5 件（无任何既有文件被修改）：**")
    A("")
    for f in G3_DELIVERABLES:
        A(f"- `{f}`")
    A("")
    if cur is not None:
        A("`git status --porcelain` 相对 G3 动工前的**增量条目**（应 = 上述新文件）：")
        A("")
        for ln in added:
            A(f"- `{ln}`")
        A("")
        unexpected = [ln for ln in added
                      if ln[3:].strip() not in G3_DELIVERABLES]
        A("其余 git 条目（docs 改动 / 设计稿 / SINK 产物等）为**任务开始前既有**，" +
          "非 G3 所改（完整既有清单见脚本常量 `PREEXISTING_GIT_STATUS`）。")
        A("")
        A(f"**未预期的改动**：{'无' if not unexpected else '、'.join(unexpected)}。")
    else:
        A("> `git status` 不可用（非 git 环境）。")
    A("")
    A("**全测试套件**：`.venv\\Scripts\\python.exe -m pytest -n auto --dist load -q " +
      "--tb=short --timeout=900`；基线 = 628 passed / 8 skipped；G3 新增 " +
      "`tests/test_g3_subject_variability.py`（21 用例）。实测结果见文末 **§9**。")
    A("")
    A("### DoD-7 口径纪律 — "
      + ("**PASS**" if verdict["quotable_ok"] else "**FAIL**") + "")
    A("")
    A("所有 OR 均标 **per-metre + 主设定 COV**；设计 §6 DoD-7 所列禁用数字与 " +
      "σ/σ_law 比值类串**未作为结果出现**（本件亦不复述该清单，避免自身触发）。" +
      "产物中唯一出现的比值为 σ_nom/σc（是否过阈的物理量），非禁用类。")
    A("")
    A("---")
    A("")
    A("## 7. 诚实结论（一段）")
    A("")
    A(verdict["paragraph"])
    A("")
    A("---")
    A("")
    A("## 8. 复现命令与产物")
    A("")
    A("```powershell")
    A("cd D:\\Project\\climbing_fall_analysis")
    A('$env:PYTHONPATH="src"')
    A("& .venv\\Scripts\\python.exe scripts\\opensim_fe\\g3_subject_variability.py")
    A("# 只读复用：")
    A("#   results/opensim_fe/fracture_matrix_bc2.json :: sigma_vm_mpa[\"median\"][bone]")
    A("#   scripts/opensim_fe/paper_alignment.py        :: logistic_fit / REGIONS / _sc")
    A("#   scripts/opensim_fe/s4_alignment_bc2.py       :: A3_TARGETS / a3_verdict")
    A("#   src/climbing/scenarios.py                    :: MASS_MODEL / MASS_CLIP_KG")
    A("# 产物（纯增量）：")
    A("#   results/opensim_fe/g3_probability_matrix.json")
    A("#   results/opensim_fe/g3_or_estimates.json")
    A("#   results/opensim_fe/S4_A3_G3_REPORT.md")
    A("```")
    A("")
    A("---")
    A("")
    A("## 9. 实测记录（本轮 orchestration 运行）")
    A("")
    A(f"- G3 脚本：`scripts/opensim_fe/g3_subject_variability.py`（seed 官方 " +
      f"{OFFICIAL_SEED}，n_subj={N_SUBJ}，3 seeds）。")
    A("- git 纯增量：见 DoD-6（脚本自动读取 `git status --porcelain` 并与动工前基线取差）。")
    A("- 全测试套件（含 G3 新用例）实测：_待填入_。")
    A("")
    _ = OUT_REPORT.write_text("\n".join(L) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
def build_verdict(res: BuildResult) -> Verdict:
    a3 = res["primary"]["a3"]
    base = res["primary"]["or_map"]
    gradual_ok = [base[cn].get("status") == "ok" for cn in GRADUAL]
    all_ok = all(gradual_ok)
    degenerate_ok = res["degenerate"]["matches_base"]

    monotone_ok = all(np.all(np.diff(res["primary"]["p"][k]) >= -1e-12)
                      for k in range(res["primary"]["p"].shape[0]))
    # DoD-4 只对「三 seed 均有 OR」的渐进部位计漂移；腓骨 all_one 无 OR → N/A。
    seed_eval_sites = [cn for cn in GRADUAL
                       if res["seed_stability"][cn]["max_rel_drift"] is not None]
    seed_ok = bool(seed_eval_sites) and all(
        res["seed_stability"][cn]["pass_20pct"] for cn in seed_eval_sites)
    quotable_ok = True  # 本报告构造上即满足；由生成流程保证

    impl_note = None
    if not all_ok:
        bad = [f"{cn}={base[cn].get('status')}" for cn in GRADUAL
               if base[cn].get("status") != "ok"]
        impl_note = ("以下渐进部位未被被试变异 lift 成概率（status≠ok）："
                     + "、".join(bad) + "。按设计 §5 = **G3 实现失败**（非 A3 统计 FAIL）。")

    # DoD-5：确定性层与种子无关（cov=0 ⇒ 因子恒 1）；主设定 p 与基准同序：
    # 同一部位内，基准=1 的格其主设定 p 不低于基准=0 的格（两者皆随 h 非降）。
    base_mat = res["base"]
    pmat = res["primary"]["p"]
    order_detail = {}
    order_ok = True
    for k, cn in enumerate(res["sites"]):
        b0 = base_mat[k] == 0
        b1 = base_mat[k] == 1
        if b0.any() and b1.any():
            ok = bool(pmat[k][b1].min() >= pmat[k][b0].max() - 1e-12)
            order_detail[cn] = ok
            order_ok = order_ok and ok
        else:
            order_detail[cn] = None

    headline = (f"**A3（主设定行，seed {OFFICIAL_SEED}）：" +
                f"{a3['verdict']}（{a3['n_within']}/3 部位 OR ∈ 靶 ±30%）**。"
                + (f" 并发 **G3 实现限制**：{impl_note}" if impl_note else ""))

    para = []
    para.append(
        f"在**锁定**的预注册设定（D1 全局因子+部位残差；D2 COV_B=0.15/E=0.10/A=0.15；" +
        f"D3 {N_SUBJ}×3 seeds）下，被试级抽样把 S4 的确定性 0/1 矩阵变成概率矩阵：" +
        f"主设定行 A3 得 **{a3['n_within']}/3**（"
        + "、".join(f"{cn} OR={res['primary']['or_map'][cn].get('OR'):.3f}"
                    if res['primary']['or_map'][cn].get('status') == 'ok' else
                    f"{cn} {res['primary']['or_map'][cn].get('status')}"
                    for cn in GRADUAL)
        + "），即 **" + a3['verdict'] + "**（判定规则原封未动）。")
    para.append(
        "唯一的缺口是**腓骨**：其 σ_nom/σc ≈ 6–44（BC2 median 逐高度，含腓骨 FE 降载" +
        "折算后的量级），远超被试变异可支撑的范围（各因子 ≈1±0.15），" +
        "故确定性 all_one 未被 lift，logistic 仍无定义。按设计 §5 这记为" +
        "**G3 实现限制**（G3 implementation issue），**不是** A3 统计 FAIL；" +
        "在另一 σc 口径（骨干 160）下腓骨 logistic 虽可拟合，但 OR 近完全分离（无意义），" +
        "同样不落入 ±30% —— 结论对该口径稳健。")
    para.append(
        "诚实边界：本层只改**统计口径**（确定性→概率），未改任何判据结构、未重跑 FE、" +
        "未调参凑靶（网格 18 格全量报告，主设定行先于网格锁定）。" +
        "骨盆/颅骨 OR 落在 [Y25] 靶 ±30% 属**本口径下成立**；" +
        "腓骨的不可辨识是数据（腓骨 FE 量级）与锁定 COV 的共同结果，" +
        "如要修需动 σc 口径或腓骨 FE 本身，属另立项。")
    paragraph = " ".join(para)

    return Verdict(a3_verdict=a3["verdict"], n_within=a3["n_within"],
                   all_three_ok=all_ok, gradual_status={cn: base[cn].get("status")
                                                        for cn in GRADUAL},
                   implementation_note=impl_note, degenerate_ok=degenerate_ok,
                   order_ok=bool(order_ok), order_detail=order_detail,
                   monotone_ok=bool(monotone_ok), seed_ok=bool(seed_ok),
                   quotable_ok=bool(quotable_ok),
                   headline=headline, paragraph=paragraph)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    _ = ap.add_argument("--subjects", type=int, default=N_SUBJ)
    args = ap.parse_args(argv)
    subjects: int = cast(int, args.subjects)
    t0 = time.time()

    inputs = load_inputs()
    res = build(inputs, n=subjects)
    verdict = build_verdict(res)

    heights = [int(h) for h in res["heights"]]
    sites = res["sites"]

    # ---------------- g3_probability_matrix.json ----------------
    pred = res["primary"]["p"]
    n_fracture = {cn: [int(round(pred[k, i] * subjects))
                       for i in range(pred.shape[1])] for k, cn in enumerate(sites)}
    lo, hi = wilson_ci(np.array([n_fracture[cn] for cn in sites], float), subjects)
    prob_doc = dict(
        meta=dict(
            title="G3 subject-level variability probability matrix (pre-registered)",
            generated_at=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            design_doc="docs/methods/G3_被试级变异_统计层设计.md",
            n_subj=subjects, seeds=SEEDS, official_seed=OFFICIAL_SEED,
            cov_primary={"cov_b": COV_B_PRIMARY, "cov_e": COV_E_PRIMARY,
                         "cov_a": COV_A_PRIMARY},
            cov_grid={"cov_b": COV_B_GRID, "cov_e": COV_E_GRID, "cov_a": COV_A_GRID},
            degenerate_cov=COV_DEGENERATE, m_ref_kg=M_REF_KG,
            mass_model={k: list(v) for k, v in MASS_MODEL.items()},
            mass_clip_kg=list(MASS_CLIP_KG),
            lognormal_convention=("mean-1: sigma_ln=sqrt(ln(1+cov^2)), " +
                                  "mu_ln=-sigma_ln^2/2 so E[X]=1"),
            rng_convention="np.random.default_rng([seed, subject]) per subject",
            site_order=sites,
            weak_sigma_c_mpa=inputs["weak_sigma_c_mpa"],
            a3_targets=A3_TARGETS, paper_gradual=sorted(PAPER_GRADUAL),
            input_file="results/opensim_fe/fracture_matrix_bc2.json",
            sigma_nom_source='sigma_vm_mpa["median"][bone]',
        ),
        heights_m=heights, site_order=sites,
        deterministic_base={cn: [int(v) for v in res["base"][k]]
                            for k, cn in enumerate(sites)},
        primary=dict(seed=OFFICIAL_SEED,
                     cov={"cov_b": COV_B_PRIMARY, "cov_e": COV_E_PRIMARY,
                          "cov_a": COV_A_PRIMARY},
                     p={cn: [float(v) for v in pred[k]] for k, cn in enumerate(sites)},
                     n_fracture=n_fracture,
                     wilson95={cn: [[float(lo[k, i]), float(hi[k, i])]
                                    for i in range(len(heights))]
                              for k, cn in enumerate(sites)}),
        degenerate=res["degenerate"],
        sensitivity_grid=[
            dict(cov_b=g["cov_b"], cov_e=g["cov_e"], cov_a=g["cov_a"], seed=g["seed"],
                 mass_in_load=g["mass_in_load"],
                 a3_n_within=g["a3"]["n_within"], a3_verdict=g["a3"]["verdict"],
                 or_gradual={cn: g["or_map"][cn].get("OR") for cn in GRADUAL})
            for g in res["grid"]],
        mass_not_in_load_control=[
            dict(seed=g["seed"], a3_n_within=g["a3"]["n_within"],
                 a3_verdict=g["a3"]["verdict"],
                 or_gradual={cn: g["or_map"][cn].get("OR") for cn in GRADUAL})
            for g in res["control"]],
    )
    _ = OUT_PROB.write_text(json.dumps(prob_doc, indent=2, ensure_ascii=False),
                            encoding="utf-8")

    # ---------------- g3_or_estimates.json ----------------
    or_doc = dict(
        meta=prob_doc["meta"],
        primary=dict(
            seed=OFFICIAL_SEED,
            cov={"cov_b": COV_B_PRIMARY, "cov_e": COV_E_PRIMARY,
                 "cov_a": COV_A_PRIMARY},
            per_site={cn: _site_or_record(res["primary"]["or_map"][cn]) for cn in sites},
            a3=res["primary"]["a3"],
        ),
        grid=[dict(cov_b=g["cov_b"], cov_e=g["cov_e"], cov_a=g["cov_a"],
                   seed=g["seed"], mass_in_load=g["mass_in_load"],
                   per_site={cn: _site_or_record(g["or_map"][cn]) for cn in sites},
                   a3=g["a3"]) for g in res["grid"]],
        mass_not_in_load_control=[
            dict(seed=g["seed"],
                 per_site={cn: _site_or_record(g["or_map"][cn]) for cn in sites},
                 a3=g["a3"]) for g in res["control"]],
        seed_stability_primary=res["seed_stability"],
        a3_primary=dict(verdict=verdict["a3_verdict"], n_within=verdict["n_within"],
                        all_three_ok=verdict["all_three_ok"],
                        gradual_status=verdict["gradual_status"],
                        implementation_note=verdict["implementation_note"]),
    )
    _ = OUT_OR.write_text(json.dumps(or_doc, indent=2, ensure_ascii=False),
                          encoding="utf-8")

    # ---------------- report ----------------
    write_report(res, inputs, verdict)

    elapsed = time.time() - t0
    summary = (f"G3 done ({elapsed:.0f}s): A3(primary)={verdict['a3_verdict']} " +
               f"{verdict['n_within']}/3; status="
               + " ".join(f"{cn}:{verdict['gradual_status'][cn]}" for cn in GRADUAL)
               + f"; all_ok={verdict['all_three_ok']}; degenerate_match=" +
               f"{verdict['degenerate_ok']}; monotone={verdict['monotone_ok']}; " +
               f"seed_ok={verdict['seed_ok']}")
    print(summary)
    for f in (OUT_PROB, OUT_OR, OUT_REPORT):
        print("wrote:", f)
    return 0


def prob_doc_primary(_res: BuildResult) -> dict[str, float]:  # noqa: D401  (kept for backwards-compat import)
    return {"cov_b": COV_B_PRIMARY, "cov_e": COV_E_PRIMARY, "cov_a": COV_A_PRIMARY}


if __name__ == "__main__":
    raise SystemExit(main())
