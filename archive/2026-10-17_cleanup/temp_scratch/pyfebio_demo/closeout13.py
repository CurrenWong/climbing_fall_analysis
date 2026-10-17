"""压缩前收口 ①：session note §十三 + daily 索引 + MEMORY 等量替换。"""
from pathlib import Path

WS = Path(r"C:\Users\curre\.copaw\workspaces\default")
SESS = WS / "memory/2026-10-04/ankle-fe-solve-and-handover.md"
DAILY = WS / "memory/2026-10-04.md"
M = WS / "MEMORY.md"
LIM = 10240

SEC = """

---

## 十三、软骨重分网：打通「可建可解」，tie 成新瓶颈（2026-10-04 末）

承接 `[LRN-052]`（三条捷径排除：改偏移方向 / Taubin 平滑 / 纯细分）。

### 13.1 ✅ 成果：细网格可建、可解
新增 `scripts/ankle_fe/cart_patch.py`（**2D 参数平面重心网格 + 精确抬回**）：
- 不做边界裁剪/点定位：在**每个原三角内部**打重心网格 ⇒ 域 = 原补丁，覆盖 100%
- 用**同一组重心坐标**精确抬回 ⇒ 落在原分片线性曲面上 ⇒ **翘曲恒 0**

| m | 底面三角 | 软骨tet | 接触面(tib/tal) | 耗时 | rc |
|---|---|---|---|---|---|
| 基线 | 1119 | — | 152/200 | 18s | 0 |
| 1 | 1808 | 1074 | 152/200 | 59s | **0** |
| 2 | 5424 | 3222 | 456/600 | 90s | **0** |
| 3 | 10848 | 6444 | 915/1206 | 156s | **0** |
| 4 | 14464 | — | — | 218s | **0** |

质量**不比基线差**（nv=V/lmax³）：m=2 nv 中位 1.75e-02 vs 基线 1.31e-02；
nv<1e-4 占 **0.65%** vs 基线 **1.29%**。

### 13.2 🔴 三个真 bug（全已修）
**① `build_feb.py` DET_MIN 逻辑漏洞**（基线受益）
`if d > DET_MIN: keep` 否则无条件翻转 ⇒ **正但小**的单元被翻成负 ⇒ 当坏单元剔除 = 留洞。
改成 `if d > 0: keep` ⇒ **基线剔除 9 → 0**，`rc=0/17.1s` 不变 ✓

**② 底面三角对角线必须与骨面一致（0-2）** ← 最初「怎么都不收敛」的真凶
为优化最小内角换 1-3 对角 ⇒ primary 与 secondary 不重合 ⇒ **tie 拉错 ⇒ 虚假刚度 ⇒ 翻转**。
强制 `CP_DIAG=02` 后接触面/tie primary **152/200 与基线逐字一致**。
⇒ **判据：tied-elastic 的 primary 三角化必须与 secondary（骨面）同源同对角。**

**③ 细网格必须 `LAUGON=AUGLAG`**，PENALTY 必失败
**决定性判据**：失败点恒为 `1/(11·TIME_STEPS)`，与步长无关（扫 50/100/200 全是 /11）
⇒ **近乎零载荷就翻转** ⇒ 不是几何/载荷，**是接触求解器**。
且 TIME_STEPS 扫描 50→200 **完全无效** ⇒ 排除「步长太大」。

### 13.3 🚧 剩余瓶颈：tie 界面
- **`|u|max = 1.35~2.08 mm`，而指定顶面位移只有 0.1472 mm** ⇒ 超出 **10×以上**，
  物理上不可能
- 支反力**不收敛**：penalty=1 → **308.5 / 205.1 / 227.4 N**（m=1/2/3）
- `TIE_PENALTY=1` 太软 ⇒ 软骨底面滑移 ~1mm；**加硬到 4/16/64 全部 rc=1**
  （与 `[LRN-046~050]` 粗网格 penalty 单调收敛的结论**相反**）
- 线索：m=2 时 **primary 456 面 vs secondary 仅 117 面**（`TIE_SEC_TOL=0.05` 筛出），
  比例从 1.3:1 恶化到 3.9:1

### 13.4 ➡️ 下一步
1. 查 tie 发散：`TIE_SEC_TOL` 放大 / `search_radius` 调整 / secondary 覆盖核对
2. 或绕开 tie：**关节区局部**骨面细分到与软骨共形（整根骨加密曾失败，局部可能可行）
3. **加 sanity 判据**：`|u|max ≈ 指定位移`（现在 10× 异常应当场报警）

### 13.5 产物
`scripts/ankle_fe/cart_patch.py`（新）· `build_feb.py`（修 DET_MIN）·
`cart_hex.py`（v1，K=1 复现基线）· `temp/thums/cartilage_mesh_patch_m{1,2,3,4}.npz` ·
`temp/pyfebio_demo/{feas_patch2,diag_degen,diag_nv_m,conv_final,pen_conv2,conv_study}.py` ·
`[LRN-20261004-052]` `[LRN-20261004-053]`
"""
if SESS.exists():
    t = SESS.read_text(encoding="utf-8")
    if "十三、软骨重分网" not in t:
        SESS.write_text(t.rstrip() + SEC, encoding="utf-8")
    print(f"session note -> {SESS.stat().st_size:,} B")

if DAILY.exists():
    t = DAILY.read_text(encoding="utf-8")
    add = ("\n- `2026-10-04/ankle-fe-solve-and-handover.md` → §十三 软骨重分网："
           "`cart_patch.py` 可建可解(m1..4 全 rc=0，底面格 5424@m2)、"
           "3 个 bug 已修(DET_MIN 漏洞/对角线必须 0-2/AUGLAG 必需)、"
           "tie 界面成新瓶颈(|u|max 10× 超标、F 不收敛)\n")
    if "十三 软骨重分网" not in t:
        DAILY.write_text(t.rstrip() + "\n" + add, encoding="utf-8")
    print(f"daily note -> {DAILY.stat().st_size:,} B")

# MEMORY.md：等量替换，加入本轮状态
s = M.read_text(encoding="utf-8")
OLD = ("- ✅ **tied-elastic 验证通过**（`cart_refine.py`+`TIE=1`）：L=0 罚法**单调收敛**到共形解\n"
       "  436.817 N → penalty 1/2/**4**=404.7/420.9/**432.5**（**−0.99%**）。坑：tie 面须是\n"
       "  **真实 facet**（tet4 只三角）；**secondary 有洞 ⇒ 伪装成「硬 2.2×」**；\n"
       "  ⚠️ 面内加密仍被**母单元扭曲**挡住（hex8 恒 28% 翻负，`[LRN-046~051]`）")
NEW = ("- ✅ tied-elastic（`TIE=1`）L=0 单调收敛 436.817→432.5 N（−0.99%）；\n"
       "  🟡 **软骨已可重分网**（`cart_patch.py`，m=1..4 全 rc=0，底面格 5424@m2）\n"
       "  ⚠️ 铁律：**底面三角对角必须=骨面(0-2)**，否则 tie 拉错⇒假刚度⇒翻转；\n"
       "  **细网格必须 `LAUGON=AUGLAG`**。🚧 瓶颈=tie 界面（`|u|max` 超指定 10×、F 不收敛）\n"
       "  `[LRN-052/053]`")
if OLD in s:
    s = s.replace(OLD, NEW, 1)
    M.write_text(s, encoding="utf-8")
sz = M.stat().st_size
print(f"MEMORY.md disk = {sz:,} / {LIM}  余量 {LIM - sz:+,}  "
      f"[{'OK' if 'cart_patch.py' in M.read_text(encoding='utf-8') else 'XX'}]")
