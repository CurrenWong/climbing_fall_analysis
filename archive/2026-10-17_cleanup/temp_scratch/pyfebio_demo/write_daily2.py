"""写 2026-10-06 日记 + 最终校验（口径统一为解码字节数）。"""
import io
import os
import re

R = r"D:\Project\climbing_fall_analysis"
W = r"C:\Users\curre\.copaw\workspaces\default"
DAY = os.path.join(W, "memory", "2026-10-06.md")

TXT = """# 2026-10-06

## 🎯 踝 FE：C2 重定义完成（①），并解开"标签/应变"之谜

接压缩前收口的 4 项待办，本轮做完 **① 重定义 C2**，并**意外解开了两轮悬案**。

### ✅ C2 重定义（`[LRN-20261006-074]`）

**原 C2 为什么废**：它测"距骨 t=0→末态转动"，而 `POSE_IE` 把足部组**烘焙进参考构型**
⇒ 结构上只能测到 0.0014~0.0033°。**架构错位**（Anderson 的姿势是输出，我的是输入）。

**新 C2 = 验"姿势作为输入被忠实实现 + 语义方向正确"**，三项全 PASS：

| 子判据 | 实测 |
|---|---|
| C2a-1 幅值（\\|实测相对转角 − 请求角\\|<0.5°）| **误差 0.000°** |
| C2a-2 轴向（`\\|轴·AP\\|→1`）| **1.0000** |
| C2a-3 语义（三锚定一致）| **全一致、单调线性** |

**三个锚定判据**：① 足底法向 ② 足外侧缘高度 ③ 外侧韧带长变。
**"外侧"由腓骨定义**（`lat = 腓骨顶 − 胫骨顶`，⊥加载轴，实测 `|lat·axis|=0.0000`）
—— 与模型左右/镜像无关 ⇒ 不被模型自身符号约定污染。

**⇒ 结论：`+IE = 内翻`，代码注释本来就对。**

### 🚨 撤回我自己的错结论：`[ERR-20261006-002]`

上一轮写进 MEMORY/使用说明/TODO 的「`POSE_IE` 标签反了、+IE 实际是外翻」**是错的**。
**产生机制是循环论证**：
1. 看到「CFL 应变在 `−3` 最高」→ 2. 反推「所以 −3 必是内翻」→ 3. 用这个标签宣布「拮抗关系物理正确」
`temp/pyfebio_demo/save_atfl.py` 里 "现象：内翻（`POSE_IE=−3`）…" 是这条链的化石。

### 🔑 真正的错：**引用了错的口径**（`[LRN-20261006-075]`）

韧带应变有**三个口径**，旧结论用的是 (b) 却当成总应变量解读：

| 口径 | 定义 | 用途 |
|---|---|---|
| (a) 预应变 | (L_姿势 − L_中性)/L_中性 | 姿势本身拉了多少 ← **扭伤机制** |
| (b) 载荷增量 | (L_载荷 − L_姿势)/L_姿势 | = FEBio t=1 元素应变 |
| **(c) 总应变** | (L_载荷 − L_中性)/L_中性 | ✅ **阈值设在它** |

**⚠️ (b) 与损伤风险**反相关**：松弛韧带被压缩"拉直"时增量大，已绷紧的增量≈0
⇒ 增量最大的恰是**松弛（无风险）**那根。**

**`POSE_IE=+5`（内翻，`NORMAL` 7.6 s / 0 负 J）实测**：

| 韧带 | 侧 | (a)预应变% | (b)增量% | **(c)总应变%** |
|---|---|---|---|---|
| PTFL | 外 | +4.58 | −0.17 | **+4.42** |
| **CFL** | 外 | +8.98 | +2.58 | **+11.80** |
| ATFL | 外 | +5.12 | −0.04 | **+5.08** |
| MD1~4 | 内 | −3.4~−4.3 | +0.4~+3.7 | **−0.7~−3.1** |

⇒ 外侧预应变 **+6.22%** / 内侧 **−3.86%** ⇒ **内翻拉外侧，解剖学正确 ✓**
⇒ **旧"ATFL 内翻时几乎不动(0.17~0.20%)"之谜解开**：那是 (b)，ATFL 本就已绷住 (+5.12% 预应变)。

### ⚠️ 新发现：内外翻**不对称**（`[LRN-20261006-076]`）
- `+5`（内翻）`NORMAL`；**`−5`（外翻）`rc=1` / 1 个负 J / t=0 即失败**
- 外翻落位平移需 **+0.186 mm**（内翻只需 −0.092）⇒ 该侧几何更"顶"
- **⇒ 目前只有内翻侧可用**，"判姿势"覆盖范围受限

### 🐛 工具类 bug：`[ERR-20261006-001]`
`extreme_faces()` 选出的 NodeSet（`calcaneus_bottom_ns` 等）**跨姿势不稳定**
（中性 194 → 某姿势 227 节点）⇒ 跨姿势比较**必须固定用一侧的 ID 列表**。

### 🛠️ 新工具（进 git）
| 脚本 | 作用 |
|---|---|
| `scripts/ankle_fe/check_c2.py` | C2a 全检（`NO_RUN=1` 只建不解，零求解成本）|
| `scripts/ankle_fe/lig_strain.py` | **权威应变读数**（三口径，从几何直算）|
| `ankle_lig_feb.py` 新增 `NO_RUN=1` | 只写 `.feb` 不调 FEBio（几何检查用）|

### 落盘
- LRN-20261006-074/075/076、ERR-20261006-001/002（LEARNINGS 89 条 / ERRORS 20 条，**重号 0**）
- 新 `references/ankle-fe-criteria.md`（C2 三锚定 / 三口径实测表 / 限制 / 工具）
- MEMORY 更正 + 压缩到 **10,226 / 10,240 B** ✅
- 使用说明 §4 判据表换成 C2 结果 + 三口径说明；§5 加两条限制；下一步重排

### 下一步
1. **查/修外翻侧不收敛**（`POSE_IE<0`）—— 双侧可用的前提
2. **修韧带插入点几何**（ATFL 36.99 → 15~20 mm 量级）—— 修完 (c) 才有资格当**定量**指标
3. 给 C3 定量阈值（用 (c)）
"""

os.makedirs(os.path.dirname(DAY), exist_ok=True)
old = io.open(DAY, encoding="utf-8").read() if os.path.exists(DAY) else ""
if "C2 重定义完成（①）" in old:
    print("· 日记已存在该节")
else:
    io.open(DAY, "w", encoding="utf-8").write(TXT + ("\n\n" + old if old else ""))
    print(f"✓ {DAY}  {len(TXT.encode('utf-8'))} B")

print("\n" + "=" * 72)
print("最终校验（口径：解码后字节数）")
print("=" * 72)
ok = True
CHK = [
    (os.path.join(W, "MEMORY.md"), 10240),
    (os.path.join(W, "memory", "2026-10-06.md"), None),
    (os.path.join(W, "references", "ankle-fe-criteria.md"), None),
    (os.path.join(R, "docs", "踝FE可用模型使用说明.md"), None),
    (os.path.join(R, "ANKLE_FIX_TODO.md"), None),
    (os.path.join(R, ".learnings", "LEARNINGS.md"), None),
    (os.path.join(R, ".learnings", "ERRORS.md"), None),
    (os.path.join(R, "scripts", "ankle_fe", "check_c2.py"), None),
    (os.path.join(R, "scripts", "ankle_fe", "lig_strain.py"), None),
]
for p, cap in CHK:
    if not os.path.exists(p):
        print(f"  ✗ 缺失 {p}")
        ok = False
        continue
    t = io.open(p, encoding="utf-8").read()
    n = len(t.encode("utf-8"))
    flag = ""
    if cap:
        flag = " ✅" if n <= cap else f" ❌ 超 {n-cap}"
        ok &= n <= cap
    nm = os.path.relpath(p, R) if p.startswith(R) else os.path.relpath(p, W)
    print(f"  {nm:<46} {n:>8} B{flag}")

# 关键内容存在性
D = os.path.join(R, "docs", "踝FE可用模型使用说明.md")
d = io.open(D, encoding="utf-8").read()
need = ["C2 重定义", "+11.80%", "lig_strain.py", "三个口径", "内外翻不对称",
        "36.99 mm", "check_c2.py", "0.000°", "目前只有内翻侧可用"]
miss = [k for k in need if k not in d]
print(f"\n  使用说明关键内容缺失：{miss or '无 ✓'}")
ok &= not miss

m = io.open(os.path.join(W, "MEMORY.md"), encoding="utf-8").read()
need2 = ["C2 已重定义并 PASS", "ankle-fe-criteria.md", "反相关", "外翻 `−5` 不收敛"]
miss2 = [k for k in need2 if k not in m]
print(f"  MEMORY 关键内容缺失：{miss2 or '无 ✓'}")
ok &= not miss2

print("\n结论：", "PASS ✅" if ok else "有问题 ⚠")
