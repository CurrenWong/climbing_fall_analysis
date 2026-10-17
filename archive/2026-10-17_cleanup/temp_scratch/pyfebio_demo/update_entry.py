"""更新 FE_PIPELINE_TODO.md 接手入口的「下一步」，反映韧带工作的真实状态。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
t = P.read_text(encoding="utf-8")
b = len(t.encode("utf-8"))

old = """1. **加 8 组踝韧带** ⭐ **建议先做这个**
   数据已在 `temp/thums/thums_lowerlimb.npz`（`81100101~81101401_shell`）。
   理由：**不经过接触这条路径** ⇒ 不会被下面的 2 卡住；
   而且能复现吴恺的**定性排序**（**胫跟＞胫舟＞跟腓＞下胫腓前**），
   这比 3.74 MPa 的绝对值**验证价值更高**（项目计划里验证目标就是「量级 + 定性排序」，
   量级已拿到，排序还是空白）"""
assert t.count(old) == 1, f"hit {t.count(old)}"

new = """1. **加 8 组踝韧带** —— ⚠️ **解剖/材料/模型都做完，但卡在"韧带单元表示"**
   **进度与全部细节见 `ANKLE_LIG_TODO.md`**（含三条路各自的硬墙与两条下一步）。
   **一句话状态**：吴恺 4 条韧带**已全部定位且有几何证据**（胫跟=`81101202`、
   胫舟=`81101203`、跟腓=`81100801`、下胫腓=`81101101`），
   材料已解出（EA 510/510/510/**75.4** MPa、厚 2.0 mm），
   模型已组装到"读文件 + 跑首个载荷步"（32,488 实体 + 胫距接触 + 600 N + 固定跟骨/约束距骨）。
   **卡点**：韧带 FE 表示 —— ① `elastic-truss` 真实模型**读文件即崩**
   （实验性特性，最小实验能过）② THUMS 原始 shell 引入 **137 个只属韧带的节点**
   ⇒ 转动自由度无约束 ⇒ 刚体模态 ⇒ 首步 22 个单元翻负 ③ 只用骨节点的 tri3 网读文件崩。
   **下一步最可能成功**：把韧带专属节点**焊到最近骨表面节点**后用 shell
   （或重做骨网格保留全部 THUMS 表面节点）。
   **吴恺的目标值**（务必照抄）：**胫跟 5.24 > 胫舟 5.04 > 跟腓 2.05 > 下胫腓前 1.92 mm**；
   边界=胫腓骨上截面 600 N + 固定跟骨 + 约束距骨。
   ✅ **已确认**：下胫腓**必然垫底**（材料软 6.8 倍，几何先验也最低）——
   两边一致，这是目前唯一能下的排序结论。"""
t = t.replace(old, new)

# 顺手在关键文件表里加两行
old2 = "| `docs/踝关节FE建模交接.md` | 全程历史（§1–§13，46 KB）|"
assert t.count(old2) == 1
t = t.replace(old2, "| `ANKLE_LIG_TODO.md` | **韧带工作全部细节**（3 条路的墙 + 下一步）|\n"
                    "| `scripts/ankle_fe/ligaments.py` | 从 THUMS 壳导出韧带几何 + 附着 |\n"
                    "| `scripts/ankle_fe/ankle_lig_feb.py` | 韧带模型组装（9 个环境开关）|\n"
                    "| `docs/踝关节FE建模交接.md` | 全程历史（§1–§13，46 KB）|")

P.write_text(t, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md {b:,} -> {len(t.encode('utf-8')):,} bytes")
