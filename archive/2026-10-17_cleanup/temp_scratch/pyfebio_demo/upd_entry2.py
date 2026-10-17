"""更新 FE_PIPELINE_TODO.md 接手入口的韧带条目（模型已跑通，排序为负结果）。"""
from pathlib import Path

P = Path(r"D:\Project\climbing_fall_analysis\FE_PIPELINE_TODO.md")
t = P.read_text(encoding="utf-8")
b = len(t.encode("utf-8"))

i0 = t.index("1. **加 8 组踝韧带**")
i1 = t.index("2. **修 subtalar 接触面 2-流形**")

new = """1. **加 8 组踝韧带** —— ✅ **模型已跑通（5 秒）** / ❌ **排序复现为负结果**
   **全部细节见 `ANKLE_LIG_TODO.md`**。
   **可用配置**（5.4 s 收敛，20 步）：
   `LIG_MODE=shell DRIVER=disp TRUNC_MM=60 DISP_MM=0.20 FIX_TALUS=all`
   - **解剖/材料**：吴恺 4 条韧带全部定位（胫跟=`81101202`、胫舟=`81101203`、
     跟腓=`81100801`、下胫腓=`81101101`），EA = 510/510/510/**75.4** MPa、厚 2.0 mm
   - **提速关键**：`TRUNC_MM=60` 截断胫腓骨 → 32,488→**9,109 单元**，
     **900+ s → 5.4 s（>160 倍）**。踝关节模型不需要整根胫骨（LRN-027）
   - **另三个必要条件**：不写空 Surface/NodeSet + 手册 §3.3 的 Control
     （含 `qn_method BFGS` + 显式容差）+ 位移控制（力控制病态）
   - ❌ **负结果**：本模型 ΔL = **0.001~0.08 mm**，且 **`胫跟` 缩短最多（最后一名）**；
     吴恺是 **`胫跟` 拉伸最大（5.24 mm，第一名）** —— **量级差 2 个数量级 + 名次颠倒**。
     纯轴向压缩下胫骨朝固定足骨**靠近**⇒胫-足韧带**都在缩短**；
     要让 `胫跟` 拉伸必须关节分离 ⇒ **纯压缩做不到**。
   - **下一步最高价值**：**拿到吴恺「图 4」的原始数值**（现文本抽取不含）
     + 确认「位移（拉伸长度）」的定义（长度变化 / 节点位移模 / 相对位移向量模）
   """
t = t[:i0] + new + t[i1:]

P.write_text(t, encoding="utf-8")
print(f"FE_PIPELINE_TODO.md {b:,} -> {len(t.encode('utf-8')):,}")
