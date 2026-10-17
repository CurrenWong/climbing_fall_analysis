"""把两条精确性结论补进 ANKLE_FIX_TODO（不动横幅）。"""
import io

P = r"D:\Project\climbing_fall_analysis\ANKLE_FIX_TODO.md"
ADD = """

### 📌 压缩前补记的两条精确性结论（2026-10-05 末）

1. **压缩标定曲线有配置依赖性**：§2 的 −132.56/−292.22/−479.60/**−608.82**/−741.55/−1017.36 N
   **只在 `CART_NEW=1`+`TIE=1`+`CART_DIR=temp/thums/cart_m1`（新软骨）下成立**。
   换回**旧软骨格式**（不设 `CART_DIR`）同一 `DISP_MM=-1.0` 给 **−465.88 N**（差 24%）
   ⇒ **换软骨必须重标定**。
2. **旧格式 + `TRUNC_MM=60` 路径本就不 NORMAL**（`DISP=-0.1472` 给 236.14 N、非 NORMAL）
   —— **与本轮改动无关**（本轮前后逐字相同），已确认非回归，但**未查明根因**，列为待办。
3. **可支持并已复验的零回归口径**：`TIE=1`+`CART_NEW=1`+`cart_m1`+`DISP_MM=1.6`
   ⇒ **NORMAL / 592.6020 N，本轮多次逐字一致** ✅
   （⚠️ 该 592.6 N 是**拉伸侧**值，仅作"输出可复现"的回归基准，**不是物理结论**。）
"""
s = io.open(P, encoding="utf-8").read()
io.open(P, "a", encoding="utf-8").write(ADD)
print(f"ANKLE_FIX_TODO.md {len(s)} -> {len(s)+len(ADD)}")
