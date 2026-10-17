"""使用说明 §2：把"压缩标定表"换成真·压缩曲线 + 显式符号警告。"""
import io
import os

D = r"D:\Project\climbing_fall_analysis\docs\踝FE可用模型使用说明.md"
d = io.open(D, encoding="utf-8").read()
n0 = len(d.encode("utf-8"))
ls = d.splitlines()
i0 = [i for i, l in enumerate(ls) if "压缩侧标定曲线" in l]
print(f"锚点行：{i0}")
if i0:
    k = i0[0]
    # 找到该小节结束（下一个 "## " 或 "### " 开头）
    j = k + 1
    while j < len(ls) and not ls[j].startswith(("## ", "### ")):
        j += 1
    print(f"替换 [{k+1}..{j}]  共 {j-k} 行")
    for l in ls[k:j]:
        print("   " + l[:110])
    NEW = """### ⭐ 真·**压缩**曲线（`DISP_MM > 0`，全部 `NORMAL`）

> 🚨 **符号陷阱（`[ERR-20261006-003]`）**：`DRIVER=dispx` 下 **`DISP_MM > 0` 才是压实**，
> `DISP_MM < 0` 是**把关节拉开**（拉伸）。三证：`.feb` 写入的 z 值 = `−DISP_MM`；
> 顶面 dz 朝距骨（距骨在下方）；距骨全固定时接触力 `+1.0→149.6 N` vs `−1.0→0.0000 N`。
> **⚠️ `LOAD_SIGN` 对 `dispx` 分支无效，不要靠它"纠符号"。**
> **⇒ 本文件旧版那张 `DISP_MM=−0.25…−1.6` 的表是**拉伸**值（由韧带撑），已作废。**

| `DISP_MM` | +0.05 | +0.10 | +0.25 | +0.50 | +1.00 | **+1.60** | +2.50 |
|---|---|---|---|---|---|---|---|
| 两顶面合力 (N) | 36.08 | 58.25 | 124.04 | 230.71 | 418.10 | **592.60** | 615.35 |

**⇒ 600 N 在 `DISP_MM ≈ +1.89 mm` 达成**（在 +1.60/+2.50 之间插值）。
**⇒ 历史零回归值 `TIE=1 CART_NEW=1 DISP_MM=+1.6` = 592.6020 N 逐位复现** —— 它本来就是压实。

### 🧭 载荷路径（`[LRN-20261006-080]`，做对照前必读）

距骨全固定 + 压缩 +1.0 mm，四组对照（全 `NORMAL`）：

| 配置 | 两顶面合力 |
|---|---|
| 韧带开 · 接触开 | **541.61 N** |
| **韧带关** · 接触开 | **149.59 N** |
| 韧带开 · **接触关** | **414.58 N** |
| 两者都关 | 0.0000 N |

- **⇒ 压缩下韧带占 ~77%、接触 ~28%**
- `FIX_TALUS=transverse` 时**几乎全走韧带**（`LIG_MODE=none` ⇒ 0.0042 N）
- **⚠️ 路径对照必须用 `FIX_TALUS=all`**：否则去掉韧带会让模型**近奇异**（实测 `|u|max` 冲到 10⁶）
"""
    ls[k:j] = NEW.splitlines()
    io.open(D, "w", encoding="utf-8").write("\n".join(ls) + "\n")
    print(f"\n使用说明 {n0} -> "
          f"{len(io.open(D, encoding='utf-8').read().encode('utf-8'))} B")
else:
    print("⚠ 未找到锚点")
