# UPB-CLI 恒星演化计算器

**UPB-PROFESSIONAL RELEASE 1.0** · 版本号 1.0 · 发布日期 2026-10-14

UPB-CLI 是一个用 C 语言实现的单星演化计算命令行程序。它把
[Hurley, Pols & Tout (2000)](https://doi.org/10.1046/j.1365-8711.2000.03426.x)
的 SSE 解析演化公式（恒星演化拟合公式）完整移植为 C 代码，并以命令行界面
（REPL + 一次性命令）的形式提供：输入初始质量与金属度，输出恒星从零龄主序
（ZAMS）一路演化到白矮星 / 中子星 / 黑洞等残骸的完整轨迹。

零龄主序的光度与半径采用 Tout et al. (1996) 的拟合公式，质量损失按论文第 7 节
叠加 Nieuwenhuijzen & de Jager (1990)、Reimers (1975)、Vassiliadis & Wood (1993)
三种星风处方。

---

## 目录

- [特性](#特性)
- [目录结构](#目录结构)
- [快速开始](#快速开始)
- [命令行参数](#命令行参数)
- [指令参考](#指令参考)
- [使用示例](#使用示例)
- [单位系统与配置文件](#单位系统与配置文件)
- [输出文件格式](#输出文件格式)
- [物理模型](#物理模型)
- [适用范围与限制](#适用范围与限制)
- [源码模块说明](#源码模块说明)
- [从源码编译](#从源码编译)
- [已知问题](#已知问题)
- [参考文献](#参考文献)

---

## 特性

- **纯解析模型，无需积分恒星结构方程** —— 毫秒级完成一次完整演化
  （例如 1 M<sub>☉</sub> 恒星演化到碳氧白矮星约 1136 步，记录 1136 个点）。
- **完整的 SSE 阶段覆盖**：主序 → 赫氏空隙 → 第一巨星支 → 核心氦燃烧 →
  早期 AGB → 热脉动 AGB → 裸氦星 → 各类致密残骸（共 16 种恒星类型）。
- **叠加质量损失**：三套星风处方 + 小包层 Wolf-Rayet 星风，演化过程中实时
  扣减总质量与包层。
- **交互式与批处理两种用法**：既可作为 REPL 交互使用，也可一行命令完成计算，
  还可批量演化最多 100 颗恒星。
- **可切换输出单位**：质量、光度、半径、温度各自支持多种单位（含木星质量、
  地球质量、太阳质量、L<sub>R136a1</sub>、瓦特、AU、km、°C、°F 等），
  设置会持久化到 `upb_config.json`。
- **带自描述头部的文本输出**：结果文件为纯文本表格，头部记录版本、模型、
  初始参数、阶段边界时间、最终状态，便于归档与二次绘图。
- **零依赖**：只依赖 C 标准库与 `libm`，无需任何第三方库。

---

## 目录结构

```
UPB_PROF_RELEASE_1.0/
├── upb_cli.exe              # 已编译好的 Windows 可执行文件（MinGW，64 位）
├── upb_config.json          # 输出单位设置的持久化文件
├── README.md                # 本文件
├── UPB_PROF_TRAILOR.png     # 宣传图（注意：文件内容实为 JPEG，见「已知问题」）
├── results/                 # 演化结果输出目录（自动创建）
│   ├── Sun.txt
│   ├── Trappist-1.txt
│   ├── 比邻星.txt
│   ├── star_M0.33_Z0.02.txt
│   ├── star_M0.7_Z0.02.txt
│   ├── star_M50_Z0.02.txt
│   └── star_M100_Z0.02.txt
└── src/                     # 全部 C 源码
    ├── main.c               # 命令行入口、REPL、指令分发
    ├── cli.c / cli.h        # 单位换算、路径与文件名工具、阶段名
    ├── engine.c / engine.h  # 面向用户的演化引擎（输入校验、结果容器）
    ├── evolve.c / evolve.h  # 演化主循环（步进、阶段切换、质量损失叠加）
    ├── phases.c / phases.h  # 赫罗图诊断 hrdiag 与包层估计 mrenv（论文第 5 节）
    ├── timescales.c / .h    # 时标与阶段边界量（论文附录系数驱动）
    ├── massloss.c / .h      # 质量损失处方与观测换算
    ├── zcoeffs.c / .h       # 随金属度 Z 变化的附录系数表
    ├── zams_coeffs.c / .h   # Tout et al. (1996) 零龄主序拟合
    ├── icon.ico / icon.rc   # 程序图标资源
    └── 本程序采用C语言开发     # 空占位文件
```

---

## 快速开始

### 1. 直接运行（Windows）

双击 `upb_cli.exe`，或从命令行启动：

```bat
upb_cli.exe
```

程序会打印封面并进入交互模式：

```
==============================================================================
  ★  U P B _ C L I  恒星演化计算器  ★
  UPB-PROFESSIONAL RELEASE 1.0
  发布日期: 2026-10-14
==============================================================================
  模型: Hurley, Pols & Tout (2000) MNRAS 315, 543 单星解析演化公式
  指令: 输入 help 查看全部指令, exit / quit / q 退出
==============================================================================

输入 help 查看指令列表; exit / quit / q 退出程序

UPB-CLI>
```

输入指令，例如：

```
UPB-CLI> evolution_single 1.0 0.02 0 0
UPB-CLI> phases 5.0 0.02
UPB-CLI> exit
```

### 2. 一次性执行（不进交互模式）

把指令直接写在命令行参数里，程序执行完自动退出：

```bat
upb_cli.exe --no-banner evolution_single 1.0 0.02 0 0
upb_cli.exe --no-banner phases 5.0 0.02
```

### 3. 一条命令串起多条指令

用 `;` 分隔（括号内的分号不会被切分）：

```bat
upb_cli.exe "unit mass 2; evolution_single 1.0 0.02 0 0"
```

### 4. 计算结果在哪里

所有结果写入**当前工作目录**下的 `results/`（不存在会自动创建），
文件名默认按 `star_M{质量}_Z{金属度}.txt` 生成，例如
`results/star_M1_Z0.02.txt`。

---

## 命令行参数

| 参数 | 说明 |
| --- | --- |
| 无参数 | 打印封面并进入交互模式（REPL） |
| `--no-banner` / `-n` / `--quiet` | 不打印封面 |
| `--interactive` / `-i` | 强制进入交互模式，忽略后续指令 |
| 其它参数 | 视为待执行指令，多个参数用空格连接后按 `;` 拆分逐条执行 |

退出码：全部成功为 `0`；某条指令失败则返回该指令的非零返回值。

---

## 指令参考

| 指令 | 说明 |
| --- | --- |
| `evolution_single [mass] [metallicity] [step_size] [end_time] [filename]` | 演化单颗恒星并写出结果文件 |
| `evolution_multiple (m,z,step,end[,name]) ...` | 批量演化，最多 100 组 |
| `phases [mass] [metallicity]` | 只显示各演化阶段的时间与特征光度，不写文件 |
| `results` | 列出结果目录中的 `.txt` 文件 |
| `results open` | 在文件管理器中打开结果目录 |
| `remove [filename] [--yes]` | 删除某个结果文件，`--yes` 跳过确认 |
| `unit` | 查看当前输出单位 |
| `unit [quantity] [code]` | 设置输出单位 |
| `unit reset` / `unit clear` | 恢复所有单位代码为 1 |
| `help` / `?` | 显示指令列表 |
| `version` | 显示版本、模型与结果目录 |
| `about` | 显示物理模型、演化阶段与质量损失说明 |
| `exit` / `quit` / `q` | 退出程序 |

### `evolution_single` 参数说明

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `mass` | 1.0 | 初始质量，单位 M<sub>☉</sub>，允许 0.077 – 100 |
| `metallicity` | 0.02 | 金属度 Z，允许 0.0001 – 0.03 |
| `step_size` | 0 | 固定步长（Myr）。`0` 表示默认步长 `t_MS / 1000` |
| `end_time` | 0 | 演化终点年龄（Myr）。`0` 表示一直演化到残骸阶段 |
| `filename` | 自动命名 | 输出文件名（自动补 `.txt`，非法字符替换为 `_`） |

交互模式下若不带参数调用 `evolution_single`，程序会逐项询问（直接回车即用
默认值）。

两条内置规则：

- 当 `step_size > 0` 且 `end_time` 超过主序寿命 `t_MS` 时，`end_time` 会被
  自动限制为 `t_MS`，并打印提示。
- 默认步长下，如果某个阶段内步数超过 60000 步（核心质量/光度停止变化），
  演化会被判定为停滞并截断，结果中给出 `[注意] 演化被截断 (未到残骸)`。

### 批量演化语法

每组参数写在括号里，字段用逗号分隔，组与组之间用空格分隔：

```
evolution_multiple (1.0,0.02,0,0) (2.0,0.02,0,0) (5.0,0.02,0.1,120,my_star)
```

每组字段为 `质量, 金属度, 步长, 终点[, 文件名]`，字段数不足 4 组会被跳过。

---

## 使用示例

以下输出均为 `upb_cli.exe` 的真实运行结果。

### 查看 1 M<sub>☉</sub>、Z = 0.02 的各阶段时间

```
> phases 1.0 0.02
恒星参数: M = 1 Msun, Z = 0.02  (X=0.7000, Y=0.2800)
主序寿命 t_MS = 11003.13 Myr

演化阶段 (时间单位: Myr)
--------------------------------------------------------------------------
阶段                            开始          结束      持续时间
--------------------------------------------------------------------------
主序 (MS)                       0.00      11003.13      11003.13
赫氏空隙 (HG)               11003.13      11582.24        579.11
第一巨星支 (GB)             11582.24      12311.14        728.89
核心氦燃烧 (CHeB)           12311.14      12463.95        152.81
早期 AGB (EAGB)             12463.95      12469.63          5.68
热脉动 AGB (TPAGB)          12469.63      12474.24          4.61
--------------------------------------------------------------------------
核演化终点 t_nuc = 12474.24 Myr

特征光度 (单位: Lsun)
  L_ZAMS  = 0.697717
  L_TMS   = 2.1197
  L_BGB   = 2.51543
  L_HeI   = 249.547
  L_BAGB  = 117.851
  临界质量: Mhook=1.0185  MHeF=1.9950  MFGB=13.0324 Msun
```

### 演化一颗 5 M<sub>☉</sub> 恒星，指定步长、终点与文件名

```bat
upb_cli.exe --no-banner evolution_single 5.0 0.02 0.1 120 my_star
```

程序结束时会打印摘要（最终年龄、阶段、质量、光度、半径、温度、光谱型、
输出文件路径）。

### 切换输出单位后演化

```bat
upb_cli.exe "unit radius 3; unit temperature 2; evolution_single 1 0.02 0 0"
```

此后半径以 AU、温度以 °C 写入结果文件，设置同时保存到 `upb_config.json`。

### 批量演化并列出结果

```bat
upb_cli.exe --no-banner "evolution_multiple (0.5,0.02,0,0) (1.0,0.02,0,0) (20,0.02,0,0); results"
```

---

## 单位系统与配置文件

程序内部一律使用 M<sub>☉</sub>、L<sub>☉</sub>、R<sub>☉</sub>、K、Myr
计算，仅在**写出结果文件与打印摘要**时按单位代码换算。

| 类别 | 代码 | 单位 | 备注 |
| --- | --- | --- | --- |
| `mass` | 1 / 2 / 3 / 4 | M<sub>☉</sub> / M<sub>Jupiter</sub> / M<sub>Earth</sub> / kg | 默认 1 |
| `luminosity` | 1 / 2 / 3 | L<sub>☉</sub> / L<sub>R136a1</sub> / W | L<sub>R136a1</sub> = 8×10⁶ L<sub>☉</sub> |
| `radius` | 1 / 2 / 3 / 4 / 5 | R<sub>☉</sub> / R<sub>Stephenson_218</sub> / AU / R<sub>Earth</sub> / km | R<sub>Stephenson_218</sub> = 2150 R<sub>☉</sub> |
| `temperature` | 1 / 2 / 3 | K / °C / °F | 仿射换算 |

别名：`lum` 等价于 `luminosity`，`temp` 等价于 `temperature`。

单位设置保存在工作目录的 `upb_config.json`：

```json
{
  "mass": 1,
  "luminosity": 1,
  "radius": 1,
  "temperature": 1
}
```

程序启动时读取该文件（文件不存在或损坏则保持默认值 1），每次 `unit` 设置后
立即回写。因此该文件是「当前运行目录」的语言环境状态，删掉即恢复默认。

---

## 输出文件格式

结果文件是带 `#` 注释头部的纯文本表格，可直接用 `numpy.loadtxt`、`gnuplot`、
Excel 等读取。以 `results/Sun.txt` 为例：

```
# ============================================================================
# UPB-CLI 恒星演化输出文件
# 程序版本   : UPB-PROFESSIONAL RELEASE 1.0  发布日期: 2026-10-14
# 物理模型   : Hurley, Pols & Tout (2000) MNRAS 315, 543 SSE
# ----------------------------------------------------------------------------
# 初始质量   : 1 Msun
# 金属度 Z   : 0.02   (X = 0.7000, Y = 0.2800)
# 步长       : 11.0031 Myr
# 主序寿命   : 11003.1302 Myr
# 红巨星底端 : 11582.2424 Myr
# 氦点火     : 12311.1362 Myr
# 演化终点   : 残骸阶段
# 演化点数   : 1136   步数: 1136
# 最终状态   : t = 12488.8644 Myr, 阶段 = 碳氧白矮星 (type 11)
# 最终质量   : 0.512227 Msun   核心质量: 0.512227 Msun
# 生成时间   : 2026-10-14 12:52:52
# 单位代码   : mass=1 lum=1 radius=1 temp=1
# ----------------------------------------------------------------------------
# 数据列: 1=年龄(Myr) 2=亮度 3=半径 4=温度 5=质量
# ============================================================================
       年龄(Myr)                 亮度                 半径             温度                 质量
          0.0000         0.6977165691         0.8882494503          5598.42                    1
         11.0031         0.6981025936         0.8884707058          5598.49                    1
         22.0063          0.698489523         0.8886920174          5598.57                    1
...
```

数据列固定为 5 列：

| 列 | 含义 | 单位 |
| --- | --- | --- |
| 1 | 年龄 | Myr |
| 2 | 光度 | 由 `unit luminosity` 决定 |
| 3 | 半径 | 由 `unit radius` 决定 |
| 4 | 有效温度 | 由 `unit temperature` 决定 |
| 5 | 当前总质量 | 由 `unit mass` 决定 |

---

## 物理模型

### 演化阶段（论文第 4 节的恒星类型 `kw`）

| kw | 阶段 | kw | 阶段 |
| --- | --- | --- | --- |
| 0 | 主序（全对流，M ≲ 0.7 M<sub>☉</sub>） | 8 | 裸氦星赫氏空隙 |
| 1 | 主序 | 9 | 裸氦星巨星支 |
| 2 | 赫氏空隙 | 10 | 氦白矮星 |
| 3 | 第一巨星支 | 11 | 碳氧白矮星 |
| 4 | 核心氦燃烧 | 12 | 氧氖白矮星 |
| 5 | 早期渐近巨星支（EAGB） | 13 | 中子星 |
| 6 | 热脉动 AGB（TPAGB） | 14 | 黑洞 |
| 7 | 裸氦星主序 | 15 | 无质量残骸 |

kw ≥ 10 视为残骸，演化在主循环中到达残骸即停止。

### 质量损失（论文第 7 节）

四套处方叠加，单位为 M<sub>☉</sub>/yr，程序内部换算为 M<sub>☉</sub>/Myr：

1. **Nieuwenhuijzen & de Jager (1990)** —— 全赫罗图的高光度星风；
2. **Reimers (1975)** —— 巨星阶段（kw 2–9），η = 0.5；
3. **Vassiliadis & Wood (1993)** —— 脉动 AGB（kw 5、6）；
4. **减弱的 Wolf-Rayet 星风** —— 包层被剥离、暴露裸核的恒星。

主序与赫氏空隙上，序列参数会随当前质量重新锚定（论文第 7.1 节）；半径公式与
质量损失率始终使用当前质量。

### 数值方法

- 默认步长为 `t_MS / 1000`，各阶段统一使用；
- 使用默认步长时，步进会自动缩短以避免跨过阶段边界；
- 步数上限 10⁵；单阶段步数超过 60000 判定为停滞并终止；
- 到达残骸（kw 10–15）或质量归零即结束；
- 每次演化记录全部轨迹点（`record = 1`），并在结果头中给出阶段边界时间。

---

## 适用范围与限制

- **质量范围**：0.077 – 100 M<sub>☉</sub>（0.077 M<sub>☉</sub> 约为氢燃烧下限）。
- **金属度范围**：Z = 0.0001 – 0.03，太阳取 Z = 0.02。
- **元素丰度**：按 `Y = 0.24 + 2Z`、`X = 0.76 − 3Z` 给定
  （Z = 0.02 时 X = 0.70、Y = 0.28）。Z > 0.0267 时 X 会变为负值，
  属于模型外推，请勿将其视为物理结果。
- **单星模型**：不考虑双星相互作用、潮汐、自转、磁场与化学演化。
  `mass_loss_rate()` 保留 Roche 瓣半径参数，但单星演化固定传 0。
- **解析拟合的精度**：SSE 公式本身是对恒星演化计算的拟合，
  质量损失与晚期阶段（尤其是 TPAGB 与残骸质量）的不确定度较大。
- **结果目录**：按实现，`results/` 与 `upb_config.json` 都相对于
  **当前工作目录**解析，而不是可执行文件所在目录。
  想固定输出位置，请先 `cd` 到目标目录再运行程序。
- **控制台编码**：程序输出为 GBK（代码页 936）中文。
  在默认的简体中文 Windows 控制台可直接正常显示；
  若终端被设为 UTF-8（65001），中文可能显示为乱码。

---

## 源码模块说明

| 文件 | 职责 |
| --- | --- |
| `src/main.c` | 封面、指令分词（`;` 拆分、括号分组）、REPL、全部指令实现、结果文件写出 |
| `src/cli.c` / `cli.h` | 单位换算表与符号、结果目录与文件名校验、`upb_config.json` 读写、阶段名表 |
| `src/engine.c` / `engine.h` | 对外引擎 API：`engine_evolve`、`engine_get_main_sequence_lifetime`、`engine_get_phase_times`、批量演化；`StellarState` / `EvolutionResult` 容器；光谱型判定（O/B/A/F/G/K/M 阈值 30000/10000/7500/6000/5200/3700 K） |
| `src/evolve.c` / `evolve.h` | 演化主循环：步进、阶段切换、边界不跨越、质量损失叠加、停滞保护 |
| `src/phases.c` / `phases.h` | `hrdiag`（论文第 5 节：给定年龄返回 L、R、核心质量、类型）与 `mrenv`（第 7.2 节对流包层与回转半径） |
| `src/timescales.c` / `timescales.h` | `Sequence`：与时间无关的 (M, Z) 派生量 —— 临界质量、特征光度、阶段边界时间、巨星支参数 |
| `src/massloss.c` / `massloss.h` | 质量损失率、由 L/R 求 T<sub>eff</sub>、光谱型字符串、热改正（仅报告用） |
| `src/zcoeffs.c` / `zcoeffs.h` | 论文附录中随 Z 变化的系数 a1–a81、b1–b57 及与 Z 无关常数 |
| `src/zams_coeffs.c` / `zams_coeffs.h` | Tout et al. (1996) 的 L<sub>ZAMS</sub>、R<sub>ZAMS</sub> 拟合多项式 |
| `src/icon.rc` / `src/icon.ico` | Windows 可执行文件图标资源 |

模块划分与原 Python 项目一一对应：`main.py` → `main.c` + `cli.c`，
`engine.py` → `engine.c`，`sse/evolve.py` → `evolve.c`，
`sse/phases.py` → `phases.c`，`sse/timescales.py` → `timescales.c`，
`sse/massloss.py` → `massloss.c`，`sse/zcoeffs.py` → `zcoeffs.c`，
`sse/zams_coeffs.py` → `zams_coeffs.c`。

---

## 从源码编译

发行包中已包含可直接运行的 `upb_cli.exe`，一般无需自行编译。
若需重新构建，使用 MinGW-w64 的 gcc：

```bat
cd src
gcc -O2 -std=c99 -o upb_cli.exe main.c cli.c engine.c evolve.c phases.c ^
    timescales.c massloss.c zcoeffs.c zams_coeffs.c -lm
```

带图标资源的完整构建（需要 `windres`）：

```bat
cd src
windres icon.rc -O coff -o icon.o
gcc -O2 -std=c99 -o upb_cli.exe main.c cli.c engine.c evolve.c phases.c ^
    timescales.c massloss.c zcoeffs.c zams_coeffs.c icon.o -lm
```

说明：

- 源码使用 C89/C99 兼容写法，只依赖 C 标准库与 `-lm`；
- Windows 下 `results open` 依赖 `start`，类 Unix 下依赖 `xdg-open` / `ls`，
  这些外部命令通过 `system()` 调用，属于平台相关分支；
- 加上 `-Wall` 会产生若干「定义了但未使用」的函数/变量告警
  （如 `count_rows`、`is_remnant_phase`），不影响构建与运行。

---

## 已知问题

- `UPB_PROF_TRAILOR.png` 的实际内容是 JPEG（文件头 `FF D8 FF E0 ... JFIF`），
  扩展名与格式不符。部分严格校验格式的查看器/工具可能拒绝打开。
- `src/本程序采用C语言开发` 是 0 字节的占位文件，不含内容。
- 结果目录与配置文件按当前工作目录解析（见「适用范围与限制」），
  从其它目录启动程序时，`results/` 会出现在那个目录下。
- `results` 指令依赖系统 `dir` / `ls` 输出，在受限环境下可能无输出。

---

## 参考文献

1. Hurley, J. R., Pols, O. R., & Tout, C. A. (2000).
   *Comprehensive analytic formulae for stellar evolution as a function of mass
   and metallicity.* MNRAS, 315, 543.
2. Tout, C. A., Pols, O. R., Eggleton, P. P., & Han, Z. (1996).
   *Zero-age main-sequence radii and luminosities as analytic functions of mass
   and metallicity.* MNRAS, 281, 257.
3. Nieuwenhuijzen, H., & de Jager, C. (1990). A&A, 231, 134.
4. Reimers, D. (1975). Mém. Soc. R. Sci. Liège, 8, 369.
5. Vassiliadis, E., & Wood, P. R. (1993). ApJ, 413, 641.
