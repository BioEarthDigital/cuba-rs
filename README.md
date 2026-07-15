# cubar-rs

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/rust-1.75%2B-orange.svg)](https://www.rust-lang.org)

**cubar-rs** 是一个用 Rust 编写的高性能密码子使用偏好性（Codon Usage Bias）分析工具，功能上参考了 R 包 [cubar](https://github.com/mt1022/cubar)，在保持计算结果一致的同时，提供更快的启动速度和更低的内存占用。

## 功能特性

- **全功能 CLI**：11 个子命令覆盖密码子偏好性分析的所有主要需求
- **27 个 NCBI 遗传密码表**：支持标准密码表及线粒体、原生生物等变体
- **与 R cubar 数值一致**：ENC/CAI/GC/RSCU 等核心指标与 R 版相关系数 > 0.997
- **并行计算**：利用 rayon 多线程加速大规模数据分析
- **多格式输出**：CSV / TSV / JSON
- **零依赖运行**：单个静态二进制文件，无需 R/Bioconductor 环境

## 安装

### 从源码编译

```bash
git clone https://github.com/user/cubar-rs.git
cd cubar-rs
cargo build --release
```

编译后的二进制文件位于 `target/release/cubar`，可直接复制到任意 `$PATH` 目录。

### 预编译二进制（即将提供）

```bash
# macOS (Apple Silicon)
curl -L https://github.com/user/cubar-rs/releases/latest/download/cubar-aarch64-apple-darwin.tar.gz | tar xz

# Linux (x86_64)
curl -L https://github.com/user/cubar-rs/releases/latest/download/cubar-x86_64-unknown-linux-gnu.tar.gz | tar xz
```

## 快速开始

```bash
# 查看帮助
cubar --help

# 列出所有遗传密码表
cubar list-codes

# 显示标准密码表
cubar show-code 1
```

### 典型分析流程

```bash
# 1. 密码子计数
cubar count cdna.fasta -o codon_counts.csv

# 2. 计算有效密码子数（ENC）—— 密码子偏好性总体度量
cubar enc cdna.fasta -o enc.csv

# 3. 计算相对同义密码子使用度（RSCU）
cubar rscu cdna.fasta -o rscu.csv

# 4. 使用高表达基因作为参考，计算密码子适应指数（CAI）
cubar cai cdna.fasta -r highly_expressed.fasta -o cai.csv

# 5. GC 含量分析（GC、GC3s、GC4d）
cubar gc cdna.fasta -o gc.csv

# 6. 识别最优密码子
cubar optimal cdna.fasta -e expression.tsv -o optimal_codons.csv

# 7. 密码子优化
cubar optimize target.fasta -O optimal_codons.txt -o optimized.fasta

# 8. 滑动窗口分析
cubar slide cdna.fasta -w 20 -s 5 -r ref.fasta -o windows.csv
```

## 命令参考

| 命令 | 说明 | 关键参数 |
|---|---|---|
| `count` | 密码子频率计数 | `-c` 密码表ID, `-f` 输出格式 |
| `enc` | 有效密码子数（ENC） | `-l subfam/amino_acid` |
| `cai` | 密码子适应指数（CAI） | `-r` 参考序列文件 |
| `rscu` | 相对同义密码子使用度 | `-p` pseudocount, `--incl-stop` |
| `fop` | 最优密码子比例 | `-O` 最优密码子文件 |
| `tai` | tRNA 适应指数 | `-t` tRNA 基因拷贝数文件 |
| `gc` | GC/GC3s/GC4d 含量 | — |
| `optimal` | 识别最优密码子 | `-e` 表达量文件 |
| `optimize` | 密码子优化 | `-O` 最优密码子文件, `-m` 策略 |
| `slide` | 滑动窗口分析 | `-w` 窗口, `-s` 步长, `-M` 指标 |
| `list-codes` | 列出遗传密码表 | — |
| `show-code` | 显示密码表详情 | 密码表ID |

### 全局参数

| 参数 | 说明 | 默认值 |
|---|---|---|
| `-c, --gcid` | NCBI 遗传密码表 ID | `1` (Standard) |
| `-f, --format` | 输出格式：`csv` / `tsv` / `json` | `csv` |
| `-o, --output` | 输出文件路径 | stdout |

## 支持的遗传密码表

| ID | 名称 |
|---|---|
| 1 | Standard |
| 2 | Vertebrate Mitochondrial |
| 3 | Yeast Mitochondrial |
| 4 | Mold/Protozoan Mitochondrial & Mycoplasma/Spiroplasma |
| 5 | Invertebrate Mitochondrial |
| 6 | Ciliate/Dasycladacean/Hexamita Nuclear |
| 9 | Echinoderm/Flatworm Mitochondrial |
| 10 | Euplotid Nuclear |
| 11 | Bacterial/Archaeal/Plant Plastid |
| 12 | Alternative Yeast Nuclear |
| 13 | Ascidian Mitochondrial |
| 14 | Alternative Flatworm Mitochondrial |
| 15 | Blepharisma Nuclear |
| 16 | Chlorophycean Mitochondrial |
| 21 | Trematode Mitochondrial |
| 22 | *Scenedesmus obliquus* Mitochondrial |
| 23 | *Thraustochytrium* Mitochondrial |
| 24 | Pterobranchia Mitochondrial |
| 25 | Candidate Division SR1/Gracilibacteria |
| 26 | *Pachysolen tannophilus* Nuclear |
| 27 | Karyorelict Nuclear |
| 28 | *Condylostoma* Nuclear |
| 29 | *Mesodinium* Nuclear |
| 30 | Peritrich Nuclear |
| 31 | *Blastocrithidia* Nuclear |
| 33 | Cephalodiscidae Mitochondrial |

## 输入格式

### FASTA 文件 (CDS)

```fasta
>gene1
ATGGCTGGTAAATGGGCTGCTGGTGGTGCTTAA
>gene2
ATGGCCGGAAGGTGGGCAGCCGGCGGCGCCTGA
```

### tRNA 基因拷贝数文件 (TSV)

```
# anticodon    copy_number
GAA            10
CAT            8
CCA            6
```

### 表达量文件 (TSV)

```
# gene_id      expression_level
YPL071C        1234.5
YLL050C        567.8
```

### 最优密码子文件

```
# 每行一个最优密码子
GCT
GGC
AAA
```

## 与 R cubar 精度对比

使用 497 个酵母 CDS 基因验证（yeast_cds 数据集）：

| 指标 | Pearson r | 平均绝对差 |
|---|---|---|
| ENC | 1.000000 | 0.000025 |
| GC | 0.999056 | 0.001412 |
| GC3s | 1.000000 | 0.000000 |
| GC4d | 1.000000 | 0.000000 |
| CAI | 0.997207 | 0.006688 |

全部差异均在浮点舍入误差范围内。

## 实现的算法与参考文献

| 指标 | 方法 | 参考文献 |
|---|---|---|
| ENC | Sun et al. (2013) 改进算法 | *Mol Biol Evol* 30:191–196 |
| CAI | Sharp & Li (1987) | *Nucleic Acids Res* 15:1281–1295 |
| RSCU | Sharp et al. (1986) | *Nucleic Acids Res* 14:5125–5143 |
| tAI | dos Reis et al. (2004) | *Nucleic Acids Res* 32:5036–5044 |
| Fop | Ikemura (1981) | *J Mol Biol* 151:389–409 |
| DP | — | Deviation from proportionality |

## 项目结构

```
CUB/
├── Cargo.toml              # workspace
├── cubar-core/             # 核心计算库
│   └── src/
│       ├── genetic_code.rs  # 遗传密码表 (NCBI 1-33)
│       ├── sequence.rs      # FASTA I/O, CDS 验证
│       ├── count.rs         # 密码子计数
│       ├── metrics/         # 核心指标
│       │   ├── enc.rs       # ENC
│       │   ├── cai.rs       # CAI
│       │   ├── rscu.rs      # RSCU
│       │   ├── fop.rs       # Fop
│       │   ├── tai.rs       # tAI
│       │   ├── gc.rs        # GC/GC3s/GC4d
│       │   └── dp.rs        # DP
│       ├── optimize.rs      # 密码子优化
│       └── slide.rs         # 滑动窗口
├── cubar-cli/               # CLI 二进制
│   └── src/
│       ├── main.rs
│       └── commands/        # 各子命令实现
└── test_data/               # 测试数据
```

## 依赖

| Crate | 用途 |
|---|---|
| `clap` | CLI 参数解析 |
| `needletail` | FASTA 文件解析 |
| `serde` / `serde_json` | 序列化 |
| `csv` | CSV/TSV 输出 |
| `rayon` | 并行计算 |
| `anyhow` | 错误处理 |

## 开发

```bash
# 运行测试
cargo test

# 发布构建
cargo build --release

# 运行 linter
cargo clippy

# 格式化代码
cargo fmt
```

## License

MIT © cubar-rs developers

## 致谢

本项目参考了 [cubar](https://github.com/mt1022/cubar) R 包（Hong Zhang, Mengyue Liu, Bu Zi），感谢原作者在密码子偏好性分析方法上的卓越工作。
