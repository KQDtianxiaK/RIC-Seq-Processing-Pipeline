# RIC-seq 数据处理流程

[English](README.md) | **中文**

本流程将原始 RIC-seq（RNA 原位构象测序）FASTQ 数据经过质量控制、接头修剪、PCR 重复去除、rRNA 过滤、基因组比对、配对末端标签收集等一系列步骤处理，最终鉴定分子内和分子间的 RNA-RNA 相互作用。整个工作流程完全容器化（使用 Singularity）并由 Snakemake 调度，确保可重复性，并严格遵循 **Nature Protocols 2021** RIC-seq 标准。

---

## 目录

- [RIC-seq 简介](#ric-seq-简介)
- [工作流程概述](#工作流程概述)
- [系统要求](#系统要求)
- [安装 Singularity](#安装-singularity)
- [安装 Snakemake](#安装-snakemake)
- [参考数据准备](#参考数据准备)
- [构建 Singularity 容器](#构建-singularity-容器)
- [运行流程](#运行流程)
- [输出目录结构](#输出目录结构)
- [结果解读](#结果解读)
- [版本更新记录](#版本更新记录)
- [常见问题排查](#常见问题排查)
- [重要注意事项](#重要注意事项)
- [引用](#引用)
- [许可证](#许可证)

---

## RIC-seq 简介

RIC-seq（RNA In situ Conformation sequencing，RNA 原位构象测序）是一种基于邻近连接（proximity ligation）的技术，用于在活细胞中捕获 RNA-RNA 相互作用，同时保留 RNA 分子的空间信息。该技术可以鉴定：

- **分子内 RNA 相互作用（intra-molecular）**：同一 RNA 分子内的二级/三级结构
- **分子间 RNA 相互作用（inter-molecular）**：不同 RNA 分子之间的相互作用

---

## 工作流程概述

本流程严格遵循 **Nature Protocols 2021** 计算分析工作流（第 175–186 步），主要分析步骤如下：

### 预处理与质量控制（第 175–179 步）

| 步骤 | 工具 | 说明 |
|------|------|------|
| 175–176：原始数据 QC | FastQC | 评估原始测序数据质量 |
| 177：接头修剪 | Trim Galore | 去除接头序列，过滤低质量碱基 |
| 178：PCR 去重 | 自定义 Perl 脚本 | 去除 PCR 重复，避免分析偏差 |
| 179：低复杂度过滤 | cutadapt | 去除 poly-N 末端和同聚物序列 |
| 最终 QC | FastQC | 验证预处理后的数据质量 |
| rRNA 过滤 | STAR | 比对至 rRNA 索引，去除核糖体 RNA 污染 |

### 核心分析流程（第 180–186 步）

| 步骤 | 工具/脚本 | 说明 |
|------|-----------|------|
| 180–181：基因组比对 | STAR | 对 R1 和 R2 分别进行比对，检测嵌合读取（chimeric reads） |
| 182：配对末端标签收集 | collect_pair_tags.pl | 收集并合并配对比对信息 |
| 183：分子内/间分离 | BEDtools + 自定义脚本 | 区分分子内和分子间相互作用 |
| 184：分子内分类 | category_intra_reads.pl | 将分子内读取分类为嵌合读取和单体读取 |
| 185：分子内聚类 | cluster_intra_reads.pl | 对嵌合读取进行聚类（连接评分过滤 ≥0.01） |
| 186：分子间网络分析 | Monte Carlo 模拟 | 蒙特卡罗模拟 + 局部多重检验校正 |

整个流程使用 Singularity 容器化，所有生物信息学工具和 Perl 模块均已预装，只需一条 Snakemake 命令即可执行完整流程。

---

## 系统要求

### 推荐系统配置

- **CPU**：16 核处理器（最低 8 核）
- **内存**：64 GB RAM（最低 32 GB）
- **存储**：每个样品至少 200 GB 可用空间

### 软件依赖

| 软件 | 版本 | 用途 |
|------|------|------|
| Singularity | ≥4.0 | 容器运行时 |
| Snakemake | ≥7.0 | 流程调度 |
| Python | ≥3.8 | Snakemake 依赖 |

容器内已预装工具（无需手动安装）：

| 工具 | 版本 | 用途 |
|------|------|------|
| FastQC | 0.12.1 | 质量控制 |
| STAR | 2.7.11b | 基因组比对 |
| SAMtools | 1.22.1 | SAM 文件处理 |
| BEDtools | 2.28.0 | 基因组区间操作 |
| Trim Galore | 0.6.10 | 接头修剪 |
| cutadapt | 3.4 | 低复杂度过滤 |
| MultiQC | 1.19 | QC 报告汇总 |

---

## 安装 Singularity

以下为 Ubuntu 22.04 的详细安装步骤。其他操作系统请参考 [Singularity 官方安装指南](https://docs.sylabs.io/guides/latest/user-guide/)。

**第 1 步：安装系统依赖**

```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    libfuse3-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    curl wget git
```

**第 2 步：安装 Go 语言**

```bash
wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
rm go1.21.3.linux-amd64.tar.gz

echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc
```

**第 3 步：下载、编译并安装 Singularity**

```bash
cd /path/to/software
wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz
tar -xvzf singularity-ce-4.0.1.tar.gz
cd singularity-ce-4.0.1
./mconfig
cd builddir
make
sudo make install
```

**第 4 步：验证安装**

```bash
singularity --version
singularity -h
```

---

## 安装 Snakemake

Snakemake 需要 Python 3，可通过 pip 安装：

```bash
pip install snakemake
```

---

## 参考数据准备

流程需要以下人类基因组参考文件（以 hg38 为例）：

### STAR 基因组索引

```bash
mkdir -p References/hg38
cd References

# 下载基因组 FASTA 和 GTF 注释文件
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

# 解压
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz

# 构建 STAR 索引（需要 Singularity 容器）
singularity exec --cleanenv RIC-seq.sif STAR \
    --runMode genomeGenerate \
    --genomeDir ./hg38/STAR_index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v38.primary_assembly.annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 16
```

### rRNA 索引

```bash
# 下载 rRNA 序列（示例：NCBI 45S pre-rRNA）
wget -O rRNA.fa "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=555853"

# 构建 rRNA STAR 索引
singularity exec --cleanenv RIC-seq.sif STAR \
    --runMode genomeGenerate \
    --genomeDir ./hg38/rRNA_STAR_index \
    --genomeFastaFiles rRNA.fa \
    --genomeSAindexNbases 8 \
    --runThreadN 16
```

### 基因注释文件

流程需要两个特殊的 BED 文件：

- **Junction BED**：包含剪接位点配对信息
- **基因注释 BED**：包含全基因区域信息和基因 ID（已提供示例文件：`References/hg38/whole_gene_region.bed`）

生成 Junction BED 文件：

```bash
# GTF 转 BED
perl scripts/gtf_to_bed.pl gencode.v38.primary_assembly.annotation.gtf > gencode.v38.annotation.bed

# 生成 Junction BED
perl scripts/creat_junction_bed.pl gencode.v38.annotation.bed > gencode.v38.all_exon_junction.bed
```

生成全基因区域 BED 文件（第 4 列必须包含基因 ID）：

```bash
awk 'BEGIN{OFS="\t"} $3=="gene" {print $1, $4-1, $5, $10, "60", $7}' \
    gencode.v38.primary_assembly.annotation.gtf | \
    sed 's/"//g' | sed 's/;//g' > whole_gene_region.bed
```

### 接头序列文件

准备包含 Illumina 接头序列的 FASTA 文件：

```
cat > adapter.fa << 'EOF'
>PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
EOF
```

---

## 构建 Singularity 容器

使用仓库提供的定义文件（`Containers/RIC-seq.def`）构建容器：

```bash
# 构建容器（需要 sudo 权限）
sudo singularity build RIC-seq.sif Containers/RIC-seq.def
```

容器内包含：
- **生物信息学工具**：FastQC v0.12.1、STAR v2.7.11b、SAMtools v1.22.1、BEDtools v2.28.0
- **质控工具**：Trim Galore v0.6.10、cutadapt v3.4、MultiQC v1.19
- **Perl 模块**：List::Util、File::Spec、Getopt::Long、Graph::Undirected、Statistics::Distributions、Math::CDF
- **分析脚本**：完整的 RICpipe 脚本，预装于容器内 `/mnt/RICseq_scripts_root/`

---

## 运行流程

### 必要的文件目录结构

```
project_directory/
├── Scripts/
│   ├── config.yaml        # 配置文件
│   └── RICseq.smk         # Snakemake 工作流脚本
├── Containers/
│   └── RIC-seq.sif        # Singularity 容器
├── References/
│   ├── hg38/
│   │   ├── STAR_index/           # STAR 基因组索引
│   │   └── rRNA_STAR_index/      # rRNA STAR 索引
│   ├── gencode.v38.all_exon_junction.bed  # 剪接位点注释
│   └── whole_gene_region.bed              # 全基因区域注释
└── rawdata/
    ├── sample1_R1.fastq.gz
    └── sample1_R2.fastq.gz
```

### 配置文件（config.yaml）

编辑 `Scripts/config.yaml`，指定样品路径和分析参数：

```yaml
# 样品配置（支持多样品）
samples:
  sample1:
    R1: "/path/to/rawdata/sample1_R1.fastq.gz"
    R2: "/path/to/rawdata/sample1_R2.fastq.gz"
  sample2:
    R1: "/path/to/rawdata/sample2_R1.fastq.gz"
    R2: "/path/to/rawdata/sample2_R2.fastq.gz"

# 输出目录
outputdir: "/path/to/output"

# Singularity 容器路径
sif: "/path/to/Containers/RIC-seq.sif"

# 并发线程数
threads: 16

# 参考文件
rRNA_index: "/path/to/References/hg38/rRNA_STAR_index"
star_index: "/path/to/References/hg38/STAR_index"
junction_bed: "/path/to/References/hg38/gencode.v38.all_exon_junction.bed"
gene_annotation_bed: "/path/to/References/hg38/whole_gene_region.bed"

# 脚本路径（容器内路径，无需修改）
scripts_root: "/mnt/RICseq_scripts_root"

# 分析参数
# 第 3/4 步：分子内相互作用参数
fragment_len_cutoff: 1000            # 分类用片段长度阈值（bp）
step4_fragment_cutoff: 2             # 聚类所需最小唯一片段数
step4_connection_score_cutoff: 0.01  # 连接评分阈值（Nature Protocols 默认值）

# 第 5 步：分子间相互作用参数
pvalue_cutoff: 0.05       # 显著相互作用 p 值阈值
monte_carlo_iters: 100000 # Monte Carlo 迭代次数（Nature Protocols：100,000）
monte_carlo_threads: 16   # Monte Carlo 模拟并行线程数
```

### 执行流程

**第 1 步：预运行（可选，推荐）**

预览工作流但不实际执行：

```bash
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --cores 16 -n
```

**第 2 步：执行流程**

```bash
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --singularity-args "-B /path/to/project_directory" \
    --cores 16 -p --rerun-incomplete \
    --latency-wait 30
```

**第 3 步：后台执行（推荐）**

```bash
nohup snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --singularity-args "-B /path/to/project_directory" \
    --cores 16 -p --rerun-incomplete \
    --latency-wait 30 > pipeline.log 2>&1 &
```

### 常用 Snakemake 参数说明

| 参数 | 说明 |
|------|------|
| `--configfile` | 配置 YAML 文件路径 |
| `--use-singularity` | 在 Singularity 容器内执行规则 |
| `--singularity-args "-B"` | 挂载宿主机目录到容器（多路径用逗号分隔） |
| `--cores` | 最大 CPU 核心数 |
| `-p` | 打印执行的 shell 命令 |
| `-n` | 预运行（预览工作流，不执行） |
| `--rerun-incomplete` | 重新运行上次未完成的任务 |
| `--latency-wait` | 等待输出文件出现的时间（秒，处理文件系统延迟） |

**解锁目录（如果流程中断）：**

```bash
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml --unlock
```

---

## 输出目录结构

```
output/
├── qc/                          # 质量控制目录
│   ├── raw/                     # 原始读取 FastQC（第 175-176 步）
│   │   ├── {sample}_R1_fastqc.{html,zip}
│   │   └── {sample}_R2_fastqc.{html,zip}
│   └── final/                   # 处理后读取 FastQC（验证）
│       ├── {sample}_read1_fastqc.{html,zip}
│       └── {sample}_read2_fastqc.{html,zip}
├── trim/                        # 修剪和低复杂度过滤
│   ├── {sample}_R{1,2}_trimming_report.txt      # 修剪报告
│   ├── {sample}_R{1,2}_trimmed.fq.gz            # Trim Galore 后结果
│   └── {sample}_read{1,2}.clean.rmDup.rmPoly.fq # cutadapt 后结果
├── dedup/                       # PCR 去重（第 178 步）
│   ├── {sample}_read{1,2}.clean.rmDup.fq
│   └── {sample}_dedup_stats.txt
├── rRNA_filter/                 # rRNA 过滤结果（第 180 步）
│   ├── {sample}.no_rRNA.{1,2}.fq
│   └── {sample}.Log.final.out
├── alignment/                   # STAR 比对（第 181 步）
│   ├── {sample}.read{1,2}.Aligned.out.sam   # 普通比对结果
│   ├── {sample}.read{1,2}.Chimeric.out.sam  # 嵌合读取结果
│   └── {sample}.read{1,2}.Log.final.out     # 比对统计日志
├── step1_pair_tags/             # 配对末端标签收集（第 182 步）
│   ├── {sample}.interaction.sam
│   └── {sample}_num_of_interactions.list
├── step2_separate/              # 分子内/间分离（第 183 步）
│   ├── {sample}.intraMolecular.sam
│   ├── {sample}.interMolecular.sam
│   └── {sample}.pets_in_same_gene.list
├── step3_category/              # 分子内分类（第 184 步）
│   ├── {sample}.intraMolecular.Chimeric.sam
│   └── {sample}.intraMolecular.Singleton.sam
├── step4_intra_cluster/         # 分子内聚类（第 185 步）
│   └── {sample}.cluster.withScore.highQuality.list  ⭐ 关键结果
├── step5_inter_network/         # 分子间网络分析（第 186 步）
│   ├── {sample}.merged.network
│   ├── {sample}_sim/            # Monte Carlo 模拟结果
│   │   ├── 1.base_on_observed/  # 基于观测数据的模拟
│   │   └── 2.base_on_random/    # 基于随机数据的模拟
│   ├── {sample}_post/           # 后处理
│   │   ├── 2.pre-process/
│   │   ├── 3.calculate_pvalue/
│   │   └── 4.recalibrate_pvalue/
│   └── {sample}.significant.interMolecular.interaction.list  ⭐ 关键结果
├── multiqc/                     # 汇总质量控制报告
│   ├── multiqc_report.html
│   └── multiqc_data/
└── logs/                        # 各步骤详细日志
    └── {sample}_{step}.log
```

⭐ **关键结果文件：**
- **分子内**：`{sample}.cluster.withScore.highQuality.list`
- **分子间**：`{sample}.significant.interMolecular.interaction.list`

---

## 结果解读

### 关键结果文件说明

**分子内相互作用：**

- **`{sample}.cluster.withScore.highQuality.list`** — 连接评分 ≥0.01 的高置信度分子内 RNA 结构聚类，代表单个转录本内的 RNA 二级/三级结构。

**分子间相互作用：**

- **`{sample}.significant.interMolecular.interaction.list`** — 经 Monte Carlo 模拟和局部多重检验校正后，统计学显著（p<0.05）的分子间 RNA-RNA 相互作用，代表不同转录本之间的 RNA-RNA 相互作用。

### 质量控制指标解读

**FastQC 报告（两个质控检查点）：**

1. **原始数据质控（`qc/raw/`）** — 任何处理前的质量评估：
   - 逐碱基质量评分（基础质量水平）
   - 接头含量检测
   - 过度代表序列
   - GC 含量分布
   - **用途**：判断是否需要重新测序

2. **最终质控（`qc/final/`）** — 所有预处理后的质量评估：
   - 验证修剪后质量是否提升
   - 确认接头已去除
   - 检查是否有残余问题
   - **用途**：验证预处理效果

**修剪报告：**

- **`*_trimming_report.txt`** — Trim Galore 统计信息：处理读取数量、检测并去除的接头序列、质量修剪统计

**比对统计：**

- **`*.Log.final.out`** — STAR 比对日志，包含：
  - 总读取数
  - 唯一比对率（期望 >70%）
  - 多位点比对读取
  - 嵌合读取（RIC-seq 关键指标，期望 5–15%）

### 期望质量指标

成功的 RIC-seq 实验期望结果：

| 质量指标 | 期望范围 | 说明 |
|---------|---------|------|
| 原始读取质量 | >Q30（大多数碱基） | 良好测序质量 |
| 原始数据接头含量 | <1% | 低接头污染 |
| 修剪后接头含量 | 0% | 接头成功去除 |
| 基因组比对率 | >70% | 高比对成功率 |
| 嵌合读取率 | 5–15% | RIC-seq 关键指标 |
| PCR 重复率 | <30% | 良好文库复杂度 |
| rRNA 污染率 | 10–40% | 取决于样品制备 |
| 分子间相互作用数 | 10³–10⁵ | 取决于测序深度和 p 值阈值 |
| 分子内结构聚类数 | 10²–10⁴ | 结构多样性 |

---

## 版本更新记录

### 版本 3.0 — Nature Protocols 合规版

**主要改进：**

1. **✅ 严格遵循 Nature Protocols** — 完整实现 Nature Protocols 2021 第 175–186 步
2. **✅ 原始数据质控（新增）** — 添加原始读取 FastQC 作为第一步（第 175–176 步）
3. **✅ 双质控检查点** — 原始数据质控 + 处理后质控，全面追踪数据质量
4. **✅ 正确处理顺序** — 修正工作流顺序：
   ```
   FastQC（原始）→ Trim Galore → PCR 去重 → cutadapt → FastQC（最终）→ rRNA 过滤 → 比对
   ```
5. **✅ 增强故障排除** — 早期质量检测，避免在低质量数据上浪费计算资源

**关键工作流变化：**

| 方面 | 旧版本 | v3.0（当前） |
|------|--------|-------------|
| 初始质控 | 仅修剪后 | **任何处理前**（第 175-176 步） |
| 处理顺序 | Trim → cutadapt → FastQC → 去重 | **FastQC → Trim → 去重 → cutadapt → FastQC** |
| 质控检查点 | 1 个（修剪后） | **2 个（处理前 + 处理后）** |
| 协议符合性 | 部分 | **完整（第 175-186 步）** |
| 输出结构 | `qc/` | `qc/raw/` + `qc/final/` |

### 工作流演变历史

| 步骤 | v1 | v2 | v3（当前） |
|------|----|----|-----------|
| 原始 QC | ❌ 无 | ❌ 无 | ✅ **FastQC（第 175-176 步）** |
| 修剪工具 | Trimmomatic | Trim Galore | Trim Galore（第 177 步） |
| 去重时机 | 手动 | FastQC 后 | **修剪后（第 178 步）** |
| Poly-N 去除 | ❌ 无 | cutadapt | cutadapt（第 179 步） |
| 最终 QC | ❌ 无 | FastQC | ✅ **FastQC（验证）** |
| rRNA 过滤 | ❌ 无 | STAR | STAR（第 180 步） |
| Monte Carlo | 1,000 次迭代 | 100,000 次 | 100,000 次（第 186 步） |
| 协议符合性 | 部分 | 部分 | ✅ **100%（第 175-186 步）** |

---

## 常见问题排查

**问题：找不到 FastQC 输出 / MissingOutputException**

解决方案：FastQC 根据输入文件名生成输出。如遇到文件命名问题：

```bash
# 增加文件系统延迟等待时间
snakemake --use-singularity --cores 16 --latency-wait 30

# 检查日志中实际生成的文件名
cat output/logs/{sample}_fastqc_raw.log
```

**问题：原始读取质量差（<Q20）**

解决方案：查看 `qc/raw/` 中的原始 FastQC 报告：
- 若 >50% 碱基中位质量 <Q20，考虑重新测序
- 若仅有接头污染，流程会自动处理
- 若检测到过度代表序列，调查其来源

**问题：比对率低（<50%）**

解决方案：
- 检查原始和最终 FastQC 报告，比较质量变化
- 验证接头是否已成功去除
- 确认 STAR 索引与样品物种匹配
- 检查 rRNA 污染是否过高

**问题：PCR 重复率高（>40%）**

解决方案：
- 检查原始 FastQC 中的重复序列
- 说明文库复杂度低
- 建议在文库制备中使用更多起始材料
- 探索性分析可继续，但需注明此限制

**问题：找不到 Math::CDF 模块**

解决方案：Singularity 容器已包含 Math::CDF。请确保使用正确的容器版本：

```bash
# 验证 Math::CDF 已安装
singularity exec RIC-seq.sif perl -MMath::CDF -e 'print "Math::CDF: OK\n"'
```

**问题：未找到显著相互作用**

解决方案：
- 检查测序深度（建议每个样品 >5000 万配对端读取）
- 降低 `pvalue_cutoff`（探索性分析可尝试 0.1）
- 确保处理了生物学重复以提高统计效力
- 确认 Monte Carlo 模拟已成功完成（应有 16 个线程输出）
- 检查 `qc/raw/` 中的原始数据质量——数据质量差会影响下游结果

**问题：第 3 步输出为空**

解决方案：
- 查看 `output/logs/{sample}_step3.log` 日志
- 验证基因注释 BED 和 Junction BED 文件格式是否正确
- 确保第 2 步产生了足够的分子内相互作用
- 检查最终 QC 报告，验证输入数据质量

**问题：流程卡住或运行缓慢**

解决方案：
- 检查系统资源（CPU、RAM、磁盘 I/O）
- 第 5 步 Monte Carlo 模拟计算量大（默认参数下需 2–6 小时）
- 增加 `monte_carlo_threads` 以充分利用可用 CPU 核心
- 测试时可减少 `monte_carlo_iters`（最小 10,000 次）

**问题：Singularity 绑定错误**

解决方案：
- 使用 `--singularity-args "-B /path1,/path2"` 确保挂载所有必要路径
- 路径必须在宿主机上存在
- 配置文件中使用绝对路径

---

## 重要注意事项

### 参考文件格式要求

- **Junction BED**：必须使用 `creat_junction_bed.pl` 生成的配对剪接位点格式（共 12 列）
- **基因注释 BED**：第 4 列必须包含基因 ID（而非转录本 ID）
- **rRNA 索引**：必须专门针对 rRNA 序列构建（小基因组使用较小的 `genomeSAindexNbases`）

### 脚本路径配置

`scripts_root` 参数指向容器内 `/mnt/RICseq_scripts_root`，所有 Perl 分析脚本均位于此处：

```
/mnt/RICseq_scripts_root/
├── step0.remove_PCR_duplicates/     # PCR 去重脚本
├── step1.collect_pair_tags/         # 配对标签收集脚本
├── step2.separate_intra_inter_molecular/  # 分子内/间分离脚本
├── step3.category_intra_reads/      # 分子内分类脚本
├── step4.cluster_intramolecular/    # 分子内聚类脚本
└── step5.screen_high-confidence_intermolecular/  # 高置信度分子间筛选脚本
```

### 计算资源估算

| 步骤 | 预计时间（50M 读取对） |
|------|----------------------|
| 第 175–179 步（预处理） | 1–2 小时 |
| 第 180–184 步（比对与分离） | 2–3 小时 |
| 第 185 步（聚类） | 0.5–1 小时 |
| 第 186 步（Monte Carlo，默认参数） | 2–6 小时 |
| **总计** | **约 6–12 小时** |

### 存储空间估算

| 文件类型 | 每样品大小 |
|---------|-----------|
| 原始 FASTQ | 5–20 GB |
| SAM 文件 | 10–50 GB |
| STAR 比对文件 | 20–100 GB |
| **总计（含中间文件）** | **100–300 GB** |

> 💡 建议使用 Snakemake 的 `--delete-temp-output` 标志，在成功完成后删除中间文件。

### 质控最佳实践

1. **始终先检查原始 FastQC** — 数据质量不足时可节省计算时间
2. **比对原始与最终 QC** — 验证预处理效果
3. **监控嵌合读取率** — 5–15% 是 RIC-seq 的理想范围
4. **检查 rRNA 过滤效果** — rRNA 过高（>50%）可能表明样品制备存在问题
5. **查看 MultiQC 报告** — 所有样品的综合质量概览

---

## 引用

如果您在研究中使用了本流程，请引用以下文献：

**RIC-seq 实验方案：**

Cao C, Cai Z, Ye R, et al. Global in situ profiling of RNA-RNA spatial interactions with RIC-seq. *Nature Protocols*. 2021;16:2916-2946. doi:10.1038/s41596-021-00524-2

**RIC-seq 原始方法：**

Lu Z, Zhang QC, Lee B, et al. RNA Duplex Map in Living Cells Reveals Higher-Order Transcriptome Structure. *Cell*. 2016;165(5):1267-1279. doi:10.1016/j.cell.2016.04.028

**更新协议：**

Cai Z, Cao C, Ji L, et al. RIC-seq for global in situ profiling of RNA-RNA spatial interactions. *Nature*. 2020;582(7812):432-437. doi:10.1038/s41586-020-2249-1

**原始仓库：**

https://github.com/caochch/RICpipe

---

## 许可证

本流程依据 MIT 许可证分发。详见 LICENSE 文件。

---

**流程版本**：3.0  
**最后更新**：2026 年 1 月  
**协议符合性**：Nature Protocols 2021 第 175-186 步 ✅
