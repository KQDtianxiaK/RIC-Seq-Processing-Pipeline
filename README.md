# RIC-seq Processing Pipeline

This pipeline processes raw RIC-seq (RNA In situ Conformation sequencing) FASTQ data through a comprehensive workflow including quality control, adapter trimming, PCR duplicate removal, genome alignment, paired-end tag collection, and downstream analysis to identify both intra-molecular and inter-molecular RNA-RNA interactions. The entire workflow is fully containerized using Singularity and orchestrated by Snakemake, ensuring reproducibility and ease of deployment.

## Part I Overview

## What is RIC-seq?

RIC-seq (RNA In situ Conformation sequencing) is a proximity ligation-based technique designed to capture RNA-RNA interactions in living cells. The method preserves spatial information about RNA molecules that are in close proximity, enabling the identification of both intra-molecular RNA structures (within the same RNA molecule) and inter-molecular RNA interactions (between different RNA molecules).

## Workflow Summary

The RIC-seq pipeline consists of five major analytical steps:

1. **Step 0: Preprocessing** - Quality control, adapter trimming, and PCR duplicate removal
2. **Step 1: Alignment & Pairing** - STAR alignment and paired-end tag collection
3. **Step 2: Separation** - Distinguish intra-molecular from inter-molecular interactions
4. **Step 3: Categorization** - Classify intra-molecular reads into normal transcripts vs. chimeric structures
5. **Step 4: Clustering** - Cluster intra-molecular chimeric reads into high-confidence interactions
6. **Step 5: Network Analysis** - Screen significant inter-molecular RNA-RNA interactions using Monte Carlo simulation

The workflow is fully containerized using Singularity with all required bioinformatics tools and Perl modules pre-installed, and the entire pipeline can be executed with a single Snakemake command.

## Part II Requirements

## 1. Recommended System Configuration

* **CPU** : 16-core processor (minimum 8 cores)
* **RAM** : 64 GB (minimum 32 GB)
* **Storage** : At least 200 GB free space for intermediate files

## 2. Singularity

Singularity must be installed on your system. Below are detailed installation steps for Ubuntu 22.04. For other operating systems, refer to the [official Singularity installation guide](https://docs.sylabs.io/guides/3.0/user-guide/installation.html).

**Step 1: Install System Dependencies**

```
# Update package lists and install dependencies
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

**Step 2: Install Go Language**

```
# Download and install Go
wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
rm go1.21.3.linux-amd64.tar.gz

# Configure Go environment variables
echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
source ~/.bashrc
```

**Step 3: Download, Build, and Install Singularity**

```
# Navigate to your preferred directory for source code
cd /mnt/share/software

# Download Singularity CE source
wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

# Extract and build
tar -xvzf singularity-ce-4.0.1.tar.gz
cd singularity-ce-4.0.1
./mconfig
cd builddir
make
sudo make install
```

**Step 4: Verify Installation**

```
# Check installed version
singularity --version

# Display help
singularity -h
```

## 3. Snakemake

Snakemake must be installed and requires Python 3.

```
pip install snakemake
```

## 4. Reference Data

The pipeline requires several reference files for the human genome (hg38 is used as an example):

**4.1 STAR Genome Index**

```
mkdir -p References/hg38
cd References

# Download genome FASTA and GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

# Unzip files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz

# Build STAR index
singularity exec --cleanenv RIC-seq.sif STAR \
    --runMode genomeGenerate \
    --genomeDir ./hg38/STAR_index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v38.primary_assembly.annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 16
```

**4.2 Gene Annotation Files**

The pipeline requires two specialized BED files:

* **Junction BED** : Contains pairwise splicing site information
* **Gene Annotation BED** : Contains whole gene region information with Gene IDs

Generate junction BED file:

```
# Convert GTF to BED
perl scripts/gtf_to_bed.pl gencode.v38.primary_assembly.annotation.gtf > gencode.v38.annotation.bed

# Create junction BED
perl scripts/creat_junction_bed.pl gencode.v38.annotation.bed > gencode.v38.all_exon_junction.bed
```

Create whole gene region BED (must include Gene ID and detailed information):

```
# Example format: chr, start, end, gene_name, score, strand
awk 'BEGIN{OFS="\t"} $3=="gene" {print $1, $4-1, $5, $10, "60", $7}' \
    gencode.v38.primary_assembly.annotation.gtf | \
    sed 's/"//g' | sed 's/;//g' > whole_gene_region.bed
```

**4.3 Adapter Sequences**

Prepare an adapter FASTA file containing Illumina adapter sequences for trimming.

```
# Example adapter.fa
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

## 5. Required File Structure

```
project_directory/
├── Scripts/
│   ├── config.yaml
│   └── RICseq.smk
├── Containers/
│   └── RIC-seq.sif
├── References/
    ├── hg38/
    │   └── STAR_index/
    ├── gencode.v38.all_exon_junction.bed
    ├── whole_gene_region.bed
    └── adapter.fa
```

**File Descriptions:**

* **RICseq.smk** — Main Snakemake workflow script
* **config.yaml** — Configuration file with paths, parameters, and sample information (⚠️ must be in the same directory as RICseq.smk)
* **RIC-seq.sif** — Singularity container with all required tools and Perl modules
* **STAR_index/** — STAR genome index directory
* **gencode.v38.all_exon_junction.bed** — Pairwise splicing site annotation (generated from GTF)
* **whole_gene_region.bed** — Whole gene regions with Gene IDs
* **adapter.fa** — Illumina adapter sequences
* **RICseq_scripts_root/** — Directory containing all Perl analysis scripts

## Part III Running the Pipeline

## Example Execution

**Step 1: Edit `config.yaml`**

```
# Sample information
prefix: "RIC-seq_HeLa_Total_rep1"

# Input files
fastq:
  R1: "/path/to/rawdata/RIC-seq_HeLa_Total_rep1_1.fastq.gz"
  R2: "/path/to/rawdata/RIC-seq_HeLa_Total_rep1_2.fastq.gz"

# Output directory
outputdir: "/path/to/output"

# Singularity container
sif: "/path/to/Containers/RIC-seq-v1.2.sif"

# Resources
threads: 16

# Reference files
star_index: "/path/to/References/hg38/STAR_index"
junction_bed: "/path/to/References/hg38/gencode.v38.all_exon_junction.bed"
gene_annotation_bed: "/path/to/References/hg38/whole_gene_region.bed"
adapterFa: "/path/to/References/adapter.fa"
trimmomatic_jar: "/mnt/software/Trimmomatic-0.36/trimmomatic-0.36.jar"

# Script paths (container internal path)
scripts_root: "/mnt/RICseq_scripts_root"

# Trimmomatic parameters
trimmomatic_params: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

# Analysis parameters
fragment_len_cutoff: 1000      # For Step 3/4
pvalue_cutoff: 0.05             # For Step 5
monte_carlo_iters: 1000         # For Step 5
monte_carlo_threads: 16         # For Step 5
```

**Step 2: Run Snakemake**

```
# Dry run to check workflow
snakemake -s RICseq.smk --use-singularity --cores 16 -n

# Execute pipeline
snakemake -s RICseq.smk --use-singularity --cores 16 \
    --singularity-args "--bind /path/to/project_directory:/path/to/project_directory"
```

## Command Parameters

**Configuration File (`config.yaml`):**

* `prefix`: Sample name prefix for output files (required)
* `fastq`: Paths to paired-end FASTQ files (R1 and R2) (required)
* `outputdir`: Output directory path (required)
* `sif`: Path to Singularity container image (required)
* `threads`: Number of CPU threads to use (default: 16)
* `star_index`: Path to STAR genome index directory (required)
* `junction_bed`: Path to pairwise splicing junction BED file (required)
* `gene_annotation_bed`: Path to whole gene region BED file with Gene IDs (required)
* `adapterFa`: Path to adapter FASTA file (required)
* `trimmomatic_jar`: Path to Trimmomatic JAR file (required)
* `scripts_root`: Root directory containing all Perl scripts (container internal path) (required)
* `trimmomatic_params`: Trimmomatic parameter string (optional)
* `fragment_len_cutoff`: Minimum fragment length cutoff for categorization (default: 1000)
* `pvalue_cutoff`: P-value threshold for significant interactions (default: 0.05)
* `monte_carlo_iters`: Number of Monte Carlo simulation iterations (default: 1000)
* `monte_carlo_threads`: Number of threads for Monte Carlo simulations (default: 16)

**Snakemake Execution:**

* `--use-singularity`: Execute rules within Singularity container
* `--cores`: Maximum number of CPU cores for parallel execution
* `--singularity-args "--bind"`: Mount directories into the container (format: `/host/path:/container/path`, separate multiple paths with commas)
* `-n`: Dry run (preview workflow without execution)

## Part IV Output

## Output Structure

```
output/
├── rawdata.qc/
│   ├── {prefix}_R1_fastqc.html
│   ├── {prefix}_R1_fastqc.zip
│   ├── {prefix}_R2_fastqc.html
│   └── {prefix}_R2_fastqc.zip
├── trim/
│   ├── {prefix}_R1.paired.fq.gz
│   ├── {prefix}_R2.paired.fq.gz
│   ├── {prefix}_R1.clean.fq
│   └── {prefix}_R2.clean.fq
├── dedup/
│   ├── read1.clean.rmDup.fq
│   └── read2.clean.rmDup.fq
├── alignment/
│   ├── {prefix}.read1.Aligned.out.sam
│   ├── {prefix}.read1.Chimeric.out.sam
│   ├── {prefix}.read2.Aligned.out.sam
│   └── {prefix}.read2.Chimeric.out.sam
├── step1_pair_tags/
│   └── {prefix}.interaction.sam
├── step2_separate/
│   ├── {prefix}.interMolecular.sam
│   ├── {prefix}.intraMolecular.sam
│   └── {prefix}.interMolecular.withGene.sam
├── step3_category/
│   ├── {prefix}.intraMolecular.Singleton.sam
│   └── {prefix}.intraMolecular.Chimeric.sam
├── step4_intra_cluster/
│   └── {prefix}.cluster.withScore.highQuality.list
├── step5_inter_network/
│   └── {prefix}.significant.interMolecular.interaction.list
├── multiqc/
│   └── multiqc_report.html
└── logs/
    └── (various log files)
```

## Output Interpretation

**Quality Control Files:**

* **`*_fastqc.html/zip`** : FastQC reports for raw sequencing data quality assessment. Check for per-base quality scores, adapter content, and overrepresented sequences.

**Trimmed and Deduplicated Reads:**

* **`*.clean.fq`** : Adapter-trimmed and quality-filtered reads ready for alignment
* **`*.clean.rmDup.fq`** : PCR duplicate-removed reads, ensuring unique molecular observations

**Alignment Files:**

* **`*.Aligned.out.sam`** : Primary STAR alignments for each read
* **`*.Chimeric.out.sam`** : Chimeric alignments capturing potential RNA-RNA interaction junctions

**Key Results:**

* **`{prefix}.interMolecular.sam`** : Paired-end tags representing inter-molecular RNA-RNA interactions (different RNA molecules)
* **`{prefix}.intraMolecular.sam`** : Paired-end tags representing intra-molecular interactions (same RNA molecule)
* **`{prefix}.cluster.withScore.highQuality.list`** : High-confidence clustered intra-molecular RNA structures with connection scores
* **`{prefix}.significant.interMolecular.interaction.list`** : Statistically significant inter-molecular RNA-RNA interactions after Monte Carlo simulation and multiple testing correction

**MultiQC Report:**

* **`multiqc_report.html`** : Comprehensive quality control summary integrating FastQC, alignment statistics, and processing metrics across all pipeline steps. Open in a web browser for interactive exploration.

## Expected Results

For a successful RIC-seq experiment:

1. **Mapping Rate** : >70% uniquely mapped reads
2. **PCR Duplication Rate** : <30% (higher rates suggest low library complexity)
3. **Inter-molecular Interactions** : Typically 10³–10⁵ significant RNA-RNA interactions depending on sequencing depth
4. **Intra-molecular Clusters** : Typically 10²–10⁴ high-quality structural clusters

## Part V Important Notes

## 1. Reference File Format Requirements

* **Junction BED** : Must be in pairwise splicing site format generated by `creat_junction_bed.pl` (6 columns for each junction pair)
* **Gene Annotation BED** : Must contain Gene ID information in the 4th column (not just transcript IDs)

## 2. Script Path Configuration

The `scripts_root` parameter must point to the directory **inside the container** where all Perl scripts are located. This directory should contain:

* `remove_PCR_duplicates.pl` at the root
* A `scripts/` subdirectory with all supporting Perl scripts

## 3. Computational Resources

* **Step 5 (Monte Carlo simulation)** is computationally intensive. With default parameters (1000 iterations, 16 threads), expect 2-6 hours runtime.
* Increase `monte_carlo_threads` to match available CPU cores for faster execution.

## 4. Troubleshooting

**Issue: "Math::CDF module not found"**

* Solution: Rebuild the Singularity container ensuring Math::CDF and its dependencies are installed

**Issue: Low mapping rate (<50%)**

* Check adapter contamination in FastQC reports
* Verify STAR index matches the species of your sample

**Issue: "No significant interactions found"**

* Check sequencing depth (recommend >50M paired-end reads)
* Verify `pvalue_cutoff` is not too stringent (try 0.1 for exploratory analysis)
* Ensure biological replicates are processed together for better statistical power

## Reference

For more information about the RIC-seq method and its applications:

**Original Publication:**

* Lu Z, Zhang QC, Lee B, et al. RNA Duplex Map in Living Cells Reveals Higher-Order Transcriptome Structure.  *Cell* . 2016;165(5):1267-1279. doi:10.1016/j.cell.2016.04.028

**Method Details:**

* Zhang QC, Ma C, Chu L, et al. RNA In situ Conformation sequencing (RIC-seq) enables genome-wide analysis of RNA structure.  *Methods* . 2018;143:78-86. doi:10.1016/j.ymeth.2018.04.022

**Original Repositories:**

* https://github.com/caochch/RICpipe

## Citation

If you use this pipeline in your research, please cite both the original RIC-seq method paper and this pipeline repository.

---

**Pipeline Version:** v1.2

**Last Updated:** November 2025

**Maintainer:** [Your contact information]
