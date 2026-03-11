# RIC-seq Processing Pipeline

**English** | [中文](README_zh.md)

This pipeline processes raw RIC-seq (RNA In situ Conformation sequencing) FASTQ data through a comprehensive workflow including quality control, adapter trimming, PCR duplicate removal, rRNA filtration, genome alignment, paired-end tag collection, and downstream analysis to identify both intra-molecular and inter-molecular RNA-RNA interactions. The entire workflow is fully containerized using Singularity and orchestrated by Snakemake, ensuring reproducibility and compliance with the **Nature Protocols 2021** RIC-seq standard.

## Overview

### What is RIC-seq?

RIC-seq (RNA In situ Conformation sequencing) is a proximity ligation-based technique designed to capture RNA-RNA interactions in living cells. The method preserves spatial information about RNA molecules that are in close proximity, enabling the identification of both intra-molecular RNA structures (within the same RNA molecule) and inter-molecular RNA interactions (between different RNA molecules).

### Workflow Summary

The RIC-seq pipeline follows the **Nature Protocols 2021** computational workflow (Steps 175-186) and consists of the following major analytical steps:

#### Preprocessing & Quality Control (Steps 175-179)
- **Step 175-176: Raw Data QC** - FastQC analysis on raw sequencing data to assess initial quality
- **Step 177: Adapter Trimming** - Trim Galore for adapter removal and quality filtering
- **Step 178: PCR Deduplication** - Remove PCR duplicates to avoid bias
- **Step 179: Low-Complexity Filtering** - cutadapt removes poly-N tails and homopolymer sequences
- **Final QC** - FastQC on processed reads to verify quality improvement
- **rRNA Filtration** - STAR alignment to rRNA index to remove ribosomal contamination

#### Core Analysis Pipeline (Steps 180-186)
- **Step 180-181: Genome Alignment** - STAR alignment with chimeric read detection for both R1 and R2
- **Step 182: Pair Tag Collection** - Collect and merge paired-end alignment information
- **Step 183: Intra/Inter Separation** - Distinguish intra-molecular from inter-molecular interactions
- **Step 184: Intra-molecular Categorization** - Classify reads into chimeric vs. singleton structures
- **Step 185: Intra-molecular Clustering** - Cluster chimeric reads with connection score filtering (≥0.01)
- **Step 186: Inter-molecular Network** - Monte Carlo simulation with local multiple testing correction

The workflow is fully containerized using Singularity with all required bioinformatics tools and Perl modules pre-installed, and the entire pipeline can be executed with a single Snakemake command.

## Requirements

### 1. Recommended System Configuration

- **CPU**: 16-core processor (minimum 8 cores)
- **RAM**: 64 GB (minimum 32 GB)  
- **Storage**: At least 200 GB free space for intermediate files

### 2. Singularity

Singularity must be installed on your system. Below are detailed installation steps for Ubuntu 22.04. For other operating systems, refer to the [official Singularity installation guide](https://docs.sylabs.io/guides/latest/user-guide/).

**Step 1: Install System Dependencies**

```bash
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

```bash
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

```bash
# Navigate to your preferred directory for source code
cd /path/to/software

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

```bash
# Check installed version
singularity --version

# Display help
singularity -h
```

### 3. Snakemake

Snakemake must be installed and requires Python 3.

```bash
pip install snakemake
```

### 4. Reference Data

The pipeline requires several reference files for the human genome (hg38 is used as an example):

**4.1 STAR Genome Index**

```bash
mkdir -p References/hg38
cd References

# Download genome FASTA and GTF annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz

# Unzip files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v38.primary_assembly.annotation.gtf.gz

# Build STAR index (requires Singularity container)
singularity exec --cleanenv RIC-seq.sif STAR \
    --runMode genomeGenerate \
    --genomeDir ./hg38/STAR_index \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v38.primary_assembly.annotation.gtf \
    --sjdbOverhang 100 \
    --runThreadN 16
```

**4.2 rRNA Index**

Build a separate STAR index for rRNA sequences to filter out ribosomal RNA contamination:

```bash
# Download rRNA sequences (example: 45S pre-rRNA from NCBI)
wget -O rRNA.fa "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=555853"

# Build rRNA STAR index
singularity exec --cleanenv RIC-seq.sif STAR \
    --runMode genomeGenerate \
    --genomeDir ./hg38/rRNA_STAR_index \
    --genomeFastaFiles rRNA.fa \
    --genomeSAindexNbases 8 \
    --runThreadN 16
```

**4.3 Gene Annotation Files**

The pipeline requires two specialized BED files:

- **Junction BED**: Contains pairwise splicing site information
- **Gene Annotation BED**: Contains whole gene region information with Gene IDs

Generate junction BED file:

```bash
# Convert GTF to BED
perl scripts/gtf_to_bed.pl gencode.v38.primary_assembly.annotation.gtf > gencode.v38.annotation.bed

# Create junction BED  
perl scripts/creat_junction_bed.pl gencode.v38.annotation.bed > gencode.v38.all_exon_junction.bed
```

Create whole gene region BED (must include Gene ID and detailed information):

```bash
# Example format: chr, start, end, gene_name, score, strand
awk 'BEGIN{OFS="\t"} $3=="gene" {print $1, $4-1, $5, $10, "60", $7}' \
    gencode.v38.primary_assembly.annotation.gtf | \
    sed 's/"//g' | sed 's/;//g' > whole_gene_region.bed
```

**4.4 Adapter Sequences**

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

### 5. Building the Singularity Container

The Singularity container includes all required tools and Perl modules. Build it using the provided definition file:

```bash
# Build the container (requires sudo)
sudo singularity build RIC-seq.sif RIC-seq.def
```

The container includes:
- **Bioinformatics Tools**: FastQC v0.12.1, STAR v2.7.11b, SAMtools v1.22.1, BEDtools v2.28.0
- **Quality Control**: Trim Galore v0.6.10, cutadapt v3.4, MultiQC v1.19
- **Perl Modules**: List::Util, File::Spec, Getopt::Long, Graph::Undirected, Statistics::Distributions, Math::CDF
- **Analysis Scripts**: Complete RICpipe scripts pre-installed at `/mnt/RICseq_scripts_root/`

### 6. Required File Structure

```
project_directory/
├── Scripts/
│   ├── config.yaml
│   └── RICseq.smk
├── Containers/
│   └── RIC-seq.sif
├── References/
│   ├── hg38/
│   │   ├── STAR_index/
│   │   └── rRNA_STAR_index/
│   ├── gencode.v38.all_exon_junction.bed
│   └── whole_gene_region.bed
└── rawdata/
    ├── sample1_R1.fastq.gz
    └── sample1_R2.fastq.gz
```

**File Descriptions:**

- **RICseq.smk** — Main Snakemake workflow script following Nature Protocols 2021
- **config.yaml** — Configuration file supporting multiple samples
- **RIC-seq.sif** — Singularity container with all dependencies
- **STAR_index/** — STAR genome index directory
- **rRNA_STAR_index/** — STAR rRNA index for pre-filtering
- **gencode.v38.all_exon_junction.bed** — Pairwise splicing site annotation
- **whole_gene_region.bed** — Whole gene regions with Gene IDs

## Running the Pipeline

### Configuration

**Edit `config.yaml` to specify your samples and parameters:**

```yaml
# Samples configuration (supports multiple samples)
samples:
  sample1:
    R1: "/path/to/rawdata/sample1_R1.fastq.gz"
    R2: "/path/to/rawdata/sample1_R2.fastq.gz"
  sample2:
    R1: "/path/to/rawdata/sample2_R1.fastq.gz"
    R2: "/path/to/rawdata/sample2_R2.fastq.gz"

# Output directory
outputdir: "/path/to/output"

# Singularity container
sif: "/path/to/Containers/RIC-seq.sif"

# Resources
threads: 16

# Reference files
rRNA_index: "/path/to/References/hg38/rRNA_STAR_index"
star_index: "/path/to/References/hg38/STAR_index"
junction_bed: "/path/to/References/hg38/gencode.v38.all_exon_junction.bed"
gene_annotation_bed: "/path/to/References/hg38/whole_gene_region.bed"

# Script paths (container internal path)
scripts_root: "/mnt/RICseq_scripts_root"

# Analysis parameters
# Step 3/4: Intramolecular interaction parameters
fragment_len_cutoff: 1000              # Fragment length cutoff for categorization
step4_fragment_cutoff: 2                # Minimum unique fragment number for clustering
step4_connection_score_cutoff: 0.01    # Connection score cutoff (Nature Protocols default)

# Step 5: Intermolecular interaction parameters
pvalue_cutoff: 0.05                    # P-value threshold for significant interactions
monte_carlo_iters: 100000              # Monte Carlo iterations (Nature Protocols: 100,000)
monte_carlo_threads: 16                # Parallel threads for Monte Carlo simulation
```

### Execution

**Step 1: Dry Run (Optional)**

Verify the workflow without executing:

```bash
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --cores 16 -n
```

**Step 2: Execute Pipeline**

```bash
# Run with Singularity support
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --singularity-args "-B /path/to/project_directory" \
    --cores 16 -p --rerun-incomplete \
    --latency-wait 30
```

**Step 3: Background Execution (Recommended)**

```bash
# Run in background with nohup
nohup snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml \
    --use-singularity \
    --singularity-args "-B /path/to/project_directory" \
    --cores 16 -p --rerun-incomplete \
    --latency-wait 30 > pipeline.log 2>&1 &
```

### Command Parameters

**Snakemake Execution:**

- `--configfile` — Path to configuration YAML file
- `--use-singularity` — Execute rules within Singularity container
- `--singularity-args "-B"` — Mount directories into the container (format: `/host/path:/container/path`, separate multiple paths with commas)
- `--cores` — Maximum number of CPU cores for parallel execution
- `-p` — Print shell commands for transparency
- `-n` — Dry run (preview workflow without execution)
- `--rerun-incomplete` — Rerun incomplete jobs from previous runs
- `--latency-wait` — Wait time (seconds) for output files to appear (filesystem latency)

**Unlock Directory (If Pipeline Interrupted):**

```bash
snakemake -s Scripts/RICseq.smk \
    --configfile Scripts/config.yaml --unlock
```

## Output Structure

```
output/
├── qc/                          # Quality control directories
│   ├── raw/                     # FastQC on raw reads (Step 175-176)
│   │   ├── {sample}_R1_fastqc.{html,zip}
│   │   └── {sample}_R2_fastqc.{html,zip}
│   └── final/                   # FastQC on processed reads
│       ├── {sample}_read1_fastqc.{html,zip}
│       └── {sample}_read2_fastqc.{html,zip}
├── trim/                        # Trimming and low-complexity filtering
│   ├── {sample}_R{1,2}_trimming_report.txt
│   ├── {sample}_R{1,2}_trimmed.fq.gz       # After Trim Galore
│   └── {sample}_read{1,2}.clean.rmDup.rmPoly.fq  # After cutadapt
├── dedup/                       # PCR duplicate removal (Step 178)
│   ├── {sample}_read{1,2}.clean.rmDup.fq
│   └── {sample}_dedup_stats.txt
├── rRNA_filter/                 # rRNA filtration results (Step 180)
│   ├── {sample}.no_rRNA.{1,2}.fq
│   └── {sample}.Log.final.out
├── alignment/                   # STAR alignment (Step 181)
│   ├── {sample}.read{1,2}.Aligned.out.sam
│   ├── {sample}.read{1,2}.Chimeric.out.sam
│   └── {sample}.read{1,2}.Log.final.out
├── step1_pair_tags/             # Paired-end tag collection (Step 182)
│   ├── {sample}.interaction.sam
│   └── {sample}_num_of_interactions.list
├── step2_separate/              # Intra/inter-molecular separation (Step 183)
│   ├── {sample}.intraMolecular.sam
│   ├── {sample}.interMolecular.sam
│   └── {sample}.pets_in_same_gene.list
├── step3_category/              # Intra-molecular categorization (Step 184)
│   ├── {sample}.intraMolecular.Chimeric.sam
│   └── {sample}.intraMolecular.Singleton.sam
├── step4_intra_cluster/         # Intra-molecular clustering (Step 185)
│   └── {sample}.cluster.withScore.highQuality.list
├── step5_inter_network/         # Inter-molecular network analysis (Step 186)
│   ├── {sample}.merged.network
│   ├── {sample}_sim/            # Monte Carlo simulation results
│   │   ├── 1.base_on_observed/
│   │   └── 2.base_on_random/
│   ├── {sample}_post/           # Post-processing
│   │   ├── 2.pre-process/
│   │   ├── 3.calculate_pvalue/
│   │   └── 4.recalibrate_pvalue/
│   └── {sample}.significant.interMolecular.interaction.list
├── multiqc/                     # Aggregated quality control report
│   ├── multiqc_report.html
│   └── multiqc_data/
└── logs/                        # Detailed logs for all steps
    └── {sample}_{step}.log
```

## Output Interpretation

### Key Result Files

**Intra-molecular Interactions:**

- **`{sample}.cluster.withScore.highQuality.list`** — High-confidence clustered intra-molecular RNA structures with connection scores ≥0.01. These represent RNA secondary/tertiary structures within individual transcripts.

**Inter-molecular Interactions:**

- **`{sample}.significant.interMolecular.interaction.list`** — Statistically significant inter-molecular RNA-RNA interactions (p < 0.05) after Monte Carlo simulation and local multiple testing correction. These represent RNA-RNA interactions between different transcripts.

### Quality Control Metrics

**FastQC Reports (Two Checkpoints):**

1. **Raw Data QC (`qc/raw/`)** — Quality assessment BEFORE any processing:
   - Per-base quality scores (baseline quality)
   - Adapter content detection
   - Overrepresented sequences
   - GC content distribution
   - **Purpose**: Decide if re-sequencing is needed

2. **Final QC (`qc/final/`)** — Quality assessment AFTER all preprocessing:
   - Verify quality improvement after trimming
   - Confirm adapter removal
   - Check for remaining artifacts
   - **Purpose**: Validate preprocessing effectiveness

**Trimming Reports:**

- **`*_trimming_report.txt`** — Trim Galore statistics:
  - Number of reads processed
  - Adapter sequences detected and removed
  - Quality trimming statistics

**Alignment Statistics:**

- **`*.Log.final.out`** — STAR alignment logs containing:
  - Total reads
  - Uniquely mapped reads percentage (expect >70%)
  - Multi-mapped reads
  - Chimeric reads (important for RIC-seq, expect 5-15%)

**PCR Duplication:**

- **`*_dedup_stats.txt`** — PCR duplication removal statistics
  - Expect <30% duplication rate for good library complexity

**rRNA Filtration:**

- **rRNA mapping rate** — Percentage of reads mapped to rRNA (varies by sample prep, typically 10-40%)

**MultiQC Report:**

- **`multiqc_report.html`** — Comprehensive quality control summary integrating all QC metrics across samples, including:
  - Raw vs. final quality comparison
  - Trimming statistics
  - Alignment metrics
  - Duplication rates
  - Open in a web browser for interactive exploration

### Expected Results

For a successful RIC-seq experiment:

1. **Raw Read Quality**: >Q30 for most bases (check in `qc/raw/`)
2. **Adapter Content**: <1% in raw data, 0% after trimming
3. **Mapping Rate**: >70% uniquely mapped reads to genome
4. **Chimeric Read Rate**: 5-15% (critical for RIC-seq)
5. **PCR Duplication Rate**: <30% (higher rates suggest low library complexity)
6. **rRNA Contamination**: 10-40% (depends on sample preparation)
7. **Inter-molecular Interactions**: Typically 10³–10⁵ significant RNA-RNA interactions depending on sequencing depth and p-value cutoff
8. **Intra-molecular Clusters**: Typically 10²–10⁴ high-quality structural clusters

## Pipeline Updates

### Version 3.0 - Nature Protocols Compliance

**Major Improvements:**

1. **✅ Strict Nature Protocols Adherence** — Complete implementation of Steps 175-186 from Nature Protocols 2021
2. **✅ Raw Data Quality Control (NEW)** — Added FastQC on raw reads as first step (Steps 175-176)
3. **✅ Dual QC Checkpoints** — Raw data QC + final processed data QC for comprehensive quality tracking
4. **✅ Correct Processing Order** — Fixed workflow sequence:
   ```
   FastQC (raw) → Trim Galore → PCR Duplication → cutadapt → FastQC (final) → rRNA filter → Alignment
   ```
5. **✅ Enhanced Troubleshooting** — Early quality detection prevents wasting compute on poor-quality data

**Key Workflow Changes:**

| Aspect | Previous Version | v3.0 (Current) |
|--------|-----------------|----------------|
| **Initial QC** | After trimming only | **Before any processing** (Step 175-176) |
| **Processing Order** | Trim → cutadapt → FastQC → Dedup | **FastQC → Trim → Dedup → cutadapt → FastQC** |
| **QC Checkpoints** | 1 (post-trim) | **2 (pre + post)** |
| **Protocol Compliance** | Partial | **Full (Steps 175-186)** |
| **Output Structure** | `qc/` | `qc/raw/` + `qc/final/` |

**Why This Matters:**

- **Early Detection**: Identifies sequencing quality issues before wasting computational resources
- **Decision Point**: Raw FastQC allows early termination if re-sequencing is needed
- **Reproducibility**: Matches published Nature Protocols workflow exactly
- **Validation**: Two QC checkpoints verify each processing step's effectiveness

### Previous Updates

1. **Trim Galore Integration** — Replaced Trimmomatic for streamlined adapter trimming
2. **Poly-N Tail Removal** — Added cutadapt homopolymer filtering
3. **rRNA Pre-filtration** — Dedicated STAR rRNA filtering step
4. **Multi-sample Support** — Native YAML configuration for batch processing
5. **Improved Error Handling** — Robust file handling in Steps 3 and 5
6. **Updated Parameters** — Nature 2020 alignment (100,000 Monte Carlo iterations)
7. **Complete Containerization** — All dependencies in Singularity container

### Workflow Evolution

| Step | v1 | v2 | v3 (Current) |
|------|----|----|--------------|
| **Raw QC** | ❌ Not included | ❌ Not included | ✅ **FastQC (Step 175-176)** |
| **Trimming** | Trimmomatic | Trim Galore | Trim Galore (Step 177) |
| **Deduplication** | Manual | After FastQC | **After Trimming (Step 178)** |
| **Poly-N removal** | ❌ Not included | cutadapt | cutadapt (Step 179) |
| **Final QC** | ❌ Not included | FastQC | ✅ **FastQC (validation)** |
| **rRNA filtering** | ❌ Not included | STAR | STAR (Step 180) |
| **Monte Carlo** | 1,000 iterations | 100,000 | 100,000 (Step 186) |
| **Protocol Match** | Partial | Partial | ✅ **100% (Steps 175-186)** |

## Troubleshooting

### Common Issues

**Issue: FastQC outputs not found / MissingOutputException**

Solution: FastQC generates files based on input filenames. If you encounter file naming issues:

```bash
# Increase filesystem latency wait time
snakemake --use-singularity --cores 16 --latency-wait 30

# Check log files for actual generated filenames
cat output/logs/{sample}_fastqc_raw.log
```

**Issue: Poor raw read quality (<Q20)**

Solution: Check raw FastQC reports in `qc/raw/`:
- If median quality <Q20 for >50% of bases, consider re-sequencing
- If only adapter contamination, pipeline will handle it
- If overrepresented sequences detected, investigate source

**Issue: Low mapping rate (<50%)**

Solutions:
- Check both raw and final FastQC reports to compare quality
- Verify adapter removal was successful
- Confirm STAR index matches sample species
- Check for high rRNA contamination

**Issue: High PCR duplication (>40%)**

Solutions:
- Check raw FastQC for duplicate sequences
- Indicates low library complexity
- Consider using more input material in library prep
- For exploratory analysis, proceed but note limitation

**Issue: "Math::CDF module not found"**

Solution: The Singularity container includes Math::CDF. Ensure you're using the correct container version.

```bash
# Verify Math::CDF is installed
singularity exec RIC-seq.sif perl -MMath::CDF -e 'print "Math::CDF: OK\n"'
```

**Issue: "No significant interactions found"**

Solutions:
- Check sequencing depth (recommend >50M paired-end reads per sample)
- Verify `pvalue_cutoff` is not too stringent (try 0.1 for exploratory analysis)
- Ensure biological replicates are processed for better statistical power
- Check Monte Carlo simulation completed successfully (16 thread outputs should exist)
- Review raw data quality in `qc/raw/` — poor quality impacts downstream results

**Issue: Step 3 outputs are empty**

Solutions:
- Check logs in `output/logs/{sample}_step3.log`
- Verify gene annotation BED and junction BED files are correct format
- Ensure sufficient intra-molecular interactions from Step 2
- Review final QC reports to verify input quality

**Issue: Pipeline stuck or slow**

Solutions:
- Check system resources (CPU, RAM, disk I/O)
- Step 5 Monte Carlo simulation is computationally intensive (2-6 hours with default parameters)
- Increase `monte_carlo_threads` to match available CPU cores
- Consider reducing `monte_carlo_iters` for faster testing (min 10,000)

**Issue: Singularity binding errors**

Solutions:
- Ensure all required paths are bound using `--singularity-args "-B /path1,/path2"`
- Paths must exist on host system
- Use absolute paths in configuration file

**Issue: Workflow order confusion**

Solution: The current pipeline strictly follows Nature Protocols 2021 Steps 175-186:
1. Step 175-176: FastQC on **raw reads** first
2. Step 177: Adapter trimming
3. Step 178: PCR deduplication
4. Step 179: Low-complexity filtering
5. FastQC on **processed reads** (optional validation)
6. Step 180-186: Core analysis pipeline

## Important Notes

### Reference File Format Requirements

- **Junction BED**: Must be in pairwise splicing site format generated by `creat_junction_bed.pl` (12 columns total)
- **Gene Annotation BED**: Must contain Gene ID information in the 4th column (not just transcript IDs)
- **rRNA Index**: Must be built specifically for rRNA sequences (use lower `genomeSAindexNbases` for small genomes)

### Script Path Configuration

The `scripts_root` parameter points to `/mnt/RICseq_scripts_root` inside the container where all Perl scripts are located. This directory structure is:

```
/mnt/RICseq_scripts_root/
├── step0.remove_PCR_duplicates/
├── step1.collect_pair_tags/
├── step2.separate_intra_inter_molecular/
├── step3.category_intra_reads/
├── step4.cluster_intramolecular/
└── step5.screen_high-confidence_intermolecular/
```

### Computational Resources

- **Steps 175-179 (Preprocessing)**: 1-2 hours for 50M read pairs
- **Steps 180-184 (Alignment & Separation)**: 2-3 hours
- **Step 185 (Clustering)**: 0.5-1 hour
- **Step 186 (Monte Carlo simulation)**: Computationally intensive
  - With default parameters (100,000 iterations, 16 threads): 2-6 hours
  - Increase `monte_carlo_threads` to match available CPU cores for faster execution
  - Memory usage scales with network complexity (~4-8 GB per thread)

### Data Storage

Intermediate files can be large:
- Raw FASTQ: 5-20 GB per sample
- SAM files: 10-50 GB per sample
- STAR alignment: 20-100 GB per sample  
- Total pipeline: 100-300 GB per sample (including temp files)

Consider using `--delete-temp-output` flag in Snakemake to remove intermediate files after successful completion.

### Quality Control Best Practices

1. **Always check raw FastQC first** — Saves time if data quality is insufficient
2. **Compare raw vs. final QC** — Validates preprocessing effectiveness
3. **Monitor chimeric read rate** — 5-15% is ideal for RIC-seq
4. **Check rRNA filtration** — High rRNA (>50%) may indicate sample prep issues
5. **Review MultiQC report** — Comprehensive overview of all samples

## Citation

If you use this pipeline in your research, please cite:

**RIC-seq Protocol:**

Cao C, Cai Z, Ye R, et al. Global in situ profiling of RNA-RNA spatial interactions with RIC-seq. *Nature Protocols*. 2021;16:2916-2946. doi:10.1038/s41596-021-00524-2

**Original RIC-seq Method:**

Lu Z, Zhang QC, Lee B, et al. RNA Duplex Map in Living Cells Reveals Higher-Order Transcriptome Structure. *Cell*. 2016;165(5):1267-1279. doi:10.1016/j.cell.2016.04.028

**Updated Protocol:**

Cai Z, Cao C, Ji L, et al. RIC-seq for global in situ profiling of RNA-RNA spatial interactions. *Nature*. 2020;582(7812):432-437. doi:10.1038/s41586-020-2249-1

**Original Repository:**

https://github.com/caochch/RICpipe

## License

This pipeline is distributed under the MIT License. See LICENSE file for details.

---
**Pipeline Version**: 3.0  
**Last Updated**: January 2026  
**Protocol Compliance**: Nature Protocols 2021 Steps 175-186 ✅
```