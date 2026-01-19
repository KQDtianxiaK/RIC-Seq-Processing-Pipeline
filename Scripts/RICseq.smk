import os
from pathlib import Path

configfile: "Scripts/config.yaml"

# Directory paths
OUT_DIR = Path(config["outputdir"])
QC_DIR = OUT_DIR / "qc"
QC_RAW_DIR = QC_DIR / "raw"
QC_FINAL_DIR = QC_DIR / "final"
TRIM_DIR = OUT_DIR / "trim"
DEDUP_DIR = OUT_DIR / "dedup"
RRNA_DIR = OUT_DIR / "rRNA_filter"
ALIGN_DIR = OUT_DIR / "alignment"
STEP1_DIR = OUT_DIR / "step1_pair_tags"
STEP2_DIR = OUT_DIR / "step2_separate"
STEP3_DIR = OUT_DIR / "step3_category"
STEP4_DIR = OUT_DIR / "step4_intra_cluster"
STEP5_DIR = OUT_DIR / "step5_inter_network"
LOG_DIR = OUT_DIR / "logs"
MULTIQC_DIR = OUT_DIR / "multiqc"

SCRIPTS = config["scripts_root"]
SAMPLES = list(config["samples"].keys())

# Step5 parameters
THREADS = list(range(1, config["monte_carlo_threads"] + 1))

rule all:
    input:
        expand(str(MULTIQC_DIR / "multiqc_report.html")),
        expand(str(STEP5_DIR / "{sample}.significant.interMolecular.interaction.list"), sample=SAMPLES),
        expand(str(STEP4_DIR / "{sample}.cluster.withScore.highQuality.list"), sample=SAMPLES)

# ==========================================
# Step 175-176: FastQC on Raw Reads
# ==========================================

rule fastqc_raw_reads:
    """
    Step 175-176: Quality control of raw sequencing data.
    """
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["R2"]
    output:
        html1 = QC_RAW_DIR / "{sample}_R1_fastqc.html",
        html2 = QC_RAW_DIR / "{sample}_R2_fastqc.html",
        zip1 = QC_RAW_DIR / "{sample}_R1_fastqc.zip",
        zip2 = QC_RAW_DIR / "{sample}_R2_fastqc.zip"
    container: config["sif"]
    threads: 2
    log: LOG_DIR / "{sample}_fastqc_raw.log"
    shell:
        """
        mkdir -p {QC_RAW_DIR}
        
        # Run FastQC
        fastqc -t {threads} -o {QC_RAW_DIR} {input.r1} {input.r2} > {log} 2>&1
        
        # Get base names from input files (remove .fastq.gz extension)
        BASE1=$(basename {input.r1} .fastq.gz)
        BASE2=$(basename {input.r2} .fastq.gz)
        
        # Rename output files to match expected names
        echo "Renaming FastQC outputs..." >> {log}
        echo "BASE1=$BASE1, BASE2=$BASE2" >> {log}
        
        mv {QC_RAW_DIR}/${{BASE1}}_fastqc.html {output.html1}
        mv {QC_RAW_DIR}/${{BASE1}}_fastqc.zip {output.zip1}
        mv {QC_RAW_DIR}/${{BASE2}}_fastqc.html {output.html2}
        mv {QC_RAW_DIR}/${{BASE2}}_fastqc.zip {output.zip2}
        
        echo "FastQC completed successfully for {wildcards.sample}" >> {log}
        """

# ==========================================
# Step 177: Adapter Trimming with Trim Galore
# ==========================================

rule trim_galore:
    """
    Step 177: Trim adapters and low-quality bases using Trim Galore.
    """
    input:
        r1 = lambda wildcards: config["samples"][wildcards.sample]["R1"],
        r2 = lambda wildcards: config["samples"][wildcards.sample]["R2"],
        
        qc_check = QC_RAW_DIR / "{sample}_R1_fastqc.zip"
    output:
        r1_trim = temp(TRIM_DIR / "{sample}_R1_trimmed.fq.gz"),
        r2_trim = temp(TRIM_DIR / "{sample}_R2_trimmed.fq.gz"),
        report1 = TRIM_DIR / "{sample}_R1_trimming_report.txt",
        report2 = TRIM_DIR / "{sample}_R2_trimming_report.txt"
    container: config["sif"]
    threads: 8
    log: LOG_DIR / "{sample}_trim_galore.log"
    shell:
        """
        mkdir -p {TRIM_DIR}
        
        trim_galore --paired --phred33 --stringency 3 --length 30 \
            --cores {threads} --gzip \
            -o {TRIM_DIR} \
            {input.r1} {input.r2} > {log} 2>&1
        
        # Get base names
        BASE1=$(basename {input.r1} .fastq.gz)
        BASE2=$(basename {input.r2} .fastq.gz)
        
        # Rename trimmed files
        mv {TRIM_DIR}/${{BASE1}}_val_1.fq.gz {output.r1_trim}
        mv {TRIM_DIR}/${{BASE2}}_val_2.fq.gz {output.r2_trim}
        
        # Save trimming reports
        mv {TRIM_DIR}/${{BASE1}}.fastq.gz_trimming_report.txt {output.report1} 2>/dev/null || true
        mv {TRIM_DIR}/${{BASE2}}.fastq.gz_trimming_report.txt {output.report2} 2>/dev/null || true
        """

# ==========================================
# Step 178: Remove PCR Duplicates
# ==========================================

rule remove_pcr_duplicates:
    """
    Step 178: Remove PCR duplicates.
    """
    input:
        r1_trim = TRIM_DIR / "{sample}_R1_trimmed.fq.gz",
        r2_trim = TRIM_DIR / "{sample}_R2_trimmed.fq.gz"
    output:
        r1_dedup = DEDUP_DIR / "{sample}_read1.clean.rmDup.fq",
        r2_dedup = DEDUP_DIR / "{sample}_read2.clean.rmDup.fq",
        stats = DEDUP_DIR / "{sample}_dedup_stats.txt"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step0.remove_PCR_duplicates/remove_PCR_duplicates.pl",
        out_dir = DEDUP_DIR
    log: LOG_DIR / "{sample}_dedup.log"
    shell:
        """
        mkdir -p {params.out_dir}
        
        gunzip -c {input.r1_trim} > {params.out_dir}/{wildcards.sample}_R1.tmp.fq
        gunzip -c {input.r2_trim} > {params.out_dir}/{wildcards.sample}_R2.tmp.fq
        touch {params.out_dir}/{wildcards.sample}_dummy.fq

        perl {params.script} \
            {params.out_dir}/{wildcards.sample}_R1.tmp.fq \
            {params.out_dir}/{wildcards.sample}_R2.tmp.fq \
            {params.out_dir}/{wildcards.sample}_dummy.fq \
            {params.out_dir}/{wildcards.sample}_dummy.fq \
            {params.out_dir} > {params.out_dir}/{wildcards.sample}_run_dedup.sh

        sh {params.out_dir}/{wildcards.sample}_run_dedup.sh > {log} 2>&1

        mv {params.out_dir}/read1.clean.rmDup.fq {output.r1_dedup}
        mv {params.out_dir}/read2.clean.rmDup.fq {output.r2_dedup}
        
        # Generate dedup statistics
        echo "Sample: {wildcards.sample}" > {output.stats}
        echo "Deduplication completed" >> {output.stats}
        wc -l {output.r1_dedup} | awk '{{print "Reads after dedup: " $1/4}}' >> {output.stats}
        
        rm {params.out_dir}/{wildcards.sample}_R1.tmp.fq {params.out_dir}/{wildcards.sample}_R2.tmp.fq {params.out_dir}/{wildcards.sample}_dummy.fq
        """

# ==========================================
# Step 179: Trim Low-Complexity Sequences with cutadapt
# ==========================================

rule cutadapt_low_complexity:
    """
    Step 179: Trim low-complexity fragments from each end of reads by cutadapt.
    """
    input:
        r1 = DEDUP_DIR / "{sample}_read1.clean.rmDup.fq",
        r2 = DEDUP_DIR / "{sample}_read2.clean.rmDup.fq"
    output:
        r1_clean = TRIM_DIR / "{sample}_read1.clean.rmDup.rmPoly.fq",
        r2_clean = TRIM_DIR / "{sample}_read2.clean.rmDup.rmPoly.fq"
    container: config["sif"]
    threads: 4
    log: LOG_DIR / "{sample}_cutadapt.log"
    shell:
        """
        cutadapt -b A{{100}} -b C{{100}} -b G{{100}} -b T{{100}} \
            -B A{{100}} -B C{{100}} -B G{{100}} -B T{{100}} \
            -n 3 -e 0.1 --minimum-length 30 \
            -o {output.r1_clean} -p {output.r2_clean} \
            {input.r1} {input.r2} > {log} 2>&1
        """

# ==========================================
# Final FastQC (Optional but Recommended)
# ==========================================

rule fastqc_final:
    """
    Optional: Run FastQC on final processed reads before mapping.
    """
    input:
        r1 = TRIM_DIR / "{sample}_read1.clean.rmDup.rmPoly.fq",
        r2 = TRIM_DIR / "{sample}_read2.clean.rmDup.rmPoly.fq"
    output:
        html1 = QC_FINAL_DIR / "{sample}_read1_fastqc.html",
        html2 = QC_FINAL_DIR / "{sample}_read2_fastqc.html",
        zip1 = QC_FINAL_DIR / "{sample}_read1_fastqc.zip",
        zip2 = QC_FINAL_DIR / "{sample}_read2_fastqc.zip"
    container: config["sif"]
    threads: 2
    log: LOG_DIR / "{sample}_fastqc_final.log"
    shell:
        """
        mkdir -p {QC_FINAL_DIR}
        
        # Run FastQC
        fastqc -t {threads} -o {QC_FINAL_DIR} {input.r1} {input.r2} > {log} 2>&1
        
        # Get base names (remove .fq extension)
        BASE1=$(basename {input.r1} .fq)
        BASE2=$(basename {input.r2} .fq)
        
        echo "Renaming final FastQC outputs..." >> {log}
        echo "BASE1=$BASE1, BASE2=$BASE2" >> {log}
        
        # Rename to standardized names
        mv {QC_FINAL_DIR}/${{BASE1}}_fastqc.html {output.html1}
        mv {QC_FINAL_DIR}/${{BASE1}}_fastqc.zip {output.zip1}
        mv {QC_FINAL_DIR}/${{BASE2}}_fastqc.html {output.html2}
        mv {QC_FINAL_DIR}/${{BASE2}}_fastqc.zip {output.zip2}
        
        echo "Final FastQC completed successfully for {wildcards.sample}" >> {log}
        """

# ==========================================
# Step 180: rRNA Filtration
# ==========================================

rule rrna_filter:
    """
    Step 180: Map reads to rRNA index and keep unmapped reads.
    """
    input:
        r1 = TRIM_DIR / "{sample}_read1.clean.rmDup.rmPoly.fq",
        r2 = TRIM_DIR / "{sample}_read2.clean.rmDup.rmPoly.fq"
    output:
        r1_clean = RRNA_DIR / "{sample}.no_rRNA.1.fq",
        r2_clean = RRNA_DIR / "{sample}.no_rRNA.2.fq",
        log_final = RRNA_DIR / "{sample}.Log.final.out"
    container: config["sif"]
    threads: config["threads"]
    params:
        index = config["rRNA_index"],
        prefix = str(RRNA_DIR / "{sample}.")
    log: LOG_DIR / "{sample}_rrna.log"
    shell:
        """
        mkdir -p {RRNA_DIR}
        STAR --runThreadN {threads} \
             --genomeDir {params.index} \
             --readFilesIn {input.r1} {input.r2} \
             --outFileNamePrefix {params.prefix} \
             --outReadsUnmapped Fastx \
             --outSAMmode None \
             > {log} 2>&1
        
        mv {params.prefix}Unmapped.out.mate1 {output.r1_clean}
        mv {params.prefix}Unmapped.out.mate2 {output.r2_clean}
        """

# ==========================================
# Step 181: Genome Alignment
# ==========================================

rule star_align_r1:
    """Step 181: Align read1 to reference genome"""
    input: RRNA_DIR / "{sample}.no_rRNA.1.fq"
    output:
        sam = ALIGN_DIR / "{sample}.read1.Aligned.out.sam",
        chim = ALIGN_DIR / "{sample}.read1.Chimeric.out.sam",
        log_final = ALIGN_DIR / "{sample}.read1.Log.final.out"
    container: config["sif"]
    threads: config["threads"]
    params:
        index = config["star_index"],
        prefix = str(ALIGN_DIR / "{sample}.read1.")
    log: LOG_DIR / "{sample}_star_r1.log"
    shell:
        """
        mkdir -p {ALIGN_DIR}
        STAR --runThreadN {threads} --genomeDir {params.index} \
             --readFilesIn {input} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype SAM \
             --outSAMunmapped Within \
             --chimSegmentMin 15 --chimJunctionOverhangMin 15 \
             --chimOutType SeparateSAMold \
             --outFilterMultimapNmax 100 \
             --alignIntronMin 1 \
             --scoreGapNoncan -4 --scoreGapATAC -4 \
             --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 \
             --alignSJstitchMismatchNmax 5 -1 5 5 \
             --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
             > {log} 2>&1
        """

rule star_align_r2:
    """Step 181: Align read2 to reference genome"""
    input: RRNA_DIR / "{sample}.no_rRNA.2.fq"
    output:
        sam = ALIGN_DIR / "{sample}.read2.Aligned.out.sam",
        chim = ALIGN_DIR / "{sample}.read2.Chimeric.out.sam",
        log_final = ALIGN_DIR / "{sample}.read2.Log.final.out"
    container: config["sif"]
    threads: config["threads"]
    params:
        index = config["star_index"],
        prefix = str(ALIGN_DIR / "{sample}.read2.")
    log: LOG_DIR / "{sample}_star_r2.log"
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {params.index} \
             --readFilesIn {input} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype SAM \
             --outSAMunmapped Within \
             --chimSegmentMin 15 --chimJunctionOverhangMin 15 \
             --chimOutType SeparateSAMold \
             --outFilterMultimapNmax 100 \
             --alignIntronMin 1 \
             --scoreGapNoncan -4 --scoreGapATAC -4 \
             --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 \
             --alignSJstitchMismatchNmax 5 -1 5 5 \
             --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
             > {log} 2>&1
        """

# ==========================================
# Step 182: Collect Pair Tags
# ==========================================

rule step1_collect_pairs:
    """Step 182: Collect pairwise tags from the same paired-end read"""
    input:
        sam1 = ALIGN_DIR / "{sample}.read1.Aligned.out.sam",
        sam2 = ALIGN_DIR / "{sample}.read2.Aligned.out.sam",
        chim1 = ALIGN_DIR / "{sample}.read1.Chimeric.out.sam",
        chim2 = ALIGN_DIR / "{sample}.read2.Chimeric.out.sam"
    output:
        merged_sam = STEP1_DIR / "{sample}.interaction.sam",
        stats = STEP1_DIR / "{sample}_num_of_interactions.list"
    container: config["sif"]
    threads: 10
    params:
        script = f"{SCRIPTS}/step1.collect_pair_tags/collect_pair_tags.pl",
        file_prefix = "{sample}",
        script_dir = f"{SCRIPTS}/step1.collect_pair_tags/scripts"
    log: LOG_DIR / "{sample}_step1.log"
    shell:
        """
        mkdir -p {STEP1_DIR}
        
        if [ -f {params.script_dir}/precess_Chimeric_sam.pl ] && [ ! -f {params.script_dir}/process_Chimeric_sam.pl ]; then
            ln -s {params.script_dir}/precess_Chimeric_sam.pl {params.script_dir}/process_Chimeric_sam.pl
        fi

        cd {STEP1_DIR}
        
        ln -sf {input.sam1} .
        ln -sf {input.sam2} .
        ln -sf {input.chim1} .
        ln -sf {input.chim2} .

        perl {params.script} \
            $(basename {input.sam1}) $(basename {input.sam2}) \
            $(basename {input.chim1}) $(basename {input.chim2}) \
            {threads} \
            {params.file_prefix} > {wildcards.sample}_run_step1.sh

        perl -i -pe 's/samtools sort -n -@ (\\d+) (\\S+) (\\S+)/samtools sort -n -@ $1 -o $3.bam $2/' {wildcards.sample}_run_step1.sh
        perl -i -pe 's/samtools sort -@ (\\d+) (\\S+) (\\S+)/samtools sort -@ $1 -o $3.bam $2/' {wildcards.sample}_run_step1.sh

        sh {wildcards.sample}_run_step1.sh >> {log} 2>&1
        
        if [ -f num_of_interactions.list ]; then
            mv num_of_interactions.list {output.stats}
        fi
        """

# ==========================================
# Step 183: Separate Intra/Inter
# ==========================================

rule sam_to_bed:
    """Step 183: Convert SAM to BED format"""
    input: STEP1_DIR / "{sample}.interaction.sam"
    output:
        bed1 = temp(STEP2_DIR / "{sample}.read_1.bed"),
        bed2 = temp(STEP2_DIR / "{sample}.read_2.bed")
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/from_sam_to_pair_reads_bed.pl"
    log: LOG_DIR / "{sample}_sam_to_bed.log"
    shell:
        """
        mkdir -p {STEP2_DIR}
        cd {STEP2_DIR}
        ln -sf {input} {wildcards.sample}.interaction.sam
        perl {params.script} {wildcards.sample}.interaction.sam > {log} 2>&1
        
        mv read_1.bed {output.bed1}
        mv read_2.bed {output.bed2}
        rm {wildcards.sample}.interaction.sam
        """

rule intersect_genes:
    """Step 183: Intersect reads with gene annotations"""
    input:
        bed1 = STEP2_DIR / "{sample}.read_1.bed",
        bed2 = STEP2_DIR / "{sample}.read_2.bed",
        genes = config["gene_annotation_bed"]
    output:
        ov1 = temp(STEP2_DIR / "{sample}.gene_overlap_with_read1.bed"),
        ov2 = temp(STEP2_DIR / "{sample}.gene_overlap_with_read2.bed")
    container: config["sif"]
    log: LOG_DIR / "{sample}_intersect.log"
    shell:
        """
        bedtools intersect -wa -wb -a {input.genes} -b {input.bed1} > {output.ov1} 2> {log}
        bedtools intersect -wa -wb -a {input.genes} -b {input.bed2} > {output.ov2} 2>> {log}
        """

rule find_intra_candidates:
    """Step 183: Find intramolecular interaction candidates"""
    input:
        ov1 = STEP2_DIR / "{sample}.gene_overlap_with_read1.bed",
        ov2 = STEP2_DIR / "{sample}.gene_overlap_with_read2.bed"
    output: STEP2_DIR / "{sample}.pets_in_same_gene.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/intra_pets_list.pl"
    log: LOG_DIR / "{sample}_intra_list.log"
    shell:
        "perl {params.script} {input.ov1} {input.ov2} > {output} 2> {log}"

rule separate_files:
    """Step 183: Separate intra- and inter-molecular interactions"""
    input:
        sam = STEP1_DIR / "{sample}.interaction.sam",
        list_file = STEP2_DIR / "{sample}.pets_in_same_gene.list"
    output:
        intra = STEP2_DIR / "{sample}.intraMolecular.sam",
        inter = STEP2_DIR / "{sample}.interMolecular.sam"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/separate_intra_inter_pets.pl"
    log: LOG_DIR / "{sample}_separate.log"
    shell:
        """
        cd {STEP2_DIR}
        cp {input.sam} {wildcards.sample}.temp.sam
        
        perl {params.script} {wildcards.sample}.temp.sam {input.list_file} > {log} 2>&1
        
        mv {wildcards.sample}.temp.intraMolecular.sam {output.intra}
        mv {wildcards.sample}.temp.interMolecular.sam {output.inter}
        rm {wildcards.sample}.temp.sam
        """

# ==========================================
# Step 184: Categorize Intra-molecular
# ==========================================

rule step3_category:
    """Step 184: Categorize intramolecular pairwise tags"""
    input:
        intra_sam = STEP2_DIR / "{sample}.intraMolecular.sam",
        gene_bed = config["gene_annotation_bed"],
        junction_bed = config["junction_bed"]
    output:
        chim = STEP3_DIR / "{sample}.intraMolecular.Chimeric.sam",
        sing = STEP3_DIR / "{sample}.intraMolecular.Singleton.sam"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step3.category_intra_reads/category_intra_reads.pl",
        frag_cutoff = config["fragment_len_cutoff"],
        work_dir = STEP3_DIR
    log: LOG_DIR / "{sample}_step3.log"
    shell:
        """
        mkdir -p {params.work_dir}
        cd {params.work_dir}
        
        cp {input.intra_sam} .
        BASENAME=$(basename {input.intra_sam})
        
        echo "=== Running step3 category script ===" > {log} 2>&1
        echo "Working directory: $(pwd)" >> {log} 2>&1
        echo "Input file: $BASENAME" >> {log} 2>&1
        
        perl {params.script} \
            $BASENAME \
            {params.frag_cutoff} \
            {input.gene_bed} \
            {input.junction_bed} >> {log} 2>&1
        
        echo "=== Script completed, checking outputs ===" >> {log} 2>&1
        ls -lh *.sam >> {log} 2>&1
        
        if [ -f "${{BASENAME%.sam}}.Chimeric.sam" ]; then
            mv "${{BASENAME%.sam}}.Chimeric.sam" {output.chim}
        elif [ -f "{wildcards.sample}.intraMolecular.Chimeric.sam" ]; then
            mv "{wildcards.sample}.intraMolecular.Chimeric.sam" {output.chim}
        else
            touch {output.chim}
        fi
        
        if [ -f "${{BASENAME%.sam}}.Singleton.sam" ]; then
            mv "${{BASENAME%.sam}}.Singleton.sam" {output.sing}
        elif [ -f "{wildcards.sample}.intraMolecular.Singleton.sam" ]; then
            mv "{wildcards.sample}.intraMolecular.Singleton.sam" {output.sing}
        else
            touch {output.sing}
        fi
        
        rm -f $BASENAME
        """

# ==========================================
# Step 185: Cluster Intramolecular
# ==========================================

rule cluster_intramolecular:
    """Step 185: Cluster intramolecular chimeric reads"""
    input:
        intra_chim_sam = STEP3_DIR / "{sample}.intraMolecular.Chimeric.sam"
    output:
        final_list = STEP4_DIR / "{sample}.cluster.withScore.highQuality.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step4.cluster_intramolecular/cluster_intra_reads.pl",
        frag_cutoff = config["step4_fragment_cutoff"],
        score_cutoff = config["step4_connection_score_cutoff"],
        tmp_prefix = str(STEP4_DIR / "{sample}")
    log: LOG_DIR / "{sample}_step4.log"
    shell:
        """
        mkdir -p {STEP4_DIR}
        cd {STEP4_DIR}

        perl {params.script} \
            {params.frag_cutoff} \
            {params.score_cutoff} \
            {params.tmp_prefix} \
            {input.intra_chim_sam} > {wildcards.sample}_run_step4.sh 2>> {log}

        sh {wildcards.sample}_run_step4.sh >> {log} 2>&1
        """

# ==========================================
# Step 186: Intermolecular Network (Monte Carlo)
# ==========================================

rule step5_prepare:
    """Step 186: Prepare intermolecular interaction network"""
    input:
        inter_sam = STEP2_DIR / "{sample}.interMolecular.sam",
        gene_bed = config["gene_annotation_bed"]
    output:
        network = STEP5_DIR / "{sample}.merged.network",
        dir_marker = touch(STEP5_DIR / "{sample}.dirs_prepared")
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/0.prepare/from_sam_to_pair_reads_bed.pl",
        count_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/0.prepare/count_RRI_multiple_details_addMultiple.pl",
        out_prefix = "{sample}"
    log: LOG_DIR / "{sample}_step5_prepare.log"
    shell:
        """
        mkdir -p {STEP5_DIR}
        mkdir -p {STEP5_DIR}/{wildcards.sample}_sim/1.base_on_observed
        mkdir -p {STEP5_DIR}/{wildcards.sample}_sim/2.base_on_random
        
        WORKDIR={STEP5_DIR}/{wildcards.sample}_work
        mkdir -p $WORKDIR
        cd $WORKDIR
        
        perl {params.script} {input.inter_sam} 2>> {log}
        
        bedtools intersect -wa -wb -a {input.gene_bed} -b read_1.bed > gene_overlap_with_read1.bed 2>> {log}
        bedtools intersect -wa -wb -a {input.gene_bed} -b read_2.bed > gene_overlap_with_read2.bed 2>> {log}
        
        perl {params.count_script} gene_overlap_with_read1.bed gene_overlap_with_read2.bed > {output.network} 2>> {log}
        
        rm -rf $WORKDIR
        """

rule step5_sim_observed:
    """Step 186: Monte Carlo simulation based on observed data"""
    input:
        network = STEP5_DIR / "{sample}.merged.network",
        marker = STEP5_DIR / "{sample}.dirs_prepared"
    output:
        res = STEP5_DIR / "{sample}_sim/1.base_on_observed/result1.thread{tid}.simulation.list",
        pairs = STEP5_DIR / "{sample}_sim/1.base_on_observed/all_pairs.thread{tid}.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/1.MonteCarlo_simulation.pl",
        iters_per_thread = int(config["monte_carlo_iters"] / config["monte_carlo_threads"])
    wildcard_constraints:
        tid = "\\d+"
    log: LOG_DIR / "{sample}_sim_obs_{tid}.log"
    shell:
        """
        cd {STEP5_DIR}/{wildcards.sample}_sim/1.base_on_observed
        
        perl {params.script} \
            ../../{wildcards.sample}.merged.network \
            {params.iters_per_thread} \
            thread{wildcards.tid} \
            > {output.res} 2>> {log}
        """

rule step5_sim_random:
    """Step 186: Monte Carlo simulation based on random data"""
    input:
        network = STEP5_DIR / "{sample}.merged.network",
        marker = STEP5_DIR / "{sample}.dirs_prepared"
    output:
        res = STEP5_DIR / "{sample}_sim/2.base_on_random/result1.thread{tid}.simulation_on_random.list",
        pairs = STEP5_DIR / "{sample}_sim/2.base_on_random/all_pairs.thread{tid}.list"
    container: config["sif"]
    params:
        sim_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/1.MonteCarlo_simulation.pl",
        rand_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/0.creat_random_interaction.pl",
        iters_per_thread = int(config["monte_carlo_iters"] / config["monte_carlo_threads"])
    wildcard_constraints:
        tid = "\\d+"
    log: LOG_DIR / "{sample}_sim_rnd_{tid}.log"
    shell:
        """
        cd {STEP5_DIR}/{wildcards.sample}_sim/2.base_on_random
        
        perl {params.rand_script} ../../{wildcards.sample}.merged.network > random_net_{wildcards.tid}.network 2>> {log}
        
        perl {params.sim_script} \
            random_net_{wildcards.tid}.network \
            {params.iters_per_thread} \
            thread{wildcards.tid} \
            > {output.res} 2>> {log}
        """

rule step5_post_process:
    """Step 186: Post-process Monte Carlo results and calculate p-values"""
    input:
        obs_res = expand(STEP5_DIR / "{{sample}}_sim/1.base_on_observed/result1.thread{tid}.simulation.list", tid=THREADS),
        rnd_res = expand(STEP5_DIR / "{{sample}}_sim/2.base_on_random/result1.thread{tid}.simulation_on_random.list", tid=THREADS)
    output:
        final_list = STEP5_DIR / "{sample}.significant.interMolecular.interaction.list"
    container: config["sif"]
    params:
        script_split = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/3.pre_process/split_and_transpose_large_matrix.pl",
        script_pval = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/4.calculate_pvalue/pvalue_calculator.pl",
        script_merge = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/5.recalibrate_pvalue/merge_pvalue.pl",
        script_fmt = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/5.recalibrate_pvalue/1.convert_format_needed_by_CloseCall.pl",
        script_corr = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/5.recalibrate_pvalue/2.multiple_testing_correction.pl",
        script_repl = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/5.recalibrate_pvalue/4.replace_rawPvalue_by_LocalCorrectPvalue.pl",
        script_filt = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/5.recalibrate_pvalue/5.filter_interaction.pl",
        pval_cutoff = config["pvalue_cutoff"],
        threads_list = " ".join([str(t) for t in THREADS]),
        step5_dir = STEP5_DIR
    threads: 16
    log: LOG_DIR / "{sample}_step5_post.log"
    shell:
        """
        cd {params.step5_dir}
        
        WORKBASE={wildcards.sample}_post
        mkdir -p $WORKBASE/2.pre-process/1.base_on_observed
        mkdir -p $WORKBASE/2.pre-process/2.base_on_random
        mkdir -p $WORKBASE/3.calculate_pvalue/1.base_on_observed
        mkdir -p $WORKBASE/3.calculate_pvalue/2.base_on_random
        mkdir -p $WORKBASE/4.recalibrate_pvalue

        SIM_DIR={wildcards.sample}_sim

        # 1. Split matrices
        for i in {params.threads_list}; do
            perl {params.script_split} \
                $SIM_DIR/1.base_on_observed/result1.thread$i.simulation.list \
                100 thread$i &
        done
        wait
        mv *.transposed.matrix $WORKBASE/2.pre-process/1.base_on_observed/ 2>/dev/null || true

        for i in {params.threads_list}; do
            perl {params.script_split} \
                $SIM_DIR/2.base_on_random/result1.thread$i.simulation_on_random.list \
                100 thread$i &
        done
        wait
        mv *.transposed.matrix $WORKBASE/2.pre-process/2.base_on_random/ 2>/dev/null || true

        # 2. Calculate P-values
        for i in {params.threads_list}; do
            perl {params.script_pval} \
                $SIM_DIR/1.base_on_observed/all_pairs.thread$i.list \
                $WORKBASE/2.pre-process/1.base_on_observed/thread$i.transposed.matrix \
                thread$i &
        done
        wait
        mv comparison.observed_to_simulated.*.xls $WORKBASE/3.calculate_pvalue/1.base_on_observed/ 2>/dev/null || true

        for i in {params.threads_list}; do
            perl {params.script_pval} \
                $SIM_DIR/2.base_on_random/all_pairs.thread$i.list \
                $WORKBASE/2.pre-process/2.base_on_random/thread$i.transposed.matrix \
                thread$i &
        done
        wait
        mv comparison.observed_to_simulated.*.xls $WORKBASE/3.calculate_pvalue/2.base_on_random/ 2>/dev/null || true

        # 3. Recalibrate p-values
        cd $WORKBASE/4.recalibrate_pvalue
        
        perl {params.script_merge} \
            ../3.calculate_pvalue/1.base_on_observed/comparison.*.xls \
            > comparison.observed_to_simulated.basedOnObserved.finalMerge.xls 2>> {log}

        perl {params.script_merge} \
            ../3.calculate_pvalue/2.base_on_random/comparison.*.xls \
            > comparison.observed_to_simulated.basedOnRandom.finalMerge.xls 2>> {log}
            
        perl {params.script_fmt} comparison.observed_to_simulated.basedOnObserved.finalMerge.xls > result1.results_file.for_Test.list
        perl {params.script_fmt} comparison.observed_to_simulated.basedOnRandom.finalMerge.xls > result1.control_file.for_Test.list

        perl {params.script_corr} --control result1.control_file.for_Test.list --results result1.results_file.for_Test.list --window 500 2>> {log}
        
        gzip -d -f result1.results_file.for_Test.list.window_500.qval.txt.gz

        perl {params.script_repl} comparison.observed_to_simulated.basedOnObserved.finalMerge.xls result1.results_file.for_Test.list.window_500.qval.txt > result2.comparison.observed_to_simulated.Add_Local_Corrected_pvalue.list
        
        perl {params.script_filt} result2.comparison.observed_to_simulated.Add_Local_Corrected_pvalue.list {params.pval_cutoff} > {output.final_list}
        """

# ==========================================
# MultiQC: Aggregate all QC reports
# ==========================================

rule multiqc:
    input:
        # Raw FastQC
        expand(str(QC_RAW_DIR / "{sample}_R1_fastqc.zip"), sample=SAMPLES),
        expand(str(QC_RAW_DIR / "{sample}_R2_fastqc.zip"), sample=SAMPLES),
        # Final FastQC
        expand(str(QC_FINAL_DIR / "{sample}_read1_fastqc.zip"), sample=SAMPLES),
        expand(str(QC_FINAL_DIR / "{sample}_read2_fastqc.zip"), sample=SAMPLES),
        # Trimming reports
        expand(str(TRIM_DIR / "{sample}_R1_trimming_report.txt"), sample=SAMPLES),
        # Dedup stats
        expand(str(DEDUP_DIR / "{sample}_dedup_stats.txt"), sample=SAMPLES),
        # rRNA filter logs
        expand(str(RRNA_DIR / "{sample}.Log.final.out"), sample=SAMPLES),
        # STAR alignment logs
        expand(str(ALIGN_DIR / "{sample}.read1.Log.final.out"), sample=SAMPLES),
        expand(str(ALIGN_DIR / "{sample}.read2.Log.final.out"), sample=SAMPLES),
        # Step1 stats
        expand(str(STEP1_DIR / "{sample}_num_of_interactions.list"), sample=SAMPLES),
        # Final outputs
        expand(str(STEP5_DIR / "{sample}.significant.interMolecular.interaction.list"), sample=SAMPLES)
    output:
        report = MULTIQC_DIR / "multiqc_report.html",
        data = directory(MULTIQC_DIR / "multiqc_data")
    container: config["sif"]
    log: LOG_DIR / "multiqc.log"
    shell:
        """
        multiqc {OUT_DIR} -o {MULTIQC_DIR} -f --verbose 2> {log}
        """
