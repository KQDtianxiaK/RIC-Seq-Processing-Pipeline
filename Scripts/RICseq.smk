import os
from pathlib import Path

configfile: "config.yaml"

# Directory paths
OUT_DIR = Path(config["outputdir"])
QC_DIR = OUT_DIR / "rawdata.qc"  # Quality control output directory
TRIM_DIR = OUT_DIR / "trim"  # Trimmed reads directory
DEDUP_DIR = OUT_DIR / "dedup"  # Deduplicated reads directory
ALIGN_DIR = OUT_DIR / "alignment"  # Alignment results directory
STEP1_DIR = OUT_DIR / "step1_pair_tags"  # Step1: Collect paired tags
STEP2_DIR = OUT_DIR / "step2_separate"  # Step2: Separate intra/inter-molecular interactions
STEP3_DIR = OUT_DIR / "step3_category"  # Step3: Categorize interactions
STEP4_DIR = OUT_DIR / "step4_intra_cluster"  # Step4: Cluster intramolecular interactions
STEP5_DIR = OUT_DIR / "step5_inter_network"  # Step5: Analyze intermolecular interaction network
LOG_DIR = OUT_DIR / "logs"  # Log files directory
MULTIQC_DIR = OUT_DIR / "multiqc"  # MultiQC report directory

PREFIX = config["prefix"]
SCRIPTS = config["scripts_root"]

# Step5 parameters
THREADS = list(range(1, config["monte_carlo_threads"] + 1))

rule all:
    input:
        # Step 2 outputs
        str(STEP2_DIR / f"{PREFIX}.interMolecular.sam"),
        str(STEP2_DIR / f"{PREFIX}.intraMolecular.sam"),
        # Step 5 results
        str(STEP5_DIR / f"{PREFIX}.significant.interMolecular.interaction.list"),
        # Step 4 (optional, high memory usage)
        str(STEP4_DIR / f"{PREFIX}.cluster.withScore.highQuality.list"),
        str(MULTIQC_DIR / "multiqc_report.html"), 
        str(QC_DIR / f"{PREFIX}_R1_fastqc.zip"),

# ==========================================
# Quality Control Rules
# ==========================================

rule fastqc_raw:
    input:
        r1 = config["fastq"]["R1"],
        r2 = config["fastq"]["R2"]
    output:
        zip1 = QC_DIR / f"{PREFIX}_R1_fastqc.zip",
        html1 = QC_DIR / f"{PREFIX}_R1_fastqc.html",
        zip2 = QC_DIR / f"{PREFIX}_R2_fastqc.zip",
        html2 = QC_DIR / f"{PREFIX}_R2_fastqc.html"
    container: config["sif"]
    threads: 4
    log: LOG_DIR / "fastqc_raw.log"
    shell:
        """
        mkdir -p {QC_DIR}
        fastqc -t {threads} -o {QC_DIR} {input.r1} {input.r2} 2> {log}

        BASE_R1=$(basename {input.r1} .fastq.gz | sed 's/.fq.gz//')
        BASE_R2=$(basename {input.r2} .fastq.gz | sed 's/.fq.gz//')

        if [ -f {QC_DIR}/${{BASE_R1}}_fastqc.zip ]; then
            mv {QC_DIR}/${{BASE_R1}}_fastqc.zip {output.zip1}
            mv {QC_DIR}/${{BASE_R1}}_fastqc.html {output.html1}
        fi

        if [ -f {QC_DIR}/${{BASE_R2}}_fastqc.zip ]; then
            mv {QC_DIR}/${{BASE_R2}}_fastqc.zip {output.zip2}
            mv {QC_DIR}/${{BASE_R2}}_fastqc.html {output.html2}
        fi
        """

# ==========================================
# Read Preprocessing and Deduplication
# ==========================================

rule trimmomatic:
    input:
        r1 = config["fastq"]["R1"],
        r2 = config["fastq"]["R2"]
    output:
        r1_p = TRIM_DIR / f"{PREFIX}_R1.paired.fq.gz",
        r1_u = TRIM_DIR / f"{PREFIX}_R1.unpaired.fq.gz",
        r2_p = TRIM_DIR / f"{PREFIX}_R2.paired.fq.gz",
        r2_u = TRIM_DIR / f"{PREFIX}_R2.unpaired.fq.gz"
    container: config["sif"]
    threads: config["threads"]
    params:
        jar = config["trimmomatic_jar"],
        trim_cmd = config["trimmomatic_params"].replace("{adapter}", config["adapterFa"])
    log: LOG_DIR / "trimmomatic.log"
    shell:
        """
        mkdir -p {TRIM_DIR}
        java -jar {params.jar} PE -threads {threads} -phred33 \
            {input.r1} {input.r2} \
            {output.r1_p} {output.r1_u} \
            {output.r2_p} {output.r2_u} \
            {params.trim_cmd} 2> {log}
        """

rule cutadapt:
    input:
        r1 = TRIM_DIR / f"{PREFIX}_R1.paired.fq.gz",
        r2 = TRIM_DIR / f"{PREFIX}_R2.paired.fq.gz"
    output:
        r1 = TRIM_DIR / f"{PREFIX}_R1.clean.fq",
        r2 = TRIM_DIR / f"{PREFIX}_R2.clean.fq"
    container: config["sif"]
    threads: 8
    log: LOG_DIR / "cutadapt.log"
    shell:
        """
        cutadapt -q 25 -m 36 -j {threads} \
            -o {output.r1} -p {output.r2} \
            {input.r1} {input.r2} > {log}
        """

rule unzip_unpaired:
    input:
        r1_u = TRIM_DIR / f"{PREFIX}_R1.unpaired.fq.gz",
        r2_u = TRIM_DIR / f"{PREFIX}_R2.unpaired.fq.gz"
    output:
        r1_u = temp(TRIM_DIR / f"{PREFIX}_R1.unpaired.fq"),
        r2_u = temp(TRIM_DIR / f"{PREFIX}_R2.unpaired.fq")
    shell:
        "gunzip -c {input.r1_u} > {output.r1_u} && gunzip -c {input.r2_u} > {output.r2_u}"

rule remove_pcr_duplicates:
    input:
        r1_c = TRIM_DIR / f"{PREFIX}_R1.clean.fq",
        r2_c = TRIM_DIR / f"{PREFIX}_R2.clean.fq",
        r1_u = TRIM_DIR / f"{PREFIX}_R1.unpaired.fq",
        r2_u = TRIM_DIR / f"{PREFIX}_R2.unpaired.fq"
    output:
        r1_dedup = DEDUP_DIR / "read1.clean.rmDup.fq",
        r2_dedup = DEDUP_DIR / "read2.clean.rmDup.fq"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step0.remove_PCR_duplicates/remove_PCR_duplicates.pl",
        out_dir = DEDUP_DIR
    log: LOG_DIR / "dedup.log"
    shell:
        """
        mkdir -p {params.out_dir}
        perl {params.script} \
            {input.r1_c} {input.r2_c} \
            {input.r1_u} {input.r2_u} \
            {params.out_dir} > {params.out_dir}/run_dedup_commands.sh

        sh {params.out_dir}/run_dedup_commands.sh > {log} 2>&1
        """

# ==========================================
# Step 1: Alignment & Collecting Pairs
# ==========================================

rule star_align_r1:
    input: DEDUP_DIR / "read1.clean.rmDup.fq"
    output:
        sam = ALIGN_DIR / f"{PREFIX}.read1.Aligned.out.sam",
        chim = ALIGN_DIR / f"{PREFIX}.read1.Chimeric.out.sam"
    container: config["sif"]
    threads: config["threads"]
    params:
        index = config["star_index"],
        prefix = str(ALIGN_DIR / f"{PREFIX}.read1.")
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
             --alignIntronMax 1000000 --alignMatesGapMax 1000000
        """

rule star_align_r2:
    input: DEDUP_DIR / "read2.clean.rmDup.fq"
    output:
        sam = ALIGN_DIR / f"{PREFIX}.read2.Aligned.out.sam",
        chim = ALIGN_DIR / f"{PREFIX}.read2.Chimeric.out.sam"
    container: config["sif"]
    threads: config["threads"]
    params:
        index = config["star_index"],
        prefix = str(ALIGN_DIR / f"{PREFIX}.read2.")
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {params.index} \
             --readFilesIn {input} \
             --outFileNamePrefix {params.prefix} \
             --outSAMtype SAM \
             --outSAMunmapped Within \
             --chimSegmentMin 15 --chimJunctionOverhangMin 15 \
             --chimOutType SeparateSAMold \
             --alignIntronMax 1000000 --alignMatesGapMax 1000000
        """

rule step1_collect_pairs:
    input:
        sam1 = ALIGN_DIR / f"{PREFIX}.read1.Aligned.out.sam",
        sam2 = ALIGN_DIR / f"{PREFIX}.read2.Aligned.out.sam",
        chim1 = ALIGN_DIR / f"{PREFIX}.read1.Chimeric.out.sam",
        chim2 = ALIGN_DIR / f"{PREFIX}.read2.Chimeric.out.sam"
    output:
        merged_sam = STEP1_DIR / f"{PREFIX}.interaction.sam",
        stats = STEP1_DIR / "num_of_interactions.list"
    container: config["sif"]
    threads: 10
    params:
        script = f"{SCRIPTS}/step1.collect_pair_tags/collect_pair_tags.pl",
        file_prefix = f"{PREFIX}",
        script_dir = f"{SCRIPTS}/step1.collect_pair_tags/scripts"
    log: LOG_DIR / "step1_collect.log"
    shell:
        """
        mkdir -p {STEP1_DIR}

        # Fix typo in original script filename if needed
        if [ -f {params.script_dir}/precess_Chimeric_sam.pl ] && [ ! -f {params.script_dir}/process_Chimeric_sam.pl ]; then
            ln -s {params.script_dir}/precess_Chimeric_sam.pl {params.script_dir}/process_Chimeric_sam.pl
        fi

        cd {STEP1_DIR}

        perl {params.script} \
            {input.sam1} {input.sam2} \
            {input.chim1} {input.chim2} \
            {threads} \
            {params.file_prefix} > run_step1.sh

        # Fix samtools sort syntax
        perl -i -pe 's/samtools sort -n -@ (\d+) (\S+) (\S+)/samtools sort -n -@ $1 -o $3.bam $2/' run_step1.sh
        perl -i -pe 's/samtools sort -@ (\d+) (\S+) (\S+)/samtools sort -@ $1 -o $3.bam $2/' run_step1.sh

        echo "=== FINAL COMMANDS ===" >> {log}
        cat run_step1.sh >> {log}
        echo "======================" >> {log}

        sh run_step1.sh >> {log} 2>&1
        """

# ==========================================
# Step 2: Separate Intra/Inter
# ==========================================

rule sam_to_bed:
    """Convert SAM format to BED format for downstream analysis"""
    input: STEP1_DIR / f"{PREFIX}.interaction.sam"
    output:
        bed1 = STEP2_DIR / "read_1.bed",
        bed2 = STEP2_DIR / "read_2.bed"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/from_sam_to_pair_reads_bed.pl"
    shell:
        """
        mkdir -p {STEP2_DIR}
        cd {STEP2_DIR}
        perl {params.script} {input}
        """

rule intersect_genes:
    """Intersect reads with gene annotations to identify overlapping genes"""
    input:
        bed1 = STEP2_DIR / "read_1.bed",
        bed2 = STEP2_DIR / "read_2.bed",
        genes = config["gene_annotation_bed"]
    output:
        ov1 = STEP2_DIR / "gene_overlap_with_read1.bed",
        ov2 = STEP2_DIR / "gene_overlap_with_read2.bed"
    container: config["sif"]
    log: LOG_DIR / "step2_intersect.log"
    shell:
        """
        bedtools intersect -wa -wb -a {input.genes} -b {input.bed1} > {output.ov1} 2> {log}
        bedtools intersect -wa -wb -a {input.genes} -b {input.bed2} > {output.ov2} 2>> {log}
        """

rule find_intra_candidates:
    """Identify read pairs within the same gene (intramolecular)"""
    input:
        ov1 = STEP2_DIR / "gene_overlap_with_read1.bed",
        ov2 = STEP2_DIR / "gene_overlap_with_read2.bed"
    output: STEP2_DIR / "pets_in_same_gene.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/intra_pets_list.pl"
    shell:
        "perl {params.script} {input.ov1} {input.ov2} > {output}"

rule separate_files:
    """Separate intramolecular and intermolecular interactions into different files"""
    input:
        sam = STEP1_DIR / f"{PREFIX}.interaction.sam",
        list_file = STEP2_DIR / "pets_in_same_gene.list"
    output:
        intra = STEP2_DIR / f"{PREFIX}.intraMolecular.sam",
        inter = STEP2_DIR / f"{PREFIX}.interMolecular.sam"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step2.separate_intra_inter_molecular/scripts/separate_intra_inter_pets.pl"
    shell:
        """
        cd {STEP2_DIR}
        cp {input.sam} temp_for_split.sam
        perl {params.script} temp_for_split.sam {input.list_file}
        mv temp_for_split.intraMolecular.sam {output.intra}
        mv temp_for_split.interMolecular.sam {output.inter}
        rm -f temp_for_split.sam
        """

# ==========================================
# Step 4: Cluster Intramolecular (可选)
# ==========================================

rule cluster_intramolecular:
    """Cluster intramolecular interactions and score clusters (high memory usage)"""
    input:
        intra_sam = STEP2_DIR / f"{PREFIX}.intraMolecular.sam"
    output:
        final_list = STEP4_DIR / f"{PREFIX}.cluster.withScore.highQuality.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step4.cluster_intramolecular/cluster_intra_reads.pl",
        frag_cutoff = 2,
        score_cutoff = 10,
        tmp_prefix = str(STEP4_DIR / f"{PREFIX}")
    log: LOG_DIR / "step4_cluster.log"
    shell:
        """
        mkdir -p {STEP4_DIR}
        cd {STEP4_DIR}

        perl {params.script} \
            {params.frag_cutoff} \
            {params.score_cutoff} \
            {params.tmp_prefix} \
            {input.intra_sam} > run_step4.sh

        sh run_step4.sh > {log} 2>&1
        """

# ==========================================
# Step 5: Intermolecular Network
# ==========================================

rule step5_prepare:
    input:
        inter_sam = STEP2_DIR / f"{PREFIX}.interMolecular.sam",
        gene_bed = config["gene_annotation_bed"]
    output:
        network = STEP5_DIR / f"{PREFIX}.merged.network",
        dir_marker = touch(STEP5_DIR / ".dirs_prepared")
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/0.prepare/from_sam_to_pair_reads_bed.pl",
        count_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/0.prepare/count_RRI_multiple_details_addMultiple.pl",
        out_prefix = f"{PREFIX}"
    log: LOG_DIR / "step5_prepare.log"
    shell:
        """
        mkdir -p {STEP5_DIR}
        mkdir -p {STEP5_DIR}/1.run_simulation/1.base_on_observed
        mkdir -p {STEP5_DIR}/1.run_simulation/2.base_on_random

        cd {STEP5_DIR}

        # Generate BED files
        perl {params.script} {input.inter_sam} 2>> {log}

        # Intersect reads with gene annotations
        bedtools intersect -wa -wb -a {input.gene_bed} -b read_1.bed > gene_overlap_with_read1.bed 2>> {log}
        bedtools intersect -wa -wb -a {input.gene_bed} -b read_2.bed > gene_overlap_with_read2.bed 2>> {log}

        # Count RNA-RNA interactions
        perl {params.count_script} gene_overlap_with_read1.bed gene_overlap_with_read2.bed > {output.network} 2>> {log}

        # Clean up temporary files
        rm -f read_1.bed read_2.bed gene_overlap_with_read1.bed gene_overlap_with_read2.bed
        """

# Generate observed simulation results
rule step5_sim_observed:
    """Run Monte Carlo simulation based on observed interaction network"""
    input:
        network = STEP5_DIR / f"{PREFIX}.merged.network",
        marker = STEP5_DIR / ".dirs_prepared"
    output:
        res = STEP5_DIR / "1.run_simulation/1.base_on_observed/result1.thread{tid}.simulation.list",
        pairs = STEP5_DIR / "1.run_simulation/1.base_on_observed/all_pairs.thread{tid}.list"
    container: config["sif"]
    params:
        script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/1.MonteCarlo_simulation.pl",
        iters_per_thread = int(config["monte_carlo_iters"] / config["monte_carlo_threads"])
    wildcard_constraints:
        tid = "\d+"
    threads: 1
    log: LOG_DIR / "step5_sim_observed_{tid}.log"
    shell:
        """
        cd {STEP5_DIR}/1.run_simulation/1.base_on_observed

        perl {params.script} \
            ../../{PREFIX}.merged.network \
            {params.iters_per_thread} \
            thread{wildcards.tid} \
            > {output.res} 2>> {log}
        """

# Generate random simulation results
rule step5_sim_random:
    """Run Monte Carlo simulation based on randomized interaction network"""
    input:
        network = STEP5_DIR / f"{PREFIX}.merged.network",
        marker = STEP5_DIR / ".dirs_prepared"
    output:
        res = STEP5_DIR / "1.run_simulation/2.base_on_random/result1.thread{tid}.simulation_on_random.list",
        pairs = STEP5_DIR / "1.run_simulation/2.base_on_random/all_pairs.thread{tid}.list"
    container: config["sif"]
    params:
        sim_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/1.MonteCarlo_simulation.pl",
        rand_script = f"{SCRIPTS}/step5.screen_high-confidence_intermolecular/scripts/2.run_simulation/0.creat_random_interaction.pl",
        iters_per_thread = int(config["monte_carlo_iters"] / config["monte_carlo_threads"])
    wildcard_constraints:
        tid = "\d+"
    threads: 1
    log: LOG_DIR / "step5_sim_random_{tid}.log"
    shell:
        """
        cd {STEP5_DIR}/1.run_simulation/2.base_on_random

        # Create random network (independent for each thread)
        perl {params.rand_script} ../../{PREFIX}.merged.network > random_net_{wildcards.tid}.network 2>> {log}

        # Run simulation
        perl {params.sim_script} \
            random_net_{wildcards.tid}.network \
            {params.iters_per_thread} \
            thread{wildcards.tid} \
            > {output.res} 2>> {log}
        """

# Post-processing: Aggregate simulation results and calculate p-values
rule step5_post_process:
    """Merge multi-threaded simulation results, calculate p-values, and apply multiple testing correction"""
    input:
        obs_res = expand(STEP5_DIR / "1.run_simulation/1.base_on_observed/result1.thread{tid}.simulation.list", tid=THREADS),
        obs_pairs = expand(STEP5_DIR / "1.run_simulation/1.base_on_observed/all_pairs.thread{tid}.list", tid=THREADS),
        rnd_res = expand(STEP5_DIR / "1.run_simulation/2.base_on_random/result1.thread{tid}.simulation_on_random.list", tid=THREADS),
        rnd_pairs = expand(STEP5_DIR / "1.run_simulation/2.base_on_random/all_pairs.thread{tid}.list", tid=THREADS)
    output:
        final_list = STEP5_DIR / f"{PREFIX}.significant.interMolecular.interaction.list"
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
        prefix = f"{PREFIX}",
        threads_list = " ".join([str(t) for t in THREADS])
    threads: 16
    log: LOG_DIR / "step5_post.log"
    shell:
        """
        cd {STEP5_DIR}

        mkdir -p 2.pre-process/1.base_on_observed
        mkdir -p 2.pre-process/2.base_on_random
        mkdir -p 3.calculate_pvalue/1.base_on_observed
        mkdir -p 3.calculate_pvalue/2.base_on_random
        mkdir -p 4.recalibrate_pvalue

        echo "Phase 2: Pre-process..." >> {log}

        # Split and transpose observed simulation results
        for i in {params.threads_list}; do
            perl {params.script_split} \
                1.run_simulation/1.base_on_observed/result1.thread$i.simulation.list \
                100 thread$i &
        done
        wait

        mv *.transposed.matrix 2.pre-process/1.base_on_observed/ 2>/dev/null || true

        # Split and transpose random simulation results
        for i in {params.threads_list}; do
            perl {params.script_split} \
                1.run_simulation/2.base_on_random/result1.thread$i.simulation_on_random.list \
                100 thread$i &
        done
        wait

        mv *.transposed.matrix 2.pre-process/2.base_on_random/ 2>/dev/null || true

        echo "Phase 3: Calculate P-value..." >> {log}

        # Calculate p-value for observed interactions
        for i in {params.threads_list}; do
            perl {params.script_pval} \
                1.run_simulation/1.base_on_observed/all_pairs.thread$i.list \
                2.pre-process/1.base_on_observed/thread$i.transposed.matrix \
                thread$i &
        done
        wait

        mv comparison.observed_to_simulated.*.xls 3.calculate_pvalue/1.base_on_observed/ 2>/dev/null || true

        # Calculate p-value for random interactions
        for i in {params.threads_list}; do
            perl {params.script_pval} \
                1.run_simulation/2.base_on_random/all_pairs.thread$i.list \
                2.pre-process/2.base_on_random/thread$i.transposed.matrix \
                thread$i &
        done
        wait

        mv comparison.observed_to_simulated.*.xls 3.calculate_pvalue/2.base_on_random/ 2>/dev/null || true

        echo "Phase 4: Recalibrate..." >> {log}

        cd 4.recalibrate_pvalue

        # Merge p-values from multiple threads
        perl {params.script_merge} \
            ../3.calculate_pvalue/1.base_on_observed/comparison.*.xls \
            > comparison.observed_to_simulated.basedOnObserved.finalMerge.xls 2>> {log}

        perl {params.script_merge} \
            ../3.calculate_pvalue/2.base_on_random/comparison.*.xls \
            > comparison.observed_to_simulated.basedOnRandom.finalMerge.xls 2>> {log}

        # Convert format for multiple testing correction
        perl {params.script_fmt} \
            comparison.observed_to_simulated.basedOnObserved.finalMerge.xls \
            > result1.results_file.for_Test.list

        perl {params.script_fmt} \
            comparison.observed_to_simulated.basedOnRandom.finalMerge.xls \
            > result1.control_file.for_Test.list

        # Apply multiple testing correction (CloseCall algorithm)
        perl {params.script_corr} \
            --control result1.control_file.for_Test.list \
            --results result1.results_file.for_Test.list \
            --window 500 2>> {log}

        gzip -d -f result1.results_file.for_Test.list.window_500.qval.txt.gz

        # Replace raw p-values with corrected p-values and filter significant interactions
        perl {params.script_repl} \
            comparison.observed_to_simulated.basedOnObserved.finalMerge.xls \
            result1.results_file.for_Test.list.window_500.qval.txt \
            > result2.comparison.observed_to_simulated.Add_Local_Corrected_pvalue.list

        perl {params.script_filt} \
            result2.comparison.observed_to_simulated.Add_Local_Corrected_pvalue.list \
            {params.pval_cutoff} \
            > {output.final_list}
        """

rule multiqc_report:
    """Generate MultiQC report summarizing QC metrics from all analysis steps"""
    input:
        STEP2_DIR / f"{PREFIX}.interMolecular.sam"
    output:
        MULTIQC_DIR / "multiqc_report.html"
    container: config["sif"]
    threads: 1
    log: LOG_DIR / "multiqc.log"
    shell:
        """
        mkdir -p {MULTIQC_DIR}

        echo "=== MultiQC Analysis ===" > {log}
        echo "Prefix: {PREFIX}" >> {log}  # Fixed: Use PREFIX instead of wildcards.sample
        echo "Date: $(date)" >> {log}
        echo "" >> {log}

        # Derive parent directory from MULTIQC_DIR
        OUTPUT_ROOT=$(dirname {MULTIQC_DIR})
        echo "Scanning directory: $OUTPUT_ROOT" >> {log}

        # MultiQC recursive search
        multiqc "$OUTPUT_ROOT" \
            -o {MULTIQC_DIR} \
            --force \
            --filename multiqc_report.html \
            --ignore "*.sam" \
            --ignore "*.bam" \
            --ignore "*.fastq*" \
            --ignore "*.fq*" \
            --dirs-depth 3 \
            --verbose \
            2>> {log} || (echo "MultiQC completed with warnings" >> {log}; exit 0)

        # Verify output
        echo "" >> {log}
        if [ -f "{output}" ]; then
            echo "✅ MultiQC report generated successfully" >> {log}
            echo "Report location: {output}" >> {log}
        else
            echo "⚠️  Warning: MultiQC report may not have been generated" >> {log}
        fi

        """
