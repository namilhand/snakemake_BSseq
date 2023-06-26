# Snakemake workflow for Bismark

# Usage (snakemake --cores should reflect available cores):
# conda env create --file environment.yaml --name bsseq
# conda activate bsseq
# snakemake -p --cores 30 
# conda deactivate

# Specify config file parameters
configfile: "config.yaml"
## library info
SAMPLE = config["SAMPLES"]
REFGENOME = config["GENOME"]["refgenome"]

# Specify the desired end target file(s)
rule all:
    input:
        # expand("results/bams/{sample}",
        #         sample = SAMPLE),
        expand("results/bams/{sample}_pe.bam",
                sample = SAMPLE),
        expand("results/bams/{sample}_PE_report.txt",
                sample = SAMPLE),
        expand("results/bams/{sample}.dedup.bam",
                sample = SAMPLE),
        expand("results/bams/{sample}.dedup.report.txt",
                sample = SAMPLE),
        expand("results/methylext/{sample}",
                sample = SAMPLE)

rule bismark_pe:
    """Align BS-seq reads using Bismark"""
    output:
        # bam="results/bams/{sample}_pe.bam",
        # report="results/bams/{sample}_PE_report.txt"
        bam="results/bams/{sample}_pe.bam",
        report="results/bams/{sample}_PE_report.txt",
    input:
        fq_1="raw/{sample}/{sample}_1.fq.gz",
        fq_2="raw/{sample}/{sample}_2.fq.gz"
        # fq_1="exercise_raw/subsample/{sample}_1.fq",
        # fq_2="exercise_raw/subsample/{sample}_2.fq"
    log:
        "log/bams/{sample}.log"
    resources:
        tmpdir = "tmp"
    threads: config["THREADS"]["alignment"]
    params:
        N=1,
        temp_dir="tmp/{sample}",
        basename="{sample}"
    shell:
        r"""
        mkdir -p "results/bams"
        bismark --genome {REFGENOME} \
            -o results/bams \
            -N {params.N} \
            -1 {input.fq_1} \
            -2 {input.fq_2} \
            --temp_dir {params.temp_dir} \
            --multicore {threads} \
            &> {log}
        mv results/bams/{wildcards.sample}_1_bismark_bt2_PE_report.txt {output.report}
        mv results/bams/{wildcards.sample}_1_bismark_bt2_pe.bam {output.bam}
        """
        
rule deduplicate_bismark:
    input: "results/bams/{sample}_pe.bam"
    output:
        bam="results/bams/{sample}.dedup.bam",
        report="results/bams/{sample}.dedup.report.txt"
    log:
        "log/dedup/{sample}.dedup.log"
    resources:
        tmpdir="tmp"
    wrapper:
        "v2.0.0/bio/bismark/deduplicate_bismark"

rule methylext:
    output:
        directory("results/methylext/{sample}")
    input: "results/bams/{sample}.dedup.bam"
    log:
        "log/methylext/{sample}.methylext.log"
    resources:
        tmpdir = "tmp"
    threads: config["THREADS"]["methylext"]
    params:
        refgenome=REFGENOME,
        CG_methylext="CpG_context_{sample}.dedup.txt",
        CHG_methylext="CHG_context_{sample}.dedup.txt",
        CHH_methylext="CHH_context_{sample}.dedup.txt"
    shell:
        r"""
        mkdir -p {output}
        bismark_methylation_extractor \
            -p \
            --comprehensive \
            --bedGraph \
            --CX \
            --cytosine_report \
            --genome_folder {params.refgenome} \
            --multicore {threads} \
            -o {output} \
            {input} &> {log}
        rm results/methylext/{wildcards.sample}/{params.CG_methylext}
        rm results/methylext/{wildcards.sample}/{params.CHG_methylext}
        rm results/methylext/{wildcards.sample}/{params.CHH_methylext}
        """
    # --bedGraph: The output of the methylation extractor will be transformed into a bedgraph and coverage file.
    # Coverage file contains two additional column to the bedGraph format: [count methylated] [count unmethylated]
    # By default, --bedGraph only consider CpG, but it can be extended to cytosines in any context by using the option --CX.
    # --cytosine_report: Every cytosines will be considered irrespective of whethere they were actually covered by any reads in the experiment or not. As for the bedGraph mode, --CX enables cytosines in any context to be considered.
