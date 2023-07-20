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
REFGENOME_SIZE = config["GENOME"]["refgenome_size"]
REFGENOME_NAME = config["GENOME"]["name"]

## tiled genome
WINDOW_SIZE_NAME = config["GENOMETILE"]["window_size_name"]
STEP_SIZE_NAME = config["GENOMETILE"]["step_size_name"]


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
                sample = SAMPLE),
        expand("qc/bisulfite_conversion_rate/{sample}_bs-conversion-rate.txt",
                sample = SAMPLE),
        expand("results/bedg/{sample}.dedup.{context}.bedg",
                sample = SAMPLE,
                context = ["C", "CG", "CHG", "CHH", "nonCG"]),
        expand("results/bw/{sample}.dedup.{context}.bw",
                sample = SAMPLE,
                context = ["C", "CG", "CHG", "CHH", "nonCG"]),
        expand("data/tiled_genome_bed/{refgenome}_window_{window_size_name}_step_{step_size_name}.bed",
                refgenome = REFGENOME_NAME,
                window_size_name = WINDOW_SIZE_NAME,
                step_size_name = STEP_SIZE_NAME),
        expand("results/bedg/tiled/window_{window_size_name}_step_{step_size_name}/{sample}_{context}_{refgenome}_window{window_size_name}_step{step_size_name}.bedg",
                sample = SAMPLE,
                window_size_name = WINDOW_SIZE_NAME,
                step_size_name = STEP_SIZE_NAME,
                refgenome = REFGENOME_NAME,
                context = ["C", "CG", "CHG", "CHH", "nonCG"]),
        expand("results/tsv/tiled/window_{window_size_name}_step_{step_size_name}/{sample}_{context}_{refgenome}_window{window_size_name}_step{step_size_name}.tsv",
                sample = SAMPLE,
                window_size_name = WINDOW_SIZE_NAME,
                step_size_name = STEP_SIZE_NAME,
                refgenome = REFGENOME_NAME,
                context = ["C", "CG", "CHG", "CHH", "nonCG"]),
        expand("results/metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.gz",
                sample = SAMPLE,
                region = ["longTE"],
                refpoint = ["TSS"],
                context = ["CG", "CHG"]),
        expand("results/metaprofile/{region}/tsv/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.tsv",
                sample = SAMPLE,
                region = ["longTE"],
                refpoint = ["TSS"],
                context = ["CG", "CHG"])
        #context = ["C", "CG", "CHG", "CHH"])

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
        directory = directory("results/methylext/{sample}"),
        cx_report = "results/methylext/{sample}/{sample}.dedup.CX_report.txt"
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
        mkdir -p {output.directory}
        bismark_methylation_extractor \
            -p \
            --comprehensive \
            --bedGraph \
            --CX \
            --cytosine_report \
            --genome_folder {params.refgenome} \
            --multicore {threads} \
            -o {output.directory} \
            {input} &> {log}
        rm results/methylext/{wildcards.sample}/{params.CG_methylext}
        rm results/methylext/{wildcards.sample}/{params.CHG_methylext}
        rm results/methylext/{wildcards.sample}/{params.CHH_methylext}
        """
    # --bedGraph: The output of the methylation extractor will be transformed into a bedgraph and coverage file.
    # Coverage file contains two additional column to the bedGraph format: [count methylated] [count unmethylated]
    # By default, --bedGraph only consider CpG, but it can be extended to cytosines in any context by using the option --CX.
    # --cytosine_report: Every cytosines will be considered irrespective of whethere they were actually covered by any reads in the experiment or not. As for the bedGraph mode, --CX enables cytosines in any context to be considered.

rule bsconversion_rate:
    output:
        "qc/bisulfite_conversion_rate/{sample}_bs-conversion-rate.txt",
    input:
        "results/methylext/{sample}/{sample}.dedup.CX_report.txt"
    shell:
        r"""
        Rscript src/bisulfite_conversion_rate.R {input} {output}
        """

rule cx2bedg:
    output:
        "results/bedg/{sample}.dedup.{context}.bedg"
    input:
        "results/methylext/{sample}/{sample}.dedup.CX_report.txt"
    params:
        context="{context}",
        temp_dir="tmp"
    shell:
        r"""
        bash src/cx2bedGraph.sh {input} results/bedg {params.context} {params.temp_dir}
        """

rule bedg2bw:
    output:
        "results/bw/{sample}.dedup.{context}.bw"
    input:
        "results/bedg/{sample}.dedup.{context}.bedg"
    params:
        refgenome_size=REFGENOME_SIZE
    shell:
        r"""
        bedGraphToBigWig {input} {params.refgenome_size} {output}
        """
rule tiling_genome:
    output:
        "data/tiled_genome_bed/{refgenome}_window_{window_size_name}_step_{step_size_name}.bed"
    input:
        config["GENOME"]["refgenome"]
    params:
        refgenome = REFGENOME_SIZE,
        window = config["GENOMETILE"]["window_size"],
        step = config["GENOMETILE"]["step_size"]
    shell:
        r"""
        bash src/tiling_refgenome_to_bed.sh {params.refgenome} {params.window} {params.step} {output}
        """

rule tiling_cx_bedg:
    output:
        "results/bedg/tiled/window_{window_size_name}_step_{step_size_name}/{sample}_{context}_{refgenome}_window{window_size_name}_step{step_size_name}.bedg"
    input:
        bedg = "results/bedg/{sample}.dedup.{context}.bedg",
        bed = "data/tiled_genome_bed/{refgenome}_window_{window_size_name}_step_{step_size_name}.bed"
    shell:
        # -c: the column from the B file to map onto intervals in A.
        r"""
        bedtools map -c 4 -o mean \
                -a {input.bed} \
                -b {input.bedg} \
                > {output}
        """

rule bedg2tsv:
    # Caution: It only works when window size = step size
    output:
        "results/tsv/tiled/window_{window_size_name}_step_{step_size_name}/{sample}_{context}_{refgenome}_window{window_size_name}_step{step_size_name}.tsv"
    input:
        "results/bedg/tiled/window_{window_size_name}_step_{step_size_name}/{sample}_{context}_{refgenome}_window{window_size_name}_step{step_size_name}.bedg"
    params:
        refgenome = config["GENOME"]["refgenome_fasta"],
        window_size = config["GENOMETILE"]["window_size"]
    shell:
        r"""
        Rscript src/genomeBin_bedgraphToTSV.R {input} {params.refgenome} {params.window_size} {output}
        """

rule metaprofile:
    output:
        gzip = "results/metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.gz",
        tab = "results/metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.tab",
        bed = "results/metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.bed"
    input:
        "results/bw/{sample}.dedup.{context}.bw"
    params:
        refpoint="{refpoint}",
        region = config["METAPROFILE"]["region"],
        binSize = config["METAPROFILE"]["binSize"]
    threads:
        config["METAPROFILE"]["threads"]
    log:
        "log/metaprofile/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.log"
    shell:
        r"""
        computeMatrix reference-point \
                --referencePoint {params.refpoint} \
                -b 2000 -a 2000 \
                -R {params.region} \
                -S {input} \
                --binSize {params.binSize} \
                --sortRegions "keep" \
                --sortUsing "mean" \
                --samplesLabel {wildcards.sample} \
                --nanAfterEnd \
                -p {threads} \
                --outFileName {output.gzip} \
                --outFileNameMatrix {output.tab} \
                --outFileSortedRegions {output.bed} &> {log}
        """
rule avgBins:
    output: "results/metaprofile/{region}/tsv/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.tsv"
    input: "results/metaprofile/{region}/{sample}_metaprofile_region-{region}_referencePoint-{refpoint}_{context}.tab"
    params:
        upstream = config["METAPROFILE"]["upstream"],
        downstream = config["METAPROFILE"]["downstream"],
        binSize = config["METAPROFILE"]["binSize"]
    shell:
        r"""
        Rscript src/averageComputeMatrixTab.R {input} \
                {params.upstream} {params.downstream} \
                {params.binSize} \
                {output}
        """
