#!/bin/bash

# input bw
bw_arp6="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/03_bamCoverage/bw/mnase-arp6_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"
bw_wt="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/03_bamCoverage/bw/mnase-col0_MappedOn_tair10_nuclear_sort_md_norm_1bp-resolution.bw"
bw_control="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_MNase-control/results/03_bamCoverage/bw/MNase-ctrl-merged_1bp-resolution.bw"

# output directory
dirout="/datasets/data_4/nison/NGS_libraries/MNase/Choi18_WT_suf3-1/results/04_log2ChIP"


# misc
threads=30
toTSV="genomeBin_bedgraphToTSV.R"
refgenome_fasta="/home/nison/work/refgenome/TAIR10/TAIR10.fasta"
#======================#

mkdir -p $dirout

function bwcompare {
    chip=$1
    control=$2
    prefix=$3
	binSize=$4
	binName=$5
	dirout=$6
	toBedg=${7:-1}
	toTsv=${8:-1}

    bigwigCompare -b1 $chip -b2 $control -of bigwig \
        --binSize $binSize \
        -p $threads \
		--pseudocount 1 \
		--operation log2 \
        --skipZeroOverZero \
        -o $dirout/${prefix}_log2ChIP_binSize${binName}.bw
	if [[ $toBedg == 1 ]]; then
    bigwigCompare -b1 $chip -b2 $control -of bedgraph \
        --binSize $binSize \
        -p $threads \
		--pseudocount 1 \
		--operation log2 \
        --skipZeroOverZero \
        -o $dirout/${prefix}_log2ChIP_binSize${binName}.bg
    # --skipZeroOverZero: Skip bins where BOTH BAM files lack coverage.
	fi
	
	if [[ $toTsv == 1 ]]; then
		Rscript $toTSV $dirout/${prefix}_log2ChIP_binSize${binName}.bg $refgenome_fasta ${binSize}
	fi
}

#bwcompare $bw_arp6 $bw_control MNase_suf3-bud 1 1bp $dirout 0 0
#bwcompare $bw_wt $bw_control MNase_WT-bud 1 1bp $dirout 0 0
bwcompare $bw_arp6 $bw_control MNase_suf3-bud 100000 100kb $dirout 1 1
bwcompare $bw_wt $bw_control MNase_WT-bud 100000 100kb $dirout 1 1

