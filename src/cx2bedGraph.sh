#!/bin/bash

########
# Date: 2023/07/06
# Author: namilhand
# Description: Convert CX_report produced by Bismark into bedgraph format. Bedgraph files for each context will be generated as a result.
########

cx_input=$1
dir_output=$2
context=$3
temp_dir=$4
# temporary directory for bash sort.
prefix=$(basename ${cx_input})
prefix=${prefix%\.CX_report\.txt}

CNN_context=("C" "CG" "CHG" "CHH")

mkdir -p $dir_output

if [[ $context = "C" ]]; then
    awk 'BEGIN{FS="\t"; OFS="\t";} (( $4 + $5 ) >= 4) {chrom=$1; start=$2-1; end=$2; score=( $4 / ( $4 + $5 )); print chrom, start, end, score}' $cx_input | sort -T $temp_dir -k1,1 -k2,2n | sed '1s/^/track type=bedGraph\n/' > ${dir_output}/${prefix}.${context}.bedg
elif [[ ${CNN_context[@]} =~ $context ]]; then
    awk -v context=$context 'BEGIN{FS="\t"; OFS="\t";} (( $6 == context )) && (( $4 + $5 ) >= 4) {chrom=$1; start=$2-1; end=$2; score=( $4 / ( $4 + $5 )); print chrom, start, end, score}' $cx_input | sort -T $temp_dir -k1,1 -k2,2n | sed '1s/^/track type=bedGraph\n/' > ${dir_output}/${prefix}.${context}.bedg 
elif [[ $context == "nonCG" ]]; then
    awk 'BEGIN{FS="\t"; OFS="\t";} (( $6 == "CHG" ) || ( $6 == "CHH" )) && (( $4 + $5 ) >= 4) {chrom=$1; start=$2-1; end=$2; score=( $4 / ( $4 + $5 )); print chrom, start, end, score}' $cx_input | sort -T $temp_dir -k1,1 -k2,2n | sed '1s/^/track type=bedGraph\n/' > ${dir_output}/${prefix}.${context}.bedg 
else
    echo context should be one of C, CG, CHG, CHH, or nonCG
fi
