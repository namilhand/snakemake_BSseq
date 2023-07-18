#!/bin/bash

dirin="/datasets/data_4/nison/clpt/wgbs/cx_report/2022_12/bedGraph"
region="/home/nison/work/various_regions/tair10_chrom15_window10000_step5000.bed"
dirout="/datasets/data_4/nison/clpt/wgbs/mC_landscape/2022_12/10kb_mean_meth"
mkdir -p $dirout

#bedtools map : apply a function to a column from B intervals that overlap A
#-c : columns from the B file to map onto intervals in A
#-o : operation that should be applied to -c

n=0
for f in $dirin/*.bedg; do
		((n=n%32)); ((n++==0)) && wait;
        base=$(basename $f)
        outname=${base%.bedg}.map_tair10_chrom15.bedg
        bedtools map -c 4 -o mean \
                -a $region \
                -b $f \
                -null 0 \
                > $dirout/$outname &
        done
