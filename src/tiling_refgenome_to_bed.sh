#!/bin/bash

#TAIR10_genome="../TAIR10_chromosome_size.txt"
refgenome=$1
window=$2
step=$3
output=$4

# window = 100kb, step = 50kb
# bedtools makewindows -g $TAIR10_genome -w 100000 -s 50000 > ../TAIR10_window_100kb_step_50kb.bed
bedtools makewindows -g $refgenome -w $window -s $step > $output
