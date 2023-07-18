library(tidyverse)

########### NOTE #############################
# Yu Zhang et al., 2018, PNAS
# Filtered public BS-seq data with low bisulfite conversion rate.
# They estimated bisulfite conversion rate by using chloroplast genome.
# In the mathod part of paper, they described they used chloroplast mC < 1% as a rule of thumb
# for good bisulfite conversion rate.
# However, in the supplemental data they provided (Supplemental data 2), I found they considered both Chloroplast genome and nucleus genome. If chloroplast genome shows > 1% mC & nucleus genome shows relatively low mC level, they considered it as inefficiently bisulfite converted library.
##############################################


args <- commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

#' read input
report <- read_tsv(input, col_names=FALSE)
colnames(report) <- c("chr", "pos", "strand", "mC", "C", "context", "remove")


#' filter chloroplast
chl <- report %>%
		filter(chr == "ChrC")
chr15 <- report %>%
		filter(chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"))
#' binning by context

chl.bin <- chl %>%
		group_by(context) %>%
		summarise(mC=sum(mC), C=sum(C)) %>%
		add_row(context="C", mC=sum(.$mC), C=sum(.$C)) %>%
		mutate(mC_percent=round(100 * ( mC / ( mC + C ) ), digits=1))

chr15.bin <- chr15 %>%
		group_by(context) %>%
		summarise(mC=sum(mC), C=sum(C)) %>%
		add_row(context="C", mC=sum(.$mC), C=sum(.$C)) %>%
		mutate(mC_percent=round(100 * ( mC / ( mC + C ) ), digits=1))

all <- list(chr15=chr15.bin, chloroplast=chl.bin)
all.bind <- bind_rows(all, .id="chr")

write_tsv(all.bind, output)
