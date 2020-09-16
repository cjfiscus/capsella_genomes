#!/usr/bin/env Rscript
# plot .het from vcftools
# cjfiscus
# 2020-09-16
# usage
# Rscript plot_het.R het.txt out_prefix
# args[1] is file produced by vcftools --gzvcf "$RESULTS"/"$SPP"_filter3.vcf.gz --het --stdout
# args[2] is output prefix

args = commandArgs(trailingOnly=TRUE)
library(pacman)
p_load(ggplot2)

# read in data
df<-read.delim(args[1])

# calc prop het
df$PROP_HET<-(df$N_SITES-df$O.HOM.)/df$N_SITES

# sort by het
df<-df[order(df$PROP_HET, decreasing=F),]
df$INDV<-factor(df$INDV, levels=df$INDV)

# plot
p1<-ggplot(df, aes(x=INDV, y=PROP_HET)) + geom_bar(stat="identity") +
  theme_classic() + coord_flip()
out<-paste0(args[2], "_het.jpeg")
ggsave(out, p1, height=10, width=8)

