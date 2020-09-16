#!/usr/bin/env Rscript
# plot depth from bcftools
# cjfiscus
# 2020-09-16
# usage
# Rscript plot_dp.R depth.txt out_prefix
# args[1] is file produced by bcftools query -f '%CHROM %POS %DP\n' test.vcf
# args[2] is output prefix
args = commandArgs(trailingOnly=TRUE)

library(pacman)
p_load(ggplot2)

# read in data
df<-read.table(args[1])

# hard cut threshold
thres<-quantile(df$V3)[4] + 1.5*IQR(df$V3)

# plot with threshold 
p1<-ggplot(df, aes(x=V3)) + geom_density() + 
  geom_vline(aes(xintercept=thres), color="red") +
  theme_classic() +
  xlab("INFO/DP")
out<-paste0(args[2], "_depth.jpeg")
ggsave(out, p1, height=5, width=5)

# write out depth value
out<-paste0(args[2], "_depth_cutval.txt")
write.table(as.numeric(thres), out, quote=F, row.names=F, col.names=F)
