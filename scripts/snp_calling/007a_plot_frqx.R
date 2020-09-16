# plot plink.frqx
# cjfiscus
# 2020-09-12
#
# usage
# Rscript plot_frqx.R plink.frqx sites.txt 0.1 prefix
# args[1] is plink frqx file
# args[2] is list of CHR POS from bcftools query -f'%CHROM %POS\n' test.vcf > sites.txt
# args[3] is min het threshold for filtered list
# args[4] is output prefix
args = commandArgs(trailingOnly=TRUE)

library(pacman)
p_load(data.table,ggplot2, reshape2)

# read in .frqx from plink
df<-fread(args[1])
df$TOTAL<-df$`C(HOM A1)`+df$`C(HET)`+df$`C(HOM A2)`
df$FRAC_HET<-df$`C(HET)`/df$TOTAL
df$FRAC_A1<-df$`C(HOM A1)`/df$TOTAL
df$FRAC_A2<-df$`C(HOM A2)`/df$TOTAL

# read in sites from bcftools query
df1<-fread(args[2])
names(df1)<-c("CHR", "POS")
sub<-df[,c("FRAC_HET", "FRAC_A1", "FRAC_A2")]

# combine data and write out
sub<-as.data.frame(cbind(df1, sub))
out<-paste0(args[4], "_frqx.txt")
fwrite(sub, out, quote=F, sep="\t")

# list exceeding filter
sub1<-sub[sub$FRAC_HET > as.numeric(args[3]),]
out<-paste0(args[4], "_hetmin",".txt")
fwrite(sub1, out, quote=F, sep="\t")

# prepare for plotting
long<-melt(sub, id.vars=c("CHR", "POS"))

# plotdis
p1<-ggplot(long, aes(x=value, color=variable), alpha=0.5) + 
  geom_density() +
  theme_classic() +
  theme(legend.position="top")
out<-paste0(args[4], "_frqx.jpeg")
ggsave(out, p1, height=3, width=5)
