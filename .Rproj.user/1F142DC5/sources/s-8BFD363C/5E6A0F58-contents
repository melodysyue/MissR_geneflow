library(tidyverse)
library(ggsignif) #geom_signif
library(scales)
library(patchwork)
library(ggpubr)
library(grid)


rm(list=ls())

aln <- read.table("./data/alignment/ersn.bwa.aln.stats.txt",header=T)
dxy <- read.table("./data/dxy/ersn.dxy.txt",header=T)
ol <- read.table("./data/outlier_neutral/ersn.fst_gea_ol.list",header=T)
ol <- ol %>% 
  mutate(snp = paste0("tig_", chrom, ":", pos))


aln.stat <- left_join(aln, dxy, by="locus")
c=9
ld <- read.table(paste0("./data/ld/ersn.chr9.ld"),header=T)

###filter ild.stat to only the specified chromosome, and order by SNP position
aln.stat <- aln.stat %>% 
  mutate(tag_id=paste0("tig_", tag, ":", tag.pos)) %>% 
  mutate(ol = tag %in% ol$chrom ) %>% 
  filter(chr==c) %>% 
  arrange(chr.pos) %>% 
  rownames_to_column("rowID") %>% 
  select(rowID, locus, tag, tag.pos, tag_id, chr, chr.pos, chr.cumpos, ol, Fstp, Ho, dxy_avg)

###Get LD information
###format ld data. 
#Get position information
aln <- aln %>% 
  select(.,locus:chr.cumpos, maf) %>% 
  mutate(snp=paste0("tig_",tag,":",tag.pos))

#Fill in positions in the LD table
ld.pos <- left_join(ld, aln, by=c("SNP_A"="snp")) %>% 
  dplyr::select(chr_A=chr, chr.pos_A=chr.pos, SNP_A, CHR_B, BP_B, SNP_B, R2)
ld.pos <- left_join(ld.pos, aln, by=c("SNP_B"="snp")) %>% 
  dplyr::select(chr_A, chr.pos_A, SNP_A, chr_B=chr, chr.pos_B=chr.pos, SNP_B,R2)

#filter loci with maf >=0.01 for LD calculation
maf <- aln %>% 
  filter(maf>=0.01) %>% 
  pull(snp)

ld.pos <- ld.pos %>% 
  mutate(dist_kb=abs(chr.pos_A-chr.pos_B)/1000) %>% 
  filter(SNP_A %in% maf & SNP_B %in%maf)


###perform permutation tests

list.files("./scripts/permutation_functions/", full.names = TRUE) %>% sapply(source) %>% invisible #load all the related functions
num_permutations <- 10000
num.outliers <- sum(aln.stat$ol)

# draw 10000 samples of num.outliers random loci, take the mean, and return the ecdf and mean
null <- replicate(num_permutations, calculate.null.metrics(aln.stat, ld.pos, num.outliers))

null <- t(null) %>% 
  as.data.frame()
colnames(null) <- c("fstp", "ho", "dxy", "ld")

head(null)


# calculate the estimate mean null 
null.mean <- apply(null, 2, mean, na.rm=TRUE)
null.ecdf <- apply(null, 2, ecdf) #empirical cumulative distribution function

# calculate the empirical nndist for real outliers
empirical.mean <- calculate.emp.metrics(aln.stat, ld.pos)
names(empirical.mean) <- c("dxy", "ho", "dxy", "ld")

#number of total loci
n.sites <- aln.stat %>% filter(!is.na(chr.pos)) %>% select(chr.pos) %>% unlist %>% length

perm.stat <- data.frame(
  chr = unique(aln.stat$chr),
  n.sites = n.sites,
  num.outliers = num.outliers,
  #fstp
  fstp.mean.emp = empirical.mean[1],
  fstp.mean.null = null.mean[1], 
  fstp.pvalue = two_side_p(null$fstp, empirical.mean[1]),
  #ho
  ho.mean.emp = empirical.mean[2],
  ho.mean.null = null.mean[2], 
  ho.pvalue = two_side_p(null$ho, empirical.mean[2]),
  #dxy
  dxy.mean.emp = empirical.mean[3],
  dxy.mean.null = null.mean[3], 
  dxy.pvalue = two_side_p(null$dxy, empirical.mean[3]),
  #ld
  ld.mean.emp = empirical.mean[4],
  ld.mean.null = null.mean[4], 
  ld.pvalue = two_side_p(null$ld, empirical.mean[4])
)


perm.stat



