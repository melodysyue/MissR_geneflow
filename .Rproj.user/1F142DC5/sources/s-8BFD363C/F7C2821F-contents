library(tidyverse)

rm(list=ls())

###########
###input###
###########
aln <- read.table("./data/alignment/ersn.bwa.uni.map20.aln.txt", header=T)
aln <- aln %>% 
  select(tag=snp, chr, chr.pos=pos)

#order chr
aln$chr <- as.factor(aln$chr)
nCHR=length(unique(aln$chr))
levels(aln$chr)
#comment the following out for gizzard shad
aln$chr <- factor(aln$chr, levels=c(seq(1,nCHR-1,1),"Unplaced"))
levels(aln$chr)


ol <- read.table("./data/outlier_neutral/ersn.fst_gea_ol.list", header=T)
colnames(ol) <- c("tag", "tag.pos")

aln.ol <- aln %>% 
  mutate(ol = tag %in% ol$tag)

############################################################################
###filter chromosomes with low data, e.g., n.sites >= 30, n.outliers >= 3###
############################################################################

chr.target <- aln.ol %>% 
  group_by(chr) %>% 
  summarize(n.site=n(), n.ol=sum(ol)) %>% 
  arrange(chr) %>% 
  filter(n.site>=30 & n.ol>=3) %>% 
  pull(chr) %>% 
  as.character()

#######################
###Permutation tests###
#######################
##load all the related functions
list.files("./scripts/permutation_functions/", full.names = TRUE) %>% sapply(source) %>% invisible #load all the related functions
num_permutations <- 10000

nnd.stat <- calculate_nndist_all_chr(aln.ol, chr.target, num_permutations) #it takes a few mins
nnd.stat




