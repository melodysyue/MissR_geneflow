library(tidyverse)

rm(list=ls())

##############
###Table S4###
##############

loc.sum <- read.table("./data/outlier_neutral/ersn_outlier_summary.txt", header=T)
head(loc.sum)

table(loc.sum$bs)
table(loc.sum$arl)
table(loc.sum$ofk)
table(loc.sum$pcadapt)
table(loc.sum$rda)
table(loc.sum$lfmm)
table(loc.sum$bf)


###################
###FST outliers ###
###################

#definition:by at least 2 Fst methods

fst <- loc.sum %>% 
  mutate(fst.tot=rowSums(select(.,bs:pcadapt))) %>% 
  filter(fst.tot>1)


##################
###GEA outliers###
##################

#definition:by at least 2 GEA methods

gea <- loc.sum %>% 
  mutate(gea.tot=rowSums(select(.,rda:bf))) %>% 
  filter(gea.tot>1)


##################
###Outlier SNPs###
##################
#definition:union of Fst outlier and GEA outlier

fst_gea <- loc.sum %>% 
  filter(locus %in% fst$locus | locus %in% gea$locus) %>% 
  mutate(fst.tot=rowSums(select(.,bs:pcadapt))) %>% 
  mutate(gea.tot=rowSums(select(.,rda:bf)))


##################
###Neutral SNPs###
##################

#definition:not by any methods
neu <- loc.sum %>% 
  mutate(ol.tot=rowSums(select(.,bs:bf))) %>% 
  filter(ol.tot==0) 




#########################
###prepare whitelist ####
#########################

fst.list <- fst %>% 
  select(chrom,pos)

gea.list <- gea %>% 
  select(chrom, pos)

fst_gea.list <- fst_gea %>% 
  select(chrom,pos)

neu.list <- neu %>% 
  select(chrom, pos)

write.table(fst.list, "./data/outlier_neutral/ersn.fstol.list", quote=F, row.names = F, sep="\t")
write.table(gea.list, "./data/outlier_neutral/ersn.geaol.list", quote=F, row.names = F, sep="\t")
write.table(fst_gea.list,"./data/outlier_neutral/ersn.fst_gea_ol.list", quote=F, row.names = F, sep="\t")
write.table(neu.list, "./data/outlier_neutral/ersn.neutral.list", quote=F, row.names = F, sep="\t")




