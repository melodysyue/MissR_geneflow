library(adegenet) #dapc
library(here)
library(tidyverse)
library(hierfstat)#basic.stats
library(poppr) #popsub, aboot
library(RColorBrewer)
library(ggtree)
library(ape) #write.tree
library(pegas)
library(zvau)


rm(list=ls())

##############
###Table S5###
##############

########################
## Make Genepop Object##
########################
#genepop files can be found on Dryad. 

#double check if your .gen file has a title, if it is blank, add a title. 
#Otherwise it will give you an error. Read.genepop file can take 10 min. 
#read.genepop reads .gen file and convert it into a genind object. 

gen <- read.genepop("ersn.neutral_thin.gen", ncode=3L)
gen

s.keep <- indNames(gen) #to get the list of individuals

#add pop information
p.keep =read.table("/data/popmap/popmap_filtered_ersn.txt",header=TRUE)
colnames(p.keep)=c("sample","pop")
table(s.keep == p.keep$sample)
pop(gen)=p.keep$pop
gen

#double check
nPop(gen)
popNames(gen)
table(p.keep$pop)

#################
###    Fst    ###
#################

bs <- basic.stats(gen)

##overall
names(bs)
bs$overall

##pairwise fst
ppfst <- genet.dist(gen, method="WC84")
ppfst <- as.data.frame(as.matrix(ppfst))
ppfst


##pairwise fst significance
hf <- genind2hierfstat(gen) 
pop.levels <- levels(gen$pop)
pop.levels

##for boot.ppfst to run properly, the pop column needs to be in numeric form and sorted.
pop <- sort(as.numeric(hf[,1])) 
pop
ppfst.ci<- boot.ppfst(data.frame(pop, hf[,-1]), nboot=1000, quant=c(0.025, 0.975))

pop.levels
ppfst.ci$ll
ppfst.ci$ul

