library(tidyverse)

rm(list=ls())
#data for HMM states, usually 1 for low, 2 for intermediate, 3 for high. But, do check, cause sometimes state 1 is for high.
state <- read.table("./data/hmm_inversions/ersn.aln.Fstp_3state_HMMstates.txt",header=T)
#alignment
aln <- read.table("./data/alignment/ersn.bwa.aln.stats.txt", header=T)
#ol
ol <- read.table("./data/outlier_neutral/ersn.fst_gea_ol.list", header=T)
#ol summary
ol.s <- read.table("./data/outlier_neutral/ersn_outlier_summary.txt", header=T)

#############################################
###decide on which state is the high state###
#############################################
aln.stat <- cbind(aln, state)
aln.stat <- left_join(aln.stat, ol.s, by=c("locus", "tag"="chrom", "tag.pos"="pos"))
##make state as factor
aln.stat$x <- as.factor(aln.stat$x)
levels(aln.stat$x)

aln.stat %>% 
  group_by(x) %>% 
  summarize(mean=mean(Fstp))

hs <- aln.stat %>% 
  group_by(x) %>% 
  summarize(avg=mean(Fstp)) %>% 
  top_n(1) %>% 
  pull(x) %>% 
  as.character()

#######################################
###Sig Island SNPs excluding PCAdapt###
#######################################
#PCAdapt seems to have a lot of false positive
#Only bayescan, arlequin and outflank.

##identify islands of differentiation as high state & contains at least ONE fst.ol
aln.stat <- aln.stat %>% 
  rownames_to_column("pos.index") %>% 
  mutate(ol=tag %in% ol$chrom) %>% 
  mutate(id=paste0(tag,"_", tag.pos)) %>% 
  mutate(fst.ol=rowSums(select(.,bs:ofk))) %>%  #no pcadapt
  mutate(sig.state=ifelse(x==hs & fst.ol > 0,"high", "non-high")) 

aln.stat$pos.index=as.numeric(aln.stat$pos.index)

###########################################
###Extract Islands along with their SNPs###
###########################################

#create a column called "window" based on change of states
aln.stat$win <- cumsum(c(1,as.numeric(diff(aln.stat$x))!=0))

#generate info for windows
aln.stat <- aln.stat %>% 
  group_by(win) %>% 
  mutate(win.size=n()) 

#order chr
aln.stat$chr <- as.factor(aln.stat$chr)
nCHR=length(unique(aln.stat$chr))
levels(aln.stat$chr)
#comment for gizzard shad
aln.stat$chr <- factor(aln.stat$chr, levels=c(seq(1,nCHR-1,1),"Unplaced"))
levels(aln.stat$chr)

write.table(aln.stat, "./data/alignment/ersn.aln.hmm.islands.stats.txt",quote=F, row.names = F)


###################
##filter islands###
###################
#Remove islands if it is a single-SNP island and it is the only island on that chromosomes.
#Remove islands located on the unplaced chromosomes

# extract SNPs in the island windows
win.target <- aln.stat %>% 
  filter(sig.state=="high") %>% 
  distinct(win) %>% 
  pull(win)

wl <- aln.stat %>% 
  filter(win %in% win.target) %>% 
  dplyr::select(win, win.size, fst.ol,bs:bf, ol, x, sig.state, id, tag, tag.pos, chr, chr.pos) 


#how many islands after filtering?
win.target.f <- wl %>% 
  group_by(chr) %>% 
  mutate(n.win.chr=n()) %>% 
  filter(n.win.chr>1) %>% 
  filter(chr!="Unplaced") %>% 
  pull(win) %>% 
  unique()
length(win.target.f)

#how many chromosomes were these islands located?
chr.target <- wl %>% 
  group_by(chr) %>% 
  mutate(n.win.chr=n()) %>% 
  filter(n.win.chr>1) %>% 
  filter(chr!="Unplaced") %>% 
  pull(chr) %>% 
  unique()
length(chr.target)


wl %>% 
  filter(win %in% win.target.f) %>% 
  select(win, chr) %>% 
  unique() %>% 
  group_by(chr) %>% 
  summarise(n=n())


#distribution of outlier SNPs among these HMM islands
win.ol <- wl %>% 
  filter(win %in% win.target.f & ol==TRUE) %>% 
  select(win, win.size, chr, chr.pos)

nrow(win.ol)

length(unique(win.ol$win))
length(unique(win.ol$chr))

win.ol %>% 
  group_by(win, chr) %>% 
  summarize(n=n()) %>% 
  arrange(desc(n)) %>% 
  group_by(chr) %>% 
  summarize(sum=sum(n))


write.table(win.target.f, "./data/hmm_inversions/ersn.aln.hmm.islands.target.list", quote=F, row.names = F)


