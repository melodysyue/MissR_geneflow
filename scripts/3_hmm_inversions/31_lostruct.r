library(lostruct)
library(tidyverse)
library(gridExtra)
library(grid)
library(spdep)
library(gghighlight)
library(ggsignif)


rm(list=ls())

################
###read input###
################

ws=20 #ws=50 if # of aligned loci > 10,000, else ws=20.
mds=40
gen <- read.table("./data/geno_input/ersn.rda.geno",header=FALSE)
loc <- read.table("./data/rowID_tag_tagpos/ersn.chrom.pos.rowID.list",header=FALSE)
popmap <- read.table("data/popmap/popmap_filtered_ersn.txt",header=TRUE)
colnames(loc) <- c("id","tag","pos")
colnames(gen)=loc$id
rownames(gen)=popmap$sample
gen <- apply(gen, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))


##reorder the SNPs based on their positions and associated chromosomes
aln <- read.table("./data/alignment/ersn.bwa.uni.map20.aln.txt",header=T)
aln <- left_join(aln, loc, by=c("tag"))
aln <- aln %>% 
  dplyr::select(id, tag, tag.pos=pos, reference, chr, chr.pos)  

#order chromosomes  
aln$chr <- as.factor(aln$chr)
levels(aln$chr)
nCHR <- length(unique(aln$chr))
nCHR
#comment out for gizzard shad
aln$chr <- factor(aln$chr, levels=c(seq(1,nCHR-1,1),"unplaced"))
levels(aln$chr)=c(seq(1,nCHR-1,1),"Unplaced")
levels(aln$chr)

aln.chr <- aln %>% 
  filter(chr != "Unplaced")

##order rows of gen file by position of aligned SNPs
gen[1:6, 1:6]
snp <- t(gen)
snp[1:6, 1:6]
snp <- snp %>% 
  as.data.frame() %>% 
  rownames_to_column("id")

aln.snp <- left_join(aln.chr, snp, by="id")

mat <- aln.snp %>% 
  dplyr::select(id, starts_with("19")) %>% 
  column_to_rownames(var="id") %>% 
  as.matrix()
dim(mat)

##############
###lostruct###
##############
pcs <- eigen_windows(mat, k=2, win=ws) #window size,Omits the last (short) window
pcdist <- pc_dist(pcs,npc=2)
mds <- cmdscale(pcdist, eig=TRUE, k=mds)
mds.coords <- mds$points
colnames(mds.coords) <- paste0("MDS", 1:ncol(mds.coords))

yy <- mds.coords %>% 
  as.data.frame() %>% 
  mutate(window=seq(1,nrow(mds.coords),1))


dd <- aln.snp %>%
  rownames_to_column("rowID") %>% 
  dplyr::select(-starts_with("19"))

dd <- dd %>% 
  mutate(rowID=as.numeric(rowID)) %>% 
  mutate(window=floor((rowID-1)/ws + 1))


xx <- left_join(yy, dd, by="window")


#pdf(paste0("./Lostruct/ersn/ersn.Ch.lostruct.ws",ws,".mds1_40.pdf"), width = 12, height = 4)

for (i in 1:40){
  cutoff <- mean(xx[,i])+c(-1,1)*4*sd(xx[,i])
  p <- xx %>% 
    ggplot(aes(x=window, y= xx[,i])) + 
    geom_point()+
    geom_hline(yintercept = cutoff, color="red") +
    theme_bw(base_size = 15)+
    xlab("Position on Chromosomes")+
    ylab(paste0("MDS",i))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    facet_wrap(~chr, nrow=1, scales="free_x")
  print(p)
}
#dev.off()

#write.table(xx,"./Lostruct/ersn/ersn.Ch.lostruct.output.txt",quote=F, row.names = F)
#write.table(mat, "./Lostruct/ersn/ersn.Ch.lostruct.input.matrix.txt", quote=F, row.names = T)



