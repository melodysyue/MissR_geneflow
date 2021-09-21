library(tidyverse)
library(scales)
library(ggsignif) #geom_sig
library(ggpubr) #ggarrange
library(patchwork)

rm(list=ls())

###############
###Figure S2###
###############

#make LD decay plots

myld_chr <- function(aln, ld_files, name) {
  ld <- ld_files[1]
  for(i in 2: length(ld_files)){
    ld <- bind_rows(ld, ld_files[i])
  }
  
  ###Get position information
  aln <- aln %>% 
    mutate(snp=paste0("tig_",tag,":",tag.pos))
  
  ###Fill in positions in the LD table
  ld.pos <- left_join(ld, aln, by=c("SNP_A"="snp")) %>% 
    dplyr::select(chr_A=chr, chr.pos_A=chr.pos, SNP_A, CHR_B, BP_B, SNP_B, R2)
  ld.pos <- left_join(ld.pos, aln, by=c("SNP_B"="snp")) %>% 
    dplyr::select(chr_A, chr.pos_A, SNP_A, chr_B=chr, chr.pos_B=chr.pos, SNP_B, R2)
  
  ld.pos <- ld.pos %>% 
    mutate(dist_kb=abs(chr.pos_A-chr.pos_B)/1000)
  
  ###Order chromosomes###
  ld.pos$chr_A <- as.factor(ld.pos$chr_A)
  nCHR=length(levels(ld.pos$chr_A))
  #ld.pos$chr_A <- factor(ld.pos$chr_A, levels=c(seq(1,nCHR-1,1),"Unplaced")) #comment this out for gizzard shad
  
  #plot
  
  ld_chr <- ld.pos %>% 
    ggplot(aes(x=dist_kb, y=R2), alpha=0.5)+
    geom_point(col="lightgray")+
    facet_wrap(.~chr_A)+
    theme_classic(base_size = 15)+
    xlab("Distance between loci (Kbp)")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          strip.background = element_rect(size=0.5))+
    scale_x_continuous(expand=c(0,0), breaks=scales::pretty_breaks(n=3))+
    scale_y_continuous(expand=c(0,0), breaks=scales::pretty_breaks(n=3))+
    ylab(expression(r^2))+
    xlab("Distance between loci (Kbp)")+
    ggtitle(name)
  
  return(ld_chr)
  
}


#example
aln <- read.table("./data/alignment/ersn.bwa.uni.map20.aln.txt",header=T)
ld_files <- list.files(path = "./ld/", pattern = "ersn.chr*", full.names=TRUE) %>% 
  lapply(read.table,header=TRUE)
p.ersn <- myld_chr(aln,ld_files, "Emerald Shiner")

