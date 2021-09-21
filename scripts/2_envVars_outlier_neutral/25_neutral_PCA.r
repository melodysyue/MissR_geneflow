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
library(ggpubr)

rm(list=ls())

#genepopfile (neutral SNPs after thinning) can be downloaded from Dryad. 
gen.files <- list.files(path="./genepop_files/neutralSNPs_thin/", pattern = ".gen", full.names = TRUE)
p.keep.files <- list.files(path="./data/popmap", pattern="popmap_filtered", full.names = TRUE)
species.name <- c("Bullhead Minnow", "Bluegill", "Freshwater Drum", "Channel Catfish", "Gizzard Shad", "Emerald Shiner")

#make sure they are in the same order
gen.files
names(gen.files)=c("./genepop_files/bhmw.neutral.gen", "./genepop_files/blgl.neutral.gen", "./genepop_files/cncf.neutral.gen", 
                   "./genepop_files/ersn.neutral.gen", "./genepop_files/fwdm.neutral.gen", "./genepop_files/gzsd.neutral.gen")

gen.files.sorted <- gen.files[c("./genepop_files/bhmw.neutral.gen", "./genepop_files/blgl.neutral.gen", "./genepop_files/fwdm.neutral.gen",
                                "./genepop_files/cncf.neutral.gen", "./genepop_files/gzsd.neutral.gen", "./genepop_files/ersn.neutral.gen")]

gen.files.sorted


p.keep.files
names(p.keep.files)=c("./pop_info/popmap_filtered_bhmw.txt", "./pop_info/popmap_filtered_blgl.txt", "./pop_info/popmap_filtered_cncf.txt", 
                      "./pop_info/popmap_filtered_ersn.txt","./pop_info/popmap_filtered_fwdm.txt", "./pop_info/popmap_filtered_gzsd.txt")

p.keep.files.sorted <- p.keep.files[c("./pop_info/popmap_filtered_bhmw.txt", "./pop_info/popmap_filtered_blgl.txt", "./pop_info/popmap_filtered_fwdm.txt",
                                      "./pop_info/popmap_filtered_cncf.txt", "./pop_info/popmap_filtered_gzsd.txt", "./pop_info/popmap_filtered_ersn.txt")]
p.keep.files.sorted

species.name


l <- length(gen.files.sorted)

#create a list to store plots
pcaList <- list()

for (i in 1:l){
  gen <- read.genepop(gen.files.sorted[i], ncode=3L)
  p.keep <- read.table(p.keep.files.sorted[i], header=TRUE)
  
  s.keep <- indNames(gen) #to get the list of individuals
  colnames(p.keep)=c("sample","pop")
  levels(p.keep$pop)
  p.keep <- p.keep[p.keep$sample %in% s.keep,]
  pop(gen)=p.keep$pop
  popNames(gen)=c("Pool 4", "Pool 8", "Pool 13", "Pool 26", "Open River", "La Grange")
  print(gen)
  
  #PCA
  trans <- scaleGen(gen, NA.method="mean")
  pca_obj <- dudi.pca(trans,cent=FALSE,scale=FALSE,scannf=FALSE, nf=3)
  eigenvects_to_plot=pca_obj$li
  variance_explained=(pca_obj$eig)/sum(pca_obj$eig)
  xlabel=paste0("PC",1," (",round(variance_explained[1],3)*100,"%)")
  ylabel=paste0("PC",2," (",round(variance_explained[2],3)*100,"%)")
  eigenvects_to_plot$pop <- gen$pop 
  pca_p <- 
    ggplot(eigenvects_to_plot, aes(x=Axis1, y=Axis2, color=pop)) + #don't use data$column inside aes
    geom_point(size=2, alpha=0.8)+
    labs(x=xlabel, y=ylabel,
         col="Study Reaches")+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(xintercept=0, linetype="dashed")+
    theme_bw(base_size=12)+
    ggtitle(label=species.name[i])+ #label each plot with species name
    scale_color_brewer(palette="Dark2")+
    theme(legend.position = "bottom")+
    guides(colour=guide_legend(nrow=1)) #legend as one row
  pcaList[[i]]=pca_p
}

pdf("Fig3_PCA_neutralSNPs.pdf", width = 12, height=9)
do.call("ggarrange", c(pcaList, ncol=3, nrow=2, common.legend=TRUE, legend="bottom"))
dev.off()


