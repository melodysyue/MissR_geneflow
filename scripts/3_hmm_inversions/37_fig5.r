library(adegenet) #dapc
library(here)
library(tidyverse)
library(poppr) #popsub, aboot
library(RColorBrewer)
library(genepopedit)


rm(list=ls())

#vcf file can be downloaded from Dryad. 
path1 <- "./data/geno_input/ersn.noclone.G50.l90.G60.l70.G70.l50.mac3.recode.singleton.KH.gen"
outdir <- "./data/geno_input/"

genepop_detective(path, "Loci")


#loci list
aln <- read.table("./data/alignment/ersn.bwa.uni.map20.aln.txt", header=T)
loc <- read.table("./data/rowID_tag_tagpos/ersn.chrom.pos.rowID.list",header=F)
colnames(loc) <- c("id","tag","tag.pos")
aln.loc <- left_join(aln,loc,by=c("tag")) %>% 
  dplyr::select(id, tag, tag.pos, ref.chr=chr, chr.pos, reference)

chr6 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds3.chr6.ol.txt", header=T)
chr9 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds1.chr9.ol.txt", header=T)
chr19 <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds2.chr19.ol.txt", header=T)

all <- aln.loc %>% 
  pull(id) %>% 
  as.vector()

chr6 <- chr6 %>% 
  pull(id) %>% 
  as.vector()

chr9 <- chr9 %>% 
  pull(id) %>% 
  as.vector()

chr19 <- chr19 %>% 
  pull(id) %>% 
  as.vector()

three <- c(chr6, chr9,chr19)


#subset genepop
subset_genepop(genepo=path1, keep=TRUE, subs=all, path=paste0(outdir, "ersn.all_aln.gen"))

path2 <- "./genetic_input/ersn.all_aln.gen"
subset_genepop(genepo=path2, keep=TRUE, subs=three,path=paste0(outdir, "ersn.3invt.gen"))
subset_genepop(genepo=path2, keep=FALSE, subs=three,path=paste0(outdir, "ersn.no3Invt.gen"))



###PCA plot
popmap <- read.table("./data/popmap/popmap_filtered_ersn.txt", header=T)

mypca <- function(gen, popmap, label){
  s.keep <- indNames(gen) #to get the list of individuals
  pop(gen)=popmap$pop
  popNames(gen)=c("Pool 4", "Pool 8", "Pool 13", "Pool 26", "Open River", "La Grange")
  
  #PCA
  trans <- scaleGen(gen, NA.method="mean")
  pca_obj <- dudi.pca(trans,cent=FALSE,scale=FALSE,scannf=FALSE, nf=3)
  eigenvects_to_plot=pca_obj$li
  variance_explained=(pca_obj$eig)/sum(pca_obj$eig)
  xlabel=paste0("PC",1," (",round(variance_explained[1],3)*100,"%)")
  ylabel=paste0("PC",2," (",round(variance_explained[2],3)*100,"%)")
  eigenvects_to_plot$pop <- gen$pop 
  pca_p <- 
    #ggplot(eigenvects_to_plot, aes(x=Axis1, y=Axis2, color=pop)) + #don't use data$column inside aes
    #ggplot(eigenvects_to_plot, aes(x=Axis1, y=-Axis2, color=pop)) + #flip y axis for p.all
    ggplot(eigenvects_to_plot, aes(x=-Axis1, y=Axis2, color=pop)) + #flip x axis for p.3invt
    geom_point(size=2, alpha=0.8)+
    labs(x=xlabel, y=ylabel,
         col="Study Reaches")+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_vline(xintercept=0, linetype="dashed")+
    theme_bw(base_size=15)+
    ggtitle(label=label)+ #label each plot with species name
    scale_color_brewer(palette="Dark2")+
    theme(legend.position = "bottom")+
    guides(colour=guide_legend(nrow=1))
  return(pca_p)
}


gen.all <- read.genepop("./data/geno_input/ersn.all_aln.gen", ncode=3L)
gen.3invt <- read.genepop("./data/geno_input/ersn.3invt.gen", ncode=3L)
gen.no3invt <- read.genepop("./data/geno_input/ersn.no3Invt.gen", ncode=3L)


gen.all
gen.3invt
gen.no3invt

p.all <- mypca(gen.all, popmap, "All loci")
p.3invt <- mypca(gen.3invt, popmap, "Three inversions")
p.no3invt <- mypca(gen.no3invt, popmap, "No inversions")



pdf("./Fig5.pdf", width = 12, height = 6)
p.all+p.3invt + p.no3invt+
  plot_layout(guides="collect") & theme(legend.position = "bottom")
dev.off()

