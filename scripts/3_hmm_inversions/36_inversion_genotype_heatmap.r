library(RColorBrewer)
library(vcfR) #read vcf
library(pheatmap)
library(gplots)
library(tidyverse)
library(patchwork)


rm(list=ls())

aln <- read.table("./data/alignment/ersn.bwa.uni.map20.aln.txt", header=T)
inv <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds1.chr9.ol.txt", header=T)
vcf<-read.vcfR("./data/hmm_inversions/ersn.chr19.inv.recode.vcf")
cluster <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds1.chr9.ol.cluster.genotype.txt", check.names = F, header=T)


#####################################
####Get genotype and make heatmap####
#####################################

#extract genotypes in numeric format
genotypes_numeric<-extract.gt(vcf, element = "GT", mask = FALSE, as.numeric = FALSE,
                              return.alleles = FALSE, IDtoRowNames = TRUE, extract = TRUE,
                              convertNA = TRUE)


#function to convert genotypes to 0,1,2 format used for heatmap
convertGenos<-function(genotype){
  alleles<-as.numeric(str_split_fixed(genotype,"/",2))
  genoCode<-alleles[1]+alleles[2]
  return(genoCode)
}

#convert to 0,1,2 format
genosCode<-apply(genotypes_numeric,1:2,convertGenos)
genosCode<-t(genosCode)

#format for heatmap2
heatmap2<-genosCode
heatmap2[is.na(heatmap2)] <- -1 #replace msising info with -1
heatmap2<-apply(heatmap2,1:2,function(x) as.numeric(x))
dim(heatmap2)

head(heatmap2)


###fancy heatmap
###by cluster
cluster$name <- str_remove_all(cluster$name, "X")
cluster$name <- str_replace_all(cluster$name, "\\.", "-")
cluster <- cluster %>% 
  select(sample=name, cluster)
cluster$cluster <- as.factor(cluster$cluster)


cluster_colors <- list(cluster=c("red", "purple", "blue"))
cluster_colors
names(cluster_colors$cluster)=c("0","1", "2")

#order rows by cluster
row.order <- cluster %>% 
  arrange(cluster) %>% 
  pull(sample)

heatmap2 <- heatmap2 %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")

heatmap2_sorted <- heatmap2[match(row.order, heatmap2$sample),]
heatmap2_sorted[1:5,1:5]
rownames(heatmap2_sorted) <- heatmap2_sorted$sample
heatmap2_sorted <- heatmap2_sorted[,-1]


cluster_sorted <- cluster[match(row.order, cluster$sample),]
rownames(cluster_sorted) <- cluster_sorted$sample
cluster_sorted <- cluster_sorted %>% 
  select(-sample)
head(cluster_sorted)
str(cluster_sorted)

#order columns by snp position
col.order <- aln %>% 
  filter(tag %in% inv$tag) %>% 
  select(tag, tag.pos, chr, chr.pos) %>% 
  mutate(loc=paste0(tag, "_", tag.pos)) %>% 
  arrange(chr.pos) %>% 
  select(loc, chr.pos, tag)

heatmap2_sorted <- heatmap2_sorted[, match(col.order$loc, colnames(heatmap2_sorted))]
table(colnames(heatmap2_sorted)==col.order$loc)
colnames(heatmap2_sorted)=col.order$chr.pos

heat_p <- pheatmap(heatmap2_sorted,
                   annotation_row = cluster_sorted,
                   annotation_colors  = cluster_colors,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_colnames = TRUE,
                   show_rownames = FALSE,
                   color=c("snow", "lightblue","steelblue", "darkblue"),
                   legend_breaks = c(-1,0,1,2),
                   legend_labels =c("NA","homo1","het","homo2")
)



pdf("./FigS6_ersn_chr9_ivt_genotype_heatmap.pdf", width=18, height = 6 )
heat_p
dev.off()

