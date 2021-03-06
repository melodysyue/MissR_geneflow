library(tidyverse)
library(lostruct)
library(ggsignif)
library(patchwork)
library(ggpubr)
library(forcats) #fct_rev
rm(list=ls())

ol <- read.table("./data/hmm_inversions/ersn.danRer11.ws20.mds3.chr9.ol.txt", header=T)

mat <- read.table("./Lostruct/ersn/ersn.danRer11.lostruct.input.matrix.txt",header=T)
popmap <- read.table("./data/popmap/popmap_filtered_ersn.txt", header=T)

popmap$pop=factor(popmap$pop, levels=c("Pool04_MN","Pool08_WI","Pool13_IA","Pool26_IL","OR_MO","LG_IL"))
levels(popmap$pop) <- c("Pool 4", "Pool 8", "Pool 13", "Pool 26", "Open River", "La Grange")
popmap <- popmap %>% 
  mutate(name=paste0("X", sample))
popmap$name <- gsub("-",".", popmap$name)

#subset mat according to loci in the outlier windows and make PCA plot
mat_sub <- mat[rownames(mat)%in%ol$id,]
out <- cov_pca(mat_sub, k=2) # out is a numeric vector
#the first entry gives the total sum of eigenvectors, and the next k entries give eivgenvaluesr for each PC;
#The following values are eigenvectors for each sample in each PCA. 
PC1_perc <- out[2]/out[1]
PC2_perc <- out[3]/out[1]
samples <- colnames(mat_sub)
matrix.out <- t(matrix(out[4:length(out)], ncol=length(samples), byrow=T))
out <- as_tibble(cbind(samples, matrix.out)) %>% 
  rename(name=samples, PC1="V2", PC2="V3") %>% 
  mutate(PC1=as.double(PC1), PC2=as.double(PC2))

#kmeans to cluster along PC1 and PCA
kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
out$cluster <- kmeans_cluster$cluster - 1 #0, 1, 2
out$cluster <- as.character(out$cluster)
betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss


#PCA plot
p.pca <- out %>%
  ggplot(aes(x=PC1, y=PC2, col=cluster)) + 
  geom_point() + 
  theme_bw(base_size = 15) +
  scale_color_manual(name="Cluster", values=c("red","purple","blue")) +
  xlab(paste0("PC",1," (",round(100*PC1_perc,2),"%)")) +
  ylab(paste0("PC",2," (",round(100*PC2_perc,2),"%)")) +
  annotate(geom="text", -Inf, Inf,
           label=paste0("Discreteness=",round(betweenSS_perc, digits=4)), hjust=-0.2, vjust=2, size=4)+
  theme(legend.position = "none")


#Heterozygosity plot
mat_sub %>% as.tibble() -> snps
dim(snps)
snps_geno <- snps %>% gather("name","genotype",1:ncol(snps)) %>% 
  group_by(name, genotype) %>%
  summarize(count=n()) %>%
  spread(genotype, count)

snps_geno[is.na(snps_geno)]=0
heterozygosity <- snps_geno %>%   
  summarize(het=`1`/(`0` + `1` + `2`)) 

p.het <- inner_join(out, heterozygosity) %>% 
  ggplot(., aes(x=as.character(cluster), y=het, fill=as.character(cluster))) + 
  geom_boxplot() + 
  theme_bw(base_size = 15) +
  scale_fill_manual(name="Cluster", values=c("red","purple","blue")) +
  xlab("Cluster") + ylab("Heterozygosity") +
  geom_signif(comparisons = list(c("1", "0"), c("1","2"), c("0","2")), 
              map_signif_level=TRUE)+
  theme(legend.position = "none")

### Genetype frequency
head(popmap)

popcluster <- left_join(popmap,out, by="name")
head(popcluster)
popcluster <- popcluster %>% 
  group_by(pop, cluster) %>% 
  summarise(n=n()) %>% 
  group_by(pop) %>% 
  mutate(total=sum(n)) %>% 
  mutate(freq=n / total)
popcluster$cluster <- as.character(popcluster$cluster)

p.freq <- popcluster %>% 
  ggplot(aes(x=fct_rev(pop), y=freq, fill=cluster)) +
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(name="cluster", values=c("red","purple","blue")) +
  xlab("")+
  ylab("Frequency")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic(base_size=15)+
  theme(legend.position = "none")+
  coord_flip()


pdf("./Fig4_chr9_ABC.pdf", width = 12, height=6)
p.pca + p.het + p.freq + 
  plot_annotation(tag_levels="A")
dev.off()


popcluster %>% 
  select(-freq) %>% 
  spread(key=cluster,value=n)

