library(tidyverse)
library(lostruct)


rm(list=ls()
   
#ws=50 if # of aligned loci > 10,000, else ws=20.
#ws is 20 bp for bullhead minnow, emerald shiner and gizard shad, 
#and 50 bp for  bluegill, channel catfish and freshwater drum. 
)
ws=50
#xx <- read.table("./Lostruct/ersn/ersn.Lc.lostruct.output.txt",header=T)
#mat <- read.table("./Lostruct/ersn/ersn.Lc.lostruct.input.matrix.txt",header=T)

xx.long <- xx %>% 
  gather(MDS_axis, MDS_loading, MDS1:MDS40, factor_key = TRUE)


for (i in 1:40){
  #outlier cutoff
  cutoff <- mean(xx[,i])+c(-1,1)*4*sd(xx[,i])
  xx.olwd <- xx.long %>% 
    filter(MDS_axis==paste0("MDS",i)) %>% 
    filter(MDS_loading < cutoff[1] | MDS_loading > cutoff[2]) 
  
  #extract outlier windows
  if(nrow(xx.olwd)>0){
    nchr <- unique(xx.olwd$chr)
    nchr
    
    #plot per chr
    for (c in nchr){
      l <- xx.olwd %>% 
        filter(chr==c) %>% 
        dplyr::select(id, tag, tag.pos)
      
      #subset mat according to loci in the outlier windows and make PCA plot
      mat_sub <- mat[rownames(mat)%in%l$id,]
      out <- cov_pca(mat_sub, k=2) # out is a numeric vector
      PC1_perc <- out[2]/out[1]
      PC2_perc <- out[3]/out[1]
      samples <- colnames(mat_sub)
      matrix.out <- t(matrix(out[4:length(out)], ncol=length(samples), byrow=T))
      out <- as_tibble(cbind(samples, matrix.out)) %>% 
        rename(name=samples, PC1="V2", PC2="V3") %>% 
        mutate(PC1=as.double(PC1), PC2=as.double(PC2))
      
      #kmeans to cluster along PC1 and PCA
      try_3_clusters <-try(kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1]))))
      if("try-error" %in% class(try_3_clusters)){
        kmeans_cluster <-kmeans(matrix.out[,1], 2, centers=c(min(matrix.out[,1]), max(matrix.out[,1])))
      }else{
        kmeans_cluster <- kmeans(matrix.out[,1], 3, centers=c(min(matrix.out[,1]), (min(matrix.out[,1])+max(matrix.out[,1]))/2, max(matrix.out[,1])))
      }
      out$cluster <- kmeans_cluster$cluster - 1
      out$cluster <- as.character(out$cluster)
      betweenSS_perc <- kmeans_cluster$betweenss / kmeans_cluster$totss
      
      
      #filter outlier windows if 3 clusters are highly discrete
      if(betweenSS_perc>0.89){
        
        ###PCA plot
        pca_plot <- out %>%
          ggplot(., aes(x=PC1, y=PC2, col=cluster)) + 
          geom_point() + 
          theme_bw(base_size = 15) +
          scale_color_manual(name="Cluster", values=c("red","purple","blue")) +
          xlab(paste0("PC",1," (",round(100*PC1_perc,2),"%)")) +
          ylab(paste0("PC",2," (",round(100*PC2_perc,2),"%)")) +
          annotate(geom="text", -Inf, Inf,
                   label=paste0("Discreteness=",round(betweenSS_perc, digits=4)), hjust=-0.2, vjust=2, size=5)+
          theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches"))+
          labs(tag = "(a)")
        
        #boxplot to compare heterozygosity among clusters
        mat_sub %>% as.tibble() -> snps
        dim(snps)
        colnames(snps)
        snps_geno <- snps %>% gather("name","genotype",1:ncol(snps)) %>% 
          group_by(name, genotype) %>%
          summarize(count=n()) %>%
          spread(genotype, count)
        
        snps_geno[is.na(snps_geno)]=0
        heterozygosity <- snps_geno %>%   
          summarize(het=`1`/(`0` + `1` + `2`)) 
        het_plot <- inner_join(out, heterozygosity) %>% 
          ggplot(., aes(x=as.character(cluster), y=het, fill=as.character(cluster))) + 
          geom_boxplot() + 
          theme_bw(base_size = 15) +
          scale_fill_manual(name="Cluster", values=c("red","purple","blue")) +
          xlab("Cluster") + ylab("Heterozygosity") +
          geom_signif(comparisons = list(c("1", "0"), c("1","2"), c("0","2")), 
                      map_signif_level=TRUE)+
          theme(plot.margin = unit(c(0.3,0.3,0.3,0.3), "inches")) + labs(tag = "(b)")
        
        #output loci info in the outlier windows
        write.table(l,
                    paste0("./Lostruct/ersn/ersn.Lc.ws", ws,".mds",i,".chr",c,".ol.txt"), quote=F, row.names = F)
        
        write.table(out, 
                    paste0("./Lostruct/ersn/ersn.Lc.ws", ws,".mds",i,".chr",c, ".ol.cluster.genotype.txt"),
                    quote=F, row.names = F)
        #output figures
        pdf(paste0("./Lostruct/ersn/ersn.Lc.ws", ws,".mds",i,".chr",c,".ol.pdf"), width = 16, height = 8)
        grid.arrange(pca_plot, het_plot, 
                     nrow=1)
        dev.off()
      }
    }
  }
}





