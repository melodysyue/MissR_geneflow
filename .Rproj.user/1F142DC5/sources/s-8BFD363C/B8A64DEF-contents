####Adapted from: https://github.com/ksamuk/gene_flow_linkage/tree/master/shared_functions
#### sample N SNPs (number of outliers) and extract Fstp, Ho, Dxy, LD for these SNPs, and return the mean among these SNPs.

calculate.null.metrics <- function(stats.file.chr, ld.pos, num.outliers){
  
  site.sample <- stats.file.chr %>%
    sample_n(num.outliers) %>% 
    select(tag_id, Fstp, Ho, dxy_avg)
  
  ld.pos.sample <- ld.pos %>% 
    filter(SNP_A %in% site.sample$tag_id & SNP_B %in% site.sample$tag_id)
  
  
  #Fstp, Ho, dxy, ld
  mean <- c(mean(site.sample$Fstp),
                 mean(site.sample$Ho),
                 mean(site.sample$dxy_avg),
                 mean(ld.pos.sample$R2))
  
  return(mean)
}