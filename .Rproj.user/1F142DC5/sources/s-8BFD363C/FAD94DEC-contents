####Adapted from: https://github.com/ksamuk/gene_flow_linkage/tree/master/shared_functions
#### sample N SNPs (number of outliers) and calculate NND, return mean of NND among these SNPs;  

calculate.null.nnd <- function(stats.file.lg, num.outliers){
  
  site.sample <- stats.file.lg %>%
    filter(!is.na(chr.pos)) %>%
    select(chr.pos) %>%
    sample_n(num.outliers) %>%
    arrange(chr.pos) %>%
    mutate(dist.1 = c(NA,diff(chr.pos))) %>%
    mutate(dist.2 = c(diff(sort(chr.pos)),NA))
  
  nn.dist <- rep(NA, length(site.sample$chr.pos))
  for (k in 1:length(site.sample$chr.pos)){
    
    if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
      nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k]))
    }else if(is.na(site.sample$dist.1[k])){
      nn.dist[k] <- site.sample$dist.2[k]
    } else if(is.na(site.sample$dist.2[k])){
      nn.dist[k] <- site.sample$dist.1[k]
    }
  }
  return(mean(nn.dist, na.rm = TRUE))
}