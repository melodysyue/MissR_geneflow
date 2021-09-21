####Adapted from: https://github.com/ksamuk/gene_flow_linkage/tree/master/shared_functions
#### calculate empirical nnds from a subsetted stats file

calc_emp_nnd_dist <- function(stats.file.lg){
  
  stats.file.lg <- stats.file.lg %>%
    filter(!is.na(chr.pos))
  
  site.sample <- stats.file.lg %>%
    filter(!is.na(chr.pos)) %>%
    select(chr.pos) %>%
    arrange(chr.pos) %>%
    mutate(dist.1 = c(NA,diff(chr.pos))) %>% #distance with the following SNP
    mutate(dist.2 = c(diff(sort(chr.pos)),NA)) #distance with the prior SNP
  
  nn.dist <- rep(NA, length(site.sample$chr.pos))
  for (k in 1:length(site.sample$chr.pos)){
    
    if(!is.na(site.sample$dist.1[k]) & !is.na(site.sample$dist.2[k])){
      nn.dist[k] <- min(c(site.sample$dist.1[k],site.sample$dist.2[k])) #min of (before and after)
    }else if(is.na(site.sample$dist.1[k])){
      nn.dist[k] <- site.sample$dist.2[k] #for the 1st SNP, use the dist.2
    } else if(is.na(site.sample$dist.2[k])){
      nn.dist[k] <- site.sample$dist.1[k] #for the last SNP, use the dist.1
    }
  }
  
  return(mean(nn.dist, na.rm = TRUE))
  
}