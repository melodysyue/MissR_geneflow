#calculate empirical NND between outliers on each chromosomes
calculate_nndist_all_chr <- function (stats.file, chr.target, num_permutations) {
  
  ## first, build null distributions of nndists for each linkage group:
  nnd.stats <- data.frame()
  
  for (j in chr.target){
    
    print(paste0("Processing Chr ", j, "..."))
    
    # subset for chr j
    stats.file.chr <- stats.file %>%
      filter(chr == j)
    
    # the number of outliers on that chr
    num.outliers <- stats.file.chr %>%
      filter(!is.na(chr.pos)) %>%
      select(ol) %>%
      unlist %>%
      sum(na.rm = TRUE)
    
    print(paste0("There are ", num.outliers, " outliers.")) 
    
    if (num.outliers > 1){
      
      
      # draw 10000 samples of num.outliers random loci, take the mean, and return the ecdf and mean
      null.mean.nnds <- replicate(num_permutations, calculate.null.nnd(stats.file.chr, num.outliers))
      
      
      # calculate the estimate mean null nndist
      null.mean <- mean(null.mean.nnds, na.rm = TRUE)
      null.ecdf <- ecdf(null.mean.nnds) #empirical cumulative distribution function
      
      # calculate the empirical nndist for real outliers
      empirical.mean.nnd <- calc_emp_nnd_dist(stats.file.chr)
      
      #number of total loci
      n.sites <- stats.file.chr %>% filter(!is.na(chr.pos)) %>% select(chr.pos) %>% unlist %>% length
      
      temp <- data.frame(chr = unique(stats.file.chr$chr),
                         n.sites = n.sites,
                         num.outliers = num.outliers,
                         nnd.mean.null = null.mean, 
                         nnd.sd.null = sd(null.mean.nnds, na.rm = TRUE),
                         nnd.mean.emp = empirical.mean.nnd,
                         nnd.emp.percentile = null.ecdf(empirical.mean.nnd),
                         nnd.emp.zscore = (empirical.mean.nnd - null.mean)/sd(null.mean.nnds, na.rm = TRUE),
                         nnd.emp.pvalue = two_side_p(null.mean.nnds, empirical.mean.nnd))
      
    }
    
    nnd.stats <- rbind(nnd.stats, temp)
  }
  
  return(nnd.stats)
}


