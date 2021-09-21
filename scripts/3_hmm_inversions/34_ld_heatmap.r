library(tidyverse)

##DO NOT CHANGE THE SCRIPT##
################
###get LD(r2)###
################

ld_heatmap <- function(aln.loc,ld,c){
  aln.loc.chr <- aln.loc %>% 
    filter(ref.chr==c) %>% 
    mutate(snp=paste0("tig_",tag,":",tag.pos)) %>% 
    arrange(ref.pos) %>% 
    rownames_to_column("index")
  
  ld.pos <- left_join(ld, aln.loc.chr, by=c("SNP_A"="snp")) %>% 
    dplyr::select(CHR_A=ref.chr, BP_A=index, SNP_A, CHR_B, BP_B, SNP_B, R2)
  
  ld.pos <- left_join(ld.pos, aln.loc.chr, by=c("SNP_B"="snp")) %>% 
    dplyr::select(CHR_A, BP_A, SNP_A, CHR_B=ref.chr, BP_B=index, SNP_B, R2)
  
  ld.pos$BP_A <- as.numeric(ld.pos$BP_A)
  ld.pos$BP_B <- as.numeric(ld.pos$BP_B)
  
  ld.pos %>%   
    ggplot(aes(x=BP_A, y=BP_B))+
    geom_tile(aes(fill=R2))+
    geom_abline(slope =1, intercept = 0)+
    theme_classic(base_size=15)+
    scale_fill_gradientn(colours=c("grey95","blue","red"), values=c(0,0.5,1), name="LD") +
    scale_x_continuous(expand=c(0.02,0)) +
    scale_y_continuous(expand=c(0.02,0)) +
    coord_fixed(ratio = 1) +
    xlab("SNP1 Position Index")+
    ylab("SNP2 Position Index")
  
}
